from __future__ import annotations

from typing import Optional, List, Tuple

from logging_utils import get_conversation_logger

from .models import GeneInput, GeneOutput, WarningMessage
from .web_utils import fetch_json, fetch_text


class GeneAgent:
    """
    Deterministic gene lookup agent.

    Responsibilities:
    1. Resolve a gene symbol in the requested species.
    2. Select a transcript (prefer canonical when available).
    3. Retrieve a real CDS sequence.
    4. Clean and minimally validate the sequence.
    5. Return protein name as best-effort metadata from UniProt.

    Primary source:
    - Ensembl REST

    Fallback:
    - NCBI E-utilities
    """

    def __init__(self) -> None:
        self._logger = get_conversation_logger()

    async def _ensembl_lookup(
        self,
        inp: GeneInput,
    ) -> Tuple[Optional[str], Optional[str], List[WarningMessage]]:
        species = inp.target_species
        symbol = inp.gene_symbol

        url = f"https://rest.ensembl.org/lookup/symbol/{species}/{symbol}"
        data = await fetch_json(
            url,
            params={"expand": "1"},
            headers={"Content-Type": "application/json"},
            label="ENSEMBL_LOOKUP",
        )

        warnings: List[WarningMessage] = []

        if not data:
            warnings.append(
                WarningMessage(
                    code="ENSEMBL_LOOKUP_FAILED",
                    message=f"Could not resolve gene symbol '{symbol}' for species '{species}'.",
                )
            )
            return None, None, warnings

        gene_id = data.get("id")

        canonical_id = data.get("canonical_transcript")
        transcripts = data.get("Transcript") or data.get("transcripts") or []

        transcript_id: Optional[str] = None

        if canonical_id:
            transcript_id = canonical_id
        elif transcripts:
            first_tx = transcripts[0]
            if isinstance(first_tx, dict):
                transcript_id = first_tx.get("id")

            warnings.append(
                WarningMessage(
                    code="CANONICAL_TRANSCRIPT_MISSING",
                    message=(
                        f"No canonical transcript reported by Ensembl for '{symbol}'. "
                        "Falling back to the first available transcript."
                    ),
                )
            )

        if not gene_id:
            warnings.append(
                WarningMessage(
                    code="GENE_ID_MISSING",
                    message=f"Ensembl resolved '{symbol}' but did not return a gene id.",
                )
            )

        if not transcript_id:
            warnings.append(
                WarningMessage(
                    code="NO_TRANSCRIPT",
                    message=(
                        f"No transcript could be selected for '{symbol}' in species '{species}'. "
                        "CDS retrieval will likely fail."
                    ),
                )
            )

        return gene_id, transcript_id, warnings

    def _clean_sequence(self, raw_seq: str) -> str:
        """
        Normalize text returned by sequence endpoints.

        Handles:
        - FASTA headers
        - newlines
        - whitespace
        - lowercase
        - unexpected stray characters
        """
        if not raw_seq:
            return ""

        lines = [line.strip() for line in raw_seq.splitlines() if line.strip()]
        lines = [line for line in lines if not line.startswith(">")]

        joined = "".join(lines).upper().replace(" ", "")
        cleaned = "".join(base for base in joined if base in {"A", "T", "G", "C", "N"})
        return cleaned

    def _validate_cds(self, seq: str) -> List[WarningMessage]:
        warnings: List[WarningMessage] = []

        if not seq:
            warnings.append(
                WarningMessage(
                    code="EMPTY_CDS",
                    message="Retrieved CDS is empty after cleaning.",
                )
            )
            return warnings

        allowed = {"A", "T", "G", "C", "N"}
        if any(base not in allowed for base in seq):
            warnings.append(
                WarningMessage(
                    code="INVALID_BASES",
                    message="CDS contains invalid nucleotide characters.",
                )
            )

        if len(seq) % 3 != 0:
            warnings.append(
                WarningMessage(
                    code="CDS_NOT_MULTIPLE_OF_3",
                    message="CDS length is not divisible by 3.",
                )
            )

        if not seq.startswith("ATG"):
            warnings.append(
                WarningMessage(
                    code="CDS_NO_START_CODON",
                    message="CDS does not start with ATG.",
                )
            )

        stop_codons = {"TAA", "TAG", "TGA"}
        if len(seq) >= 3 and seq[-3:] not in stop_codons:
            warnings.append(
                WarningMessage(
                    code="CDS_NO_STOP_CODON",
                    message="CDS does not end with a standard stop codon.",
                )
            )

        return warnings

    async def _ensembl_cds(
        self,
        transcript_id: Optional[str],
        species: Optional[str] = None,
    ) -> Tuple[Optional[str], List[WarningMessage]]:
        if not transcript_id:
            return None, [
                WarningMessage(
                    code="NO_TRANSCRIPT_ID",
                    message="No transcript id available for CDS retrieval.",
                )
            ]

        warnings: List[WarningMessage] = []

        candidates = [transcript_id]
        if "." in transcript_id:
            unversioned = transcript_id.split(".")[0]
            if unversioned not in candidates:
                candidates.append(unversioned)

        for candidate in candidates:
            param_sets = [
                {"type": "cds"},
                {"type": "cds", "object_type": "transcript"},
            ]

            if species:
                param_sets.append({"type": "cds", "species": species})
                param_sets.append(
                    {"type": "cds", "species": species, "object_type": "transcript"}
                )

            for params in param_sets:
                url = f"https://rest.ensembl.org/sequence/id/{candidate}"
                text = await fetch_text(
                    url,
                    params=params,
                    headers={
                        "Content-Type": "text/plain",
                        "User-Agent": "plasmid-pipeline/0.1",
                    },
                    label="ENSEMBL_SEQUENCE",
                )

                if not text:
                    continue

                cleaned = self._clean_sequence(text)
                if cleaned:
                    if candidate != transcript_id:
                        warnings.append(
                            WarningMessage(
                                code="ENSEMBL_USED_UNVERSIONED_TRANSCRIPT",
                                message=(
                                    f"CDS retrieval succeeded using unversioned transcript id "
                                    f"'{candidate}' instead of '{transcript_id}'."
                                ),
                            )
                        )
                    warnings.extend(self._validate_cds(cleaned))
                    return cleaned, warnings

        warnings.append(
            WarningMessage(
                code="CDS_FETCH_FAILED",
                message=f"Failed to retrieve CDS sequence for transcript '{transcript_id}'.",
            )
        )
        return None, warnings

    async def _ncbi_cds_fallback(
        self,
        inp: GeneInput,
        *,
        reason: str,
    ) -> Tuple[Optional[str], List[WarningMessage]]:
        warnings: List[WarningMessage] = [
            WarningMessage(
                code="NCBI_FALLBACK_TRIGGERED",
                message=f"Using NCBI fallback because: {reason}",
            )
        ]

        species_term_map = {
            "homo_sapiens": "Homo sapiens",
            "mus_musculus": "Mus musculus",
            "saccharomyces_cerevisiae": "Saccharomyces cerevisiae",
            "escherichia_coli": "Escherichia coli",
        }
        species_name = species_term_map.get(inp.target_species, inp.target_species)

        esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

        # Step 1: search gene record
        esearch = await fetch_json(
            esearch_url,
            params={
                "db": "gene",
                "term": f'{inp.gene_symbol}[Gene Name] AND "{species_name}"[Organism]',
                "retmode": "json",
                "retmax": 1,
            },
            label="NCBI_ESEARCH_GENE",
        )

        if not esearch:
            warnings.append(
                WarningMessage(
                    code="NCBI_ESEARCH_FAILED",
                    message=f"NCBI gene search failed for {inp.gene_symbol} in {species_name}.",
                )
            )
            return None, warnings

        idlist = ((esearch.get("esearchresult") or {}).get("idlist")) or []
        if not idlist:
            warnings.append(
                WarningMessage(
                    code="NCBI_GENE_NOT_FOUND",
                    message=f"NCBI found no gene record for {inp.gene_symbol} in {species_name}.",
                )
            )
            return None, warnings

        gene_id = idlist[0]

        # Step 2: gene summary (best-effort metadata)
        esummary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        esummary = await fetch_json(
            esummary_url,
            params={
                "db": "gene",
                "id": gene_id,
                "retmode": "json",
            },
            label="NCBI_ESUMMARY_GENE",
        )

        if not esummary:
            warnings.append(
                WarningMessage(
                    code="NCBI_ESUMMARY_FAILED",
                    message=f"NCBI gene summary failed for gene id {gene_id}.",
                )
            )
        else:
            result = esummary.get("result") or {}
            gene_doc = result.get(gene_id) or {}
            nomenclature = gene_doc.get("nomenclaturesymbol")
            if nomenclature and nomenclature != inp.gene_symbol:
                warnings.append(
                    WarningMessage(
                        code="NCBI_SYMBOL_NOTE",
                        message=f"NCBI returned nomenclature symbol '{nomenclature}'.",
                    )
                )

        # Step 3: search mRNA / transcript records directly in nuccore
        nuccore_search = await fetch_json(
            esearch_url,
            params={
                "db": "nuccore",
                "term": (
                    f'{inp.gene_symbol}[Gene] AND "{species_name}"[Organism] '
                    "AND biomol_mrna[PROP]"
                ),
                "retmode": "json",
                "retmax": 5,
            },
            label="NCBI_ESEARCH_NUCCORE",
        )

        nuccore_ids = ((nuccore_search or {}).get("esearchresult") or {}).get("idlist") or []
        if not nuccore_ids:
            warnings.append(
                WarningMessage(
                    code="NCBI_NUCCORE_NOT_FOUND",
                    message=f"NCBI nuccore search found no mRNA records for {inp.gene_symbol}.",
                )
            )
            return None, warnings

        # Step 4: fetch CDS FASTA from the candidate nuccore records
        efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fasta_text = await fetch_text(
            efetch_url,
            params={
                "db": "nuccore",
                "id": ",".join(nuccore_ids),
                "rettype": "fasta_cds_na",
                "retmode": "text",
            },
            label="NCBI_EFETCH_FASTA_CDS",
        )

        if not fasta_text:
            warnings.append(
                WarningMessage(
                    code="NCBI_EFETCH_FAILED",
                    message="NCBI efetch returned no FASTA CDS text.",
                )
            )
            return None, warnings

        blocks = [b.strip() for b in fasta_text.split(">") if b.strip()]
        if not blocks:
            warnings.append(
                WarningMessage(
                    code="NCBI_EMPTY_FASTA",
                    message="NCBI efetch returned FASTA text, but no FASTA records were parsed.",
                )
            )
            return None, warnings

        # Prefer a FASTA header that explicitly mentions the requested gene symbol
        for block in blocks:
            lines = block.splitlines()
            header = lines[0] if lines else ""
            seq = self._clean_sequence("\n".join(lines[1:]))

            if inp.gene_symbol.upper() in header.upper() and seq:
                warnings.append(
                    WarningMessage(
                        code="NCBI_FALLBACK_SUCCESS",
                        message=(
                            f"Recovered CDS from NCBI using nuccore fasta_cds_na for "
                            f"{inp.gene_symbol}."
                        ),
                    )
                )
                warnings.extend(self._validate_cds(seq))
                return seq, warnings

        # Fallback: first non-empty sequence
        for block in blocks:
            lines = block.splitlines()
            seq = self._clean_sequence("\n".join(lines[1:]))
            if seq:
                warnings.append(
                    WarningMessage(
                        code="NCBI_FALLBACK_FIRST_SEQUENCE_USED",
                        message=(
                            "NCBI returned multiple CDS entries; first non-empty sequence was used."
                        ),
                    )
                )
                warnings.extend(self._validate_cds(seq))
                return seq, warnings

        warnings.append(
            WarningMessage(
                code="NCBI_NO_USABLE_CDS",
                message="NCBI returned FASTA CDS data, but no usable sequence was parsed.",
            )
        )
        return None, warnings

    async def _uniprot_name(self, symbol: str, taxid: Optional[int]) -> Optional[str]:
        if not taxid:
            return None

        url = "https://rest.uniprot.org/uniprotkb/search"
        query = f"gene:{symbol} AND organism_id:{taxid}"

        data = await fetch_json(
            url,
            params={"query": query, "format": "json", "size": 1},
            label="UNIPROT_SEARCH",
        )
        if not data:
            return None

        results = data.get("results") or []
        if not results:
            return None

        result = results[0]
        desc = result.get("proteinDescription") or {}
        rec = desc.get("recommendedName") or {}
        full = rec.get("fullName") or {}
        value = full.get("value")

        if isinstance(value, str) and value.strip():
            return value.strip()

        return None

    def _infer_taxid(self, species: str) -> Optional[int]:
        mapping = {
            "homo_sapiens": 9606,
            "mus_musculus": 10090,
            "saccharomyces_cerevisiae": 559292,
            "escherichia_coli": 562,
        }
        return mapping.get(species)

    async def run(self, inp: GeneInput) -> GeneOutput:
        self._logger.info("[GENE] INPUT %s", inp.model_dump())

        warnings: List[WarningMessage] = []

        gene_id, transcript_id, lookup_warnings = await self._ensembl_lookup(inp)
        warnings.extend(lookup_warnings)

        cds_sequence, cds_warnings = await self._ensembl_cds(
            transcript_id,
            species=inp.target_species,
        )
        warnings.extend(cds_warnings)

        if cds_sequence is None:
            ncbi_cds, ncbi_warnings = await self._ncbi_cds_fallback(
                inp,
                reason="Ensembl CDS retrieval returned no usable sequence.",
            )
            cds_sequence = ncbi_cds
            warnings.extend(ncbi_warnings)

        taxid = self._infer_taxid(inp.target_species)
        protein_name = await self._uniprot_name(inp.gene_symbol, taxid)

        out = GeneOutput(
            gene_symbol=inp.gene_symbol,
            organism=inp.target_species,
            gene_id=gene_id,
            transcript_id=transcript_id,
            cds_sequence=cds_sequence,
            protein_name=protein_name,
            warnings=warnings,
        )

        self._logger.info(
            "[GENE] OUTPUT gene_id=%s transcript_id=%s cds_len=%s warnings=%s",
            out.gene_id,
            out.transcript_id,
            len(out.cds_sequence) if out.cds_sequence else None,
            [w.model_dump() for w in out.warnings],
        )

        return out