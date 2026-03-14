from __future__ import annotations

from typing import Optional, List, Tuple, Dict, Any

from logging_utils import get_conversation_logger

from .models import GeneInput, GeneOutput, WarningMessage
from .web_utils import fetch_json, fetch_text

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:  # pragma: no cover
    SeqIO = None
    Seq = None


class GeneAgent:
    """
    Gene resolution agent.

    Responsibilities:
    1. Resolve a gene symbol in the requested species.
    2. Select a transcript (prefer canonical when available).
    3. Retrieve a real CDS sequence.
    4. If CDS is only available in a plasmid entry, extract only the CDS feature.
    5. Clean and validate the sequence.
    6. If the sequence still looks suspicious, pass it to an LLM cleanup hook.
    7. Return gene/protein metadata plus provenance.

    Primary sources:
    - Ensembl REST
    - NCBI E-utilities

    Optional fallback:
    - Addgene / plasmid-derived CDS extraction
    """

    def __init__(
        self,
        *,
        addgene_client: Optional[Any] = None,
        llm_client: Optional[Any] = None,
    ) -> None:
        self._logger = get_conversation_logger()
        self._addgene_client = addgene_client
        self._llm_client = llm_client

    # -------------------------------------------------------------------------
    # PUBLIC ENTRYPOINT
    # -------------------------------------------------------------------------

    async def run(self, inp: GeneInput) -> GeneOutput:
        self._logger.info("[GENE] INPUT %s", inp.model_dump())

        warnings: List[WarningMessage] = []
        gene_id: Optional[str] = None
        transcript_id: Optional[str] = None
        cds_sequence: Optional[str] = None
        protein_name: Optional[str] = None

        source_db: Optional[str] = None
        source_accession: Optional[str] = None
        source_method: Optional[str] = None
        was_extracted_from_plasmid = False
        extraction_coordinates: Optional[str] = None
        cleaning_method = "deterministic"

        # 1. Ensembl lookup
        gene_id, transcript_id, lookup_warnings = await self._ensembl_lookup(inp)
        warnings.extend(lookup_warnings)

        # 2. Ensembl CDS
        cds_sequence, cds_warnings = await self._ensembl_cds(
            transcript_id=transcript_id,
            species=inp.target_species,
        )
        warnings.extend(cds_warnings)

        if cds_sequence:
            source_db = "Ensembl"
            source_accession = transcript_id
            source_method = "ensembl_canonical_or_first_transcript"

        # 3. NCBI fallback
        if cds_sequence is None:
            ncbi_seq, ncbi_accession, ncbi_warnings = await self._ncbi_cds_fallback(
                inp,
                reason="Ensembl CDS retrieval returned no usable sequence.",
            )
            warnings.extend(ncbi_warnings)

            if ncbi_seq:
                cds_sequence = ncbi_seq
                source_db = "NCBI"
                source_accession = ncbi_accession
                source_method = "ncbi_nuccore_fasta_cds_na"

        # 4. Addgene / plasmid fallback
        if cds_sequence is None:
            addgene_seq, addgene_meta, addgene_warnings = await self._addgene_cds_fallback(inp)
            warnings.extend(addgene_warnings)

            if addgene_seq:
                cds_sequence = addgene_seq
                source_db = "Addgene"
                source_accession = addgene_meta.get("plasmid_id")
                source_method = "addgene_feature_extraction"
                was_extracted_from_plasmid = True
                extraction_coordinates = addgene_meta.get("coordinates")

        # 5. Deterministic cleaning + validation
        if cds_sequence is not None:
            cds_sequence = self._clean_sequence(cds_sequence)
            warnings.extend(self._validate_cds(cds_sequence))

        # 6. LLM cleanup hook if sequence still suspicious
        if cds_sequence is not None and self._sequence_suspicious(cds_sequence):
            llm_seq, llm_warnings = await self._llm_cleanup_if_needed(
                inp=inp,
                seq=cds_sequence,
                source_db=source_db,
                source_accession=source_accession,
                source_method=source_method,
            )
            warnings.extend(llm_warnings)

            if llm_seq:
                cds_sequence = llm_seq
                cleaning_method = "llm_cleanup"

        # 7. Protein metadata
        taxid = self._infer_taxid(inp.target_species)
        protein_name = await self._uniprot_name(inp.gene_symbol, taxid)

        # NOTE:
        # If your GeneOutput model does not yet support these extra fields,
        # add them there or remove them from this constructor.
        out = GeneOutput(
            gene_symbol=inp.gene_symbol,
            organism=inp.target_species,
            gene_id=gene_id,
            transcript_id=transcript_id,
            cds_sequence=cds_sequence,
            protein_name=protein_name,
            warnings=warnings,
            source_db=source_db,
            source_accession=source_accession,
            source_method=source_method,
            was_extracted_from_plasmid=was_extracted_from_plasmid,
            extraction_coordinates=extraction_coordinates,
            cleaning_method=cleaning_method,
        )

        self._logger.info(
            "[GENE] OUTPUT gene_id=%s transcript_id=%s cds_len=%s source_db=%s warnings=%s",
            out.gene_id,
            out.transcript_id,
            len(out.cds_sequence) if out.cds_sequence else None,
            getattr(out, "source_db", None),
            [w.model_dump() for w in out.warnings],
        )

        return out

    # -------------------------------------------------------------------------
    # ENSEMBL
    # -------------------------------------------------------------------------

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
                param_sets.append({"type": "cds", "species": species, "object_type": "transcript"})

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

    # -------------------------------------------------------------------------
    # NCBI
    # -------------------------------------------------------------------------

    async def _ncbi_cds_fallback(
        self,
        inp: GeneInput,
        *,
        reason: str,
    ) -> Tuple[Optional[str], Optional[str], List[WarningMessage]]:
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
        esummary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

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
            return None, None, warnings

        idlist = ((esearch.get("esearchresult") or {}).get("idlist")) or []
        if not idlist:
            warnings.append(
                WarningMessage(
                    code="NCBI_GENE_NOT_FOUND",
                    message=f"NCBI found no gene record for {inp.gene_symbol} in {species_name}.",
                )
            )
            return None, None, warnings

        gene_id = idlist[0]

        esummary = await fetch_json(
            esummary_url,
            params={
                "db": "gene",
                "id": gene_id,
                "retmode": "json",
            },
            label="NCBI_ESUMMARY_GENE",
        )

        if esummary:
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
        else:
            warnings.append(
                WarningMessage(
                    code="NCBI_ESUMMARY_FAILED",
                    message=f"NCBI gene summary failed for gene id {gene_id}.",
                )
            )

        nuccore_search = await fetch_json(
            esearch_url,
            params={
                "db": "nuccore",
                "term": (
                    f'{inp.gene_symbol}[Gene] AND "{species_name}"[Organism] '
                    "AND biomol_mrna[PROP]"
                ),
                "retmode": "json",
                "retmax": 10,
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
            return None, None, warnings

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
            return None, None, warnings

        blocks = [b.strip() for b in fasta_text.split(">") if b.strip()]
        if not blocks:
            warnings.append(
                WarningMessage(
                    code="NCBI_EMPTY_FASTA",
                    message="NCBI efetch returned FASTA text, but no FASTA records were parsed.",
                )
            )
            return None, None, warnings

        for block in blocks:
            lines = block.splitlines()
            header = lines[0] if lines else ""
            seq = self._clean_sequence("\n".join(lines[1:]))

            if inp.gene_symbol.upper() in header.upper() and seq:
                warnings.append(
                    WarningMessage(
                        code="NCBI_FALLBACK_SUCCESS",
                        message=f"Recovered CDS from NCBI for {inp.gene_symbol}.",
                    )
                )
                warnings.extend(self._validate_cds(seq))
                accession = header.split()[0] if header else None
                return seq, accession, warnings

        for block in blocks:
            lines = block.splitlines()
            header = lines[0] if lines else ""
            seq = self._clean_sequence("\n".join(lines[1:]))
            if seq:
                warnings.append(
                    WarningMessage(
                        code="NCBI_FALLBACK_FIRST_SEQUENCE_USED",
                        message="NCBI returned multiple CDS entries; first non-empty sequence was used.",
                    )
                )
                warnings.extend(self._validate_cds(seq))
                accession = header.split()[0] if header else None
                return seq, accession, warnings

        warnings.append(
            WarningMessage(
                code="NCBI_NO_USABLE_CDS",
                message="NCBI returned FASTA CDS data, but no usable sequence was parsed.",
            )
        )
        return None, None, warnings

    # -------------------------------------------------------------------------
    # ADDGENE / PLASMID EXTRACTION
    # -------------------------------------------------------------------------

    async def _addgene_cds_fallback(
        self,
        inp: GeneInput,
    ) -> Tuple[Optional[str], Dict[str, str], List[WarningMessage]]:
        warnings: List[WarningMessage] = []
        meta: Dict[str, str] = {}

        if self._addgene_client is None:
            warnings.append(
                WarningMessage(
                    code="ADDGENE_CLIENT_UNAVAILABLE",
                    message="Addgene fallback skipped because no Addgene client was configured.",
                )
            )
            return None, meta, warnings

        try:
            plasmids = await self._addgene_client.search_plasmids(inp.gene_symbol)
        except Exception as exc:
            warnings.append(
                WarningMessage(
                    code="ADDGENE_SEARCH_FAILED",
                    message=f"Addgene search failed: {exc}",
                )
            )
            return None, meta, warnings

        if not plasmids:
            warnings.append(
                WarningMessage(
                    code="ADDGENE_NO_PLASMIDS_FOUND",
                    message=f"No Addgene plasmids found for gene '{inp.gene_symbol}'.",
                )
            )
            return None, meta, warnings

        for plasmid in plasmids:
            plasmid_id = str(plasmid.get("id") or "")
            genbank_text = None

            try:
                genbank_text = await self._addgene_client.fetch_genbank(plasmid_id)
            except Exception as exc:
                warnings.append(
                    WarningMessage(
                        code="ADDGENE_GENBANK_FETCH_FAILED",
                        message=f"Failed to fetch GenBank for plasmid {plasmid_id}: {exc}",
                    )
                )
                continue

            if not genbank_text:
                continue

            seq, coords, extraction_warnings = self._extract_cds_from_genbank_text(
                genbank_text=genbank_text,
                gene_symbol=inp.gene_symbol,
            )
            warnings.extend(extraction_warnings)

            if seq:
                meta["plasmid_id"] = plasmid_id
                meta["coordinates"] = coords or ""
                warnings.append(
                    WarningMessage(
                        code="ADDGENE_FEATURE_EXTRACTION_SUCCESS",
                        message=f"Extracted CDS for {inp.gene_symbol} from Addgene plasmid {plasmid_id}.",
                    )
                )
                return seq, meta, warnings

        warnings.append(
            WarningMessage(
                code="ADDGENE_NO_USABLE_CDS",
                message=f"Addgene plasmids were found for {inp.gene_symbol}, but no usable CDS was extracted.",
            )
        )
        return None, meta, warnings

    def _extract_cds_from_genbank_text(
        self,
        *,
        genbank_text: str,
        gene_symbol: str,
    ) -> Tuple[Optional[str], Optional[str], List[WarningMessage]]:
        warnings: List[WarningMessage] = []

        if SeqIO is None:
            warnings.append(
                WarningMessage(
                    code="BIOPYTHON_NOT_INSTALLED",
                    message="Biopython is required for GenBank feature extraction.",
                )
            )
            return None, None, warnings

        from io import StringIO

        try:
            record = SeqIO.read(StringIO(genbank_text), "genbank")
        except Exception as exc:
            warnings.append(
                WarningMessage(
                    code="GENBANK_PARSE_FAILED",
                    message=f"Failed to parse GenBank text: {exc}",
                )
            )
            return None, None, warnings

        gene_symbol_upper = gene_symbol.upper()

        for feature in record.features:
            feature_type = (feature.type or "").lower()
            if feature_type not in {"cds", "gene", "misc_feature"}:
                continue

            qualifiers = feature.qualifiers or {}
            text_pool = []

            for key in ("gene", "label", "note", "product", "locus_tag"):
                values = qualifiers.get(key) or []
                if isinstance(values, str):
                    values = [values]
                text_pool.extend(values)

            joined = " | ".join(str(v) for v in text_pool).upper()

            if gene_symbol_upper not in joined:
                continue

            try:
                extracted = feature.extract(record.seq)
                seq = str(extracted).upper()
                seq = self._clean_sequence(seq)

                if not seq:
                    continue

                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                strand = feature.location.strand
                coords = f"{start}..{end} (strand={strand})"

                warnings.extend(self._validate_cds(seq))
                return seq, coords, warnings

            except Exception as exc:
                warnings.append(
                    WarningMessage(
                        code="FEATURE_EXTRACTION_FAILED",
                        message=f"Failed to extract CDS-like feature for {gene_symbol}: {exc}",
                    )
                )

        warnings.append(
            WarningMessage(
                code="NO_MATCHING_CDS_FEATURE",
                message=f"No matching CDS/gene feature found for '{gene_symbol}' in GenBank record.",
            )
        )
        return None, None, warnings

    # -------------------------------------------------------------------------
    # CLEANING / VALIDATION
    # -------------------------------------------------------------------------

    def _clean_sequence(self, raw_seq: str) -> str:
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

        if len(seq) < 60:
            warnings.append(
                WarningMessage(
                    code="CDS_TOO_SHORT",
                    message="CDS is unusually short.",
                )
            )

        n_fraction = seq.count("N") / len(seq)
        if n_fraction > 0.02:
            warnings.append(
                WarningMessage(
                    code="CDS_TOO_MANY_NS",
                    message="CDS contains too many ambiguous bases.",
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

    def _sequence_suspicious(self, seq: Optional[str]) -> bool:
        if not seq:
            return True
        if len(seq) < 60:
            return True
        if seq.count("N") / max(len(seq), 1) > 0.02:
            return True
        if len(seq) % 3 != 0:
            return True
        if not seq.startswith("ATG"):
            return True
        if seq[-3:] not in {"TAA", "TAG", "TGA"}:
            return True
        return False

    async def _llm_cleanup_if_needed(
        self,
        *,
        inp: GeneInput,
        seq: str,
        source_db: Optional[str],
        source_accession: Optional[str],
        source_method: Optional[str],
    ) -> Tuple[Optional[str], List[WarningMessage]]:
        warnings: List[WarningMessage] = []

        if not self._sequence_suspicious(seq):
            return seq, warnings

        warnings.append(
            WarningMessage(
                code="LLM_CLEANUP_TRIGGERED",
                message=(
                    f"Sequence for {inp.gene_symbol} looked suspicious; "
                    "LLM cleanup hook was triggered."
                ),
            )
        )

        if self._llm_client is None:
            warnings.append(
                WarningMessage(
                    code="LLM_CLIENT_UNAVAILABLE",
                    message="LLM cleanup was requested but no LLM client was configured.",
                )
            )
            return seq, warnings

        prompt = f"""
You are a DNA CDS cleanup assistant.

Gene symbol: {inp.gene_symbol}
Organism: {inp.target_species}
Source DB: {source_db or "unknown"}
Source accession: {source_accession or "unknown"}
Source method: {source_method or "unknown"}

Candidate sequence:
{seq}

Rules:
1. Keep only A,T,G,C,N.
2. Remove only obvious contamination.
3. Do not invent biology.
4. Do not reconstruct missing sequence.
5. Do not silently modify internal coding sequence.
6. Return JSON with:
   - cleaned_sequence
   - changed
   - issues_detected
   - confidence
   - rationale
""".strip()

        try:
            result = await self._llm_client.clean_sequence(prompt=prompt)
        except Exception as exc:
            warnings.append(
                WarningMessage(
                    code="LLM_CLEANUP_FAILED",
                    message=f"LLM cleanup failed: {exc}",
                )
            )
            return seq, warnings

        cleaned_sequence = None
        if isinstance(result, dict):
            cleaned_sequence = result.get("cleaned_sequence")

        cleaned_sequence = self._clean_sequence(cleaned_sequence or "")
        if not cleaned_sequence:
            warnings.append(
                WarningMessage(
                    code="LLM_RETURNED_EMPTY_SEQUENCE",
                    message="LLM cleanup returned no usable cleaned sequence.",
                )
            )
            return seq, warnings

        warnings.append(
            WarningMessage(
                code="LLM_CLEANUP_SUCCESS",
                message="LLM cleanup returned a cleaned sequence.",
            )
        )
        warnings.extend(self._validate_cds(cleaned_sequence))
        return cleaned_sequence, warnings

    # -------------------------------------------------------------------------
    # UNIPROT / TAXONOMY
    # -------------------------------------------------------------------------

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