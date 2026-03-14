from __future__ import annotations

import json
import os
import re
import urllib.parse
import urllib.request
from typing import Dict, Any, List, Optional, Tuple

from logging_utils import get_conversation_logger

from .addgene_client import AddgeneClient
from .models import BackboneInput, BackboneOutput, WarningMessage

try:
    from openai import OpenAI
except Exception:
    OpenAI = None


class BackboneAgent:
    """
    Backbone resolver with strict backbone-only behavior.

    Main principles:
    - Prefer true clean vector backbones.
    - Reject ambiguous "expression system" / patent-like records unless they can
      be safely cleaned.
    - If a record is loaded, try to extract ONLY the backbone by removing the
      payload/expression cassette span.
    - Use OpenAI only as a structured adjudicator for ambiguous feature
      classification. Do not let the model invent sequences.
    """

    NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    NCBI_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    SAFE_BACKBONE_MARKERS = {
        "bla", "amp", "ampr", "ampicillin", "betalactamase",
        "kan", "kana", "kanr", "kanamycin", "neo", "neor", "aph",
        "cat", "chloramphenicol", "tet", "tetracycline",
        "hygro", "hygromycin", "bsd", "blast", "zeo", "zeocin",
        "ori", "f1", "cole1", "puc", "pbr322", "rop", "sv40 ori",
        "mcs", "multiple cloning site", "polylinker",
        "lacz", "laczalpha", "rrn", "origin of replication",
        "antibiotic resistance", "selection marker"
    }

    PAYLOAD_MARKERS = {
        "gfp", "egfp", "sfgfp", "rfp", "cfp", "yfp", "bfp",
        "mcherry", "mscarlet", "mcardinal", "tdtomato", "venus",
        "reporter", "fluor", "fluorescent", "luciferase", "nanoluc",
        "fusion", "tag", "flag", "ha", "myc", "his", "strep", "v5",
        "cassette", "insert", "payload", "transgene",
        "nes", "nls", "devd", "tevp", "p2a", "t2a", "e2a", "ires",
        "cre", "cas9", "dcas9", "sensor", "biosensor",
        "tp53"
    }

    PROMOTER_MARKERS = {
        "cmv": "CMV",
        "ef1a": "EF1a",
        "ef1alpha": "EF1a",
        "cag": "CAG",
        "pgk": "PGK",
        "sv40": "SV40",
        "t7": "T7",
        "lac": "lac",
    }

    TERMINATOR_MARKERS = {
        "polya", "polyadenylation", "poly a",
        "bgh poly", "bgh polya",
        "sv40 poly", "sv40 polya",
        "terminator"
    }

    STRONG_BAD_TITLE_MARKERS = {
        "two plasmid",
        "expression system",
        "patent",
        "patent application",
        "jp ",
        "wo ",
        "us ",
        "therapeutic system",
        "gene delivery system",
    }

    CLEAN_VECTOR_TITLE_MARKERS = {
        "pcdna",
        "vector",
        "cloning vector",
        "expression vector",
        "backbone",
        "plasmid vector",
        "mcs",
    }

    def __init__(
        self,
        *,
        ncbi_api_key: Optional[str] = None,
        ncbi_email: Optional[str] = None,
        user_agent: str = "plasmid-pipeline/0.2",
        allow_clean_loaded_vectors: bool = True,
        openai_api_key: Optional[str] = None,
        openai_model: str = "gpt-5",
    ) -> None:
        self._logger = get_conversation_logger()
        self._ncbi_api_key = ncbi_api_key
        self._ncbi_email = ncbi_email
        self._user_agent = user_agent
        self._allow_clean_loaded_vectors = allow_clean_loaded_vectors
        self._openai_model = openai_model

        api_key = openai_api_key or os.getenv("OPENAI_API_KEY")
        self._openai_client = None
        if api_key and OpenAI is not None:
            try:
                self._openai_client = OpenAI(api_key=api_key)
            except Exception as e:
                self._logger.warning("[BACKBONE] OpenAI client init failed: %s", str(e))

    # ------------------------------------------------------------------
    # Basic helpers
    # ------------------------------------------------------------------

    def _clean_sequence(self, seq: Optional[str]) -> str:
        if not seq:
            return ""
        seq = seq.upper().replace("\n", "").replace("\r", "").replace(" ", "")
        return "".join(base for base in seq if base in {"A", "T", "G", "C", "N"})

    def _is_plausible_backbone_sequence(self, seq: str) -> bool:
        seq = self._clean_sequence(seq)
        return len(seq) >= 1500

    def _normalize_text(self, text: str) -> str:
        t = (text or "").lower()
        for ch in [",", ";", "(", ")", "[", "]", "{", "}", "/", "|", ":", "-", "_"]:
            t = t.replace(ch, " ")
        t = re.sub(r"\s+", " ", t).strip()
        return t

    def _http_get_json(self, url: str, params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        query = urllib.parse.urlencode(params, doseq=True)
        full_url = f"{url}?{query}"
        req = urllib.request.Request(
            full_url,
            headers={"User-Agent": self._user_agent, "Accept": "application/json"},
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as e:
            self._logger.warning("[BACKBONE] JSON GET failed url=%s error=%s", full_url, str(e))
            return None

    def _http_get_text(self, url: str, params: Dict[str, Any]) -> Optional[str]:
        query = urllib.parse.urlencode(params, doseq=True)
        full_url = f"{url}?{query}"
        req = urllib.request.Request(
            full_url,
            headers={"User-Agent": self._user_agent, "Accept": "text/plain"},
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as response:
                return response.read().decode("utf-8")
        except Exception as e:
            self._logger.warning("[BACKBONE] TEXT GET failed url=%s error=%s", full_url, str(e))
            return None

    def _ncbi_common_params(self) -> Dict[str, Any]:
        params: Dict[str, Any] = {}
        if self._ncbi_api_key:
            params["api_key"] = self._ncbi_api_key
        if self._ncbi_email:
            params["email"] = self._ncbi_email
        return params

    # ------------------------------------------------------------------
    # Query construction
    # ------------------------------------------------------------------

    def _build_search_terms(self, inp: BackboneInput) -> List[str]:
        host = inp.expression_host.lower()
        promoter = (inp.promoter or "").upper()
        vector_type = (inp.vector_type or "").lower()

        terms: List[str] = []

        if ("hek" in host or "293" in host or "mamm" in host or "cho" in host):
            terms.extend([
                'pcDNA[Title] AND vector[Title]',
                '"mammalian expression vector"[Title]',
                '"cloning vector"[Title]',
                'vector[Title] AND plasmid[Title] AND "complete sequence"',
            ])
            if promoter == "CMV":
                terms.extend([
                    '"CMV vector"[Title]',
                    '"pcDNA CMV"[Title]',
                ])
        elif ("coli" in host or "bacteria" in host or "e. coli" in host):
            terms.extend([
                'pET[Title] AND vector[Title]',
                '"bacterial expression vector"[Title]',
                '"cloning vector"[Title]',
            ])
            if promoter == "T7":
                terms.append('"T7 vector"[Title]')
        else:
            terms.extend([
                '"expression vector"[Title]',
                '"cloning vector"[Title]',
                'plasmid[Title] AND vector[Title] AND "complete sequence"',
            ])

        if vector_type:
            terms.append(f'"{vector_type}"[Title]')

        seen = set()
        out: List[str] = []
        for t in terms:
            if t not in seen:
                out.append(t)
                seen.add(t)
        return out

    # ------------------------------------------------------------------
    # NCBI retrieval
    # ------------------------------------------------------------------

    def _ncbi_esearch(self, term: str, retmax: int = 12) -> List[str]:
        data = self._http_get_json(
            self.NCBI_ESEARCH_URL,
            {
                "db": "nuccore",
                "term": term,
                "retmode": "json",
                "retmax": retmax,
                **self._ncbi_common_params(),
            },
        )
        if not data:
            return []
        return ((data.get("esearchresult") or {}).get("idlist")) or []

    def _ncbi_esummary(self, ids: List[str]) -> Dict[str, Any]:
        if not ids:
            return {}
        return self._http_get_json(
            self.NCBI_ESUMMARY_URL,
            {
                "db": "nuccore",
                "id": ",".join(ids),
                "retmode": "json",
                **self._ncbi_common_params(),
            },
        ) or {}

    def _extract_summary_docs(self, esummary: Dict[str, Any]) -> List[Tuple[str, Dict[str, Any]]]:
        result = esummary.get("result") or {}
        uids = result.get("uids") or []
        docs: List[Tuple[str, Dict[str, Any]]] = []
        for uid in uids:
            doc = result.get(uid)
            if isinstance(doc, dict):
                docs.append((uid, doc))
        return docs

    def _ncbi_fetch_genbank(self, nuccore_id: str) -> Optional[str]:
        return self._http_get_text(
            self.NCBI_EFETCH_URL,
            {
                "db": "nuccore",
                "id": nuccore_id,
                "rettype": "gbwithparts",
                "retmode": "text",
                **self._ncbi_common_params(),
            },
        )

    # ------------------------------------------------------------------
    # GenBank parsing
    # ------------------------------------------------------------------

    def _parse_genbank_sequence(self, gb_text: str) -> str:
        lines = gb_text.splitlines()
        in_origin = False
        seq_parts: List[str] = []

        for line in lines:
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin:
                if line.startswith("//"):
                    break
                seq_parts.append("".join(ch for ch in line if ch.isalpha()))

        return self._clean_sequence("".join(seq_parts))

    def _parse_genbank_features(self, gb_text: str) -> List[Dict[str, Any]]:
        lines = gb_text.splitlines()
        in_features = False
        features: List[Dict[str, Any]] = []
        current: Optional[Dict[str, Any]] = None
        current_qualifier: Optional[str] = None

        feature_key_re = re.compile(r"^\s{5}(\S+)\s+(.+)$")
        qualifier_re = re.compile(r'^\s+/([^=]+)(?:=(.*))?$')

        for line in lines:
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if in_features and line.startswith("ORIGIN"):
                break
            if not in_features:
                continue

            m = feature_key_re.match(line)
            if m:
                if current is not None:
                    features.append(current)
                current = {
                    "type": m.group(1).strip(),
                    "location": m.group(2).strip(),
                    "qualifiers": {},
                }
                current_qualifier = None
                continue

            if current is None:
                continue

            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith("/"):
                qm = qualifier_re.match(stripped)
                if not qm:
                    continue

                qkey = qm.group(1).strip()
                qval = qm.group(2)

                if qval is None:
                    value = ""
                else:
                    value = qval.strip()
                    if value.startswith('"') and value.endswith('"') and len(value) >= 2:
                        value = value[1:-1]

                current["qualifiers"].setdefault(qkey, []).append(value)
                current_qualifier = qkey
                continue

            if current_qualifier is not None:
                cont = stripped
                if cont.endswith('"'):
                    cont = cont[:-1]
                if cont.startswith('"'):
                    cont = cont[1:]

                existing = current["qualifiers"][current_qualifier]
                if existing:
                    existing[-1] = (existing[-1] + " " + cont).strip()
                else:
                    existing.append(cont)

        if current is not None:
            features.append(current)

        return features

    def _feature_text(self, feature: Dict[str, Any]) -> str:
        quals = feature.get("qualifiers", {})
        parts: List[str] = [feature.get("type", "")]
        for qvals in quals.values():
            for v in qvals:
                parts.append(v)
        return self._normalize_text(" ".join(parts))

    def _feature_hits(self, text: str, vocab: set[str]) -> List[str]:
        hits: List[str] = []
        for marker in vocab:
            if marker in text:
                hits.append(marker)
        return sorted(set(hits))

    def _location_bounds(self, location: str) -> Optional[Tuple[int, int]]:
        nums = [int(x) for x in re.findall(r"\d+", location or "")]
        if len(nums) < 2:
            return None
        return min(nums), max(nums)

    def _extract_promoters_from_features(self, features: List[Dict[str, Any]]) -> List[str]:
        promoters: List[str] = []
        for feat in features:
            text = self._feature_text(feat)
            if feat.get("type") == "promoter" or "promoter" in text:
                for k, label in self.PROMOTER_MARKERS.items():
                    if k in text and label not in promoters:
                        promoters.append(label)
        return promoters

    def _extract_terminator_like_features(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        out: List[Dict[str, Any]] = []
        for feat in features:
            text = self._feature_text(feat)
            if feat.get("type") in {"terminator", "polyA_signal", "regulatory"}:
                out.append(feat)
                continue
            if any(marker in text for marker in self.TERMINATOR_MARKERS):
                out.append(feat)
        return out

    # ------------------------------------------------------------------
    # Feature summarization for model adjudication
    # ------------------------------------------------------------------

    def _feature_to_summary(self, feat: Dict[str, Any], idx: int) -> Dict[str, Any]:
        bounds = self._location_bounds(feat.get("location", ""))
        text = self._feature_text(feat)
        return {
            "feature_index": idx,
            "type": feat.get("type"),
            "location": feat.get("location"),
            "bounds": bounds,
            "text": text[:600],
            "safe_hits": self._feature_hits(text, self.SAFE_BACKBONE_MARKERS),
            "payload_hits": self._feature_hits(text, self.PAYLOAD_MARKERS),
        }

    # ------------------------------------------------------------------
    # Annotation-based classification
    # ------------------------------------------------------------------

    def _classify_record(self, gb_text: str) -> Dict[str, Any]:
        features = self._parse_genbank_features(gb_text)
        promoters = self._extract_promoters_from_features(features)

        cds_total = 0
        non_backbone_cds = 0
        payload_markers: List[str] = []
        safe_feature_count = 0

        for feat in features:
            ftype = feat.get("type", "")
            text = self._feature_text(feat)
            safe_hits = self._feature_hits(text, self.SAFE_BACKBONE_MARKERS)
            payload_hits = self._feature_hits(text, self.PAYLOAD_MARKERS)

            if safe_hits:
                safe_feature_count += 1
            if payload_hits:
                payload_markers.extend(payload_hits)

            if ftype == "CDS":
                cds_total += 1
                if not safe_hits:
                    non_backbone_cds += 1

        payload_markers = sorted(set(payload_markers))

        is_loaded = False
        if payload_markers:
            is_loaded = True
        elif non_backbone_cds >= 1 and promoters:
            is_loaded = True
        elif non_backbone_cds >= 2:
            is_loaded = True

        return {
            "features": features,
            "promoters": promoters,
            "payload_markers": payload_markers,
            "cds_total": cds_total,
            "non_backbone_cds": non_backbone_cds,
            "safe_feature_count": safe_feature_count,
            "is_loaded": is_loaded,
        }

    # ------------------------------------------------------------------
    # OpenAI adjudication
    # ------------------------------------------------------------------

    def _llm_identify_backbone_payload_span(
        self,
        *,
        title: str,
        full_length: int,
        feature_summaries: List[Dict[str, Any]],
    ) -> Optional[Dict[str, Any]]:
        """
        Ask the model to identify which features belong to payload/expression cassette
        and which belong to the true backbone.

        The model never sees the full sequence and never invents sequence.
        It only classifies annotated features.
        """
        if self._openai_client is None:
            return None

        prompt = {
            "task": "Classify annotated GenBank features into backbone-core vs payload/expression cassette.",
            "rules": [
                "Backbone-core usually includes origin of replication, antibiotic marker, MCS/polylinker, lacZalpha, rop, bacterial maintenance features.",
                "Payload/expression cassette usually includes promoter, Kozak, CDS/transgene, tag, reporter, IRES/P2A/T2A, polyA/terminator for the transgene.",
                "Do NOT classify a feature as backbone-core just because it is in the plasmid.",
                "If the record appears to be a full mammalian expression plasmid, the promoter+gene+polyA region should usually be removable payload.",
                "Return the minimal removable contiguous span covering payload-related features.",
            ],
            "title": title,
            "full_length": full_length,
            "features": feature_summaries,
            "return_json_schema": {
                "payload_feature_indices": [0],
                "backbone_feature_indices": [1],
                "remove_span": {"start": 1, "end": 1000},
                "confidence": 0.0,
                "reason": "short explanation"
            }
        }

        try:
            resp = self._openai_client.responses.create(
                model=self._openai_model,
                input=[
                    {
                        "role": "system",
                        "content": (
                            "You are a plasmid map curator. "
                            "Classify only from provided annotations. "
                            "Never invent coordinates or sequences not implied by the data. "
                            "Return strict JSON only."
                        ),
                    },
                    {
                        "role": "user",
                        "content": json.dumps(prompt),
                    },
                ],
            )
            text = getattr(resp, "output_text", None)
            if not text:
                return None
            data = json.loads(text)
            if not isinstance(data, dict):
                return None
            return data
        except Exception as e:
            self._logger.warning("[BACKBONE] OpenAI adjudication failed: %s", str(e))
            return None

    # ------------------------------------------------------------------
    # Cleaning loaded vectors
    # ------------------------------------------------------------------

    def _find_cargo_span_heuristic(self, parsed: Dict[str, Any]) -> Optional[Tuple[int, int]]:
        features = parsed["features"]
        suspicious_bounds: List[Tuple[int, int]] = []

        terminator_like = self._extract_terminator_like_features(features)

        for feat in features:
            text = self._feature_text(feat)
            bounds = self._location_bounds(feat.get("location", ""))
            if bounds is None:
                continue

            ftype = feat.get("type", "")
            safe_hits = self._feature_hits(text, self.SAFE_BACKBONE_MARKERS)
            payload_hits = self._feature_hits(text, self.PAYLOAD_MARKERS)

            add = False
            if ftype == "promoter":
                add = True
            elif "kozak" in text:
                add = True
            elif ftype == "CDS" and not safe_hits:
                add = True
            elif payload_hits:
                add = True

            if add:
                suspicious_bounds.append(bounds)

        for feat in terminator_like:
            bounds = self._location_bounds(feat.get("location", ""))
            if bounds is not None:
                suspicious_bounds.append(bounds)

        if not suspicious_bounds:
            return None

        start = min(b[0] for b in suspicious_bounds)
        end = max(b[1] for b in suspicious_bounds)

        if end <= start:
            return None
        return start, end

    def _find_cargo_span_with_llm(
        self,
        *,
        title: str,
        seq: str,
        parsed: Dict[str, Any],
    ) -> Tuple[Optional[Tuple[int, int]], List[WarningMessage]]:
        warnings: List[WarningMessage] = []

        feature_summaries = [
            self._feature_to_summary(feat, i)
            for i, feat in enumerate(parsed["features"])
        ]

        llm = self._llm_identify_backbone_payload_span(
            title=title,
            full_length=len(seq),
            feature_summaries=feature_summaries,
        )

        if not llm:
            warnings.append(
                WarningMessage(
                    code="OPENAI_BACKBONE_ADJUDICATION_SKIPPED",
                    message="OpenAI adjudication was unavailable; falling back to heuristic cargo-span detection.",
                )
            )
            return self._find_cargo_span_heuristic(parsed), warnings

        remove_span = llm.get("remove_span")
        confidence = llm.get("confidence", 0.0)
        reason = llm.get("reason", "")

        if (
            not isinstance(remove_span, dict)
            or "start" not in remove_span
            or "end" not in remove_span
        ):
            warnings.append(
                WarningMessage(
                    code="OPENAI_BACKBONE_BAD_OUTPUT",
                    message="OpenAI adjudication returned invalid remove_span; using heuristic fallback.",
                )
            )
            return self._find_cargo_span_heuristic(parsed), warnings

        try:
            start_1 = int(remove_span["start"])
            end_1 = int(remove_span["end"])
        except Exception:
            warnings.append(
                WarningMessage(
                    code="OPENAI_BACKBONE_BAD_COORDS",
                    message="OpenAI adjudication returned non-integer coordinates; using heuristic fallback.",
                )
            )
            return self._find_cargo_span_heuristic(parsed), warnings

        if not (1 <= start_1 < end_1 <= len(seq)):
            warnings.append(
                WarningMessage(
                    code="OPENAI_BACKBONE_OUT_OF_RANGE",
                    message="OpenAI adjudication returned out-of-range cargo coordinates; using heuristic fallback.",
                )
            )
            return self._find_cargo_span_heuristic(parsed), warnings

        warnings.append(
            WarningMessage(
                code="OPENAI_BACKBONE_ADJUDICATION_USED",
                message=f"OpenAI identified removable payload span {start_1}..{end_1} with confidence={confidence}. Reason: {reason}",
            )
        )
        return (start_1, end_1), warnings

    def _clean_loaded_backbone(
        self,
        *,
        title: str,
        seq: str,
        parsed: Dict[str, Any],
    ) -> Tuple[Optional[str], List[WarningMessage], Optional[int]]:
        warnings: List[WarningMessage] = []

        cargo_span, llm_warnings = self._find_cargo_span_with_llm(
            title=title,
            seq=seq,
            parsed=parsed,
        )
        warnings.extend(llm_warnings)

        if cargo_span is None:
            warnings.append(
                WarningMessage(
                    code="LOADED_BACKBONE_CANNOT_BE_CLEANED",
                    message="The record looks loaded, but no safe removable payload span could be identified.",
                )
            )
            return None, warnings, None

        start_1, end_1 = cargo_span
        start_0 = max(0, start_1 - 1)
        end_0 = min(len(seq), end_1)

        left = seq[:start_0]
        right = seq[end_0:]
        cleaned = left + right

        if not self._is_plausible_backbone_sequence(cleaned):
            warnings.append(
                WarningMessage(
                    code="CLEANED_BACKBONE_TOO_SHORT",
                    message=f"Removing span {start_1}..{end_1} produced only {len(cleaned)} bp, too short to trust.",
                )
            )
            return None, warnings, None

        if len(cleaned) > len(seq) - 100:
            warnings.append(
                WarningMessage(
                    code="CLEANING_DID_NOT_REMOVE_ENOUGH",
                    message="The proposed cleaning removed too little sequence to plausibly strip an expression cassette.",
                )
            )
            return None, warnings, None

        insertion_site_index = len(left)
        warnings.append(
            WarningMessage(
                code="LOADED_BACKBONE_CLEANED",
                message=f"Converted loaded plasmid into putative backbone by removing payload span {start_1}..{end_1}.",
            )
        )
        return cleaned, warnings, insertion_site_index

    # ------------------------------------------------------------------
    # Candidate scoring and filtering
    # ------------------------------------------------------------------

    def _title_penalty(self, title: str) -> int:
        title_n = self._normalize_text(title)
        penalty = 0
        for marker in self.STRONG_BAD_TITLE_MARKERS:
            if marker in title_n:
                penalty -= 80
        if "system" in title_n:
            penalty -= 20
        if "patent" in title_n:
            penalty -= 50
        return penalty

    def _score_candidate_from_annotations(
        self,
        inp: BackboneInput,
        summary_doc: Dict[str, Any],
        parsed: Dict[str, Any],
        seq_len: int,
    ) -> int:
        score = 0

        title = self._normalize_text(str(summary_doc.get("title") or ""))
        host = inp.expression_host.lower()
        promoter = (inp.promoter or "").upper()

        if 2000 <= seq_len <= 12000:
            score += 30
        elif seq_len > 12000:
            score -= 20

        for marker in self.CLEAN_VECTOR_TITLE_MARKERS:
            if marker in title:
                score += 8

        score += self._title_penalty(title)

        promoters = parsed["promoters"]
        if promoter and promoter in promoters:
            # promoter on the backbone title/annotation is usually suspicious for a bare backbone
            score -= 10

        if "hek" in host or "293" in host or "mamm" in host or "cho" in host:
            if "pcdna" in title:
                score += 20
            if "mammalian" in title:
                score += 8
        elif "coli" in host or "bacteria" in host:
            if "pet" in title:
                score += 20
            if "bacterial" in title:
                score += 8

        score += 4 * parsed["safe_feature_count"]

        if parsed["is_loaded"]:
            score -= 120
        score -= 30 * parsed["non_backbone_cds"]
        score -= 8 * len(parsed["payload_markers"])

        return score

    def _is_hard_reject_title(self, title: str) -> bool:
        title_n = self._normalize_text(title)
        return any(marker in title_n for marker in self.STRONG_BAD_TITLE_MARKERS)

    # ------------------------------------------------------------------
    # Candidate retrieval
    # ------------------------------------------------------------------

    def _fetch_remote_backbone(self, inp: BackboneInput) -> Tuple[Optional[Dict[str, Any]], List[WarningMessage]]:
        warnings: List[WarningMessage] = []
        terms = self._build_search_terms(inp)

        best_clean_candidate: Optional[Dict[str, Any]] = None
        best_clean_score = -10**9

        best_cleaned_loaded_candidate: Optional[Dict[str, Any]] = None
        best_cleaned_loaded_score = -10**9

        for term in terms:
            ids = self._ncbi_esearch(term, retmax=12)
            if not ids:
                continue

            esummary = self._ncbi_esummary(ids)
            docs = self._extract_summary_docs(esummary)

            for uid, doc in docs:
                title = str(doc.get("title") or f"NCBI_nuccore_{uid}")

                gb_text = self._ncbi_fetch_genbank(uid)
                if not gb_text:
                    continue

                seq = self._parse_genbank_sequence(gb_text)
                if not self._is_plausible_backbone_sequence(seq):
                    continue

                parsed = self._classify_record(gb_text)
                score = self._score_candidate_from_annotations(
                    inp=inp,
                    summary_doc=doc,
                    parsed=parsed,
                    seq_len=len(seq),
                )

                base_notes = [
                    f"Resolved from NCBI nuccore record {uid}.",
                    f"Search term used: {term}",
                    f"Annotation-derived promoters: {parsed['promoters']}",
                    f"Annotation-derived payload markers: {parsed['payload_markers']}",
                    f"Annotation-derived non-backbone CDS count: {parsed['non_backbone_cds']}",
                    f"Annotation-derived safe backbone-feature count: {parsed['safe_feature_count']}",
                ]

                # prefer true clean vectors only if not obviously bad
                if not parsed["is_loaded"] and not self._is_hard_reject_title(title):
                    if score > best_clean_score:
                        best_clean_candidate = {
                            "name": title,
                            "sequence": seq,
                            "source": "ncbi_nuccore",
                            "accession": str(doc.get("caption") or uid),
                            "notes": base_notes,
                            "is_loaded_vector": False,
                            "backbone_promoters": parsed["promoters"],
                            "backbone_payload_markers": parsed["payload_markers"],
                            "suggested_strategy": "insert_into_clean_backbone",
                            "insertion_site_index": len(seq) // 2,
                        }
                        best_clean_score = score
                    continue

                if self._allow_clean_loaded_vectors:
                    cleaned_seq, clean_warnings, insertion_site_index = self._clean_loaded_backbone(
                        title=title,
                        seq=seq,
                        parsed=parsed,
                    )
                    if cleaned_seq and self._is_plausible_backbone_sequence(cleaned_seq):
                        cleaned_score = score + 70
                        candidate_notes = list(base_notes)
                        candidate_notes.append("Loaded record was cleaned by removing only the inferred payload/expression span.")
                        candidate = {
                            "name": title + " [backbone-only]",
                            "sequence": cleaned_seq,
                            "source": "ncbi_nuccore_cleaned",
                            "accession": str(doc.get("caption") or uid),
                            "notes": candidate_notes,
                            "is_loaded_vector": False,
                            "backbone_promoters": [],
                            "backbone_payload_markers": [],
                            "suggested_strategy": "insert_into_cleaned_backbone",
                            "insertion_site_index": insertion_site_index,
                            "extra_warnings": clean_warnings,
                        }
                        if cleaned_score > best_cleaned_loaded_score:
                            best_cleaned_loaded_candidate = candidate
                            best_cleaned_loaded_score = cleaned_score

        chosen: Optional[Dict[str, Any]] = None
        if best_clean_candidate is not None:
            chosen = best_clean_candidate
        elif best_cleaned_loaded_candidate is not None:
            chosen = best_cleaned_loaded_candidate

        if chosen is None:
            warnings.append(
                WarningMessage(
                    code="NCBI_BACKBONE_NOT_FOUND",
                    message="No suitable clean backbone could be resolved from NCBI, and no loaded candidate could be safely reduced to backbone-only.",
                )
            )
            return None, warnings

        warnings.append(
            WarningMessage(
                code="NCBI_BACKBONE_SUCCESS",
                message=f"Backbone resolved from NCBI: {chosen['name']}",
            )
        )

        for w in chosen.get("extra_warnings", []):
            warnings.append(w)

        return chosen, warnings

    # ------------------------------------------------------------------
    # Addgene-backed retrieval
    # ------------------------------------------------------------------

    def _fetch_addgene_backbone(self, inp: BackboneInput) -> Tuple[Optional[Dict[str, Any]], List[WarningMessage]]:
        warnings: List[WarningMessage] = []

        client = AddgeneClient(user_agent=self._user_agent)
        if not client.has_token():
            warnings.append(
                WarningMessage(
                    code="ADDGENE_TOKEN_MISSING",
                    message="ADDGENE_API_TOKEN is not set; skipping Addgene backbone lookup and falling back to NCBI.",
                )
            )
            return None, warnings

        promoters = inp.promoter or None
        vector_types = inp.vector_type or None

        plasmids = client.search_plasmids(
            promoters=promoters,
            vector_types=vector_types,
            page_size=10,
        )

        if not plasmids:
            warnings.append(
                WarningMessage(
                    code="ADDGENE_NO_PLASMIDS",
                    message="Addgene catalog search returned no plasmids for the requested context.",
                )
            )
            return None, warnings

        best_clean_candidate: Optional[Dict[str, Any]] = None
        best_clean_score = -10**9

        best_cleaned_loaded_candidate: Optional[Dict[str, Any]] = None
        best_cleaned_loaded_score = -10**9

        for doc in plasmids:
            plasmid_id = doc.get("id")
            if plasmid_id is None:
                continue

            plasmid = client.get_plasmid_with_sequences(int(plasmid_id))
            if not plasmid:
                continue

            gb_text = client.download_first_full_genbank(plasmid)
            if not gb_text:
                continue

            seq = self._parse_genbank_sequence(gb_text)
            if not self._is_plausible_backbone_sequence(seq):
                continue

            parsed = self._classify_record(gb_text)
            title = str(plasmid.get("name") or f"Addgene_plasmid_{plasmid_id}")

            synthetic_summary = {"title": title}
            score = self._score_candidate_from_annotations(
                inp=inp,
                summary_doc=synthetic_summary,
                parsed=parsed,
                seq_len=len(seq),
            )

            base_notes = [
                f"Resolved from Addgene plasmid {plasmid_id}.",
                f"Addgene name: {plasmid.get('name')}",
                f"Addgene description: {plasmid.get('description')}",
                f"Annotation-derived promoters: {parsed['promoters']}",
                f"Annotation-derived payload markers: {parsed['payload_markers']}",
                f"Annotation-derived non-backbone CDS count: {parsed['non_backbone_cds']}",
                f"Annotation-derived safe backbone-feature count: {parsed['safe_feature_count']}",
            ]

            if not parsed["is_loaded"] and not self._is_hard_reject_title(title):
                if score > best_clean_score:
                    best_clean_candidate = {
                        "name": title,
                        "sequence": seq,
                        "source": "addgene_catalog",
                        "accession": str(plasmid_id),
                        "notes": base_notes,
                        "is_loaded_vector": False,
                        "backbone_promoters": parsed["promoters"],
                        "backbone_payload_markers": parsed["payload_markers"],
                        "suggested_strategy": "insert_into_clean_backbone",
                        "insertion_site_index": len(seq) // 2,
                    }
                    best_clean_score = score
                continue

            if self._allow_clean_loaded_vectors:
                cleaned_seq, clean_warnings, insertion_site_index = self._clean_loaded_backbone(
                    title=title,
                    seq=seq,
                    parsed=parsed,
                )
                if cleaned_seq and self._is_plausible_backbone_sequence(cleaned_seq):
                    cleaned_score = score + 70
                    candidate_notes = list(base_notes)
                    candidate_notes.append("Loaded Addgene plasmid was reduced to backbone-only by removing the inferred payload span.")
                    candidate = {
                        "name": title + " [backbone-only]",
                        "sequence": cleaned_seq,
                        "source": "addgene_catalog_cleaned",
                        "accession": str(plasmid_id),
                        "notes": candidate_notes,
                        "is_loaded_vector": False,
                        "backbone_promoters": [],
                        "backbone_payload_markers": [],
                        "suggested_strategy": "insert_into_cleaned_backbone",
                        "insertion_site_index": insertion_site_index,
                        "extra_warnings": clean_warnings,
                    }
                    if cleaned_score > best_cleaned_loaded_score:
                        best_cleaned_loaded_candidate = candidate
                        best_cleaned_loaded_score = cleaned_score

        chosen: Optional[Dict[str, Any]] = None
        if best_clean_candidate is not None:
            chosen = best_clean_candidate
        elif best_cleaned_loaded_candidate is not None:
            chosen = best_cleaned_loaded_candidate

        if chosen is None:
            warnings.append(
                WarningMessage(
                    code="ADDGENE_BACKBONE_NOT_FOUND",
                    message="Addgene catalog did not yield a usable clean or safely cleaned backbone-only candidate.",
                )
            )
            return None, warnings

        warnings.append(
            WarningMessage(
                code="ADDGENE_BACKBONE_SUCCESS",
                message=f"Backbone resolved from Addgene: {chosen['name']}",
            )
        )

        for w in chosen.get("extra_warnings", []):
            warnings.append(w)

        return chosen, warnings

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self, inp: BackboneInput) -> BackboneOutput:
        self._logger.info("[BACKBONE] INPUT %s", inp.model_dump())
        warnings: List[WarningMessage] = []

        chosen, addgene_warnings = self._fetch_addgene_backbone(inp)
        warnings.extend(addgene_warnings)

        if chosen is None:
            ncbi_candidate, ncbi_warnings = self._fetch_remote_backbone(inp)
            warnings.extend(ncbi_warnings)
            chosen = ncbi_candidate

        if chosen is None:
            raise ValueError(
                "No validated clean backbone-only sequence could be resolved from Addgene or NCBI."
            )

        seq = self._clean_sequence(chosen.get("sequence"))
        if not self._is_plausible_backbone_sequence(seq):
            raise ValueError(
                f"Resolved backbone '{chosen.get('name', 'unknown')}' is too short ({len(seq)} bp) to be trusted."
            )

        notes = list(chosen.get("notes", []))
        if inp.promoter:
            notes.append(f"Requested promoter={inp.promoter}.")
        if inp.expression_host:
            notes.append(f"Requested expression_host={inp.expression_host}.")
        if inp.vector_type:
            notes.append(f"Requested vector_type={inp.vector_type}.")

        out = BackboneOutput(
            backbone_name=str(chosen.get("name") or "resolved_backbone"),
            backbone_sequence=seq,
            source=str(chosen.get("source") or "unknown"),
            notes=notes,
            warnings=warnings,
            is_loaded_vector=False,
            backbone_promoters=list(chosen.get("backbone_promoters", [])),
            backbone_payload_markers=list(chosen.get("backbone_payload_markers", [])),
            suggested_strategy=chosen.get("suggested_strategy"),
            insertion_site_index=chosen.get("insertion_site_index"),
        )

        self._logger.info(
            "[BACKBONE] OUTPUT name=%s source=%s length=%s warnings=%s",
            out.backbone_name,
            out.source,
            len(out.backbone_sequence),
            [w.model_dump() for w in out.warnings],
        )
        return out