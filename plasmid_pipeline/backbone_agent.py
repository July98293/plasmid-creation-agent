"""
BackboneAgent — single responsibility: return a clean, empty vector scaffold.

Strategy:
1. Build a prioritised list of specific known-good vector names based on
   expression_host / vector_type.
2. Search NCBI nuccore for each name + "complete sequence".
3. Fetch the GenBank record, parse the sequence.
4. If the record appears loaded (contains payload markers) → skip entirely.
5. Score by size fit and name match; return the best passing candidate.

No LLM calls. No promoter logic. No payload cleaning.
"""
from __future__ import annotations

import json
import os
import re
import time
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Tuple

from logging_utils import get_conversation_logger

from .models import BackboneInput, BackboneOutput, WarningMessage


# ---------------------------------------------------------------------------
# Well-known empty scaffold names, ordered by preference per context
# ---------------------------------------------------------------------------

_MAMMALIAN_VECTORS = ["pcDNA3.1", "pcDNA3.3", "pcDNA3", "pCMV-Script"]
_LENTIVIRAL_VECTORS = ["pLenti-CMV", "pLenti", "pLKO.1", "pLKO"]
_BACTERIAL_VECTORS = ["pUC19", "pBR322", "pET-28a", "pET-28"]
_GENERIC_VECTORS = ["pUC19", "pcDNA3.1"]

# Sequences outside this range are deprioritised
_IDEAL_MIN_BP = 3_000
_IDEAL_MAX_BP = 10_000

# If a GenBank record contains ≥2 of these it is treated as loaded
_PAYLOAD_MARKERS = {
    "gfp", "egfp", "rfp", "mcherry", "luciferase", "nanoluc",
    "reporter", "transgene", "insert", "payload", "cassette",
    "cas9", "cre", "p2a", "t2a", "ires",
}


class BackboneAgent:
    NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    NCBI_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    def __init__(
        self,
        *,
        ncbi_api_key: Optional[str] = None,
        ncbi_email: Optional[str] = None,
        user_agent: str = "plasmid-pipeline/0.2",
    ) -> None:
        self._logger = get_conversation_logger()
        self._ncbi_api_key = ncbi_api_key or os.getenv("NCBI_API_KEY", "")
        self._ncbi_email = ncbi_email or os.getenv("NCBI_EMAIL", "")
        self._user_agent = user_agent

    # ------------------------------------------------------------------
    # HTTP helpers
    # ------------------------------------------------------------------

    def _http_get_json(self, url: str, params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        time.sleep(0.5)
        full_url = f"{url}?{urllib.parse.urlencode(params)}"
        print(f"[BACKBONE HTTP] GET JSON → {full_url}", flush=True)
        req = urllib.request.Request(
            full_url,
            headers={"User-Agent": self._user_agent, "Accept": "application/json"},
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                print("[BACKBONE HTTP] ✓ JSON received", flush=True)
                return data
        except Exception as e:
            print(f"[BACKBONE HTTP] ✗ {e}", flush=True)
            return None

    def _http_get_text(self, url: str, params: Dict[str, Any]) -> Optional[str]:
        time.sleep(0.5)
        full_url = f"{url}?{urllib.parse.urlencode(params)}"
        print(f"[BACKBONE HTTP] GET TEXT → {full_url}", flush=True)
        req = urllib.request.Request(
            full_url,
            headers={"User-Agent": self._user_agent, "Accept": "text/plain"},
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                text = resp.read().decode("utf-8")
                print(f"[BACKBONE HTTP] ✓ TEXT received ({len(text)} bytes)", flush=True)
                return text
        except Exception as e:
            print(f"[BACKBONE HTTP] ✗ {e}", flush=True)
            return None

    # ------------------------------------------------------------------
    # Target vector selection
    # ------------------------------------------------------------------

    def _default_vectors(self, inp: BackboneInput) -> List[str]:
        host = (inp.expression_host or "").lower()
        vtype = (inp.vector_type or "").lower()

        if "lenti" in vtype:
            return _LENTIVIRAL_VECTORS
        if any(x in host for x in ["hek", "293", "cho", "cos", "hela", "mamm", "human", "mouse"]):
            return _MAMMALIAN_VECTORS
        if any(x in host for x in ["coli", "bacteria", "e.coli", "bl21"]):
            return _BACTERIAL_VECTORS
        return _GENERIC_VECTORS

    def _target_vectors(self, inp: BackboneInput) -> List[str]:
        if inp.backbone:
            return [inp.backbone] + self._default_vectors(inp)
        return self._default_vectors(inp)

    # ------------------------------------------------------------------
    # NCBI
    # ------------------------------------------------------------------

    def _ncbi_params(self) -> Dict[str, Any]:
        p: Dict[str, Any] = {}
        if self._ncbi_api_key:
            p["api_key"] = self._ncbi_api_key
        if self._ncbi_email:
            p["email"] = self._ncbi_email
        return p

    def _ncbi_search(self, vector_name: str) -> List[str]:
        data = self._http_get_json(self.NCBI_ESEARCH_URL, {
            "db": "nuccore",
            "term": f'"{vector_name}"[Title] AND "complete sequence"[Title]',
            "retmode": "json",
            "retmax": 5,
            **self._ncbi_params(),
        })
        return ((data or {}).get("esearchresult") or {}).get("idlist") or []

    def _ncbi_title(self, uid: str) -> str:
        data = self._http_get_json(self.NCBI_ESUMMARY_URL, {
            "db": "nuccore",
            "id": uid,
            "retmode": "json",
            **self._ncbi_params(),
        })
        result = (data or {}).get("result") or {}
        doc = result.get(uid) or {}
        return str(doc.get("title") or uid)

    def _ncbi_fetch_genbank(self, uid: str) -> Optional[str]:
        return self._http_get_text(self.NCBI_EFETCH_URL, {
            "db": "nuccore",
            "id": uid,
            "rettype": "gbwithparts",
            "retmode": "text",
            **self._ncbi_params(),
        })

    # ------------------------------------------------------------------
    # Sequence parsing
    # ------------------------------------------------------------------

    def _parse_sequence(self, gb_text: str) -> str:
        bases: List[str] = []
        in_origin = False
        for line in gb_text.splitlines():
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin:
                if line.startswith("//"):
                    break
                bases.append(re.sub(r"[^acgtACGT]", "", line))
        return "".join(bases).upper()

    def _is_loaded(self, gb_text: str) -> bool:
        """Return True if the record appears to carry a payload/expression insert."""
        text_lower = gb_text.lower()
        hits = sum(1 for m in _PAYLOAD_MARKERS if m in text_lower)
        return hits >= 2

    # ------------------------------------------------------------------
    # Scoring
    # ------------------------------------------------------------------

    def _score(self, seq_len: int, vector_name: str, title: str) -> int:
        score = 0
        if _IDEAL_MIN_BP <= seq_len <= _IDEAL_MAX_BP:
            score += 50
        elif seq_len < _IDEAL_MIN_BP:
            score -= 30
        else:
            score -= 10
        if vector_name.lower() in title.lower():
            score += 100
        return score

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self, inp: BackboneInput) -> BackboneOutput:
        self._logger.info("[BACKBONE] INPUT %s", inp.model_dump())
        warnings: List[WarningMessage] = []

        targets = self._target_vectors(inp)
        print(f"[BACKBONE] host={inp.expression_host!r} → target vectors: {targets}", flush=True)

        best: Optional[Dict[str, Any]] = None
        best_score = -(10 ** 9)

        for vector_name in targets:
            print(f"[BACKBONE] Searching NCBI for {vector_name!r}", flush=True)
            ids = self._ncbi_search(vector_name)
            if not ids:
                print(f"[BACKBONE] No results", flush=True)
                continue

            for uid in ids:
                title = self._ncbi_title(uid)
                print(f"[BACKBONE] uid={uid} title={title!r}", flush=True)

                gb_text = self._ncbi_fetch_genbank(uid)
                if not gb_text:
                    continue

                seq = self._parse_sequence(gb_text)
                if len(seq) < 1500:
                    print(f"[BACKBONE] uid={uid} skipped — too short ({len(seq)} bp)", flush=True)
                    continue

                if self._is_loaded(gb_text):
                    print(f"[BACKBONE] uid={uid} skipped — appears loaded", flush=True)
                    continue

                score = self._score(len(seq), vector_name, title)
                print(f"[BACKBONE] uid={uid} seq={len(seq)}bp score={score}", flush=True)

                if score > best_score:
                    best_score = score
                    best = {"name": title, "sequence": seq, "source": f"NCBI:{uid}", "genbank": gb_text}

            if best is not None:
                print(f"[BACKBONE] Candidate found after {vector_name!r} — stopping", flush=True)
                break

        if best is None:
            raise ValueError(
                f"No clean empty backbone found for expression_host={inp.expression_host!r}. "
                "All candidates were either absent, too short, or appeared loaded."
            )

        seq = best["sequence"]
        out = BackboneOutput(
            backbone_name=best["name"],
            backbone_sequence=seq,
            backbone_length=len(seq),
            source=best["source"],
            backbone_genbank=best["genbank"],
            warnings=warnings,
        )
        self._logger.info(
            "[BACKBONE] OUTPUT name=%s source=%s length=%d",
            out.backbone_name, out.source, out.backbone_length,
        )
        return out
