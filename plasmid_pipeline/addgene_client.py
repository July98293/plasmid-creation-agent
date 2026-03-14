from __future__ import annotations

import json
import os
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional

from logging_utils import get_conversation_logger


ADDGENE_BASE_URL = "https://api.developers.addgene.org"


class AddgeneClient:
    """
    Minimal synchronous Addgene client.

    Notes:
    - Reads the API token from the ADDGENE_API_TOKEN environment variable.
    - Only implements the small subset of endpoints we need:
      - /catalog/plasmid/ (search)
      - /catalog/plasmid-with-sequences/{id}/ (details + sequences)
    """

    def __init__(
        self,
        *,
        api_token: Optional[str] = None,
        user_agent: str = "plasmid-pipeline/0.1",
    ) -> None:
        self._logger = get_conversation_logger()
        self._api_token = api_token or os.getenv("ADDGENE_API_TOKEN") or ""
        self._user_agent = user_agent

    # ------------------------------------------------------------------
    # Low-level HTTP helpers
    # ------------------------------------------------------------------

    def _headers_json(self) -> Dict[str, str]:
        headers: Dict[str, str] = {
            "User-Agent": self._user_agent,
            "Accept": "application/json",
        }
        if self._api_token:
            headers["Authorization"] = f"Token {self._api_token}"
        return headers

    def _headers_text(self) -> Dict[str, str]:
        headers: Dict[str, str] = {
            "User-Agent": self._user_agent,
            "Accept": "text/plain",
        }
        if self._api_token:
            headers["Authorization"] = f"Token {self._api_token}"
        return headers

    def _http_get_json(self, path: str, params: Optional[Dict[str, Any]] = None) -> Optional[Dict[str, Any]]:
        url = f"{ADDGENE_BASE_URL}{path}"
        query = urllib.parse.urlencode(params or {}, doseq=True)
        if query:
            url = f"{url}?{query}"

        req = urllib.request.Request(url, headers=self._headers_json())

        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                text = resp.read().decode("utf-8")
        except Exception as exc:  # noqa: BLE001
            self._logger.warning("[ADDGENE] JSON GET failed url=%s error=%s", url, repr(exc))
            return None

        try:
            return json.loads(text)
        except Exception as exc:  # noqa: BLE001
            self._logger.warning("[ADDGENE] JSON decode failed url=%s error=%s", url, repr(exc))
            return None

    def _http_get_text_raw(self, url: str) -> Optional[str]:
        req = urllib.request.Request(url, headers=self._headers_text())
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read().decode("utf-8")
        except Exception as exc:  # noqa: BLE001
            self._logger.warning("[ADDGENE] TEXT GET failed url=%s error=%s", url, repr(exc))
            return None

    # ------------------------------------------------------------------
    # High-level helpers
    # ------------------------------------------------------------------

    def has_token(self) -> bool:
        return bool(self._api_token)

    def search_plasmids(
        self,
        *,
        name: Optional[str] = None,
        genes: Optional[str] = None,
        promoters: Optional[str] = None,
        backbone: Optional[str] = None,
        resistance_marker: Optional[str] = None,
        vector_types: Optional[str] = None,
        tags: Optional[str] = None,
        page_size: int = 10,
    ) -> List[Dict[str, Any]]:
        """
        Lightweight wrapper over /catalog/plasmid/.

        Only exposes the filters that are most relevant for backbone / cassette
        reasoning. The caller is responsible for applying any scoring.
        """
        params: Dict[str, Any] = {
            "page_size": page_size,
            "sort_by": "id",
        }

        if promoters:
            params["promoters"] = promoters
        if name:
            params["name"] = name
        if genes:
            params["genes"] = genes
        if backbone:
            params["backbone"] = backbone
        if resistance_marker:
            params["resistance_marker"] = resistance_marker
        if vector_types:
            params["vector_types"] = vector_types
        if tags:
            params["tags"] = tags

        data = self._http_get_json("/catalog/plasmid/", params=params)
        if not data:
            return []

        results = data.get("results") or []
        if not isinstance(results, list):
            return []
        return [r for r in results if isinstance(r, dict)]

    def get_plasmid_with_sequences(self, plasmid_id: int) -> Optional[Dict[str, Any]]:
        """
        Wrapper over /catalog/plasmid-with-sequences/{id}/.
        """
        data = self._http_get_json(f"/catalog/plasmid-with-sequences/{plasmid_id}/")
        if not data:
            return None
        return data

    def download_first_full_genbank(self, plasmid_with_sequences: Dict[str, Any]) -> Optional[str]:
        """
        Given a PlasmidWithSequences document, try to download a single
        GenBank flatfile corresponding to a full sequence.
        """
        sequences = (plasmid_with_sequences.get("sequences") or {}) if isinstance(plasmid_with_sequences, dict) else {}

        # Prefer Addgene full sequences, then user full sequences, then partials.
        priority_keys = [
            "public_addgene_full_sequences",
            "public_user_full_sequences",
            "public_addgene_partial_sequences",
            "public_user_partial_sequences",
        ]

        for key in priority_keys:
            seq_list = sequences.get(key) or []
            for seq in seq_list:
                if not isinstance(seq, dict):
                    continue
                genbank_url = seq.get("genbank_url")
                if not genbank_url:
                    continue
                gb_text = self._http_get_text_raw(genbank_url)
                if gb_text:
                    return gb_text

        return None

