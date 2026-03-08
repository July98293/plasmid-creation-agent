from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Any, List, Optional

from logging_utils import get_conversation_logger

from .models import FeatureInput, FeatureOutput, ResolvedFeature, WarningMessage


class FeatureAgent:
    """
    Resolves promoters / tags / terminators from a local JSON library.
    """

    def __init__(self, library_path: Optional[Path] = None) -> None:
        self._logger = get_conversation_logger()
        if library_path is None:
            library_path = Path(__file__).with_name("feature_library.json")
        self._library_path = library_path
        self._library: Dict[str, Any] = self._load_library()

    def _load_library(self) -> Dict[str, Any]:
        if not self._library_path.exists():
            return {}
        with self._library_path.open("r", encoding="utf-8") as f:
            return json.load(f)

    def _resolve_single(
        self,
        name: Optional[str],
        category: str,
        *,
        position_hint: Optional[str] = None,
    ) -> tuple[Optional[ResolvedFeature], Optional[WarningMessage]]:
        if not name:
            return None, None
        coll = self._library.get(category, {})
        entry = coll.get(name)
        if not entry:
            return None, WarningMessage(
                code=f"{category.upper()}_NOT_FOUND",
                message=f"Feature '{name}' not found in local {category} library.",
            )
        feat_type = "promoter" if category == "promoters" else category.rstrip("s")
        seq = (entry.get("sequence") or "").strip().upper()
        if not seq:
            return None, WarningMessage(
                code=f"{category.upper()}_SEQUENCE_EMPTY",
                message=f"Feature '{name}' in {category} library has no sequence defined.",
            )
        feat = ResolvedFeature(
            name=name,
            type=feat_type,  # type: ignore[arg-type]
            sequence=seq,
            position_hint=position_hint,
        )
        return feat, None

    def run(self, inp: FeatureInput) -> FeatureOutput:
        self._logger.info("[FEATURE] INPUT %s", inp.model_dump())

        features: List[ResolvedFeature] = []
        warnings: List[WarningMessage] = []

        promoter, w = self._resolve_single(inp.promoter, "promoters")
        if promoter:
            features.append(promoter)
        if w:
            warnings.append(w)

        n_tag, w = self._resolve_single(inp.n_terminal_tag, "tags", position_hint="N")
        if n_tag:
            features.append(n_tag)
        if w:
            warnings.append(w)

        c_tag, w = self._resolve_single(inp.c_terminal_tag, "tags", position_hint="C")
        if c_tag:
            features.append(c_tag)
        if w:
            warnings.append(w)

        term, w = self._resolve_single(inp.terminator, "terminators")
        if term:
            features.append(term)
        if w:
            warnings.append(w)

        # Always add a simple default Kozak if host looks mammalian.
        if not any(f.type == "kozak" for f in features):
            kozak_entry = self._library.get("kozak", {}).get("standard_mammalian")
            if kozak_entry:
                features.append(
                    ResolvedFeature(
                        name="standard_mammalian",
                        type="kozak",
                        sequence=kozak_entry.get("sequence", ""),
                    )
                )

        out = FeatureOutput(features=features, warnings=warnings)
        self._logger.info("[FEATURE] OUTPUT %s", out.model_dump())
        return out

