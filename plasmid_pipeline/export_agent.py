from __future__ import annotations

from typing import Any, Dict, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from logging_utils import get_conversation_logger

from .models import ExportInput, ExportOutput


class ExportAgent:
    """
    Final export layer.

    Source of truth:
    - final plasmid sequence = assembly.assembled_sequence
    - insert feature coordinates are placed relative to assembly.junction_offset
    """

    def __init__(self) -> None:
        self._logger = get_conversation_logger()

    def run(self, inp: ExportInput) -> ExportOutput:
        self._logger.info("[EXPORT] INPUT %s", inp.model_dump())

        final_seq = inp.assembly.assembled_sequence
        backbone_len = len(inp.backbone.backbone_sequence)
        insert_len = len(inp.assembly.cassette_sequence)
        total_len = len(final_seq)
        junction_offset = inp.assembly.junction_offset

        seq_record = SeqRecord(
            Seq(final_seq),
            id=inp.intent.gene_symbol,
            name=f"{inp.intent.gene_symbol}_plasmid",
            description=f"{inp.intent.gene_symbol} expression plasmid",
        )

        features: List[Dict[str, Any]] = []

        if backbone_len:
            features.append(
                {
                    "name": inp.backbone.backbone_name,
                    "start": 1,
                    "end": backbone_len,
                    "type": "backbone",
                    "source": inp.backbone.source,
                }
            )

        for ann in inp.assembly.cassette_annotations:
            features.append(
                {
                    "name": ann.name,
                    "start": junction_offset + ann.start,
                    "end": junction_offset + ann.end,
                    "type": ann.type,
                }
            )

        genbank_record: Dict[str, Any] = {
            "id": seq_record.id,
            "name": seq_record.name,
            "description": seq_record.description,
            "length": total_len,
            "sequence": str(seq_record.seq),
            "features": features,
        }

        plasmid_map_data: Dict[str, Any] = {
            "backbone_length": backbone_len,
            "insert_length": insert_len,
            "total_length": total_len,
            "junction_offset": junction_offset,
            "features": features,
        }

        primer_list: List[Dict[str, Any]] = []
        for idx, note in enumerate(inp.assembly.primer_requirements, start=1):
            primer_list.append({"name": f"primer_{idx}", "notes": note})

        bill_of_materials: List[Dict[str, Any]] = [
            {
                "type": "backbone",
                "name": inp.backbone.backbone_name,
                "source": inp.backbone.source,
            },
            {
                "type": "insert",
                "features": inp.assembly.feature_order,
                "insertion_mode": "full_cassette",
            },
        ]

        summary_json: Dict[str, Any] = {
            "intent": inp.intent.model_dump(),
            "gene": inp.gene.model_dump(),
            "features": inp.features.model_dump(),
            "expression": inp.expression.model_dump(),
            "backbone": inp.backbone.model_dump(),
            "assembly": inp.assembly.model_dump(),
        }

        out = ExportOutput(
            summary_json=summary_json,
            genbank_record=genbank_record,
            plasmid_map_data=plasmid_map_data,
            primer_list=primer_list,
            bill_of_materials=bill_of_materials,
        )

        self._logger.info("[EXPORT] OUTPUT %s", out.model_dump())
        return out