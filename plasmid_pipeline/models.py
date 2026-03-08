from __future__ import annotations

from typing import List, Optional, Literal, Any, Dict

from pydantic import BaseModel, Field


# ----------------------------
# Core shared types
# ----------------------------


class WarningMessage(BaseModel):
    code: str
    message: str


# ----------------------------
# Intent agent
# ----------------------------


class IntentInput(BaseModel):
    user_request: str


class IntentOutput(BaseModel):
    gene_symbol: str
    target_species: str
    expression_host: str
    promoter: Optional[str] = None
    assembly_method: Optional[str] = None
    n_terminal_tag: Optional[str] = None
    c_terminal_tag: Optional[str] = None
    terminator: Optional[str] = None


# ----------------------------
# Gene agent
# ----------------------------


class GeneInput(BaseModel):
    gene_symbol: str
    target_species: str


class GeneOutput(BaseModel):
    gene_symbol: str
    organism: str
    gene_id: Optional[str] = None
    transcript_id: Optional[str] = None
    cds_sequence: Optional[str] = None
    protein_name: Optional[str] = None
    warnings: List[WarningMessage] = Field(default_factory=list)


# ----------------------------
# Feature agent
# ----------------------------


class FeatureInput(BaseModel):
    promoter: Optional[str]
    n_terminal_tag: Optional[str]
    c_terminal_tag: Optional[str]
    terminator: Optional[str]


class ResolvedFeature(BaseModel):
    name: str
    type: Literal["promoter", "kozak", "tag", "terminator"]
    sequence: str
    position_hint: Optional[Literal["N", "C"]] = None


class FeatureOutput(BaseModel):
    features: List[ResolvedFeature] = Field(default_factory=list)
    warnings: List[WarningMessage] = Field(default_factory=list)


# ----------------------------
# Expression agent
# ----------------------------


class ExpressionInput(BaseModel):
    intent: IntentOutput
    gene: GeneOutput
    features: FeatureOutput


class ExpressionOutput(BaseModel):
    host_system: str
    codon_optimization_needed: bool
    compatibility_warnings: List[WarningMessage] = Field(default_factory=list)
    notes: List[str] = Field(default_factory=list)


# ----------------------------
# Backbone agent
# ----------------------------


class BackboneInput(BaseModel):
    expression_host: str
    promoter: Optional[str]
    vector_type: Optional[str] = None


class BackboneOutput(BaseModel):
    backbone_name: str
    backbone_sequence: str
    source: str
    notes: List[str] = Field(default_factory=list)
    warnings: List[WarningMessage] = Field(default_factory=list)

    # New metadata for downstream logic
    is_loaded_vector: bool = False
    backbone_promoters: List[str] = Field(default_factory=list)
    backbone_payload_markers: List[str] = Field(default_factory=list)
    suggested_strategy: Optional[str] = None
    insertion_site_index: Optional[int] = None


# ----------------------------
# Construct agent
# ----------------------------


class Annotation(BaseModel):
    name: str
    start: int
    end: int
    type: str


class ConstructInput(BaseModel):
    gene: GeneOutput
    features: FeatureOutput
    expression: ExpressionOutput
    backbone: BackboneOutput


class ConstructOutput(BaseModel):
    feature_order: List[str]
    construct_sequence: str
    annotations: List[Annotation]
    warnings: List[WarningMessage] = Field(default_factory=list)

    # New planning metadata
    promoter_included: bool = False
    terminator_included: bool = False
    insertion_mode: Literal["full_cassette", "cds_only", "cds_plus_tags"] = "full_cassette"


# ----------------------------
# Assembly agent
# ----------------------------


class AssemblyInput(BaseModel):
    construct_sequence: str
    backbone_sequence: str
    assembly_preference: Optional[str]


class AssemblyFragment(BaseModel):
    name: str
    sequence: str
    left_homology: Optional[str] = None
    right_homology: Optional[str] = None


class AssemblyOutput(BaseModel):
    assembly_method: str
    fragments: List[AssemblyFragment]
    primer_requirements: List[str]
    assembled_sequence: str
    junction_offset: int
    warnings: List[WarningMessage] = Field(default_factory=list)


# ----------------------------
# Export agent
# ----------------------------


class ExportInput(BaseModel):
    intent: IntentOutput
    gene: GeneOutput
    features: FeatureOutput
    expression: ExpressionOutput
    construct_output: ConstructOutput
    backbone: BackboneOutput
    assembly: AssemblyOutput


class ExportOutput(BaseModel):
    summary_json: Dict[str, Any]
    genbank_record: Dict[str, Any]
    plasmid_map_data: Dict[str, Any]
    primer_list: List[Dict[str, Any]]
    bill_of_materials: List[Dict[str, Any]]


# ----------------------------
# Pipeline result
# ----------------------------


class PipelineResult(BaseModel):
    intent: IntentOutput
    gene: GeneOutput
    features: FeatureOutput
    expression: ExpressionOutput
    backbone: BackboneOutput
    construct_output: ConstructOutput
    assembly: AssemblyOutput
    export_output: ExportOutput