from __future__ import annotations

from typing import List, Optional, Literal, Any, Dict

from pydantic import BaseModel, Field


# ----------------------------
# Core shared types
# ----------------------------


class WarningMessage(BaseModel):
    code: str
    message: str


class TagSpec(BaseModel):
    name: str
    position: Optional[Literal["N-terminal", "C-terminal"]] = None


# ----------------------------
# Intent agent
# ----------------------------


class IntentInput(BaseModel):
    user_request: str


class IntentOutput(BaseModel):
    # Core requested/planned cargo
    gene_symbol: Optional[str] = None
    target_species: Optional[str] = None

    # Expression context
    expression_host: Optional[str] = None
    best_expression_host: Optional[str] = None

    # Design choices
    backbone: Optional[str] = None
    promoter: Optional[str] = None
    assembly_method: Optional[str] = None

    # Features
    n_terminal_tag: Optional[str] = None
    c_terminal_tag: Optional[str] = None
    extra_tags: List[TagSpec] = Field(default_factory=list)

    selection_marker: Optional[str] = None
    origin_of_replication: Optional[str] = None
    terminator: Optional[str] = None
    polyA: Optional[str] = None
    cloning_site: Optional[str] = None

    # Metadata
    notes: List[str] = Field(default_factory=list)
    warnings: List[WarningMessage] = Field(default_factory=list)

    # Fields the LLM auto-filled (not explicitly stated by the user)
    suggested_fields: List[str] = Field(default_factory=list)


# ----------------------------
# Gene agent
# ----------------------------


class GeneInput(BaseModel):
    gene_symbol: str
    target_species: Optional[str] = None


class GeneOutput(BaseModel):
    gene_symbol: str
    organism: Optional[str] = None
    gene_id: Optional[str] = None
    transcript_id: Optional[str] = None
    cds_sequence: Optional[str] = None
    protein_name: Optional[str] = None
    warnings: List[WarningMessage] = Field(default_factory=list)


# ----------------------------
# Feature agent
# ----------------------------


class FeatureInput(BaseModel):
    promoter: Optional[str] = None
    n_terminal_tag: Optional[str] = None
    c_terminal_tag: Optional[str] = None
    extra_tags: List[TagSpec] = Field(default_factory=list)

    terminator: Optional[str] = None
    polyA: Optional[str] = None
    selection_marker: Optional[str] = None
    origin_of_replication: Optional[str] = None

    expression_host: Optional[str] = None


class ResolvedFeature(BaseModel):
    name: str
    type: Literal["promoter", "kozak", "tag", "terminator", "selection_marker", "origin_of_replication"]
    sequence: str
    source: str = "unknown"  # e.g. "NCBI:123456" | "constant"
    length: int = 0
    validated: bool = False
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

    # Recommendation layer
    requested_host: Optional[str] = None
    recommended_host: Optional[str] = None

    compatibility_warnings: List[WarningMessage] = Field(default_factory=list)
    notes: List[str] = Field(default_factory=list)


# ----------------------------
# Backbone agent
# ----------------------------


class BackboneInput(BaseModel):
    expression_host: str
    backbone: Optional[str] = None  # from IntentOutput.backbone
    vector_type: Optional[str] = None  # e.g. "plasmid", "lentiviral"


class BackboneOutput(BaseModel):
    backbone_name: str
    backbone_sequence: str
    backbone_length: int
    source: str  # e.g. "NCBI:2459837522" or "Addgene:52535"
    backbone_genbank: Optional[str] = None  # raw GenBank text for downstream feature extraction
    notes: Optional[str] = None
    warnings: List[WarningMessage] = Field(default_factory=list)


# ----------------------------
# Assembly agent
# ----------------------------


class Annotation(BaseModel):
    name: str
    start: int
    end: int
    type: str


class AssemblyInput(BaseModel):
    backbone_sequence: str
    backbone_name: str
    backbone_genbank: Optional[str] = None
    cds_sequence: str
    gene_symbol: str = ""
    features: List[ResolvedFeature]
    assembly_method: Optional[str] = None


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
    cassette_sequence: str = ""
    feature_order: List[str] = Field(default_factory=list)
    cassette_annotations: List[Annotation] = Field(default_factory=list)
    warnings: List[WarningMessage] = Field(default_factory=list)


# ----------------------------
# Export agent
# ----------------------------


class ExportInput(BaseModel):
    intent: IntentOutput
    gene: GeneOutput
    features: FeatureOutput
    expression: ExpressionOutput
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
    assembly: AssemblyOutput
    export_output: ExportOutput