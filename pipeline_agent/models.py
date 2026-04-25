from __future__ import annotations

from typing import List, Literal, Optional, Tuple, Dict
from datetime import datetime

from pydantic import BaseModel, Field


WEB_FALLBACK_WARNING = (
    "element not present in the dataset, we retrive it from web but it may be "
    "alredy with cargo - better to dubble check the plasmid!"
)

TemplateType = Literal["standard", "fusion", "multi_cassette_fusion", "synthetic_circuit"]
SeverityLevel = Literal["critical", "warning", "info"]
DifficultyLevel = Literal["easy", "moderate", "challenging"]
QualityRating = Literal["excellent", "good", "acceptable", "poor"]


# ---------------------------------------------------------------------------
# 1. BACKBONE & FEATURE DETECTION
# ---------------------------------------------------------------------------

class BackboneFeature(BaseModel):
    """One feature (promoter/marker/ori) detected in backbone."""
    type: str  # "promoter", "marker", "terminator", "ori", "MCS"
    name: str  # "CMV", "Ampicillin", "SV40", "ColE1"
    coordinates: Optional[Tuple[int, int]] = None  # (start, end) in sequence
    strength: Optional[str] = None  # "strong", "medium", "weak"


class BackboneMetadata(BaseModel):
    """Backbone resolved from Addgene with feature analysis."""
    addgene_id: Optional[int] = None
    name: str  # "pcDNA3.1"
    sequence: str
    length_bp: int
    features: List[BackboneFeature] = Field(default_factory=list)
    
    # Smart analysis (pre-computed)
    has_promoter: bool = False
    has_terminator: bool = False
    has_marker: bool = False
    has_ori: bool = False
    promoter_names: List[str] = Field(default_factory=list)
    marker_names: List[str] = Field(default_factory=list)
    
    source: str  # "addgene_backbone:12345" or "dataset:..." or "ncbi:..."
    notes: Optional[str] = None


# ---------------------------------------------------------------------------
# 2. ASSEMBLY PLANNING
# ---------------------------------------------------------------------------

class AssemblyOverlap(BaseModel):
    """Overlap between two adjacent elements for Gibson/CPEC."""
    left_element: str  # slot name
    right_element: str
    overlap_sequence: str
    overlap_length: int  # typically 15-40 bp
    gc_content: float  # 0.0-1.0


class AssemblyPlan(BaseModel):
    """Detailed assembly strategy."""
    method: str  # "Gibson", "BsaI/BbsI", "Ligation", "TOPO", "TA", "Seamless"
    elements_order: List[str]  # ordered slot names for assembly
    overlaps: List[AssemblyOverlap] = Field(default_factory=list)
    restriction_sites: Optional[Dict[str, List[int]]] = None  # e.g., {"BsaI": [123, 456]}
    oligo_design: Optional[List[Dict]] = None  # for synthetic overhangs
    
    estimated_difficulty: DifficultyLevel = "moderate"
    notes: List[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# 3. RESOLUTION TRACKING
# ---------------------------------------------------------------------------

class ResolutionAttempt(BaseModel):
    """One resolution attempt (dataset → Addgene → NCBI → placeholder)."""
    attempt_num: int
    source_type: str  # "dataset", "addgene_api", "ncbi", "placeholder"
    query: str  # what was searched
    result_found: bool
    element_id: Optional[str] = None  # Addgene ID, NCBI acc, etc
    sequence_length: int = 0
    confidence: float  # 0.0-1.0 confidence in match
    error_msg: Optional[str] = None


class ElementResolutionHistory(BaseModel):
    """Full history of resolving one slot."""
    slot_name: str
    requested_name: str
    element_type: str
    attempts: List[ResolutionAttempt] = Field(default_factory=list)
    final_result: Optional[dict] = None  # ResolvedElement dict
    warnings: List[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# 4. CODON OPTIMIZATION
# ---------------------------------------------------------------------------

class CodonUsageProfile(BaseModel):
    """Codon usage for expression host."""
    organism: str  # "E. coli", "P. pastoris", "CHO", "S. cerevisiae"
    codon_table: Dict[str, float]  # {"ATG": 0.95, "GTG": 0.04, ...}
    gc_content_target: float  # ideal %GC (e.g., 0.52 for E. coli)


class CodonOptimizationRequest(BaseModel):
    """Request codon optimization for a gene."""
    gene_sequence: str
    target_organism: str
    method: str  # "GeneDesigner", "JCAT", "IDT", "synthetic"
    preserve_motifs: List[str] = Field(default_factory=list)  # restrict sites
    
    # Output
    optimized_sequence: Optional[str] = None
    original_codon_bias: Optional[float] = None
    optimized_codon_bias: Optional[float] = None


# ---------------------------------------------------------------------------
# 5. SEQUENCE VALIDATION & QC
# ---------------------------------------------------------------------------

class ValidationRule(BaseModel):
    """One validation rule to check."""
    rule_name: str  # "no_stop_codons_in_frame", "no_homopolymers", "no_BsaI_sites"
    description: str
    severity: SeverityLevel


class SequenceValidationReport(BaseModel):
    """QC report for final construct."""
    total_length: int
    gc_content: float
    has_start_codon: bool = False
    has_stop_codon: bool = False
    frameshifts: List[Tuple[int, str]] = Field(default_factory=list)  # [(pos, description)]
    homopolymer_runs: List[Tuple[int, int, str]] = Field(default_factory=list)  # [(start, end, base)]
    restriction_sites_found: Dict[str, List[int]] = Field(default_factory=dict)  # {"BsaI": [123, 456]}
    repeat_sequences: List[Tuple[str, List[int]]] = Field(default_factory=list)  # [(seq, positions)]
    
    warnings: List[str] = Field(default_factory=list)
    errors: List[str] = Field(default_factory=list)
    is_valid: bool = True


# ---------------------------------------------------------------------------
# 6. EXPRESSION PROFILING
# ---------------------------------------------------------------------------

class ExpressionCondition(BaseModel):
    """Growth + induction condition."""
    temperature: int  # °C
    induction_agent: Optional[str] = None  # "IPTG", "Arabinose", "Doxycycline"
    induction_time: Optional[int] = None  # minutes
    expression_level: str = "medium"  # "low", "medium", "high"
    solubility: Optional[str] = None  # "insoluble", "partially soluble", "soluble"


class ExpressionProfile(BaseModel):
    """Profile for predicted expression."""
    expression_host: str
    predicted_level: str = "medium"  # "low", "medium", "high"
    predicted_solubility: str = "partially soluble"
    codon_bias_score: float = 0.5  # 0-1 (1 = optimal)
    gc_content_score: float = 0.5  # penalty if too extreme
    secondary_structure_risk: str = "low"  # "low", "moderate", "high"
    
    recommended_conditions: List[ExpressionCondition] = Field(default_factory=list)
    notes: List[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# 7. MULTI-COMPONENT SYSTEMS
# ---------------------------------------------------------------------------

class PlasmidComponent(BaseModel):
    """One plasmid in multi-plasmid system."""
    plasmid_name: str  # "pA-tetR", "pB-gfp"
    role: str  # "regulator", "target", "reporter"
    backbone: str
    elements: List[dict] = Field(default_factory=list)  # ResolvedElement dicts
    construct_sequence: str
    selection_marker: str  # must be different from other components
    copy_number: str = "medium"  # "low", "medium", "high"


class MultiPlasmidSystem(BaseModel):
    """Multi-plasmid construct (e.g., two-hybrid, split fluorescent)."""
    system_name: str
    system_type: str  # "two-plasmid", "three-plasmid", "split-function"
    components: List[PlasmidComponent] = Field(default_factory=list)
    
    # Cross-plasmid validation
    marker_conflicts: List[str] = Field(default_factory=list)
    ori_compatibility: str = "compatible"  # "compatible", "incompatible"
    notes: List[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# 8. SYNTHETIC CIRCUITS
# ---------------------------------------------------------------------------

class CircuitComponent(BaseModel):
    """One part/component in a synthetic genetic circuit."""
    name: str
    role: str  # e.g. "sensor", "actuator", "regulator", "reporter"
    element_type: str  # e.g. "promoter", "CDS", "RBS", "terminator"
    sequence_name: Optional[str] = None  # resolved element name


class CircuitDevice(BaseModel):
    """One logical device in circuit."""
    device_name: str  # "NOT gate", "AND gate"
    components: List[CircuitComponent] = Field(default_factory=list)
    logic_table: Optional[Dict] = None  # {input: output, ...}
    response_time: Optional[str] = None  # "fast", "slow"


class SyntheticCircuitSpec(BaseModel):
    """Full specification for synthetic circuit."""
    circuit_name: str
    circuit_type: str  # "logic gate", "oscillator", "toggle switch"
    devices: List[CircuitDevice] = Field(default_factory=list)
    
    # Signals
    input_signals: List[str] = Field(default_factory=list)  # ["IPTG", "aTc"]
    output_reporters: List[str] = Field(default_factory=list)  # ["GFP", "mCherry"]
    
    # Expected behavior
    expected_logic_table: Optional[Dict] = None
    expected_response_curve: Optional[str] = None
    
    # Assembly
    assembly_strategy: str = "Golden Gate"  # "Golden Gate", "Gibson", "BioBricks"
    num_parts_required: int = 0


# ---------------------------------------------------------------------------
# 9. FUSION PROTEINS
# ---------------------------------------------------------------------------

class FusionChainEntry(BaseModel):
    """One protein in a multi-part fusion chain."""
    gene: str
    tag: Optional[str] = None  # e.g. "His6", "FLAG"
    cleavage_site: Optional[str] = None  # e.g. "TEV", "3C"
    linker: Optional[str] = None  # linker AFTER this domain


# ---------------------------------------------------------------------------
# 10. QUALITY METRICS
# ---------------------------------------------------------------------------

class SequenceMetrics(BaseModel):
    """Metrics for final sequence quality."""
    total_length: int
    real_sequence_percentage: float  # 0-100
    placeholder_count: int = 0
    dataset_elements: int = 0
    web_elements: int = 0
    
    # Sequence properties
    gc_content: float = 0.0
    codon_bias_score: Optional[float] = None
    secondary_structure_complexity: Optional[str] = None


class ConstructQualityScore(BaseModel):
    """Overall quality/completeness score."""
    completeness: float = 0.0  # % non-placeholder (0-100)
    sequence_validity: float = 0.0  # passed QC checks (0-100)
    assembly_feasibility: float = 0.0  # ease of assembly (0-100)
    expression_potential: float = 0.0  # predicted expression level (0-100)
    overall_score: float = 0.0  # weighted average (0-100)
    
    rating: QualityRating = "acceptable"


# ---------------------------------------------------------------------------
# CORE INTENT & RESOLUTION (refactored with new metadata)
# ---------------------------------------------------------------------------

class IntentInput(BaseModel):
    user_request: str


class IntentOutput(BaseModel):
    gene_symbol: str
    target_species: str = "unspecified"
    expression_host: str = "unspecified"
    best_expression_host: Optional[str] = None
    backbone: Optional[str] = None
    promoter: Optional[str] = None
    promoter2: Optional[str] = None
    rbs: Optional[str] = None
    rbs2: Optional[str] = None
    assembly_method: Optional[str] = None
    n_terminal_tag: Optional[str] = None
    c_terminal_tag: Optional[str] = None
    extra_tags: List[str] = Field(default_factory=list)
    selection_marker: Optional[str] = None
    origin_of_replication: Optional[str] = None
    terminator: Optional[str] = None
    polyA: Optional[str] = None
    cloning_site: Optional[str] = None
    regulatory_element: Optional[str] = None
    construct_type: TemplateType = "standard"

    # Standard / single-gene fields
    gene_a: Optional[str] = None
    gene_b: Optional[str] = None
    linker_sequence: Optional[str] = None

    # Fusion protein fields
    fusion_chain: List[FusionChainEntry] = Field(default_factory=list)
    cleavage_sites: List[str] = Field(default_factory=list)

    # Synthetic circuit fields
    circuit_topology: Optional[str] = None
    circuit_components: List[CircuitComponent] = Field(default_factory=list)
    circuit_logic: Optional[str] = None
    input_signals: List[str] = Field(default_factory=list)
    output_reporters: List[str] = Field(default_factory=list)
    circuit_regulators: List[str] = Field(default_factory=list)

    notes: List[str] = Field(default_factory=list)
    suggested_fields: List[str] = Field(default_factory=list)


class ResolvedElement(BaseModel):
    slot: str
    requested_name: str
    element_type: str
    sequence: str
    source: str
    from_web: bool = False


class PlasmidConstructionInput(BaseModel):
    intent: IntentOutput


class PlasmidConstructionOutput(BaseModel):
    template: TemplateType
    elements: List[ResolvedElement]
    construct_sequence: str
    frontend_warnings: List[str] = Field(default_factory=list)
    
    # Enhanced with detailed metadata
    backbone_metadata: Optional[BackboneMetadata] = None
    assembly_plan: Optional[AssemblyPlan] = None
    resolution_history: List[ElementResolutionHistory] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# OUTPUT & COMPREHENSIVE DESIGN REPORT
# ---------------------------------------------------------------------------

class OutputInput(BaseModel):
    construct_output: PlasmidConstructionOutput
    name: str = "designed_plasmid"


class OutputOutput(BaseModel):
    gbk_text: str
    file_name: str
    frontend_warnings: List[str] = Field(default_factory=list)


class DesignReport(BaseModel):
    """Complete design report (JSON/PDF exportable)."""
    project_name: str
    timestamp: str = Field(default_factory=lambda: datetime.now().isoformat())
    
    # Core design
    intent: IntentOutput
    backbone_metadata: Optional[BackboneMetadata] = None
    elements: List[ResolvedElement] = Field(default_factory=list)
    
    # Detailed analysis
    assembly_plan: Optional[AssemblyPlan] = None
    expression_profile: Optional[ExpressionProfile] = None
    validation_report: Optional[SequenceValidationReport] = None
    construct_metrics: Optional[SequenceMetrics] = None
    quality_score: Optional[ConstructQualityScore] = None
    
    # Sequence outputs
    construct_sequence: str = ""
    genbank_file: Optional[str] = None
    fasta_file: Optional[str] = None
    
    # Recommendations
    next_steps: List[str] = Field(default_factory=list)
    troubleshooting_notes: List[str] = Field(default_factory=list)
    
    # Optional advanced features
    synthetic_circuit_spec: Optional[SyntheticCircuitSpec] = None
    multi_plasmid_system: Optional[MultiPlasmidSystem] = None