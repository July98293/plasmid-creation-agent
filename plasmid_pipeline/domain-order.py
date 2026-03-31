"""
DomainOrdererAgent: Suggests optimal N→C domain ordering for fusion proteins.

Uses heuristic scoring (charge balance, size ratio, flexibility, known precedents)
+ LLM refinement for biological reasoning.
"""

from dataclasses import dataclass, field
from typing import Optional
import json
from itertools import permutations
import math

from pydantic import BaseModel, Field
from anthropic import Anthropic

# ============================================================================
# Data Models
# ============================================================================

class DomainProperties(BaseModel):
    """Extracted structural properties for a single domain."""
    name: str
    length: int = Field(description="AA length")
    mw_kda: float = Field(description="Molecular weight in kDa")
    charge_net: float = Field(description="Net charge (K+R-E-D)")
    flexibility_score: float = Field(
        ge=0, le=1,
        description="Disorder prediction: 0=rigid, 1=highly flexible"
    )
    has_active_site: bool = Field(default=False, description="Has known binding site")
    n_term_ordered: bool = Field(default=True, description="Is N-terminus structured?")
    c_term_ordered: bool = Field(default=True, description="Is C-terminus structured?")
    original_context: Optional[str] = Field(
        default=None,
        description="e.g. 'originally a signal peptide' or 'transmembrane domain'"
    )


class OrderingScore(BaseModel):
    """Score breakdown for a single ordering."""
    permutation: list[str] = Field(description="Domain names in N→C order")
    total_score: float = Field(description="0-100, higher is better")
    
    # Component scores (each 0-100)
    electrostatic_score: float
    size_ratio_score: float
    flexibility_score: float
    site_exposure_score: float
    known_precedent_score: float
    
    justification: str = Field(description="Human-readable explanation")


class DomainOrdererOutput(BaseModel):
    """Final output from the agent."""
    top_orderings: list[OrderingScore] = Field(
        max_items=5,
        description="Top 5 ranked orderings with scores"
    )
    recommendation: str = Field(description="Which one to use and why")
    warnings: list[str] = Field(default_factory=list, description="Flags if any")


# ============================================================================
# Heuristic Scoring Functions
# ============================================================================

def score_electrostatic(properties_list: list[DomainProperties], order: list[int]) -> float:
    """
    Score based on charge balance between adjacent domains.
    Positive + negative adjacent = good. Like charges = penalty.
    """
    if len(order) < 2:
        return 100.0
    
    penalties = []
    for i in range(len(order) - 1):
        curr = properties_list[order[i]]
        next_ = properties_list[order[i + 1]]
        
        charge_product = curr.charge_net * next_.charge_net
        if charge_product > 0:  # Same sign → bad
            penalties.append(abs(charge_product) * 10)
        elif charge_product < 0:  # Opposite sign → good (bonus handled implicitly)
            pass
    
    penalty_sum = sum(penalties)
    return max(0, 100 - penalty_sum)


def score_size_ratio(properties_list: list[DomainProperties], order: list[int]) -> float:
    """
    Prefer balanced sizes. Avoid huge disparities.
    Ideal: no domain > 4x another. Penalise extremes at termini.
    """
    sizes = [properties_list[i].mw_kda for i in order]
    
    # Penalty for size disparities
    max_ratio = max(sizes) / (min(sizes) + 0.1)
    disparity_penalty = min(30, max(0, max_ratio - 4) * 10)
    
    # Additional penalty: heavy domain at N-term (harder to express/fold)
    if sizes[0] > 50:  # >50 kDa at N-term is awkward
        disparity_penalty += 15
    
    return max(0, 100 - disparity_penalty)


def score_flexibility(properties_list: list[DomainProperties], order: list[int]) -> float:
    """
    Prefer ordered domains at termini (better folding nucleation).
    Flexible domains in the middle (better linker accommodation).
    """
    score = 100.0
    
    # N-term should be ordered (rigid)
    if properties_list[order[0]].flexibility_score > 0.7:
        score -= 20
    
    # C-term should be ordered
    if properties_list[order[-1]].flexibility_score > 0.7:
        score -= 20
    
    # Middle domains can be flexible (good for linkers)
    for i in order[1:-1]:
        if properties_list[i].flexibility_score < 0.3:
            score -= 5  # Rigid in the middle is suboptimal
    
    return max(0, score)


def score_site_exposure(properties_list: list[DomainProperties], order: list[int]) -> float:
    """
    Domains with active sites should be at termini (exposed).
    Putting them in the middle = buried = bad.
    """
    score = 100.0
    
    for i, domain_idx in enumerate(order):
        domain = properties_list[domain_idx]
        if domain.has_active_site:
            # Active site at N or C term = good (no penalty)
            # Active site in middle = bad
            if i not in (0, len(order) - 1):
                score -= 25
    
    return max(0, score)


def score_known_precedent(
    properties_list: list[DomainProperties],
    order: list[int],
    precedent_db: Optional[dict] = None
) -> float:
    """
    Lookup known good fusion pairs in a database.
    Returns higher score if this ordering matches known successes.
    """
    if not precedent_db:
        return 50.0  # Neutral score if no DB
    
    score = 50.0  # Base score
    
    # Check consecutive pairs against known good fusions
    for i in range(len(order) - 1):
        pair = (properties_list[order[i]].name, properties_list[order[i + 1]].name)
        if pair in precedent_db:
            score += precedent_db[pair].get("boost", 10)
    
    return min(100, score)


def compute_ordering_score(
    properties_list: list[DomainProperties],
    order: list[int],
    precedent_db: Optional[dict] = None
) -> OrderingScore:
    """
    Compute composite score for a single ordering (permutation).
    
    Args:
        properties_list: List of DomainProperties (indexed)
        order: List of indices representing the N→C order
        precedent_db: Optional dict of known good pairs
    
    Returns:
        OrderingScore with breakdown
    """
    electrostatic = score_electrostatic(properties_list, order)
    size_ratio = score_size_ratio(properties_list, order)
    flexibility = score_flexibility(properties_list, order)
    site_exposure = score_site_exposure(properties_list, order)
    known = score_known_precedent(properties_list, order, precedent_db)
    
    # Weighted average (adjust weights to your priorities)
    weights = {
        "electrostatic": 0.25,
        "size_ratio": 0.20,
        "flexibility": 0.25,
        "site_exposure": 0.20,
        "known": 0.10
    }
    
    total = (
        electrostatic * weights["electrostatic"] +
        size_ratio * weights["size_ratio"] +
        flexibility * weights["flexibility"] +
        site_exposure * weights["site_exposure"] +
        known * weights["known"]
    )
    
    # Build justification
    domain_names = [properties_list[i].name for i in order]
    justification = f"Order: {' → '.join(domain_names)}. "
    justification += f"Electrostatic balance: {electrostatic:.0f}. "
    
    if size_ratio < 70:
        justification += f"(Warning: size imbalance). "
    
    if flexibility < 70:
        justification += f"(Flexibility mismatch at termini). "
    
    return OrderingScore(
        permutation=domain_names,
        total_score=total,
        electrostatic_score=electrostatic,
        size_ratio_score=size_ratio,
        flexibility_score=flexibility,
        site_exposure_score=site_exposure,
        known_precedent_score=known,
        justification=justification
    )


# ============================================================================
# DomainOrdererAgent
# ============================================================================

@dataclass
class DomainOrdererAgent:
    """
    Agent for suggesting optimal domain orderings in fusion proteins.
    
    Workflow:
    1. Extract/validate domain properties
    2. Score all permutations (or sample if >6 domains)
    3. Return top-K candidates
    4. (Optional) Use LLM to refine/explain top candidate
    """
    
    client: Anthropic
    model: str = "claude-opus-4-20250805"
    
    # Optional: pre-loaded precedent database
    precedent_db: Optional[dict] = None
    
    # Settings
    top_k: int = 3
    use_llm_refinement: bool = True
    max_domains_for_exhaustive: int = 6  # Beyond this, sample via annealing
    
    def extract_properties(
        self,
        domains: list[dict]
    ) -> list[DomainProperties]:
        """
        Convert domain objects (from your pipeline) to DomainProperties.
        Expects each domain to have: name, sequence, structure_info (optional)
        """
        properties = []
        for domain in domains:
            # Compute properties from sequence + structure
            seq = domain.get("sequence", "")
            
            # Net charge (K/R count - E/D count)
            charge = seq.count("K") + seq.count("R") - seq.count("E") - seq.count("D")
            
            # Rough MW: ~110 Da per AA
            mw_kda = len(seq) * 0.110
            
            # Flexibility: simple heuristic from seq (disorder-prone regions)
            # In production, call DSSP or IUPred API
            flexibility = self._estimate_flexibility(seq)
            
            prop = DomainProperties(
                name=domain.get("name", f"Domain_{len(properties)}"),
                length=len(seq),
                mw_kda=mw_kda,
                charge_net=charge,
                flexibility_score=flexibility,
                has_active_site=domain.get("has_active_site", False),
                n_term_ordered=domain.get("n_term_ordered", True),
                c_term_ordered=domain.get("c_term_ordered", True),
                original_context=domain.get("context", None),
            )
            properties.append(prop)
        
        return properties
    
    def _estimate_flexibility(self, sequence: str) -> float:
        """
        Quick flexibility heuristic (DSSP-like).
        Count disorder-prone residues (G, A, S, P).
        Returns 0 (rigid) to 1 (flexible).
        """
        disorder_residues = "GASPS"
        disorder_count = sum(1 for aa in sequence if aa in disorder_residues)
        return min(1.0, disorder_count / max(len(sequence), 1))
    
    def suggest_orderings(
        self,
        domains: list[dict],
        user_query: Optional[str] = None,
    ) -> DomainOrdererOutput:
        """
        Main entry point: suggest domain orderings.
        
        Args:
            domains: List of domain dicts with name, sequence, etc.
            user_query: Optional user note (e.g., "GFP should be C-terminal")
        
        Returns:
            DomainOrdererOutput with ranked orderings
        """
        # Step 1: Extract properties
        properties = self.extract_properties(domains)
        n_domains = len(properties)
        
        if n_domains < 2:
            return DomainOrdererOutput(
                top_orderings=[],
                recommendation="Need at least 2 domains for fusion.",
                warnings=["Insufficient domains"]
            )
        
        # Step 2: Score permutations
        if n_domains <= self.max_domains_for_exhaustive:
            # Exhaustive: score all N! permutations
            all_orders = list(permutations(range(n_domains)))
        else:
            # Sample: use greedy + random swaps (simulated annealing lite)
            all_orders = self._sample_orderings(n_domains, sample_size=100)
        
        scores = [
            compute_ordering_score(properties, list(order), self.precedent_db)
            for order in all_orders
        ]
        
        # Step 3: Sort and select top-K
        scores.sort(key=lambda x: x.total_score, reverse=True)
        top_orderings = scores[:self.top_k]
        
        # Step 4: Optional LLM refinement
        if self.use_llm_refinement:
            recommendation = self._refine_with_llm(
                properties, top_orderings, user_query
            )
        else:
            recommendation = (
                f"Top recommendation: {' → '.join(top_orderings[0].permutation)} "
                f"(score: {top_orderings[0].total_score:.1f}). "
                f"{top_orderings[0].justification}"
            )
        
        # Collect warnings
        warnings = self._extract_warnings(properties, top_orderings[0])
        
        return DomainOrdererOutput(
            top_orderings=top_orderings,
            recommendation=recommendation,
            warnings=warnings,
        )
    
    def _sample_orderings(self, n: int, sample_size: int = 100) -> list[tuple]:
        """
        For >6 domains, sample permutations via randomized search.
        Simple approach: start with greedy, then random swaps.
        """
        import random
        
        # Start with a greedy ordering (e.g., by size)
        current = list(range(n))
        sampled = [tuple(current)]
        
        for _ in range(sample_size - 1):
            # Random swap
            i, j = random.sample(range(n), 2)
            current[i], current[j] = current[j], current[i]
            sampled.append(tuple(current))
        
        return sampled
    
    def _refine_with_llm(
        self,
        properties: list[DomainProperties],
        top_orderings: list[OrderingScore],
        user_query: Optional[str]
    ) -> str:
        """
        Use Claude to explain and refine the top ordering.
        Catches biological edge cases and provides reasoning.
        """
        # Format properties as readable text
        domain_descriptions = []
        for prop in properties:
            desc = (
                f"- {prop.name}: {prop.length} aa, {prop.mw_kda:.1f} kDa, "
                f"net charge {prop.charge_net:+.0f}, "
                f"flexibility {prop.flexibility_score:.2f}"
            )
            if prop.has_active_site:
                desc += " [has active site]"
            if prop.original_context:
                desc += f" ({prop.original_context})"
            domain_descriptions.append(desc)
        
        domains_text = "\n".join(domain_descriptions)
        
        # Format top orderings for Claude
        orderings_text = "\n".join([
            f"{i+1}. {' → '.join(o.permutation)} (score: {o.total_score:.1f})\n"
            f"   {o.justification}"
            for i, o in enumerate(top_orderings)
        ])
        
        prompt = f"""You are a protein engineer. Given these domains and their properties:

{domains_text}

The heuristic scoring suggests these top orderings:
{orderings_text}

{"User note: " + user_query if user_query else ""}

Briefly explain:
1. Why is the top ordering the best choice?
2. Are there any biological concerns (e.g., burial of active sites, domain misfolding)?
3. Any adjustments you'd suggest?

Keep your answer concise (2-3 sentences per point)."""

        response = self.client.messages.create(
            model=self.model,
            max_tokens=500,
            messages=[{"role": "user", "content": prompt}]
        )
        
        return response.content[0].text
    
    def _extract_warnings(
        self,
        properties: list[DomainProperties],
        best_ordering: OrderingScore
    ) -> list[str]:
        """
        Extract flags from the best ordering.
        """
        warnings = []
        
        # Map back to indices
        name_to_idx = {p.name: i for i, p in enumerate(properties)}
        order_indices = [name_to_idx[name] for name in best_ordering.permutation]
        
        # Check for issues
        if best_ordering.total_score < 60:
            warnings.append("Low confidence ordering (score < 60). Consider manual review.")
        
        if best_ordering.electrostatic_score < 50:
            warnings.append("Poor electrostatic balance. Domains may repel.")
        
        if best_ordering.size_ratio_score < 60:
            warnings.append("Unbalanced domain sizes. May cause expression issues.")
        
        # Check for active site burial
        for i, domain_idx in enumerate(order_indices):
            if properties[domain_idx].has_active_site and i not in (0, len(order_indices) - 1):
                warnings.append(
                    f"Active site of {properties[domain_idx].name} may be buried. "
                    "Consider reordering or using longer linker."
                )
        
        return warnings


# ============================================================================
# Usage Example
# ============================================================================

if __name__ == "__main__":
    from anthropic import Anthropic
    
    # Example usage
    client = Anthropic()
    
    orderer = DomainOrdererAgent(
        client=client,
        model="claude-opus-4-20250805",
        use_llm_refinement=True,
    )
    
    # Mock domains (in your pipeline, these come from FetchDomainsAgent)
    example_domains = [
        {
            "name": "MBP",
            "sequence": "M" * 380,  # 380 aa
            "has_active_site": False,
            "context": "Solubility tag"
        },
        {
            "name": "GFP",
            "sequence": "M" * 238,
            "has_active_site": False,
            "context": "Fluorescent protein"
        },
        {
            "name": "His6",
            "sequence": "HHHHHH",
            "has_active_site": False,
            "context": "Purification tag"
        },
    ]
    
    result = orderer.suggest_orderings(
        example_domains,
        user_query="GFP should probably be at the C-terminus for minimal interference."
    )
    
    print("\n=== Domain Ordering Results ===")
    print(f"Recommendation:\n{result.recommendation}\n")
    print("Top orderings:")
    for i, ordering in enumerate(result.top_orderings, 1):
        print(f"{i}. {' → '.join(ordering.permutation)} (score: {ordering.total_score:.1f})")
    
    if result.warnings:
        print("\nWarnings:")
        for w in result.warnings:
            print(f"  ⚠ {w}")