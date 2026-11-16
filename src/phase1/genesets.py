"""First-neighbour gene sets (research_plan §3.3, §5; decisions D3, D23, D24).

For a compartment-free metabolite (Human-GEM `metsNoComp`, e.g. MAM02961):
  1. collapse all compartment copies (metabolite IDs = metsNoComp + one-letter compartment)
  2. take every reaction in which any copy participates
  3. exclude: boundary reactions (exchange / sink / demand: one metabolite only) and reactions in the
     subsystems 'Exchange/demand reactions', 'Artificial reactions' (biomass and biomass-pool
     reactions), 'Pool reactions' (lipid/fatty-acid pseudo-species) and 'Protein degradation'
     (lumped plasma-protein hydrolysis; D33)
  4. a reaction is a *transport* reaction for the target if, after collapsing compartments, the
     target appears among both reactants and products (X[c] -> X[e], symport/antiport of X, ...)
  5. gene set = union of all genes in the GPR rules (AND and OR ignored)
       primary         : all retained reactions (transporters included)
       enzymatic-only  : retained reactions minus transport reactions for the target
No flux analysis, solver or optimisation is used.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# D23 + D33: 'Protein degradation' holds 35 lumped plasma-protein hydrolysis reactions that attach the same
# lysosomal proteases to every amino acid; excluded. 'Aminoacyl-tRNA biosynthesis' is a separate subsystem
# and is kept.
EXCLUDED_SUBSYSTEMS = frozenset({"Exchange/demand reactions", "Artificial reactions", "Pool reactions",
                                 "Protein degradation"})


def base_id(met_id: str) -> str:
    """Compartment-free ID: Human-GEM metabolite IDs end in a single compartment letter."""
    return met_id[:-1]


def copies(model, nocomp: str) -> list:
    return [m for m in model.metabolites if base_id(m.id) == nocomp]


def exclusion_reason(rxn) -> str | None:
    if len(rxn.metabolites) == 1:
        return "boundary (exchange/sink/demand)"
    if rxn.subsystem in EXCLUDED_SUBSYSTEMS:
        return f"subsystem '{rxn.subsystem}'"
    return None


def is_transport_for(rxn, nocomp: str) -> bool:
    lhs = {base_id(m.id) for m, c in rxn.metabolites.items() if c < 0}
    rhs = {base_id(m.id) for m, c in rxn.metabolites.items() if c > 0}
    return nocomp in lhs and nocomp in rhs


@dataclass
class GeneSet:
    target: str
    n_copies: int = 0
    compartments: str = ""
    n_reactions_all: int = 0
    n_excluded: int = 0
    excluded_by: dict = field(default_factory=dict)
    n_retained: int = 0
    n_transport: int = 0
    n_enzymatic: int = 0
    n_retained_no_gpr: int = 0
    n_enzymatic_no_gpr: int = 0
    reactions_retained: list = field(default_factory=list)
    reactions_transport: list = field(default_factory=list)
    genes_primary: set = field(default_factory=set)
    genes_enzymatic: set = field(default_factory=set)


def gene_set(model, nocomp: str) -> GeneSet:
    cps = copies(model, nocomp)
    gs = GeneSet(target=nocomp, n_copies=len(cps), compartments="".join(sorted(m.compartment for m in cps)))
    rxns = {r for m in cps for r in m.reactions}
    gs.n_reactions_all = len(rxns)
    for r in sorted(rxns, key=lambda x: x.id):
        why = exclusion_reason(r)
        if why:
            gs.n_excluded += 1
            gs.excluded_by[why] = gs.excluded_by.get(why, 0) + 1
            continue
        gs.reactions_retained.append(r.id)
        genes = {g.id for g in r.genes}
        gs.genes_primary |= genes
        gs.n_retained_no_gpr += not genes
        if is_transport_for(r, nocomp):
            gs.reactions_transport.append(r.id)
        else:
            gs.genes_enzymatic |= genes
            gs.n_enzymatic_no_gpr += not genes
    gs.n_retained = len(gs.reactions_retained)
    gs.n_transport = len(gs.reactions_transport)
    gs.n_enzymatic = gs.n_retained - gs.n_transport
    return gs


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a or b) else 0.0
