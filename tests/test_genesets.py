"""Unit tests for the first-neighbour gene-set rules (src/phase1/genesets.py) on a toy COBRApy model."""
import cobra
import pytest

from src.phase1 import genesets as gsr


def _met(mid, comp):
    return cobra.Metabolite(mid, compartment=comp, name=mid)


def _rxn(model, rid, stoich, gpr="", subsystem="Metabolism"):
    r = cobra.Reaction(rid, subsystem=subsystem)
    model.add_reactions([r])
    r.add_metabolites(stoich)
    if gpr:
        r.gene_reaction_rule = gpr
    return r


@pytest.fixture
def toy():
    m = cobra.Model("toy")
    X_c, X_m, X_e = _met("MAM00001c", "c"), _met("MAM00001m", "m"), _met("MAM00001e", "e")
    Y_c, Z_c, Y_m = _met("MAM00002c", "c"), _met("MAM00003c", "c"), _met("MAM00002m", "m")
    B = _met("MAM09999c", "c")
    m.add_metabolites([X_c, X_m, X_e, Y_c, Z_c, Y_m, B])
    # enzymatic, cytosol: complex (AND) -> both subunits counted
    _rxn(m, "R_enz_c", {Y_c: -1, X_c: 1}, "G1 and G2")
    # enzymatic, mitochondria: isozymes (OR)
    _rxn(m, "R_enz_m", {X_m: -1, Z_c: 1}, "G3 or G4")
    # transport cytosol <-> mitochondria (X on both sides after collapsing)
    _rxn(m, "R_tr_cm", {X_c: -1, X_m: 1}, "G5")
    # symport of X with Y into mitochondria: transport for X
    _rxn(m, "R_symport", {X_c: -1, Y_c: -1, X_m: 1, Y_m: 1}, "G6", subsystem="Transport reactions")
    # spontaneous reaction, no GPR
    _rxn(m, "R_spont", {X_c: -1, Z_c: 1})
    # exchange (boundary: single metabolite)
    _rxn(m, "R_EX", {X_e: -1}, subsystem="Exchange/demand reactions")
    # exchange-like boundary outside the exchange subsystem
    _rxn(m, "R_sink", {X_c: -1}, "G9", subsystem="Xenobiotics metabolism")
    # biomass / pool pseudo-reactions, even with a GPR
    _rxn(m, "R_biomass", {X_c: -0.1, Y_c: -0.1, B: 1}, subsystem="Artificial reactions")
    _rxn(m, "R_pool", {X_c: -1, B: 1}, "G7", subsystem="Pool reactions")
    # extracellular <-> cytosol transport
    _rxn(m, "R_tr_ec", {X_e: -1, X_c: 1}, "G8")
    # a reaction not involving X at all
    _rxn(m, "R_other", {Y_c: -1, Z_c: 1}, "G10")
    # lumped plasma-protein hydrolysis releasing X (D33): excluded even though it has a GPR
    P = _met("MAM09998c", "c")
    m.add_metabolites([P])
    _rxn(m, "R_protdeg", {P: -1, X_c: 3, Y_c: 2}, "G11 or G12", subsystem="Protein degradation")
    # aminoacyl-tRNA synthetase consuming X: kept (separate subsystem)
    T = _met("MAM09997c", "c")
    m.add_metabolites([T])
    _rxn(m, "R_trna", {X_c: -1, T: 1}, "G13", subsystem="Aminoacyl-tRNA biosynthesis")
    return m


def test_compartment_copies_collapsed(toy):
    gs = gsr.gene_set(toy, "MAM00001")
    assert gs.n_copies == 3 and gs.compartments == "cem"


def test_exclusions(toy):
    gs = gsr.gene_set(toy, "MAM00001")
    excluded = {"R_EX", "R_sink", "R_biomass", "R_pool", "R_protdeg"}
    assert excluded.isdisjoint(gs.reactions_retained)
    assert gs.n_excluded == 5
    assert gs.excluded_by == {"boundary (exchange/sink/demand)": 2, "subsystem 'Artificial reactions'": 1,
                              "subsystem 'Pool reactions'": 1, "subsystem 'Protein degradation'": 1}
    assert {"G7", "G9", "G11", "G12"}.isdisjoint(gs.genes_primary)


def test_trna_synthetase_kept(toy):
    gs = gsr.gene_set(toy, "MAM00001")
    assert "R_trna" in gs.reactions_retained and "G13" in gs.genes_primary and "G13" in gs.genes_enzymatic


def test_reactions_not_involving_target_ignored(toy):
    gs = gsr.gene_set(toy, "MAM00001")
    assert "R_other" not in gs.reactions_retained and "G10" not in gs.genes_primary


def test_primary_union_of_gpr_genes(toy):
    gs = gsr.gene_set(toy, "MAM00001")
    assert gs.genes_primary == {"G1", "G2", "G3", "G4", "G5", "G6", "G8", "G13"}


def test_transport_detection_and_enzymatic_only(toy):
    gs = gsr.gene_set(toy, "MAM00001")
    assert set(gs.reactions_transport) == {"R_tr_cm", "R_symport", "R_tr_ec"}
    assert gs.genes_enzymatic == {"G1", "G2", "G3", "G4", "G13"}
    assert gs.n_enzymatic == gs.n_retained - gs.n_transport == 4


def test_reactions_without_gpr_counted(toy):
    gs = gsr.gene_set(toy, "MAM00001")
    assert gs.n_retained_no_gpr == 1 and gs.n_enzymatic_no_gpr == 1   # R_spont


def test_symport_is_not_transport_for_partner_unless_moved(toy):
    # Y moves c -> m in R_symport, so it is transport for Y as well
    assert gsr.is_transport_for(toy.reactions.get_by_id("R_symport"), "MAM00002")
    # Z is only produced in R_enz_m: not transport
    assert not gsr.is_transport_for(toy.reactions.get_by_id("R_enz_m"), "MAM00003")


def test_jaccard():
    assert gsr.jaccard({"a", "b"}, {"b", "c"}) == pytest.approx(1 / 3)
    assert gsr.jaccard(set(), set()) == 0.0
