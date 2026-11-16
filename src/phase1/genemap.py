"""Map Human-GEM genes (current Ensembl IDs) to CCLE 2019 expression rows (GENCODE v19).

Rule (research_plan §3.3):
  1. strip version suffixes from the GENCODE v19 IDs; drop the 47 ENSGR chrY-PAR rows (D22)
  2. direct match: Human-GEM Ensembl ID == stripped v19 gene ID
  3. fallback: Human-GEM gene symbol (genes.tsv `geneSymbols`) == v19 `gene_name` (GTF), accepted only if
     exactly one v19 gene has that name and that v19 gene is not already the direct match of another
     Human-GEM gene
  4. otherwise unmapped (reason recorded)
"""
from __future__ import annotations

import gzip
import re

import numpy as np
import pandas as pd

from ..phase0 import load
from ..phase0.sources import SOURCES
from . import humangem


def read_gtf_genes() -> pd.DataFrame:
    rows = []
    with gzip.open(SOURCES["gencode_v19_gtf"].path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            rows.append({"gene_id": a.get("gene_id"), "gene_name": a.get("gene_name"),
                         "gene_type": a.get("gene_type"), "chrom": f[0]})
    g = pd.DataFrame(rows)
    g["ensg"] = g["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    return g


def expression_genes(primary_lines: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (per-gene table for the CCLE matrix, log2(TPM+1) matrix of primary lines, ENSGR dropped)."""
    tpm, _ = load.load_rnaseq()
    tpm = load.drop_par_y(tpm)
    logx = np.log2(tpm[primary_lines].astype(np.float64) + 1.0)
    g = pd.DataFrame({"gene_id": tpm.index})
    g["ensg"] = g["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    g["mean_log2tpm1_primary"] = logx.mean(axis=1).values
    g["expressed"] = g["mean_log2tpm1_primary"] >= 1.0  # research_plan §6 (outcome-blind)
    return g, logx


def map_genes(expr_genes: pd.DataFrame, gtf: pd.DataFrame) -> pd.DataFrame:
    hg = humangem.gene_table()[["genes", "geneSymbols", "geneAliases", "geneEntrezID"]].copy()
    eg = expr_genes.merge(gtf[["gene_id", "gene_name", "gene_type", "chrom"]], on="gene_id", how="left")
    by_ensg = eg.set_index("ensg")
    name_counts = eg["gene_name"].value_counts()
    direct = set(hg["genes"]) & set(eg["ensg"])
    rows = []
    for r in hg.itertuples(index=False):
        sym = (r.geneSymbols or "").split(";")[0].strip()
        rec = {"hg_gene": r.genes, "hg_symbol": sym, "ccle_gene_id": None, "ccle_gene_name": None,
               "ccle_gene_type": None, "method": None, "reason": ""}
        if r.genes in by_ensg.index:
            e = by_ensg.loc[r.genes]
            rec.update(ccle_gene_id=e["gene_id"], ccle_gene_name=e["gene_name"], ccle_gene_type=e["gene_type"],
                       method="ensembl")
        elif sym and name_counts.get(sym, 0) == 1:
            e = eg[eg["gene_name"] == sym].iloc[0]
            if e["ensg"] in direct:
                rec.update(method=None, reason=f"symbol {sym} -> {e['gene_id']}, which is already the direct "
                                                f"Ensembl match of another Human-GEM gene")
            else:
                rec.update(ccle_gene_id=e["gene_id"], ccle_gene_name=e["gene_name"],
                           ccle_gene_type=e["gene_type"], method="symbol")
        elif sym and name_counts.get(sym, 0) > 1:
            rec["reason"] = f"symbol {sym} matches {int(name_counts[sym])} GENCODE v19 genes (ambiguous)"
        else:
            rec["reason"] = (f"Ensembl ID absent from GENCODE v19 and symbol '{sym}' not a v19 gene_name"
                             if sym else "Ensembl ID absent from GENCODE v19 and no symbol")
        rows.append(rec)
    m = pd.DataFrame(rows)
    m["mapped"] = m["ccle_gene_id"].notna()
    dup = m[m["mapped"]].groupby("ccle_gene_id")["hg_gene"].transform("count")
    m.loc[m["mapped"], "n_hg_genes_sharing_ccle_gene"] = dup
    return m
