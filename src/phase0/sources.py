"""Registry of every Phase 0 data source.

All paths are relative to the repository root. Raw files live in data/raw/<group>/
and are never committed (see .gitignore). Provenance + checksums are written to
data/raw/MANIFEST.tsv by ``python -m src.phase0.download``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw"
PROCESSED = REPO / "data" / "processed"
MANIFEST = RAW / "MANIFEST.tsv"

BROAD_CCLE = "https://data.broadinstitute.org/ccle/"

# Human-GEM pinned release (decision D7): tag -> commit resolved via the GitHub API on 2026-11-16
HUMANGEM_TAG = "v2.0.0"
HUMANGEM_COMMIT = "635f533152dc5f7290ce04d12700eaa882273c3e"
# git blob SHA-1 of each file at that commit (GitHub contents API); verified after download
HUMANGEM_BLOBS = {
    "Human-GEM.xml": "031da3f40c50c294510bff9f11b6f0393b735855",
    "metabolites.tsv": "230968f43ba8526040716b174b3de6d06edebfc2",
    "genes.tsv": "bb28a3a40fdd2dc9495d950122cd6aeadc1683d7",
    "reactions.tsv": "436ec5a5eea01344fb945d24291a137e341e0ad5",
    "README.md": "bb117f119911bd7f5aebb90234fcac0cbe418fd4",
    "dependencies.txt": "46f150c8fe5b4d8d109ea9a4356870fe20ca7f30",
}


@dataclass(frozen=True)
class Source:
    key: str
    group: str
    filename: str
    url: str | None          # None -> manual download required
    release: str
    host: str
    citation: str
    note: str = ""

    @property
    def path(self) -> Path:
        return RAW / self.group / self.filename


SOURCES: dict[str, Source] = {
    s.key: s
    for s in [
        Source(
            key="metabolomics",
            group="ccle2019",
            filename="CCLE_metabolomics_20190502.csv",
            url=BROAD_CCLE + "CCLE_metabolomics_20190502.csv",
            release="CCLE 2019 (metabolomics 2019-05-02)",
            host="Broad Institute CCLE data server (data.broadinstitute.org/ccle)",
            citation="Li H et al. Nat Med 2019;25:850-860. doi:10.1038/s41591-019-0404-8",
            note="Same file as the DepMap portal 'CCLE 2019' release; portal itself is bot-protected.",
        ),
        Source(
            key="rnaseq_tpm",
            group="ccle2019",
            filename="CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz",
            url=BROAD_CCLE + "CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz",
            release="CCLE 2019 (RNA-seq RSEM gene TPM 2018-09-29; GENCODE v19)",
            host="Broad Institute CCLE data server (data.broadinstitute.org/ccle)",
            citation="Ghandi M et al. Nature 2019;569:503-508. doi:10.1038/s41586-019-1186-3",
        ),
        Source(
            key="annotations",
            group="ccle2019",
            filename="Cell_lines_annotations_20181226.txt",
            url=BROAD_CCLE + "Cell_lines_annotations_20181226.txt",
            release="CCLE 2019 (cell line annotations 2018-12-26)",
            host="Broad Institute CCLE data server (data.broadinstitute.org/ccle)",
            citation="Ghandi M et al. Nature 2019;569:503-508. doi:10.1038/s41586-019-1186-3",
        ),
        Source(
            key="li2019_supp1",
            group="li2019",
            filename="41591_2019_404_MOESM2_ESM.xlsx",
            url=(
                "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0404-8/"
                "MediaObjects/41591_2019_404_MOESM2_ESM.xlsx"
            ),
            release="Li et al. 2019 Supplementary Tables (Supp. Table 1 sheets '1-*')",
            host="Springer Nature static content (nature.com ESM)",
            citation="Li H et al. Nat Med 2019;25:850-860. doi:10.1038/s41591-019-0404-8",
        ),
        Source(
            key="model",
            group="depmap",
            filename="Model.csv",
            url=None,
            release="(read from data/raw/depmap/RELEASE.txt)",
            host="DepMap portal (depmap.org/portal/data_page) - manual browser download",
            citation="DepMap, Broad Institute. Model.csv",
            note=(
                "depmap.org is protected by a Cloudflare Turnstile human check and asks users not to "
                "scrape it; the file must be downloaded manually."
            ),
        ),
        Source(
            key="gencode_v19_gtf",
            group="ccle2019",
            filename="gencode.v19.genes.v7_model.patched_contigs.gtf.gz",
            url=BROAD_CCLE + "gencode.v19.genes.v7_model.patched_contigs.gtf.gz",
            release="CCLE 2019 gene model (GENCODE v19) used for the RSEM expression",
            host="Broad Institute CCLE data server (data.broadinstitute.org/ccle)",
            citation="Ghandi M et al. Nature 2019;569:503-508",
            note="Annotation of the CCLE 2019 expression file (gene_name, gene_type); needed for the HGNC-symbol fallback.",
        ),
        *[
            Source(
                key=f"humangem_{fn.split('.')[0].replace('-', '_').lower()}_{fn.split('.')[-1]}",
                group="humangem",
                filename=fn,
                url=f"https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/{HUMANGEM_COMMIT}/model/{fn}",
                release=f"Human-GEM {HUMANGEM_TAG} (commit {HUMANGEM_COMMIT})",
                host="GitHub SysBioChalmers/Human-GEM (raw.githubusercontent.com, pinned commit)",
                citation="Robinson JL et al. Sci Signal 2020; Luo J et al. PNAS 2026 (v2.0.0)",
                note=f"git blob sha1 {HUMANGEM_BLOBS[fn]}",
            )
            for fn in HUMANGEM_BLOBS
        ],
        Source(
            key="model_readme",
            group="depmap",
            filename="README.txt",
            url=None,
            release="(header line of the file)",
            host="DepMap portal (depmap.org/portal/data_page) - manual browser download",
            citation="DepMap, Broad Institute. Release README.txt",
            note="Release notes downloaded with Model.csv; header names the release and documents the Model.csv schema.",
        ),
    ]
}

MODEL_RELEASE_FILE = RAW / "depmap" / "RELEASE.txt"
CELLOSAURUS_DIR = RAW / "cellosaurus"
CELLOSAURUS_API = "https://api.cellosaurus.org"
