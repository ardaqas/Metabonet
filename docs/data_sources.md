# Data sources, versions and checksums (Phase 0)

The machine-readable record is `data/raw/MANIFEST.tsv`. It is written by `src/phase0/download.py` and `src/phase0/cellosaurus.py`, and `make verify` re-checks it. Raw files are **not** committed (`.gitignore`). Every file below was retrieved on **2026-11-16**.

## Files

| Key | File (under `data/raw/`) | Source / host | Release | Bytes | SHA-256 | MD5 |
|---|---|---|---|---:|---|---|
| metabolomics | `ccle2019/CCLE_metabolomics_20190502.csv` | https://data.broadinstitute.org/ccle/CCLE_metabolomics_20190502.csv | CCLE 2019; server Last-Modified 2019-05-02 | 2,098,255 | `7c1d24aa575f4c58a29019026b5df8e6d1142a56925aba32ff3f1d1d5a7fd0ac` | `22e887c4176a42e22bafabae3f2c5cbe` |
| rnaseq_tpm | `ccle2019/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz` | https://data.broadinstitute.org/ccle/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz | CCLE 2019 RSEM gene TPM (GENCODE v19); Last-Modified 2019-01-02 | 66,221,681 | `612d2e6087364fa89d57e22676a1429f58fa3bf001dfbbe1b16addb9102c651b` | `6cf7338b532ac4fa0542619af144a7b3` |
| annotations | `ccle2019/Cell_lines_annotations_20181226.txt` | https://data.broadinstitute.org/ccle/Cell_lines_annotations_20181226.txt | CCLE 2019; Last-Modified 2019-02-11 | 335,355 | `77648d1cada2f325ba2c049a5fb5408434f73b963c8e1148fb82ba2653e469e1` | `4135a80259211b6cbb2a405be8d3fc87` |
| li2019_supp1 | `li2019/41591_2019_404_MOESM2_ESM.xlsx` | Nature Medicine ESM: https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0404-8/MediaObjects/41591_2019_404_MOESM2_ESM.xlsx | Li et al. 2019 Supplementary Tables (Supp. Table 1 = sheets `1-*`); doi:10.1038/s41591-019-0404-8 | 25,133,786 | `214a31514d813687c2c7b0786c0c57e8c6fc9aaafc3a160418c06774a73f3fea` | `a9a03f5b369502a9ec05a441ff91ff85` |
| model | `depmap/Model.csv` | DepMap portal (https://depmap.org/portal/data_page/), **manual browser download** | **DepMap Public 26Q3**, confirmed by the release README (next row) | 789,805 | `65eb1895800e71a963a3ab56e5d08c8d006cb0d9e7aa68f7e34a3810a2a3cbc7` | `a7423dbf345fcf5a822ab5fc7e854e5c` |
| model_readme | `depmap/README.txt` | DepMap portal, manual browser download (release notes) | header `# DepMap Public 26Q3` | 48,405 | `5672d9e04f18320e9ca1d25747427275dcabb44bffe45eba34d196f1b44dd3ec` | `bcd481c01f79aa280d5d2e2a6234797e` |
| cellosaurus | `cellosaurus/batch_000.txt` … `batch_018.txt`, `query_meta.json` | Cellosaurus API, https://api.cellosaurus.org/search/cell-line | **Cellosaurus 57.0** (updated 2026-10-22); queried 2026-11-16T12:58Z; 927 RRIDs | — | listed per file in MANIFEST.tsv | — |

### Added in Phase 1 (2026-11-16)

| Key | File (under `data/raw/`) | Source / host | Release | Bytes | SHA-256 |
|---|---|---|---|---:|---|
| humangem_human_gem_xml | `humangem/Human-GEM.xml` (SBML) | https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/635f533…/model/Human-GEM.xml | **Human-GEM v2.0.0**, tag commit `635f533152dc5f7290ce04d12700eaa882273c3e` | 43,115,559 | `cc5a4383c6116b0c91f4db089cc640f29aec7e840249b573b74d3792c9ca4a7a` |
| humangem_metabolites_tsv | `humangem/metabolites.tsv` | same commit | v2.0.0 | 674,780 | `9d58bb3a9cd217c47dc0d54ea9517b0a23e0603bbb60d5428551e0449af1eab3` |
| humangem_genes_tsv | `humangem/genes.tsv` | same commit | v2.0.0 | 1,132,077 | `2a6058a157b3b9f3c958ba753e71d7a702d2716742972741d2ef36535ae9dff6` |
| humangem_reactions_tsv | `humangem/reactions.tsv` | same commit | v2.0.0 | 886,012 | `1436dd90e959714366270a9904fd6a5108044a2f7f1d7f8360bf7a195fc9e28f` |
| humangem_readme_md, humangem_dependencies_txt | `humangem/README.md`, `dependencies.txt` | same commit | v2.0.0 (built with RAVEN 2.11.0, libSBML 5.20.2) | — | in MANIFEST |
| gencode_v19_gtf | `ccle2019/gencode.v19.genes.v7_model.patched_contigs.gtf.gz` | https://data.broadinstitute.org/ccle/ | CCLE 2019 gene model (GENCODE v19) | 7,509,840 | `225017a79843a51d3428810bc9656e6353dcce570c816acf02ead8f7d23e8956` |
| idlookup_kegg / idlookup_chebi | `idlookup/kegg/*.txt` (141), `idlookup/chebi/*.json` (420), `query_meta.json` | KEGG REST `rest.kegg.jp/get/`, ChEBI backend API | live APIs, queried 2026-11-16T13:59Z | — | one row per file in MANIFEST |

- **Human-GEM integrity.** The tag was resolved to its commit through the GitHub API, and the files were downloaded from that commit, not from a branch. On every download and `make verify`, each file's **git blob SHA-1** is recomputed and compared with the blob SHA that GitHub lists for the commit (`src/phase0/sources.py::HUMANGEM_BLOBS`). All six match. Human-GEM v2.0.1, v2.1.0 and v2.1.1 also exist (Sept–Oct 2026) but are not used (D30).
- **COBRApy 0.32.1** is used only to parse the SBML. No solver or FBA is called.
- **GENCODE v19 GTF.** This is the gene annotation of the CCLE 2019 expression file. It supplies the `gene_name`/`gene_type` needed for the HGNC-symbol fallback (D26). It is not a new dataset. It lacks about 2.7% of the expression matrix's gene IDs, which matters only for the symbol fallback and the gene-type labels.
- **ID verification.** KEGG and ChEBI responses are cached so the mapping can be reproduced offline. **HMDB (hmdb.ca) returns a Cloudflare challenge (HTTP 403)** to scripts. PubChem no longer resolves HMDB accessions as synonyms, so HMDB IDs cannot be checked against HMDB directly; see D28.

Li et al. 2019 Supplementary Table 1 sheets used: `1-cell line annotations` (culture medium), `1-clean data`, `1-raw data` (missingness before imputation), `1-CV`, `1-lineage effects`, `1-LC-MS parameters (polar)` and `1-LC-MS parameters (lipid)` (assay class).

## Access notes and blockers

1. **The DepMap portal is bot-protected.** Every request to `depmap.org/portal`, including its `api/download/files` listing, returns a Cloudflare Turnstile "verify you're a person" page that asks users not to scrape the portal. Phase 0 therefore does not download anything from depmap.org automatically.
   - **CCLE 2019 files.** The same files, with the same names and 2019 timestamps, are on the Broad Institute's own CCLE data server (`data.broadinstitute.org/ccle`), and they were downloaded from there. This is the original host of the same release, not a substitute dataset. Byte-identity with the DepMap portal copy was **not** checked, because the portal can't be reached by script.
   - **`Model.csv`.** This has to be a manual download (decision D10), together with the release `README.txt`. The release is taken from the README header (`# DepMap Public 26Q3`). On every `make verify`, `src/phase0/download.py` also checks two things:
     - `data/raw/depmap/RELEASE.txt` names the same release;
     - the `Model.csv` columns match the README's documented `### Model.csv` schema exactly (48 of 48 columns, same order).

     If any file is missing or any check fails, verification stops. The README gives no checksum for `Model.csv`, so this confirms the release from the schema, not from the bytes.
2. **The newest DepMap release on figshare is 24Q4** (doi:10.6084/m9.figshare.27993248). It was **not** used. The plan asks for the current `Model.csv`, and the project owner supplied 26Q3.
3. **Cellosaurus is queried live.** The cached responses are the record of what was used. Re-running with `--refresh` would query whatever the current release is at that time.

## Reproduce

```bash
make env       # pinned venv (.venv) from requirements.txt
make data      # downloads (skips existing files) + Cellosaurus cache + MANIFEST.tsv
               # (place Model.csv, README.txt and RELEASE.txt in data/raw/depmap/ first)
make test      # unit tests for medium normalisation
make audit     # verify checksums, execute notebooks/00_data_audit.ipynb (rebuilds and freezes data/processed/)
```

Frozen Phase 0 outputs and their SHA-256 values are in `data/processed/FREEZE.tsv`. Two consecutive `make audit` runs gave byte-identical outputs.
