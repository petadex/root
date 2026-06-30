# `petadex-nr` "v2" — Plastic-degrading Enzyme Search Set

This directory builds the **`petadex-nr`** search query for
the `Logan` project also called `petadex.v2` set.

Starting from a manually curated set of reference enzymes (PAZy / "plasticome v1"),
 the workflow searches the NCBI `nr` protein database for homologs,
filters them by Pfam domain models, clusters them at 90% and 30% identity,
and maps every sequence back onto the original reference components.
The final `petadex.v2.RData` are the bundles the resulting data frames for downstream
analysis and visualization.

> **Notes on naming:** 
> (1) the project was renamed during development from *plasticome* to *petadex*.
> Early files use the `plasticome` prefix and later files use `petadex`;
> they are stages of the same pipeline.
>
> (2) The output of `petadex-nr` search set are files called
> **`petadex.v2.fa`**  and **`petadex.v2.RData`**
> the former is the input for the `logan-diamond` search of `Logan v1.1` data

---

## Inputs

| File | Description |
|------|-------------|
| `plasticome.v1.fa` | Curated PAZy reference enzymes (~205 sequences). The query set. |
| `plasticome.v1.meta.csv` | Per-reference metadata: gene, component, domain, CATH, kingdom, host, substrate, PDB/UniProt accession, etc. |
| `plasticome.v1.dmnd` | DIAMOND database of the reference set (used to re-map hits to references). |
| `hmm/pfam/pfam_select.hmm` | Minimal set of ~24 Pfam HMM models that collectively span all PAZy reference domains. |
| `logan_petase_v1.fa` | Auxiliary PETase sequence set (exploratory; not part of the core index). |

---

## Workflow

The numbered scripts/notebooks run in order. The **canonical pipeline that
produces `petadex.v2.RData` is `2B_plasticome_blast_redux.Rmd`**; the earlier
`2_plasticome_blast.Rmd` is a superseded draft of the same idea.

### `0_plasticome_clustering.sh` — environment & clustering notes
Shell record of tool setup (USEARCH, BLAST+, DIAMOND, `nr` download via
`update_blastdb.pl`) and the initial all-vs-all alignment / clustering commands.
Primarily a reproducibility log.

### `1_plasticome_network.Rmd` — reference sequence-similarity network
Builds an all-vs-all similarity network of the PAZy references in Cytoscape
(via `RCy3`/`igraph`), partitions them into connected **components** (`c1`–`c43`,
each roughly an enzyme family), and writes per-component sequence lists used to
build the Pfam-based HMM profiles. Saves `net/plasticome.network.Rdata`.
These components become the top-level grouping in petadex.

### `2B_plasticome_blast_redux.Rmd` — the petadex build (→ `petadex.v2.RData`)
Four conceptual steps (BLAST steps run at the shell, parsing/clustering logic run in R):

1. **BLAST nr** — `diamond blastp` of the PAZy references against `nr`
   (`--very-sensitive`, `-k0`) → `blast/pz.nr.pro`. Reduced to best (lowest
   e-value) hit per subject → `pro.unq`.
2. **HMM filter** — `hmmscan` the unique nr hits against `pfam_select.hmm`
   (`blast/pz.nr.hmm.tsv`). `mergeHMM()` collapses domain hits per sequence;
   keep hits with **e-value ≤ 1e-5 and ≥95% model coverage** (`nr.pfam.filt`),
   then extract the matched sub-sequences → `pz.nr.hmmfilt.fa`.
3. **Merge + cluster** — concatenate references + filtered nr hits
   (`petadex.pre2.fa`) and cluster with USEARCH at:
   - **90% identity → `ez` clusters** (near-redundant variants), and
   - **30% identity → `fm` clusters** (family-level).
   `readUC()` parses the `.uc` cluster files (`ez.uc90`, `fm.uc30`).
4. **Annotate + index** — `diamond blastp` every sequence back to the references
   (`petadex.pre2.pz.pro`) to assign each cluster its **component (`c`)** and
   gene. The `petadex` data frame is assembled with, per accession: length,
   `ez`/`fm` cluster IDs, centroid flags, percent identity, component, gene, and
   PAZy membership. A search index (`petadex.search` / `petadex.kv`) is built
   from EZ centroids + full-length PAZy references, named as
   `<component>_<gene>.<fm>.<ez>:<accession>`, and emitted as
   `blast/petadex.v2.fa` (~1.05M sequences).

Final cell saves the bundle:

```r
save(ez.uc90, fm.uc30, nr.pfam, nr.pfam.filt, pazy, petadex, petadex.kv,
     petadex.pre2.pz.pro, petadex.search, pfam.meta, pro.df, pro.unq,
     file = 'petadex.v2.RData')
```

### `3_plasticome_vis.Rmd` — visualization / QC
Loads `petadex.v2.RData` and reports key metrics (PAZy count, total/unique nr
hits, post-HMM filtered counts, 90%/30% cluster counts) and per-component
summaries (`component.summary.tsv`), plus ggplot distributions of percent
identity, e-value, and per-component sequence counts.

---

## Key Output: `petadex-nr.RData` / `petadex.v2.RData`

Loading this file (`load("petadex.v2.RData")`) restores:

| Object | Contents |
|--------|----------|
| `petadex` | Master table — one row per indexed sequence with accession, length, `ez`/`fm`/`c` cluster IDs, centroid flags, percent identities, gene, and PAZy flag. |
| `petadex.search` | Subset used as the searchable index (EZ centroids + PAZy references). |
| `petadex.kv` | Key→value rename map (`accession` → `<comp>_<gene>.<fm>.<ez>:<acc>`). |
| `pazy` | Cleaned PAZy reference metadata. |
| `pro.df` / `pro.unq` | All / best-per-subject nr BLAST hits. |
| `nr.pfam` / `nr.pfam.filt` | Merged Pfam domain hits, before / after the e-5 + 95%-coverage filter. |
| `ez.uc90` / `fm.uc30` | Parsed 90% and 30% USEARCH cluster assignments. |
| `petadex.pre2.pz.pro` | Hits re-mapped to references (component/gene assignment). |
| `pfam.meta` | Statistics for the selected Pfam HMM models. |

The companion FASTA `petadex.v2.fa` (~1.05M proteins, 306 Mb of residues)
is the sequence database keyed by the `petadex.kv` naming scheme.

---

## Directory Map

| Path | Role |
|------|------|
| `blast/` | All BLAST/cluster intermediates and final `petadex.v2.fa` (large). |
| `hmm/` | Pfam model selection, per-component MSAs, and HMM profiles. |
| `net/` | Reference similarity network exports for Cytoscape. |
| `struc/`, `v0/`, `v2_hit/` | Structural analysis, prior version, and example hit screens. |
| `component.summary.tsv` | Per-component counts (references, nr hits, 90%/30% clusters). |
| `pazy.nr.pro.RData` | Intermediate R object of PAZy↔nr alignments. |
