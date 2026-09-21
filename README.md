# Salmonella Isangi

## Project

We investigated the hospital-associated emergence of extensively drug-resistant Salmonella Isangi in Malawi and an outbreak involving five hospitals in South Africa. We combined clinical surveillance, neonatal-unit and river sampling, whole-genome sequencing, antimicrobial susceptibility testing and experimental phenotyping. The analysis places these collections within a wider population of 345 Isangi genomes.

This repository accompanies **Genomic epidemiology of hospital-associated XDR Salmonella Isangi in Malawi and South Africa**. It brings together the work in the order we undertook it: isolate selection, assembly and typing, reference-based SNP analysis, resistance and plasmid analysis, and the final figures.

The publication analysis starts from the fixed trees, pairwise SNP matrices and laboratory tables in `data/`. These files are included, so the downstream analysis does not depend on access to the original sequencing server or a changing public database. The command-line stages are kept separately in `workflow/`.

## Run the publication analysis

Open `salmonella-isangi.Rproj`, or run these commands from the repository root:

```sh
Rscript environment/install.R
Rscript tests/check_R.R
Rscript scripts/reproduce.R
```

Each run writes to a new directory under `results/rebuilt/`. `results/LATEST_RUN.txt` records its location. The input files and archived reference figures are not overwritten. Figures are exported as LZW-compressed TIFFs at 600 dpi, with PNG previews and the corresponding source tables.

Stages can also be run separately:

```sh
Rscript scripts/reproduce.R cohorts,snps
Rscript scripts/reproduce.R phenotype
Rscript scripts/reproduce.R phylogenies,burden,composite
Rscript scripts/reproduce.R matrices
```

For the AMR-burden Quarto document:

```sh
quarto render scripts/Figure_3B_AMR_gene_burden.qmd
```

The Quarto project uses the repository root as its working directory. Inputs remain in `data/data_raw/`, not beside the document.

## Layout

| Directory | Contents |
|---|---|
| `data/data_raw/` | Supplied AMRFinderPlus output, the dated NCBI export, original SNP matrix and sequence-type lookup |
| `data/data_clean/` | Publication-safe metadata, recorded susceptibility results, resistance rules and figure inputs |
| `data/manifests/` | Stage-specific isolate lists, the deposition register and the publication cohort |
| `data/trees/`, `data/matrices/` | Fixed phylogenies and numeric SNP distances |
| `data/artwork/`, `data/geography/` | Original experimental/plasmid/map artwork and environmental coordinates |
| `workflow/` | Separate command-line analysis stages and PHEnix Nextflow workflow |
| `scripts/analysis/` | Cohort linkage, SNP summaries and AMR-gene burden |
| `scripts/figures/` | Main figures, labelled trees, SNP heatmaps and composite plates |
| `provenance/` | Original command records, source checksums and release validation |
| `results/reference/` | Supplied figure plates retained for comparison |

## Isolates and identifiers

The final collection contains 345 genomes. The CC25 tree contains 327: 224 ST335, 102 ST216 and one ST5028. The Malawi ST335 and South African outbreak-associated matrices contain 74 and 37 isolates, respectively. These are analytical subsets, not interchangeable recruitment denominators.

`isolate_accession_line_list_345.csv` is the main identifier register. It links analytical IDs, original sequencing aliases, run and sample accessions, metadata and figure-panel positions. `publication_cohort.tsv` provides a compact version. The 80 study-deposited isolates are indexed under ENA study **ERP189265**. Technical suffixes and the river-isolate `#` convention are retained in the source-alias fields.

The complete CC25 matrix contains 331 isolates. ERR9939684, ERR11201713, ERR9939690 and SRR3049609 occur in that matrix but not in the 327-tip tree. They remain in the complete numeric file and are excluded only from tree-aligned views. The reconciliation table records this explicitly. All between-panel comparisons remain available in the full matrices.

**CH10FJ** is the reference-assembly isolate. The original reference-building files use the label **CHF10J1**; those filenames are preserved in the upstream records. **CAAP2A** is a fixed plotting anchor for Figure 1B only. Neither the minimum spanning tree nor the SNP heatmaps is calculated from anchor distances.

## Resistance analysis

The recorded susceptibility tables contain 19 Isangi and 23 comparator isolates tested against 17 antibiotics. Measurements and recorded S/I/R interpretations are kept separately. CTX5 remains as recorded. The code does not apply new clinical breakpoints or create an AMX result. The study's operational phenotypic XDR rule requires resistant calls for AMP10, C30, SXT25 and PEF5, plus CPD10 or CTX5. It identifies 12 of the 19 Isangi isolates.

Genomic resistance profiles are calculated independently from the explicit marker dictionary. The co-trimoxazole component requires both trimethoprim and sulfonamide markers. These are genomic profiles, not measured phenotypes.

Figure 3B uses the supplied NCBI export labelled **13 January 2026** and our own AMRFinderPlus output. The NCBI count is the number of `COMPLETE` genotype entries after excluding `mdsA` and `mdsB`. The study count is the number of distinct gene symbols with both element type and subtype `AMR`; point-mutation and stress/biocide records are excluded. The external group contains records not identified as Isangi, including unresolved serovars. “Other Isangi” includes 115 isolates assigned a non-ST335 type and six with unassigned type.

These source-specific counting rules reproduce the supplied comparison; they do not constitute a rerun through a common AMRFinderPlus database. The existing zero-fill convention is retained for study records without a qualifying hit. The input audit identifies records absent from the hit-only file, rather than treating absence from that file as proof of a completed negative assay.

## Figures

Main Figure 2 is the sole 19-isolate susceptibility heatmap. Main Figure 3 combines the circular CC25 tree with the approved AMR-burden plot. Main Figure 4 contains the ST335 genomic determinant heatmap without a second phenotype panel.

The CC25 details retain isolate IDs and run accessions, with a panel index linking them to the full tree. Local SNP heatmaps use linear scales; broader overviews use square-root colour scaling with raw SNP-count labels. Row order and matrix values are exported alongside each heatmap.

Experimental microscopy, biofilm summaries, mouse survival, the plasmid comparison and plasmid network are retained as source artwork. The default Figure 1 rebuild uses the archived map for panel D; its annual surveillance, temporal and minimum-spanning-tree panels are calculated from the included data. The vector-map implementation is retained for use with the original spatial layers. Raster conversion does not add resolution to the source artwork.

`docs/figure_index.csv` maps the repository outputs to the submitted appendix. The 19 numbered source tables are in `results/tables/`; `docs/table_index.csv` records their inputs. Regenerate them with `python3 scripts/export_tables.py`. The supplied appendix repeats CC25-05 as Figures 10 and 12; both positions are recorded without counting that group twice in the analysis.

## Upstream work

`workflow/` contains the separate read-processing, assembly, typing, reference-building, PHEnix, SNP, Gubbins and phylogeny stages. `docs/upstream.md` describes their inputs and execution order. Historical commands are in `provenance/commands/`; cleaned wrappers are release code.

The fixed-output publication route and the raw-read route have different starting points. In particular, the historical PHEnix filtering YAML and the final reference FASTA are not included here. A raw-read execution requires those resources to be supplied explicitly. No replacement filtering values have been assigned to the historical run. The archived tree commands use classic RAxML; the optional RAxML-NG wrapper is a separate rerun route. Neither wrapper changes the supplied publication trees.

## Checks

```sh
python3 -m pip install -r environment/requirements.txt
python3 scripts/audit.py --burden
python3 -m unittest discover -s tests -v
```

The checks cover isolate/accession linkage, tree membership, complete and subset SNP matrices, resistance rules, susceptibility measurements, panel coverage, original-source reconciliation and shell syntax. `provenance/VALIDATION.md` records what was executed for this archive. GitHub Actions supplies an additional R parsing/data check and an optional full figure-rendering job.

`environment/capture.R` records the R environment used for a run. It can create `environment/renv.lock` when called with `--lock` after `renv` is installed. The documented upstream software versions are source records, not a reconstructed historical environment lock.

## Citation and reuse

Use `CITATION.cff` to cite this software archive. The repository is <https://github.com/PeterJohnston296/salmonella-isangi>. Version 1.0.0 is archived at <https://doi.org/10.5281/zenodo.22884463>. Code licensing and the treatment of source data and artwork are set out in `LICENSE` and `NOTICE.md`.

The public data omit patient identifiers, personal dates and free-text clinical histories. Raw reads remain in the sequence archives rather than in Git.
