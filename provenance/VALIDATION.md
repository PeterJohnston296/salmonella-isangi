# Release validation

Release: 1.0.0
Date: 21 September 2026

## Executed checks

| Check | Result | Record |
|---|---|---|
| Fixed data and identifier checks | 80 passed | `data_checks.csv` |
| Automated release tests | 9 passed | `unit_test_results.txt` |
| Shell parsing | 19 scripts passed `bash -n` | Release test suite |
| R lexical delimiter checks | 17 files passed | `R_lexical_checks.csv` |
| Source image decoding | 32 images passed | `image_inventory.csv` |
| JSON/YAML/CFF parsing | Passed | Packaging validation |
| CC25 unrooted tree topology | Identical supplied bipartitions | `tree_comparison.json` |
| File-size and restricted-binary check | Passed | Release test suite |
| R parsing and frozen-data checks | Passed | `Rscript tests/check_R.R` |
| Full R figure build | 9 stages complete | `R_stage_results.csv` |
| Quarto report render | Passed; plot embedded in self-contained HTML | `results/reports/scripts/Figure_3B_AMR_gene_burden.html` |
| Figure inventory and visual inspection | 5 main and 24 supplementary indexed TIFFs passed | 600 dpi, LZW, single-frame output inspection |

The data checks cover the 345-genome publication register, 80 study-deposited runs, 327 CC25 and 224 ST335 tree tips, all fixed matrix subsets, the four matrix-only records, panel membership, 714 zone measurements, 714 recorded susceptibility calls, the 12/19 phenotypic XDR result and the 199/224 genomic profile result. The independent Figure 3B calculation reproduces denominators of 810,319, 224 and 121 and medians of 0, 17 and 1 genes.

The NCBI reader omits 12 wholly empty tab-delimited physical rows, matching the R reader's empty-row convention. It otherwise retains the supplied row denominator, including 1,137 unresolved-serovar records and 22,540 empty genotype fields. The study audit identifies 39 isolates absent from the hit-only AMRFinder input. Their existing zero-fill convention is preserved explicitly rather than treated as evidence of a completed negative run.

The two supplied CC25 tree representations have the same 651 unrooted bipartitions. Their largest stored edge-length difference is 4.89 × 10⁻⁶ model units. Both original representations are retained; their plotting orders are indexed separately.

## Execution coverage

The Python data checks, synthetic FASTA extraction/order tests, shell parsing, structured-file checks and image decoding above were executed in the packaging environment. Runtime validation used R 4.3.3 and Quarto 1.10.18. `tests/check_R.R` passed, `scripts/reproduce.R` completed all nine stages, and the Quarto report rendered as self-contained HTML with its burden plot embedded. Visual inspection covered the five indexed main figures and 24 indexed supplementary figures; labels, legends, trees, heatmaps, retained artwork and transparency backgrounds rendered as expected. Every indexed TIFF was nonempty, single-frame, LZW-compressed and recorded at 600 dpi.

The intentional repeated CC25-05 plate at supplementary positions S10 and S12 was retained and verified byte-identical. The build can leave ignored `.partial` scratch files if a device's temporary output is interrupted; these files are excluded by `.gitignore` and the release packager and are not analytical outputs.

The exact R session, installed package inventory and `renv` lockfile captured after the successful build are stored under `environment/`. Quarto and platform versions are recorded in `environment/build_environment.txt`.

No new assembly, PHEnix run, recombination analysis, tree inference or laboratory experiment was executed while constructing this archive. Publication analyses start from the included fixed outputs. Historical source commands and new execution wrappers are distinguished in `docs/upstream.md`.

## Source preservation

`source_inventory.csv` records original source hashes separately from current archive-file hashes. `SHA256SUMS.txt` covers the public release files. The NCBI gzip is a lossless copy of the supplied TSV. No raw sequencing reads, clinical linkage keys, font binaries or unpublished correspondence are included.

The submitted appendix contains a repeated CC25-05 image at Figures 10 and 12. The archive records those two figure positions but counts the panel's isolates only once. Experimental/plasmid source artwork and the retained map are not described as regenerated raw observations.
