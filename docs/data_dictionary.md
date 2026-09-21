# Data dictionary

## Publication register

`isolate_accession_line_list_345.csv` contains one row per study genome. `isolate_id` is the analytical label; `source_isolate_id` retains the original lane/sample alias. `run_accession` and `sample_accession` refer to distinct archive entities. An ERS/SRS accession is not relabelled as a BioSample accession. `sequence_type`, `lineage` and the tree-membership flags follow the frozen analysis. `CC25_order`, `ST335_order`, `CC25_panel` and `ST335_panel` link rows to figure detail sets. `Other` denotes unavailable or grouped metadata, not a new biological category.

The tree, matrix and phenotype index tables use the same analytical IDs. `additional_hospital_metadata.csv` contains the exact archived isolate-to-hospital join; unlinked records do not receive inferred hospital labels.

## Resistance tables

`phenotype_zone_diameters_42.csv` records millimetres. `phenotype_categories_as_recorded_42.csv` records the corresponding S/I/R calls. `antibiotic_dictionary_17.csv` supplies the drug names and disk contents. `scientific_config.json` gives the operational phenotype and genomic marker definitions used by the analysis. `genomic_determinants_345.csv` contains 55 explicit binary annotations; these include resistance and other curated features, so its column total is not the Figure 3B AMR-gene count.

`AMRFinder_class_count` is a stored annotation-class summary used by tree tracks. It is neither the number of unique AMR genes in Figure 3B nor the five-category XDR rule.

The original five-strain disinfectant experiment and Kenneth Chizani's separate 42-isolate experiment are separate tables. Their source units and recorded values are retained. They are not pooled across assays.

## Pairwise matrices

Rows and columns are isolate IDs, and cells are integer SNP distances. Matrices are square, symmetric and zero on the diagonal. The canonical 331-isolate matrix reproduces `ebg_25_snpmatrix.tsv` after the explicit identifier normalisation. All smaller matrices are exact subsets, not distances inferred from tree branch lengths. Tree branch lengths remain in the supplied model units.

The older Malawi plotting matrix differs from the validated matrix in 38 cells. `Malawi_matrix_cell_reconciliation.csv` records these differences. The release uses the validated complete CC25 source when calculating Malawi distances.

## NCBI comparison export

The three fields are `#AMR genotypes`, `Computed types` and `Run`. The date in the supplied filename is 2026-01-13. The archive preserves this export rather than issuing a new database query. Gzip compression is lossless. `provenance/source_inventory.csv` records the SHA-256 of the uncompressed supplied file.

Twelve wholly empty, tab-delimited physical rows are skipped, matching the R reader.

The external comparison counts matching `=COMPLETE` entries, excluding `mdsA` and `mdsB`. It retains the original row denominator and does not deduplicate repeated run IDs. The study comparison counts distinct AMR gene symbols. Unassigned Isangi sequence types remain in “Other Isangi”. The source-specific count definitions, unresolved external serovars and zero-fill convention are exposed in the input audit.
