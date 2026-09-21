# Upstream analysis

The project was run as separate analyses rather than one uninterrupted pipeline. The wrappers below preserve those stages and write their own command, input and software records. They use Bash 4 or later, Python 3 and the named bioinformatics executables. Run them on Linux or the original analysis server; the R figure workflow can also run on macOS.

## Stage order

| Stage | Entry point | Input | Output |
|---|---|---|---|
| Public reads | `01_download_public_reads.sh` | Explicit archive-run list | Paired FASTQ files and resolved-run manifest |
| Trimming | `02_trim_bbduk.sh` | Read manifest and adapter FASTA | Trimmed read pairs |
| Assembly | `03_assemble_spades.sh` | Trimmed read manifest | SPAdes assemblies |
| Assembly assessment | `04_qc_quast.sh` | Assemblies | QUAST results and QC records |
| Serovar and ST | `05_type_sistr_mlst.sh` | Assemblies | SISTR and MLST calls |
| Resistance | `06_amrfinder.sh` | Assemblies and AMRFinder database | Per-isolate AMRFinder reports |
| Reference | `07_build_reference.sh` | Long reads and paired short reads | Flye assembly with Medaka/pypolca polishing |
| PHEnix inputs | `08_prepare_phenix_inputs.sh` | Selected reads, reference and YAML | PHEnix input manifest |
| SNP calling | `09_run_phenix.sh`, `09_phenix.nf` | PHEnix manifest | Reference-coordinate calls and FASTA files |
| Alignment | `10_build_alignments.sh` | Explicit subset list and PHEnix FASTAs | Ordered alignment with the reference record removed |
| SNP distances | `11_snp_distances.sh` | Alignment | Numeric pairwise distances |
| Recombination | `12_gubbins.sh` | Alignment | Gubbins outputs |
| Tree inference | `13_raxml.sh` or `13_raxmlng.sh` | Specified alignment | Separate tree-inference output |
| Population labels | `14_import_rhierbaps_assignments.sh` | Frozen assignment table | Identifier-matched population labels |
| Plasmids | `20_mob_typer.sh`, `21_build_hybrid_plasmid_assembly.sh` | Assemblies or matched long/short reads | Replicon and hybrid-assembly outputs |

Each script accepts its stage inputs as arguments and shared settings through `00_config.sh` or environment variables. Outputs are written under `work/`, separate from the fixed files in `data/`. Later stages can start from an existing intermediate; preceding stages do not have to be repeated.

## Historical analysis settings

The main trimming command used `ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=r trimq=20 minlength=50`. The assembly notes use paired reads and SPAdes `--careful`. QUAST was run without a reference. The later AMRFinder command specifies `--organism Salmonella` and database `2024-05-02.2`, without `--plus`.

The reference-building record names `CHF10J1_final_polished_reference.fasta`, following Flye, two Medaka passes and two pypolca passes. The analytical reference-isolate label in the current project is CH10FJ. The filename is not used to infer a sequence accession.

PHEnix used Nextflow DSL2, BWA mapping, GATK variant calling and an external YAML, followed by `vcf2fasta`. The recorded container tag is `flashton/phenix-threaded-samtools:latest`. The historical YAML, immutable image digest and final reference sequence are not part of this archive. The runnable PHEnix wrapper accepts `PHENIX_CONFIG`, `PHENIX_CONTAINER` and `REFERENCE_FASTA` from the caller and records their provenance. It is not used by the fixed-output publication analysis.

The late notebook records classic `raxmlHPC-PTHREADS`, with `GTRGAMMA`, seeds 12345 and 1,000 rapid bootstrap replicates. The separate RAxML-NG implementation uses the manuscript's specified model and bootstrap target. A tree produced by a new invocation is a new output; agreement with the publication trees must be assessed rather than assumed.

The frozen RHierBAPS labels are included in `data/data_clean/rhierbaps_assignments_327.csv`. Plasmid work is represented by original command records, curated outputs and source artwork. Pseudogene results are described in the manuscript and are not recalculated by this release. The archive does not reconstruct experimental observations from figure pixels.

## Cohort records

The historical selection lists follow assembled candidates, serovar review, assembly QC, PHEnix entry and the final tree. The publication register provides a separate current mapping to run accessions and the 224-isolate ST335 subset. Raw-input aliases and cleaned plotting IDs are distinct columns. Use `st335_224.upstream_aliases.txt` for original PHEnix directories and `st335_224.txt` for the cleaned plotting identifiers.

The 327-isolate tree, 331-isolate complete matrix and 345-isolate study register serve different purposes. Their membership is fixed in the manifests. Four matrix-only records are not retrospectively assigned a biological exclusion reason.
