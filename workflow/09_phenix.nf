nextflow.enable.dsl=2



params.reads_manifest = null
params.reference = null
params.phenix_config = null
params.output_dir = null
params.container_image = 'flashton/phenix-threaded-samtools:latest'
params.stage_name = 'phenix_bbduk'
params.max_forks = 8

process PHENIX {
    tag "${dataset_id}"
    maxForks (params.max_forks as int)
    container params.container_image
    publishDir "${params.output_dir}/${dataset_id}", mode: 'copy', overwrite: false

    input:
    tuple val(dataset_id), path(forward), path(reverse)
    path phenix_config
    path reference

    output:
    path "${params.stage_name}/${dataset_id}*"

    script:
    """
    set -euo pipefail
    mkdir -p ${params.stage_name}

    phenix.py prepare_reference \
      -r ${reference} \
      --mapper bwa \
      --variant gatk

    phenix.py run_snp_pipeline \
      -r1 ${forward} \
      -r2 ${reverse} \
      -r ${reference} \
      -c ${phenix_config} \
      --keep-temp \
      --sample-name ${dataset_id} \
      -o ${params.stage_name}

    phenix.py vcf2fasta \
      -i ${params.stage_name}/${dataset_id}.filtered.vcf \
      -o ${params.stage_name}/${dataset_id}_all.fasta \
      --reference ${reference}

    if [[ -f ${params.stage_name}/${dataset_id}.vcf ]]; then
      gzip -f ${params.stage_name}/${dataset_id}.vcf
    fi
    """
}

workflow {
    if (!params.reads_manifest || !params.reference || !params.phenix_config || !params.output_dir) {
        error 'Required parameters: --reads_manifest --reference --phenix_config --output_dir'
    }

    reads = Channel
        .fromPath(params.reads_manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row.analysis_id || !row.r1 || !row.r2) {
                error "Invalid PHEnix manifest row: ${row}"
            }
            tuple(row.analysis_id as String,
                  file(row.r1 as String, checkIfExists: true),
                  file(row.r2 as String, checkIfExists: true))
        }

    config_ch = Channel.value(file(params.phenix_config, checkIfExists: true))
    reference_ch = Channel.value(file(params.reference, checkIfExists: true))
    PHENIX(reads, config_ch, reference_ch)
}
