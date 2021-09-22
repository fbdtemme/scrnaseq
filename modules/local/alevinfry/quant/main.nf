
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ALEVINFRY_QUANT {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::alevin-fry=0.4.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/alevin-fry:0.4.1--h7d875b9_0"
    } else {
        container "quay.io/biocontainers/alevin-fry:0.4.1--h7d875b9_0"
    }

    input:
    tuple val(meta), path(input)
    path txp2gene

    output:
    tuple val(meta), path("*_quant_results")    , emit: results
    path "*.version.txt"                        , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    alevin-fry quant \\
        --threads $task.cpus \\
        --input-dir $input \\
        --tg-map $txp2gene \\
        --output-dir ${prefix}_quant_results \\
        --resolution cr-like \\
        $options.args \\
        --use-mtx
    alevin-fry --version | sed -e "s/alevin-fry //g" > ${software}.version.txt
    """
}