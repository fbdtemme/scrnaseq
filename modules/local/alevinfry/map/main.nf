
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ALEVINFRY_MAP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::salmon=1.5.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/salmon:1.5.2--h84f40af_0"
    } else {
        container "quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    }

    input:
    tuple val(meta), path(reads)
    path index
    path txp2gene
    val protocol

    output:
    path "*.version.txt"                               , emit: version
    tuple val(meta), path("*_alevin_results/*.json")   , emit: results

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    salmon alevin \\
        -l IU \\
        -p $task.cpus \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --${protocol} \\
        -i $index \\
        --tgMap $txp2gene \\
        --sketch \\
        $options.args \\
        -o ${prefix}_alevin_results
    
    salmon --version | sed -e "s/salmon //g" > ${software}.version.txt
    """
}