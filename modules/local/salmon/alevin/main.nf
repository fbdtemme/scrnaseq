
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SALMON_ALEVIN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

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
    val lib_type

    output:
    tuple val(meta), path("*_alevin_results"), emit: results
    path "*.version.txt"                     , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def strandedness_opts = [
        'A', 'U', 'SF', 'SR',
        'IS', 'IU' , 'ISF', 'ISR',
        'OS', 'OU' , 'OSF', 'OSR',
        'MS', 'MU' , 'MSF', 'MSR'
    ]
    def strandedness =  'A'
    if (lib_type) {
        if (strandedness_opts.contains(lib_type)) {
            strandedness = lib_type
        } else {
            log.info "[Salmon Alevin] Invalid library type specified '--libType=${lib_type}', defaulting to auto-detection with '--libType=A'."
        }
    } else {
        strandedness = meta.single_end ? 'U' : 'IU'
        if (meta.strandedness == 'forward') {
            strandedness = meta.single_end ? 'SF' : 'ISF'
        } else if (meta.strandedness == 'reverse') {
            strandedness = meta.single_end ? 'SR' : 'ISR'
        }
    }
    """
    salmon alevin \\
        -p $task.cpus \\
        --libType $strandedness \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --${protocol} \\
        -i $index \\
        --tgMap $txp2gene \\
        $options.args \\
        -o ${prefix}_alevin_results
    salmon --version | sed -e "s/salmon //g" > ${software}.version.txt
    """
}