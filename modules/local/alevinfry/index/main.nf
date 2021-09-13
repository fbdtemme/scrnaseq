
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ALEVINFRY_INDEX {
    tag "$splici_transcriptome"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::salmon=1.5.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/salmon:1.5.2--h84f40af_0"
    } else {
        container "quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    }

    input:
    path splici_transcriptome

    output:
    path "alevinfry"               , emit: index
    path "*.version.txt"           , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${splici_transcriptome}${options.suffix}" : "${splici_transcriptome}"

    """
     salmon \\
        index \\
        --threads $task.cpus \\
        -t $splici_transcriptome \\
        -i alevinfry \\
        ${options.args} \\

    salmon --version | sed -e "s/salmon //g" > salmon.version.txt
    """
}