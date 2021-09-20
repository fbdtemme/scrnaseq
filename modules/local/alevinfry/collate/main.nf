
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ALEVINFRY_COLLATE {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::alevin-fry=0.4.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/alevin-fry:0.4.1--h7d875b9_0"
    } else {
        container "quay.io/biocontainers/alevin-fry:0.4.1--h7d875b9_0"
    }

    input:
    tuple val(meta), path(input)
    tuple val(meta), path(raddir)

    output:
    tuple val(meta), path(input)    , emit: results
    path "*.version.txt"            , emit: version

    script:
    def software = getSoftwareName(task.process)
    
    """
    alevin-fry collate \\
        --threads $task.cpus \\
        --input-dir $input \\
        --rad-dir $raddir \\
        $options.args

    alevin-fry --version | sed -e "s/alevin-fry //g" > ${software}.version.txt
    """
}