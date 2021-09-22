
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MEANREADLENGTH {
    tag "$meta.id"
    label 'process_medium'
    publishDir "\${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.10 conda-forge::biopython=1.79" : null)
    container "registry.hub.docker.com/biopython/biopython"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), stdout    , emit: results
    path "*.version.txt"       , emit: version

    script:  // This script is bundled with the pipeline, in bin/
    def software = getSoftwareName(task.process)
    """
    mean_read_length.py \\
        $reads \\
        --threads \\
        $task.cpus \\
        $options.args
    echo "unversioned" > ${software}.version.txt
    """
}