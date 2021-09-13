
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MEAN_READ_LENGTH {
    tag "$meta.id"
    label 'process_medium'
    publishDir "\${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::python=3.8.10 conda-forge::biopython=1.79" : null)
    container "registry.hub.docker.com/biopython/biopython"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), stdout

    script:
    """
    mean_read_length.py ${reads} --threads $task.cpus
    """
}