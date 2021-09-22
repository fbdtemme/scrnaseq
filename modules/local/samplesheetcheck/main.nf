// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process SAMPLESHEETCHECK {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path samplesheet
    
    output:
    path '*.csv'

    // TODO add options.args
    // TODO generate version output with software name
    // TODO use options.suffix for the output file
    script:  // This script is bundled with the pipeline, in bin/
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}