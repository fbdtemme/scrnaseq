// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Reformat output file from the different tools to common formats
 */
process POSTPROCESS {
    tag "$name"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Matrices', publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path matrix
    path features
    path barcodes
    val  name

    output:
    path "$name", emit: outdir

    script:  // This script is bundled with the pipeline, in nf-core/scrnaseq/bin/
    """
    pip install --no-warn-script-location --user numba scipy loompy
    export PYTHONPATH=\"\$(python -m site --user-base)/bin\"
    postprocessing.py --matrix $matrix --features $features --barcodes $barcodes --output $name
    """
}