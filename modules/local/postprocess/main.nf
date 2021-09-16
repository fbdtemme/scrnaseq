// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Reformat output file from the different tools to common formats
 */
process POSTPROCESS {
    tag "$name"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    // TODO this conda recipe is probably not working
    conda (params.enable_conda ? "conda-forge::python=3.8.10 conda-forge::scanpy" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        // TODO update containers
        container "fbdtemme/scanpyfull"
    } else {
        container "fbdtemme/scanpyfull"
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
    postprocessing.py --matrix $matrix --features $features --barcodes $barcodes --output $name $options.args
    """
}