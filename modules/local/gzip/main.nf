// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GZIP {
    tag "$file"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path file

    output:
    path "$gzipped_file", emit: gzip
    path "*.version.txt", emit: version

    script:
    def software         = getSoftwareName(task.process)
    gzipped_file         = file.toString() + '.gz'
    
    realfile = "\$(readlink $file)"
    """
    gzip -c $options.args $realfile > $gzipped_file
    echo \$(gzip --version 2>&1) | head -1 | sed 's/gzip //g' > ${software}.version.txt
    """
}