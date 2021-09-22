include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CELLRANGER_MKGTF {
    tag "$gtf.baseName"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    // TODO update containers and add conda recipe (if possible)
    container "litd/docker-cellranger"   // Docker image

    input:
    path gtf

    output:
    path '*.gtf'            , emit: gtf
    path '*.version.txt'    , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    cellranger mkgtf \\
        $gtf \\
        ${gtf.baseName}_mkgtf.gtf \\
        $options.args
    echo \$(cellranger --version 2>&1) | sed 's/^.*cellranger //; s/ .*\$//' > ${software}.version.txt
    """
}