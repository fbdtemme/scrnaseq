include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CELLRANGER_MKREF {
    tag "${fasta.baseName}"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    // TODO update containers and add conda recipe (if possible)
    container "litd/docker-cellranger"   // Docker image

    input:
    path fasta
    path gtf

    output:
    path("cellranger"),   emit: reference
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    cellranger mkref \\
        --genome=cellranger \\
        --fasta=${fasta} \\
        --genes=${gtf} \\
        --memgb ${task.memory.giga} \\
        --nthreads ${task.cpus}
    cellranger --version | grep -o "[0-9\\. ]\\+" > ${software}.version.txt
    """
}
