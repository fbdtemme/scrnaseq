include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CELLRANGER_COUNT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    // TODO update containers and add conda recipe (if possible)
    container "litd/docker-cellranger"   // Docker image

    input:
    tuple val(meta), path(reads)
    path reference
    val protocol

    output:
    tuple val(meta), path("${meta.id}/outs")    , emit: results
    path "*.version.txt"                        , emit: version

    script:
    def reference_name = reference.name
    def software = getSoftwareName(task.process)

    """
    # Make sure reads are prefixed with the meta.id

    R1="${reads[0]}"
    R2="${reads[1]}"

    if [[ ! "\$R1" =~ "^${meta.id}.*" ]]; then
        mv "\$R1" "${meta.id}_\$R1"
    fi

    if [[ ! "\$R2" =~ "^${meta.id}.*" ]]; then
        mv "\$R2" "${meta.id}_\$R2"
    fi

    cellranger count --id=${meta.id} \\
        --fastqs=. \\
        --sample=${meta.id} \\
        --transcriptome=${reference_name} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --disable-ui \\
        --chemistry=${protocol} \\
        $options.args

    cellranger --version | grep -o "[0-9\\. ]\\+" > ${software}.version.txt
    """
}
