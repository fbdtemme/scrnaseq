
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BUILD_SPLICI_REF {
    tag "$genome"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "registry.hub.docker.com/combinelab/usefulaf:latest"

    input:
    path genome
    path gtf
    val read_length

    output:
    path "*.version.txt"                          , emit: version
    path "alevinfry_splici_ref/*.fa"              , emit: reference
    path "alevinfry_splici_ref/*t2g.tsv"          , emit: txp2gene  
    path "alevinfry_splici_ref/*t2g_3col.tsv"     , emit: txp2gene_3col  


    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${genome}${options.suffix}" : "${genome}"
    
    """
    build_splici_ref.R ${genome} ${gtf} ${read_length} "alevinfry_splici_ref" ${extra_args} ${options.args}

    echo "unversioned" > ${software}.version.txt
    """
}