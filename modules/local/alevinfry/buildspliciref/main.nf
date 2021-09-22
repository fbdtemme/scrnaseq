
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BUILDSPLICIREF {
    tag "$genome"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    // TODO: check if this is the right container. Support for conda, singularity and docker must be provided
    conda (params.enable_conda ? 'R-base bioconductor-eisar bioconductor-biostrings bioconductor-bsgenome r-dplyr r-stringr bioconductor-genomicfeatures r-argparser' : null)
    container "registry.hub.docker.com/combinelab/usefulaf:latest"

    input:
    path genome
    path gtf
    val read_length

    output:
    path "*.version.txt"                , emit: version
    path "*spliciref/*.fa"              , emit: reference
    path "*spliciref/*t2g.tsv"          , emit: txp2gene
    path "*spliciref/*t2g_3col.tsv"     , emit: txp2gene_3col

    script: // This script is bundled with the pipeline, in bin/
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${genome.baseName}${options.suffix}" : "${genome.baseName}"
    """
    build_splici_ref.R \\
        $genome \\
        $gtf \\
        $read_length \\
        ${prefix}_spliciref \\
        $options.args
    echo "unversioned" > ${software}.version.txt
    """
}
