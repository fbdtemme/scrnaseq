
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ALEVINFRY_GENERATE_PERMITLIST {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:meta.id) }

    conda (params.enable_conda ? 'bioconda::alevin-fry=0.4.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/alevin-fry:0.4.1--h7d875b9_0"
    } else {
        container "quay.io/biocontainers/alevin-fry:0.4.1--h7d875b9_0"
    }

    input:
    tuple val(meta), path(map)

    output:
    path "*.version.txt"                       , emit: version
    tuple val(meta), path("*_alevin_quant")  , emit: quant

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    
    // alevin-fry complains when passing empty options.args
    if (options.args && !options.args.isEmpty()) {
        """
        alevin-fry generate-permit-list \\
            --knee-distance \\
            --expected-ori fw \\
            --input ${map} \\ 
            --output-dir ${prefix}_alevin_quant \\ 
            ${options.args}

        alevin-fry --version | sed -e "s/alevin-fry //g" > ${software}.version.txt
        """
    } else {
        """
        alevin-fry generate-permit-list --knee-distance \\
            --expected-ori fw \\
            --input ${map} \\
            --output-dir ${prefix}_alevin_quant
        
        alevin-fry --version | sed -e "s/alevin-fry //g" > ${software}.version.txt
        """
    }

}