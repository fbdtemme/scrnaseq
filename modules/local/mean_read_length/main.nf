
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MEAN_READ_LENGTH {
    tag "$meta.id"
    label 'process_medium'
    publishDir "\${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::seqkit=2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:2.0.0--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/seqkit:2.0.0--h9ee0642_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), stdout

    script:
    """
    reads=()
    reads+=( ${reads} )

    count=0;
    sum=0;

    for r in \$reads; do
        r=\$( seqkit stats -T "\$r" |  sed -n 2p | cut -d\$'\t' -f7 | xargs echo -n );
        r_int=\${r%.*}
        (( sum+=r_int ))
        (( ++count ));
    done;

    echo -n \$(( sum / count ))
    """
}