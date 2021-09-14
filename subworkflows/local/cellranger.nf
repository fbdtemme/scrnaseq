


Set available_references = [ "GRCh38", "mm10" ]

include { CELLRANGER_MKREF } from '../../modules/local/cellranger/mkref/main.nf' addParams( options: [:] )
include { CELLRANGER_COUNT } from '../../modules/local/cellranger/count/main.nf' addParams( options: [:] )


/////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow CELLRANGER {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    genome_fasta        // channel: /path/to/genome.fasta
    gtf                 // channel: /path/to/annotation.gtf
    cellranger_index    // channel: /path/to/cellranger/reference
    protocol            // channel: protocol

    main:
    ch_software_versions = Channel.empty()

    // Get the protocol parameter
    (cr_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "cellranger")
    
    if (!cellranger_index && genome_fasta && gtf) {
        CELLRANGER_MKREF ( genome_fasta, gtf )
        ch_reference = CELLRANGER_MKREF.out.reference
    } else {
        ch_reference = cellranger_index
    }

    CELLRANGER_COUNT (
        reads,
        ch_reference,
        cr_protocol
    )

    ch_software_versions = ch_software_versions.mix(CELLRANGER_COUNT.out.version.ifEmpty(null))
}

// Functions needed by the workflow

def get_meta_tabs(arr) {
    def meta = [:]
    meta.gem          = arr[0]
    meta.samples      = arr[1]

    def array = []
    array = [ meta, arr[2].flatten() ]
    return array
}