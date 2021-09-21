////////////////////////////////////////////////////
/* --         CELLRANGER SUBWORKFLOW        -- */
////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def cellranger_mkref_options           = modules['cellranger_mkref']
def cellranger_mkgtf_options           = modules['cellranger_mkgtf']
def cellranger_count_options           = modules['cellranger_count']
def postprocess_options                = modules['postprocess']
postprocess_options.publish_dir        = 'cellranger'

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { CELLRANGER_MKREF }          from '../../modules/local/cellranger/mkref/main.nf'         addParams( options: cellranger_mkref_options )
include { CELLRANGER_MKGTF }          from '../../modules/local/cellranger/mkgtf/main.nf'         addParams( options: cellranger_mkgtf_options )
include { CELLRANGER_COUNT }          from '../../modules/local/cellranger/count/main.nf'         addParams( options: cellranger_count_options )
include { POSTPROCESS }               from '../../modules/local/postprocess/main'                 addParams( options: postprocess_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { UNTAR }                     from '../../modules/nf-core/modules/untar/main.nf'          addParams( options: [:] )

workflow CELLRANGER {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    genome_fasta        // channel: /path/to/genome.fasta
    gtf                 // channel: /path/to/annotation.gtf
    cellranger_index    // channel: /path/to/cellranger/reference
    protocol            // channel: protocol

    main:

    // Get the protocol parameter
    (cr_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "cellranger")
    
    // Check if a tar.gz packaged index was passed such as a url to the 10x references online
    if (!cellranger_index && genome_fasta && gtf) {
        CELLRANGER_MKGTF ( gtf )
        ch_gtf_filtered = CELLRANGER_MKGTF.out.gtf
        CELLRANGER_MKREF ( genome_fasta, ch_gtf_filtered )
        ch_reference = CELLRANGER_MKREF.out.reference
    } else {
        ch_reference = cellranger_index
    }

    // Quantify with cellranger
    CELLRANGER_COUNT (
        reads,
        ch_reference,
        cr_protocol
    )

    // Collect software versions
    ch_software_versions       = Channel.empty()
    ch_software_versions       = ch_software_versions.mix(CELLRANGER_COUNT.out.version.ifEmpty(null))

    // Reformat output and run postprocess module
    // TODO there may be a cleaner way of doing this
    // Meta and result could be set in one command and the matrix, barcodes and features could
    // be mixed into one channel
    ch_cellranger_result_files = CELLRANGER_COUNT.out.results.map{ it[1] }
    ch_meta                    = CELLRANGER_COUNT.out.results.map{ it[0] }
    ch_matrix                  = ch_cellranger_result_files.map { "${it}/raw_feature_bc_matrix/matrix.mtx.gz" }
    ch_features                = ch_cellranger_result_files.map { "${it}/raw_feature_bc_matrix/features.tsv.gz" }
    ch_barcodes                = ch_cellranger_result_files.map { "${it}/raw_feature_bc_matrix/barcodes.tsv.gz" }
    POSTPROCESS ( ch_meta, ch_matrix, ch_features, ch_barcodes )

    // Collect software versions
    ch_software_versions       = ch_software_versions.mix(POSTPROCESS.out.version.first().ifEmpty(null))

    emit:
    software_versions          = ch_software_versions
    multiqc_files              = ch_cellranger_result_files
}
