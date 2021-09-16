


Set cellranger_available_references = [ "GRCh38", "mm10" ]


////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def postprocess_options                     = modules['postprocess_transpose']
def cellranger_mkref_options                = modules['cellranger_mkref']
def cellranger_mkgtf_options                = modules['cellranger_mkgtf']
def cellranger_count_options                = modules['cellranger_count']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { POSTPROCESS }               from '../../modules/local/postprocess/main'                  addParams( options: postprocess_options )
include { CELLRANGER_GETREFERENCES }  from '../../modules/local/cellranger/get_reference/main.nf'  addParams( options: [:] )
include { CELLRANGER_MKREF }          from '../../modules/local/cellranger/mkref/main.nf'          addParams( options: cellranger_mkref_options )
include { CELLRANGER_MKGTF }          from '../../modules/local/cellranger/mkgtf/main.nf'          addParams( options: cellranger_mkgtf_options )
include { CELLRANGER_COUNT }          from '../../modules/local/cellranger/count/main.nf'          addParams( options: cellranger_count_options )

/////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow CELLRANGER {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    genome_fasta        // channel: /path/to/genome.fasta
    gtf                 // channel: /path/to/annotation.gtf
    genome              // channel: name_of_prebuild_reference_genome
    cellranger_index    // channel: /path/to/cellranger/reference
    protocol            // channel: protocol

    main:
    ch_software_versions = Channel.empty()

    // Get the protocol parameter
    (cr_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "cellranger")
    
    // Create index and references
    if (!cellranger_index && genome) {
        CELLRANGER_GETREFERENCES ( genome )
        ch_reference = CELLRANGER_GETREFERENCES.out.reference
    } else if (!cellranger_index && genome_fasta && gtf) {
        CELLRANGER_MKGTF ( gtf )
        ch_gtf_filtered = CELLRANGER_MKGTF.out.gtf
        CELLRANGER_MKREF ( genome_fasta, ch_gtf_filtered )
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

    // Reformat output
    ch_cellranger_result_files = CELLRANGER_COUNT.out.results.map{ it[1] }
    // TODO it may be better to use the filtered matrix
    ch_matrix   = ch_cellranger_result_files.map { "${it}/raw_feature_bc_matrix/matrix.mtx.gz" }
    ch_features = ch_cellranger_result_files.map { "${it}/raw_feature_bc_matrix/features.tsv.gz" }
    ch_barcodes = ch_cellranger_result_files.map { "${it}/raw_feature_bc_matrix/barcodes.tsv.gz" }
    POSTPROCESS ( ch_matrix, ch_features, ch_barcodes, "cellranger" )

    emit:
    software_versions   = ch_software_versions
    multiqc_files       = ch_cellranger_result_files
}
