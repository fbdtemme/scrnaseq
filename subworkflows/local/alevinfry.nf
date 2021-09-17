////////////////////////////////////////////////////
/* --         ALEVIN-FRY SUBWORKFLOW        -- */
////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def alevinfry_index_options                 = modules['alevinfry_index']
def alevinfry_map_options                   = modules['alevinfry_map']
def alevinfry_generate_permitlist_options   = modules['alevinfry_permitlist']
def alevinfry_collate_options               = modules['alevinfry_collate']
def alevinfry_quant_options                 = modules['alevinfry_quant']
def postprocess_options                     = modules['postprocess_transpose']
def gunzip_options                          = modules['gunzip']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { ALEVINFRY_BUILD_INDEX }           from '../../subworkflows/local/alevinfry_build_index'           addParams( options: [:] )
include { ALEVINFRY_MAP }                   from '../../modules/local/alevinfry/map/main'                   addParams( options: alevinfry_map_options )
include { ALEVINFRY_GENERATE_PERMITLIST }   from '../../modules/local/alevinfry/generate_permitlist/main'   addParams( options: alevinfry_generate_permitlist_options )
include { ALEVINFRY_COLLATE }               from '../../modules/local/alevinfry/collate/main'               addParams( options: alevinfry_collate_options )
include { ALEVINFRY_QUANT }                 from '../../modules/local/alevinfry/quant/main'                 addParams( options: alevinfry_quant_options )
include { POSTPROCESS }                     from '../../modules/local/postprocess/main'                     addParams( options: postprocess_options )
include { UNTAR }                           from '../../modules/nf-core/modules/untar/main.nf'              addParams( options: [:] )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                          from '../../modules/nf-core/modules/gunzip/main'                addParams( options: gunzip_options )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow ALEVINFRY {
    take:
    reads                   // channel: [ val(meta), [ reads ] ]
    genome_fasta            // channel: /path/to/genome.fasta
    gtf                     // channel: /path/to/annotation.gtf
    txp2gene                // channel: /path/to/txp2gene
    alevinfry_index         // channel: /path/to/alevinfry_index
    protocol                // channel: protocol
    expected_orientation    // channel:  (the options are [ fw | rc | both ]

    main:
    ch_software_versions = Channel.empty()
    ch_multiqc_files     = Channel.empty()

    // Get the protocol parameter suitable for passing to alevin and alevin-fry
    (alevin_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "alevin")

    // Build index and txp2gene mapping if no index is provided
    if (!alevinfry_index) { 
        ALEVINFRY_BUILD_INDEX ( reads, genome_fasta, gtf )
        ch_index = ALEVINFRY_BUILD_INDEX.out.index
        ch_txp2gene = ALEVINFRY_BUILD_INDEX.out.txp2gene_3col
        ch_software_versions.mix(ALEVINFRY_BUILD_INDEX.out.software_versions.ifEmpty(null))
    } else if (alevinfry_index && alevinfry_index.getName().endsWith(".tar.gz")) {
        ch_index = UNTAR ( alevinfry_index ).untar
        ch_txp2gene = txp2gene
        ch_software_versions.mix(UNTAR.out.version.ifEmpty(null))
    } else {
        ch_index = alevinfry_index
        ch_txp2gene = txp2gene
    }

    // Align reads
    ALEVINFRY_MAP(
        reads,
        ch_index,  
        ch_txp2gene,
        alevin_protocol
    )

    // Build permitlist and filter index
    ALEVINFRY_GENERATE_PERMITLIST( ALEVINFRY_MAP.out.results, expected_orientation )
    ALEVINFRY_COLLATE ( ALEVINFRY_GENERATE_PERMITLIST.out.quant, ALEVINFRY_MAP.out.results )

    // Perform quantification with alevin-fry quant
    ALEVINFRY_QUANT ( ALEVINFRY_COLLATE.out.results, ch_txp2gene )

    // Reformat output
    ch_alevin_output_dir    = ALEVINFRY_QUANT.out.results.map{ it[1] }
    ch_matrix               = ch_alevin_output_dir.map{ "${it}/alevin/quants_mat.mtx" }
    ch_features             = ch_alevin_output_dir.map{ "${it}/alevin/quants_mat_cols.txt" }
    ch_barcodes             = ch_alevin_output_dir.map{ "${it}/alevin/quants_mat_rows.txt" }
    POSTPROCESS ( ch_matrix, ch_features, ch_barcodes, "Alevinfry" )
    
    // Collect software versions
    ch_software_versions = ch_software_versions.mix(ALEVINFRY_QUANT.out.version.ifEmpty(null))

    emit:
    software_versions   = ch_software_versions
    multiqc_files       = ch_alevin_output_dir
}
