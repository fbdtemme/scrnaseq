////////////////////////////////////////////////////
/* --         SALMON ALEVIN SUBWORKFLOW        -- */
////////////////////////////////////////////////////


// Whitelist files for STARsolo and Kallisto
def whitelist_folder = "$baseDir/assets/whitelist/"


////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def gffread_txp2gene_options                = modules['gffread_tx2pgene']
def gffread_transcriptome_options           = modules['gffread_transcriptome']
def alevin_qc_options                       = modules['alevinqc']
def alevinfry_index_options                 = modules['alevinfry_index']
def alevinfry_map_options                   = modules['alevinfry_map']
def alevinfry_generate_permitlist_options   = modules['alevinfry_permitlist']
def alevinfry_collate_options               = modules['alevinfry_collate']
def alevinfry_quant_options                 = modules['alevinfry_quant']


////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { GFFREAD_TRANSCRIPTOME }           from '../../modules/local/gffread/transcriptome/main'           addParams( options: gffread_transcriptome_options )
include { ALEVINQC }                        from '../../modules/local/salmon/alevinqc/main'                 addParams( options: alevin_qc_options )
include { POSTPROCESS }                     from '../../modules/local/postprocess/main'                     addParams( options: [:] )
include { MEAN_READ_LENGTH }                from '../../modules/local/mean_read_length/main'                addParams( options: [:] )
include { BUILD_SPLICI_REF }                from '../../modules/local/alevinfry/build_splici_ref/main'      addParams( options: [:] )
include { ALEVINFRY_INDEX }                 from '../../modules/local/alevinfry/index/main'                 addParams( options: alevinfry_index_options )
include { ALEVINFRY_MAP }                   from '../../modules/local/alevinfry/map/main'                   addParams( options: alevinfry_map_options )
include { ALEVINFRY_GENERATE_PERMITLIST }   from '../../modules/local/alevinfry/generate_permitlist/main'   addParams( options: alevinfry_generate_permitlist_options )
include { ALEVINFRY_COLLATE }               from '../../modules/local/alevinfry/collate/main'               addParams( options: alevinfry_collate_options )
include { ALEVINFRY_QUANT }                 from '../../modules/local/alevinfry/quant/main'                 addParams( options: alevinfry_quant_options )



////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                      from '../../modules/nf-core/modules/gunzip/main'       addParams( options: [:] )
include { GFFREAD as GFFREAD_TXP2GENE } from '../../modules/nf-core/modules/gffread/main'      addParams( options: gffread_txp2gene_options )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow ALEVINFRY {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    genome_fasta        // channel: /path/to/genome.fasta
    gtf                 // channel: /path/to/annotation.gtf
    protocol            // channel: protocol
    barcode_whitelist   // channel: /path/to/barcode_whitelist.txt

    main:
    ch_software_versions = Channel.empty()
    ch_multiqc_files     = Channel.empty()

    // Get the protocol parameter suitable for passing to alevin and alevin-fry
    (alevin_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "alevin")

    // Flatten input for getting the mean read length
    reads
        .map { meta, reads -> [ reads[0] ] }
        .collect()
        .map { reads -> [ ["id": "FW"], reads ] }
        .set { ch_all_fw }
    reads
        .map { meta, reads -> [ reads[1] ] }
        .collect()
        .map { reads -> [ ["id": "RV"], reads ] }
        .set { ch_all_rv }

    ch_all_reads = ch_all_fw.mix ( ch_all_rv )

    // Get mean read length of fowards and reverse reads and select the second reads
    MEAN_READ_LENGTH ( ch_all_reads )
    ch_read_length =  MEAN_READ_LENGTH.out
        .toSortedList()
        .map { it[1][1] }
    
    BUILD_SPLICI_REF(
        genome_fasta,
        gtf,
        ch_read_length
    )

    ch_splici_ref    = BUILD_SPLICI_REF.out.reference
    ch_txp2gene_3col = BUILD_SPLICI_REF.out.txp2gene_3col

    // Build salmon/alevin index

    ALEVINFRY_INDEX ( ch_splici_ref )
    index = ALEVINFRY_INDEX.out.index

    ALEVINFRY_MAP(
        reads,
        index,  
        ch_txp2gene_3col,
        alevin_protocol
    )
    rad_dir = ALEVINFRY_MAP.out.results
    

    // Build permitlist and filter index
    def expected_orientation = "fw"
    ALEVINFRY_GENERATE_PERMITLIST( rad_dir, expected_orientation )
    quant_dir = ALEVINFRY_GENERATE_PERMITLIST.out.quant

    ALEVINFRY_COLLATE ( quant_dir, rad_dir )
    quant_dir = ALEVINFRY_COLLATE.out.results

    // Perform quantification with alevin-fry quant
    ALEVINFRY_QUANT ( quant_dir, ch_txp2gene_3col )

    // Reformat output
    ch_alevin_results_files = ALEVINFRY_QUANT.out.results
    ch_alevin_output_dir = ch_alevin_results_files.map{it[1]}
    
    ch_matrix   = ch_alevin_output_dir.map { "${it}/alevin/quants_mat.mtx" }
    ch_barcodes = ch_alevin_output_dir.map { "${it}/alevin/quants_mat_cols.txt" }
    ch_features = ch_alevin_output_dir.map { "${it}/alevin/quants_mat_rows.txt" }
    POSTPROCESS ( ch_matrix, ch_barcodes, ch_features, "Alevin" )
    
    // Collect software versions
    ch_software_versions = ch_software_versions.mix(ALEVINFRY_INDEX.out.version.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(ALEVINFRY_QUANT.out.version.ifEmpty(null))

    emit:
    software_versions   = ch_software_versions
    multiqc_files       = ch_multiqc_files
}
