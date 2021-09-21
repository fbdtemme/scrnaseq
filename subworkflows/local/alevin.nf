////////////////////////////////////////////////////
/* --         SALMON ALEVIN SUBWORKFLOW        -- */
////////////////////////////////////////////////////

// Whitelist files for STARsolo and Kallisto
def whitelist_folder = "$baseDir/assets/whitelist/"

////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def gffread_txp2gene_options            = modules['gffread_tx2pgene']
def gffread_transcriptome_options       = modules['gffread_transcriptome']
def salmon_index_options                = modules['salmon_index']
def salmon_alevin_options               = modules['salmon_alevin']
salmon_alevin_options.args             += ' --dumpFeatures --dumpMtx'
def salmon_alevinqc_options             = modules['salmon_alevinqc']
def postprocess_options                 = modules['postprocess']
postprocess_options.publish_dir         = 'salmon/alevin'
postprocess_options.args               += ' --transpose'
def gunzip_options                      = modules['gunzip']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { GFFREAD_TRANSCRIPTOME }         from '../../modules/local/gffread/transcriptome/main'   addParams( options: gffread_transcriptome_options )
include { SALMON_ALEVIN }                 from '../../modules/local/salmon/alevin/main'           addParams( options: salmon_alevin_options )
include { SALMON_ALEVINQC }               from '../../modules/local/salmon/alevinqc/main'         addParams( options: salmon_alevinqc_options )
include { POSTPROCESS }                   from '../../modules/local/postprocess/main'             addParams( options: postprocess_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                        from '../../modules/nf-core/modules/gunzip/main'        addParams( options: gunzip_options )
include { GFFREAD as GFFREAD_TXP2GENE }   from '../../modules/nf-core/modules/gffread/main'       addParams( options: gffread_txp2gene_options )
include { SALMON_INDEX }                  from '../../modules/nf-core/modules/salmon/index/main'  addParams( options: salmon_index_options )
include { UNTAR }                         from '../../modules/nf-core/modules/untar/main.nf'      addParams( options: [:] )

workflow ALEVIN {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    genome_fasta        // channel: /path/to/genome.fasta
    transcript_fasta    // channel: /path/to/transcriptome/fasta
    gtf                 // channel: /path/to/annotation.gtf
    txp2gene            // channel: /path/to/txp2gene.txt
    salmon_index        // channel: /path/to/salmon/index
    protocol            // channel: protocol

    main:
    ch_software_versions = Channel.empty()

    // Get the protocol parameter suitable for passing to alevin
    (alevin_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "alevin")

    // Preprocessing - Extract transcriptome fasta from genome fasta
    if (!transcript_fasta && genome_fasta && gtf) {
        GFFREAD_TRANSCRIPTOME ( genome_fasta, gtf )
        transcript_fasta = GFFREAD_TRANSCRIPTOME.out.transcriptome_extracted
        ch_software_versions = ch_software_versions.mix(GFFREAD_TRANSCRIPTOME.out.version.first().ifEmpty(null))
    }
    
    // Set up salmon index
    if (!salmon_index) {
        SALMON_INDEX ( genome_fasta, transcript_fasta )
        ch_salmon_alevin_index = SALMON_INDEX.out.index
    } else {
        ch_salmon_alevin_index = salmon_index
    }
    
    // Build txp2gene map
    if (!txp2gene){
        GFFREAD_TXP2GENE ( gtf )
        ch_txp2gene = GFFREAD_TXP2GENE.out.gtf
        // Only collect version if not already done for gffread
        if (!GFFREAD_TRANSCRIPTOME) {
            ch_software_versions = ch_software_versions.mix(GFFREAD_TXP2GENE.out.version.first().ifEmpty(null))
        }
    }
    
    // Perform quantification with salmon alevin
    // TODO: Verify the correct library type:
    // https://salmon.readthedocs.io/en/latest/library_type.html 
    // for different types of chemistry and protocol
    def lib_type = "ISR"
    SALMON_ALEVIN ( 
        reads, 
        ch_salmon_alevin_index, 
        ch_txp2gene, 
        alevin_protocol,
        lib_type
    )
    
    // Collect software versions
    ch_software_versions    = ch_software_versions.mix(SALMON_ALEVIN.out.version.first().ifEmpty(null))

    // Run alevinQC
    SALMON_ALEVINQC ( SALMON_ALEVIN.out.results )

    // Collect software versions
    ch_software_versions    = ch_software_versions.mix(SALMON_ALEVINQC.out.version.first().ifEmpty(null))

    // Reformat output and run postprocess module
    // TODO there may be a cleaner way of doing this
    // Meta and result could be set in one command and the matrix, barcodes and features could
    // be mixed into one channel
    ch_alevin_results_files = SALMON_ALEVIN.out.results.map{ it[1] }
    ch_meta                 = SALMON_ALEVIN.out.results.map{ it[0] }
    ch_matrix               = ch_alevin_results_files.map{ "${it}/alevin/quants_mat.mtx.gz" }
    ch_features             = ch_alevin_results_files.map{ "${it}/alevin/quants_mat_cols.txt" }
    ch_barcodes             = ch_alevin_results_files.map{ "${it}/alevin/quants_mat_rows.txt" }
    POSTPROCESS ( ch_meta, ch_matrix, ch_features, ch_barcodes )
    
    // Collect software versions
    ch_software_versions    = ch_software_versions.mix(POSTPROCESS.out.version.first().ifEmpty(null))
    
    emit:
    software_versions       = ch_software_versions
    multiqc_files           = ch_alevin_results_files
}
