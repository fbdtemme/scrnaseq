////////////////////////////////////////////////////
/* --         SALMON ALEVIN SUBWORKFLOW        -- */
////////////////////////////////////////////////////


// Whitelist files for STARsolo and Kallisto
def whitelist_folder = "$baseDir/assets/whitelist/"


////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def salmon_index_options                = modules['salmon_index']
def gffread_txp2gene_options            = modules['gffread_tx2pgene']
def gffread_transcriptome_options       = modules['gffread_transcriptome']
def salmon_alevin_options               = modules['salmon_alevin']
def alevin_qc_options                   = modules['alevinqc']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { GFFREAD_TRANSCRIPTOME }       from '../../modules/local/gffread/transcriptome/main'   addParams( options: gffread_transcriptome_options )
include { SALMON_ALEVIN }               from '../../modules/local/salmon/alevin/main'           addParams( options: salmon_alevin_options )
include { ALEVINQC }                    from '../../modules/local/salmon/alevinqc/main'         addParams( options: alevin_qc_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                      from '../../modules/nf-core/modules/gunzip/main'       addParams( options: [:] )
include { GFFREAD as GFFREAD_TXP2GENE } from '../../modules/nf-core/modules/gffread/main'      addParams( options: gffread_txp2gene_options )
include { SALMON_INDEX }                from '../../modules/nf-core/modules/salmon/index/main' addParams( options: salmon_index_options )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow ALEVIN {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    genome_fasta        // channel: /path/to/genome.fasta
    transcript_fasta    // channel: /path/to/transcriptome/fasta
    gtf                 // channel: /path/to/annotation.gtf
    txp2gene            // channel: /path/to/txp2gene.txt
    salmon_index        // channel: /path/to/salmon/index
    protocol            // channel: protocol
    barcode_whitelist   // channel: /path/to/barcode_whitelist.txt

    main:
    ch_software_versions = Channel.empty()

    // Get the protocol parameter suitable for passing to alevin
    (alevin_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "alevin")

    // Setup correct barcode whitelist and unzip of necessary
    if (protocol.contains("10X") && !barcode_whitelist) {
        barcode_filename = "$whitelist_folder/10x_${chemistry}_barcode_whitelist.txt.gz"
        Channel.fromPath(barcode_filename)
        .ifEmpty{ exit 1, "Cannot find ${protocol} barcode whitelist: $barcode_filename" }
        .set{ barcode_whitelist_gzipped }

        GUNZIP ( barcode_whitelist_gzipped )
        ch_barcode_whitelist = GUNZIP.out.gunzip
    } else {
         Channel.fromPath(barcode_whitelist)
        .ifEmpty{ exit 1, "Cannot find ${protocol} barcode whitelist: $barcode_filename" }
        .set{ ch_barcode_whitelist }
    }

    // Preprocessing - Extract transcriptome fasta from genome fasta
    if (!transcript_fasta && genome_fasta && gtf) {
        GFFREAD_TRANSCRIPTOME ( genome_fasta, gtf )
        transcriptome_fasta = GFFREAD_TRANSCRIPTOME.out.transcriptome_extracted
        ch_software_versions = ch_software_versions.mix(GFFREAD_TRANSCRIPTOME.out.version.ifEmpty(null))
    }
    
    // Build salmon index
    if (!salmon_index) {
        SALMON_INDEX ( genome_fasta, transcriptome_fasta )
        salmon_index_alevin = SALMON_INDEX.out.index
    } else {
        // Setup channel for salmon index if specified
        salmon_index_alevin = Channel.fromPath(salmon_index)
    }
    
    // Build txp2gene map
    if (!txp2gene){
        GFFREAD_TXP2GENE ( gtf )
        ch_txp2gene = GFFREAD_TXP2GENE.out.gtf
        // Only collect version if not already done for gffread
        if (!GFFREAD_TRANSCRIPTOME.out) {
            ch_software_versions = ch_software_versions.mix(GFFREAD_TXP2GENE.out.version.ifEmpty(null))
        }
    }
    
    // Perform quantification with salmon alevin
    SALMON_ALEVIN ( 
        reads, 
        salmon_index_alevin, 
        ch_txp2gene, 
        alevin_protocol, 
        ch_barcode_whitelist 
    )
    
    // Collect software versions
    ch_software_versions = ch_software_versions.mix(SALMON_ALEVIN.out.version.first().ifEmpty(null))

    // Run alevinQC
    ALEVINQC ( SALMON_ALEVIN.out.alevin_results )
    
    // Collect software versions
    ch_software_versions = ch_software_versions.mix(ALEVINQC.out.version.first().ifEmpty(null))

    // Collect multiqc files
    ch_salmon_multiqc   = SALMON_ALEVIN.out.alevin_results
    ch_multiqc_files    = Channel.empty()
    ch_multiqc_files    = ch_multiqc_files.mix(ch_salmon_multiqc.collect{it[1]}.ifEmpty([]))
    
    emit:
    software_versions   = ch_software_versions
    multiqc_files       = ch_multiqc_files
}
