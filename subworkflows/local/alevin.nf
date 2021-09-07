////////////////////////////////////////////////////
/* --         SALMON ALEVIN WORKFLOW           -- */
////////////////////////////////////////////////////


////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
}

// Check if GTF is supplied properly
if (params.gtf) {
    Channel
    .fromPath(params.gtf)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    .set { gtf }
}

// Check if TXP2Gene is provided for Alevin
if (!params.gtf && !params.txp2gene) {
    exit 1, "Must provide either a GTF file ('--gtf') or transcript to gene mapping ('--txp2gene') to quantify with Alevin"
}

// Setup Genome Fasta channels
if (params.genome_fasta) {
    Channel
    .fromPath(params.genome_fasta)
    .ifEmpty { exit 1, "Fasta file not found: ${params.genome_fasta}" }
    .set { genome_fasta }
} 

// Setup Transcript Fasta channels
if (params.transcript_fasta) {
  Channel
    .fromPath(params.transcript_fasta)
    .ifEmpty { exit 1, "Fasta file not found: ${params.transcript_fasta}" }
    .set { transcriptome_fasta }
} 

// Check if files for index building are given if no index is specified
if (!params.salmon_index && !params.genome_fasta) {
    exit 1, "Must provide a genome fasta file ('--genome_fasta') or a transcript fasta ('--transcript_fasta') if no index is given!"
}

// Setup channel for salmon index if specified
if (params.salmon_index) {
    Channel
    .fromPath(params.salmon_index)
    .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_index}" }
    .set { salmon_index_alevin }
}


// Check if txp2gene file has been provided and create channel
if (params.txp2gene) {
    Channel
    .fromPath(params.txp2gene)
    .set{ ch_txp2gene } 
}

// Check AWS batch settings
// TODO use the Checks.awsBatch() function instead

// Get the protocol parameter
(protocol, chemistry) = Workflow.formatProtocol(params.protocol, "alevin")

// Whitelist files for STARsolo and Kallisto
whitelist_folder = "$baseDir/assets/whitelist/"

// Automatically set up proper filepaths to the barcode whitelist files bundled with the pipeline
if (params.protocol.contains("10X") && !params.barcode_whitelist) {
    barcode_filename = "$whitelist_folder/10x_${chemistry}_barcode_whitelist.txt.gz"
    Channel.fromPath(barcode_filename)
    .ifEmpty{ exit 1, "Cannot find ${protocol} barcode whitelist: $barcode_filename" }
    .set{ barcode_whitelist_gzipped }
} else if (params.barcode_whitelist){
    Channel.fromPath(params.barcode_whitelist)
    .ifEmpty{ exit 1, "Cannot find ${protocol} barcode whitelist: $barcode_filename" }
    .set{ ch_barcode_whitelist }
}

////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def salmon_index_options                = modules['salmon_index']
def gffread_txp2gene_options            = modules['gffread_tx2pgene']
def salmon_alevin_options               = modules['salmon_alevin']
def alevin_qc_options                   = modules['alevinqc']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { GFFREAD_TRANSCRIPTOME }       from '../../modules/local/gffread/transcriptome/main'   addParams( options: [:] )
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
    reads

    main:

    ch_software_versions = Channel.empty()

  // Unzip barcodes
    if (params.protocol.contains("10X") && !params.barcode_whitelist) {
        GUNZIP ( barcode_whitelist_gzipped )
        ch_barcode_whitelist = GUNZIP.out.gunzip
    }

    // Preprocessing - Extract transcriptome fasta from genome fasta
    if (!params.transcript_fasta && params.genome_fasta && params.gtf) {
        GFFREAD_TRANSCRIPTOME ( genome_fasta, gtf )
        transcriptome_fasta = GFFREAD_TRANSCRIPTOME.out.transcriptome_extracted
        ch_software_versions = ch_software_versions.mix(GFFREAD_TRANSCRIPTOME.out.version.first().ifEmpty(null))
    }
    
    // Build salmon index
    if (!params.salmon_index) {
        SALMON_INDEX ( genome_fasta, transcriptome_fasta )
        salmon_index_alevin = SALMON_INDEX.out.index
    }

    // Build txp2gene map
    if (!params.txp2gene){
        GFFREAD_TXP2GENE ( gtf )
        ch_txp2gene = GFFREAD_TXP2GENE.out.gtf
        // Only collect version if not already done for gffread
        if (!GFFREAD_TRANSCRIPTOME.out.version) {
            ch_software_versions = ch_software_versions.mix(GFFREAD_TXP2GENE.out.version.first().ifEmpty(null))
        }
    }

    // Perform quantification with salmon alevin
    SALMON_ALEVIN ( 
        reads, 
        salmon_index_alevin, 
        ch_txp2gene, 
        protocol, 
        ch_barcode_whitelist 
    )

    // Collect software versions
    ch_software_versions = ch_software_versions.mix(SALMON_ALEVIN.out.version.first().ifEmpty(null))

    // Run alevinQC
    ALEVINQC ( SALMON_ALEVIN.out.alevin_results )
    
    // Collect software versions
    ch_software_versions = ch_software_versions.mix(ALEVINQC.out.version.first().ifEmpty(null))

    // MultiQC
    if (!params.skip_multiqc) {
        ch_salmon_multiqc   = SALMON_ALEVIN.out.alevin_results

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_salmon_multiqc.collect{it[1]}.ifEmpty([]))
    }

    emit:
    software_versions          = ch_software_versions
    multiqc_files              = ch_multiqc_files
}
