////////////////////////////////////////////////////
/* --         STARSolo WORKFLOW                -- */
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

// Setup Fasta channels
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

// Check if STAR index is supplied properly
if (params.star_index) {
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}

if (!params.star_index && (!params.gtf || !params.genome_fasta)) {
    exit 1, "STAR needs either a GTF + FASTA or a precomputed index supplied."
}

// Create a channel for input read files
if (params.input) { 
    ch_input = file(params.input)
} else { 
    exit 1, 'Input samplesheet file not specified!' 
}

// Check if txp2gene file has been provided
if (params.txp2gene) {
    Channel
    .fromPath(params.txp2gene)
    .set{ ch_txp2gene } 
}

// Check AWS batch settings
// TODO use the Checks.awsBatch() function instead

// Stage config files
ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs           = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images    = file("$projectDir/docs/images/", checkIfExists: true)

// Get the protocol parameter
(protocol, chemistry) = Workflow.formatProtocol(params.protocol, "star")

// Whitelist files for STARsolo and Kallisto
whitelist_folder = "$baseDir/assets/whitelist/"

// Automatically set up proper filepaths to the barcode whitelist files bundled with the pipeline
if (params.protocol.contains("10X") && !params.barcode_whitelist) {
    barcode_filename = "$whitelist_folder/10x_${chemistry}_barcode_whitelist.txt.gz"
    Channel.fromPath(barcode_filename)
    .ifEmpty{ exit 1, "Cannot find ${protocol} barcode whitelist: $barcode_filename" }
    .set{ barcode_whitelist_gzipped }
} else if (params.barcode_whitelist) {
    Channel.fromPath(params.barcode_whitelist)
    .ifEmpty{ exit 1, "Cannot find ${protocol} barcode whitelist: $barcode_filename" }
    .set{ ch_barcode_whitelist }
}

////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def star_genomegenerate_options     = modules['star_genomegenerate']
def star_align_options              = modules['star_align']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { STAR_ALIGN }                  from '../../modules/local/star/alignsolo/main'                 addParams( options: star_align_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                      from '../../modules/nf-core/modules/gunzip/main'               addParams( options: [:] )
include { STAR_GENOMEGENERATE }         from '../../modules/nf-core/modules/star/genomegenerate/main'  addParams( options: star_genomegenerate_options )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow STARSOLO {
    take:
    reads

    main:
    ch_software_versions = Channel.empty()

    // Unzip barcodes
    if (params.protocol.contains("10X") && !params.barcode_whitelist) {
        GUNZIP ( barcode_whitelist_gzipped )
        ch_barcode_whitelist = GUNZIP.out.gunzip
    }

    // Build STAR index if not supplied
    if (!params.star_index) {
        STAR_GENOMEGENERATE ( 
            genome_fasta, 
            gtf 
        )
        star_index = STAR_GENOMEGENERATE.out.index
    } 
  
    // Perform mapping with STAR
    STAR_ALIGN ( 
        reads,
        star_index,
        gtf,
        ch_barcode_whitelist,
        protocol
    )

    // Collect software versions
    ch_software_versions = ch_software_versions.mix(STAR_ALIGN.out.version.first().ifEmpty(null))

    ch_multiqc_files    = Channel.empty()
    ch_star_multiqc     = STAR_ALIGN.out.log_final
    ch_multiqc_files    = ch_multiqc_files.mix(ch_star_multiqc.collect{it[1]}.ifEmpty([]))


    emit:
    software_versions   = ch_software_versions
    multiqc_files       = ch_multiqc_files
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////