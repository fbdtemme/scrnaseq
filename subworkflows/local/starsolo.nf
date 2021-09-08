////////////////////////////////////////////////////
/* --         STARSolo SUBWORKFLOW             -- */
////////////////////////////////////////////////////

// Whitelist files for STARsolo and Kallisto
def whitelist_folder = "$baseDir/assets/whitelist/"

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
    genome_fasta
    gtf
    star_index
    protocol
    barcode_whitelist

    main:
    ch_software_versions = Channel.empty()

    (star_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "star")
   
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

    // Build STAR index if not supplied
    if (!star_index) {
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
        star_protocol
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
