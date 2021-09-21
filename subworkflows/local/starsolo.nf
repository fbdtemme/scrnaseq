////////////////////////////////////////////////////
/* --         STARSolo SUBWORKFLOW             -- */
////////////////////////////////////////////////////

// Whitelist files for STARsolo and Kallisto
def whitelist_folder = "$baseDir/assets/whitelist/"

////////////////////////////////////////////////////
/* --    Define command line options           -- */
////////////////////////////////////////////////////
def modules = params.modules.clone()

def star_genomegenerate_options   = modules['star_genomegenerate']
def star_align_options            = modules['star_align']
def postprocess_options           = modules['postprocess']
postprocess_options.publish_dir   = 'star'
def gunzip_options                = modules['gunzip']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { STAR_ALIGN }            from '../../modules/local/star/alignsolo/main'                 addParams( options: star_align_options )
include { POSTPROCESS }           from '../../modules/local/postprocess/main'                    addParams( options: postprocess_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                from '../../modules/nf-core/modules/gunzip/main'               addParams( options: gunzip_options )
include { STAR_GENOMEGENERATE }   from '../../modules/nf-core/modules/star/genomegenerate/main'  addParams( options: star_genomegenerate_options )
include { UNTAR }                 from '../../modules/nf-core/modules/untar/main.nf'             addParams( options: [:] )

workflow STARSOLO {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    genome_fasta       // channel: /path/to/genome.fasta
    gtf                // channel: /path/to/annotation.gtf
    star_index         // channel: /path/to/start/index
    protocol           // channel: protocol
    barcode_whitelist  // channel: /path/to/barcode_whitelist.txt

    main:
    ch_software_versions = Channel.empty()

    (star_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "star")
   
    // Setup correct barcode whitelist and unzip if necessary
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
        ch_star_index = STAR_GENOMEGENERATE.out.index
    } else {
        ch_star_index = star_index
    }
  
    // Perform mapping and quantification with STARsolo
    STAR_ALIGN ( 
        reads,
        ch_star_index,
        gtf,
        ch_barcode_whitelist,
        star_protocol
    )

    // Collect software versions
    ch_software_versions    = ch_software_versions.mix(STAR_ALIGN.out.version.first().ifEmpty(null))

    // Reformat output and run postprocess module
    // TODO there may be a cleaner way of doing this
    // Meta and result could be set in one command and the matrix, barcodes and features could
    // be mixed into one channel
    ch_star_results_files   = STAR_ALIGN.out.solo_results.map{ it[1] }
    ch_meta                 = STAR_ALIGN.out.solo_results.map{ it[0] }
    ch_matrix               = ch_star_results_files.map{ "${it}/Gene/filtered/matrix.mtx" }
    ch_features             = ch_star_results_files.map{ "${it}/Gene/filtered/features.tsv" }
    ch_barcodes             = ch_star_results_files.map{ "${it}/Gene/filtered/barcodes.tsv" }
    POSTPROCESS ( ch_meta, ch_matrix, ch_features, ch_barcodes )

    // Collect software versions
    ch_software_versions    = ch_software_versions.mix(POSTPROCESS.out.version.first().ifEmpty(null))

    // Collect multiqc files
    ch_multiqc_files        = STAR_ALIGN.out.log_final.collect{ it[1] }

    emit:
    software_versions       = ch_software_versions
    multiqc_files           = ch_multiqc_files
}
