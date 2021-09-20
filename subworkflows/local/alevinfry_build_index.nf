////////////////////////////////////////////////////
/* --         ALEVIN-FRY INDEX SUBWORKFLOW        -- */
////////////////////////////////////////////////////

def modules = params.modules.clone()

def alevinfry_index_options         = modules['salmon_alevinfry_index']
def build_splici_ref                = modules['build_splici_ref']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { MEAN_READ_LENGTH }        from '../../modules/local/mean_read_length/main'                addParams( options: [:] )
include { BUILD_SPLICI_REF }        from '../../modules/local/alevinfry/build_splici_ref/main'      addParams( options: build_splici_ref )
include { SALMON_ALEVINFRY_INDEX }  from '../../modules/local/salmon/alevinfry/index/main'          addParams( options: alevinfry_index_options )

workflow ALEVINFRY_BUILD_INDEX {
    take:
    reads
    genome_fasta
    gtf
    
    main:
      // Flatten input for getting the mean read length
    reads
        .map{ meta, reads -> [ reads[0] ] }
        .collect()
        .map{ reads -> [ ["id": "FW"], reads ] }
        .set{ ch_all_fw }
    reads
        .map{ meta, reads -> [ reads[1] ] }
        .collect()
        .map{ reads -> [ ["id": "RV"], reads ] }
        .set{ ch_all_rv }

    ch_all_reads = ch_all_fw.mix( ch_all_rv )

    // Get mean read length of fowards and reverse reads and select the second reads
    MEAN_READ_LENGTH ( ch_all_reads )
    ch_read_length = MEAN_READ_LENGTH.out.results
        .toSortedList()
        .map{ it[1][1] }
    
    // Collect software version
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(MEAN_READ_LENGTH.out.version.ifEmpty(null))

    // Build splice reference
    BUILD_SPLICI_REF(
        genome_fasta,
        gtf,
        ch_read_length
    )

    // Collect software version
    ch_software_versions = ch_software_versions.mix(BUILD_SPLICI_REF.out.version.ifEmpty(null))

    // Build alevin-fry index
    SALMON_ALEVINFRY_INDEX ( BUILD_SPLICI_REF.out.reference )

    // Collect software version
    ch_software_versions = ch_software_versions.mix(SALMON_ALEVINFRY_INDEX.out.version.ifEmpty(null))

    emit:
    index                = SALMON_ALEVINFRY_INDEX.out.index
    txp2gene_3col        = BUILD_SPLICI_REF.out.txp2gene_3col
    software_versions    = ch_software_versions
}