////////////////////////////////////////////////////
/* --         ALEVIN-FRY INDEX SUBWORKFLOW        -- */
////////////////////////////////////////////////////

def modules = params.modules.clone()

def alevinfry_index_options         = modules['salmon_alevinfry_index']
def alevinfry_buildspliciref        = modules['alevinfry_buildspliciref']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { MEANREADLENGTH }          from '../../modules/local/meanreadlength/main'                  addParams( options: [:] )
include { BUILDSPLICIREF }          from '../../modules/local/alevinfry/buildspliciref/main'        addParams( options: alevinfry_buildspliciref )
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
    MEANREADLENGTH ( ch_all_reads )
    ch_read_length = MEANREADLENGTH.out.results
        .toSortedList()
        .map{ it[1][1] }
    
    // Collect software version
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(MEANREADLENGTH.out.version.ifEmpty(null))

    // Build splice reference
    BUILDSPLICIREF(
        genome_fasta,
        gtf,
        ch_read_length
    )

    // Collect software version
    ch_software_versions = ch_software_versions.mix(BUILDSPLICIREF.out.version.ifEmpty(null))

    // Build alevin-fry index
    SALMON_ALEVINFRY_INDEX ( BUILDSPLICIREF.out.reference )

    // Collect software version
    ch_software_versions = ch_software_versions.mix(SALMON_ALEVINFRY_INDEX.out.version.ifEmpty(null))

    emit:
    index                = SALMON_ALEVINFRY_INDEX.out.index
    txp2gene_3col        = BUILDSPLICIREF.out.txp2gene_3col
    software_versions    = ch_software_versions
}