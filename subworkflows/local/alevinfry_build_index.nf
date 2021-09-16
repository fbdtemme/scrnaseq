
def modules = params.modules.clone()

def alevinfry_index_options                 = modules['alevinfry_index']
def alevinfry_map_options                   = modules['alevinfry_map']
def alevinfry_generate_permitlist_options   = modules['alevinfry_permitlist']
def alevinfry_collate_options               = modules['alevinfry_collate']
def alevinfry_quant_options                 = modules['alevinfry_quant']
def postprocess_options                     = modules['postprocess_transpose']
def gunzip_options                          = modules['gunzip']
def build_splici_ref                        = modules['build_splici_ref']
////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { MEAN_READ_LENGTH }                from '../../modules/local/mean_read_length/main'                addParams( options: [:] )
include { BUILD_SPLICI_REF }                from '../../modules/local/alevinfry/build_splici_ref/main'      addParams( options: build_splici_ref )
include { ALEVINFRY_INDEX }                 from '../../modules/local/alevinfry/index/main'                 addParams( options: alevinfry_index_options )

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
    
    // Build splice reference
    BUILD_SPLICI_REF(
        genome_fasta,
        gtf,
        ch_read_length
    )

    // Build alevin-fry index
    ALEVINFRY_INDEX ( BUILD_SPLICI_REF.out.reference )

    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(ALEVINFRY_INDEX.out.version.ifEmpty(null))


    emit:
    index               = ALEVINFRY_INDEX.out.index
    txp2gene_3col       = BUILD_SPLICI_REF.out.txp2gene_3col
    software_versions   = ch_software_versions
}