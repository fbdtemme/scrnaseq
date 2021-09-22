#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def alevinfry_generatepermitlist_options    = modules['alevinfry_generatepermitlist']
def alevinfry_collate_options               = modules['alevinfry_collate']
def alevinfry_quant_options                 = modules['alevinfry_quant']

include { UNTAR }                          from '../../../../modules/nf-core/modules/untar/main.nf'            addParams( options: [:] )
include { ALEVINFRY_GENERATEPERMITLIST }   from '../../../../modules/local/alevinfry/generatepermitlist/main'  addParams( options: alevinfry_generatepermitlist_options )
include { ALEVINFRY_COLLATE }              from '../../../../modules/local/alevinfry/collate/main'             addParams( options: alevinfry_collate_options )
include { ALEVINFRY_QUANT  }               from '../../../../modules/local/alevinfry/quant/main'               addParams( options: alevinfry_quant_options )

workflow test_alevinfry_quant {
    alevinfry_result_archive = file(params.test_data_scrnaseq["results"]["alevinfry"], checkIfExists: true)
    txp2gene_3col = params.test_data_scrnaseq["reference"]["alevinfry"]["txp2gene_3col"]

    UNTAR( alevinfry_result_archive )

    ch_input = UNTAR.out.untar
        .map { path -> [[id:"S10_L001", single_end:false], path] }

    // Build permitlist and filter index
    ALEVINFRY_GENERATEPERMITLIST( ch_input, "fw" )
    quant_dir = ALEVINFRY_GENERATEPERMITLIST.out.quant

    ALEVINFRY_COLLATE ( quant_dir, ch_input )

    ALEVINFRY_QUANT ( quant_dir, txp2gene_3col )
}
