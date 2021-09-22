#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules                               = params.modules.clone()
def alevinfry_generatepermitlist_options  = modules['alevinfry_generatepermitlist']

include { UNTAR }                         from '../../../../modules/nf-core/modules/untar/main.nf'           addParams( options: [:] )
include { ALEVINFRY_GENERATEPERMITLIST }  from '../../../../modules/local/alevinfry/generatepermitlist/main' addParams( options: alevinfry_generatepermitlist_options)

workflow test_alevinfry_generatepermitlist {
    alevinfry_result_archive = file(params.test_data_scrnaseq["results"]["alevinfry"], checkIfExists: true)

    UNTAR( alevinfry_result_archive )

    ch_input = UNTAR.out.untar
        .map { path -> [[id:"S10_L001", single_end:false], path] }

    ALEVINFRY_GENERATEPERMITLIST( ch_input, "fw" )
}
