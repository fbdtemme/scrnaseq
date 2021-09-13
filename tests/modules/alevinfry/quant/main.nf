#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR }                           from '../../../../modules/nf-core/modules/untar/main.nf'            addParams( options: [:] )
include { ALEVINFRY_GENERATE_PERMITLIST }   from '../../../../modules/local/alevinfry/generate_permitlist/main' addParams( options: [:] )
include { ALEVINFRY_COLLATE }               from '../../../../modules/local/alevinfry/collate/main'             addParams( options: [:] )
include { ALEVINFRY_QUANT  }               from '../../../../modules/local/alevinfry/quant/main'             addParams( options: [:] )

workflow test_alevinfry_quant {
    alevinfry_result_archive = file(params.test_data_scrnaseq["results"]["alevinfry"], checkIfExists: true)
    txp2gene_3col = params.test_data_scrnaseq["reference"]["alevinfry"]["txp2gene_3col"]

    UNTAR( alevinfry_result_archive )

    ch_input = UNTAR.out.untar
        .map { path -> [[id:"S10_L001", single_end:false], path] }

 // Build permitlist and filter index
    ALEVINFRY_GENERATE_PERMITLIST( ch_input, "fw" )
    quant_dir = ALEVINFRY_GENERATE_PERMITLIST.out.quant

    ALEVINFRY_COLLATE ( quant_dir, ch_input )

    ALEVINFRY_QUANT ( quant_dir, txp2gene_3col )
}
