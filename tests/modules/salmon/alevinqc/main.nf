#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules                   = params.modules.clone()
def salmon_alevinqc_options   = modules['salmon_alevinqc']

include { UNTAR }             from '../../../../modules/nf-core/modules/untar/main'   addParams( options: [:] )
include { SALMON_ALEVINQC }   from '../../../../modules/local/salmon/alevinqc/main'   addParams( options: salmon_alevinqc_options )

workflow test_salmon_alevinqc {
    alevin_result_archive = file(params.test_data_scrnaseq["results"]["alevin"], checkIfExists: true)
    UNTAR( alevin_result_archive )

    ch_input = UNTAR.out.untar
        .map { path -> [[id:"S10_L001", single_end:false], path] }

    SALMON_ALEVINQC( ch_input )
}
