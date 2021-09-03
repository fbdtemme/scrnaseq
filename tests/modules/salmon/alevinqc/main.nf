#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR }      from '../../../../modules/nf-core/modules/untar/main.nf'   addParams( options: [:] )
include { ALEVINQC }   from '../../../../modules/local/salmon/alevinqc/main.nf'   addParams( options: [:] )

workflow test_salmon_alevinqc {
    alevin_result_archive = file(params.test_data_scrnaseq["results"]["alevin"], checkIfExists: true)
    UNTAR( alevin_result_archive )

    ch_input = UNTAR.out.untar
        .map { path -> [[id:"S10_L001", single_end:false], path] }

    ALEVINQC( ch_input )
}