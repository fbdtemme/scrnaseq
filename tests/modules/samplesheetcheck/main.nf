#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMPLESHEETCHECK } from '../../../modules/local/samplesheetcheck/main' addParams( options: [:] )

workflow test_samplesheetcheck {
    samplesheet = file(params.test_data_scrnaseq["samplesheet"], checkIfExists: true)

    SAMPLESHEETCHECK ( samplesheet )
}
