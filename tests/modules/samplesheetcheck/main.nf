#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMPLESHEET_CHECK } from '../../../modules/local/samplesheetcheck/main.nf' addParams( options: [:] )

workflow test_samplesheet_check {
    samplesheet = file(params.test_data_scrnaseq["samplesheet"], checkIfExists: true)

    SAMPLESHEET_CHECK ( samplesheet )
}
