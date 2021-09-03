#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMPLESHEET_CHECK } from '../../../modules/local/samplesheetcheck/main.nf' addParams( options: [:] )

workflow test_samplesheet_check {
    gtf = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)

    GENE_MAP ( gtf )
}