#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SALMON_ALEVINFRY_INDEX } from '../../../../modules/local/salmon/alevinfry/index/main' addParams( options: [:] )

workflow test_alevinfry_index {
    splici_ref  = params.test_data_scrnaseq["reference"]["alevinfry"]["splici_ref"]

    SALMON_ALEVINFRY_INDEX ( splici_ref )
}
