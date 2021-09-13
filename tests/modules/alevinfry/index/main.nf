#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALEVINFRY_INDEX } from '../../../../modules/local/alevinfry/index/main' addParams( options: [:] )

workflow test_alevinfry_index {
    splici_ref  = params.test_data_scrnaseq["reference"]["alevinfry"]["splici_ref"]

    ALEVINFRY_INDEX ( splici_ref )
}
