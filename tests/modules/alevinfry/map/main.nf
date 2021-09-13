#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALEVINFRY_INDEX } from '../../../../modules/local/alevinfry/index/main' addParams( options: [:] )
include { ALEVINFRY_MAP } from '../../../../modules/local/alevinfry/map/main' addParams( options: [:] )

workflow test_alevinfry_map{
    splici_ref    = params.test_data_scrnaseq["reference"]["alevinfry"]["splici_ref"]
    txp2gene_3col = params.test_data_scrnaseq["reference"]["alevinfry"]["txp2gene_3col"]
    fastq         = [[id:"S10_L001", single_end:false], [
                        file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                        file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    alevin_protocol = "chromium"

    ch_index = ALEVINFRY_INDEX ( splici_ref ).index

    ALEVINFRY_MAP ( fastq, ch_index, txp2gene_3col, alevin_protocol )
}
