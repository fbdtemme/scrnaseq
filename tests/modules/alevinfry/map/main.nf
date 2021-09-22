#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules                         = params.modules.clone()
def alevinfry_index_options         = modules['alevinfry_index']
def salmon_alevin_options           = modules['salmon_alevin']
salmon_alevin_options.publish_dir   = 'alevinfry'
salmon_alevin_options.args          += ' --sketch'

include { SALMON_ALEVINFRY_INDEX }  from '../../../../modules/local/salmon/alevinfry/index/main' addParams( options: alevinfry_index_options )
include { SALMON_ALEVIN }           from '../../../../modules/local/salmon/alevin/main'          addParams( options: salmon_alevin_options )

workflow test_alevinfry_map{
    splici_ref    = params.test_data_scrnaseq["reference"]["alevinfry"]["splici_ref"]
    txp2gene_3col = params.test_data_scrnaseq["reference"]["alevinfry"]["txp2gene_3col"]
    fastq         = [[id:"S10_L001", single_end:false], [
                        file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                        file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    alevin_protocol = "chromium"

    ch_index = SALMON_ALEVINFRY_INDEX ( splici_ref ).index

    SALMON_ALEVIN ( fastq, ch_index, txp2gene_3col, alevin_protocol, "IU" )
}
