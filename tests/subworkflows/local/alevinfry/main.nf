#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules              = params.modules.clone()
def salmon_index_options = modules["salmon_index"]

include { ALEVINFRY }              from '../../../../subworkflows/local/alevinfry'
include { ALEVINFRY_BUILD_INDEX }  from '../../../../subworkflows/local/alevinfry_build_index'


workflow test_alevinfry
{
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "10XV2"

    transcript_fasta        = null
    txp2gene                = null
    alevinfry_gene_map      = null
    alevinfry_index         = null
    barcode_whitelist       = null
    expected_orientation    = "fw"

    ch_fastq = Channel.from([fastq])

    ALEVINFRY (
        ch_fastq,             
        genome_fasta,      
        gtf,
        alevinfry_gene_map,
        alevinfry_index,
        protocol,
        expected_orientation
    )
}

workflow test_alevinfry_index_and_txp2gene
{
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "10XV2"

    transcript_fasta        = null
    txp2gene                = null
    alevinfry_gene_map      = null
    alevinfry_index         = null
    barcode_whitelist       = null
    expected_orientation    = "fw"

    ch_fastq = Channel.from([fastq])

    ALEVINFRY_BUILD_INDEX ( ch_fastq, genome_fasta, gtf )
    ch_index = ALEVINFRY_BUILD_INDEX.out.index
    ch_txp2gene = ALEVINFRY_BUILD_INDEX.out.txp2gene_3col
    
    ALEVINFRY (
        ch_fastq,             
        genome_fasta,      
        gtf,
        ch_txp2gene,
        ch_index,
        protocol,
        expected_orientation
    )
}