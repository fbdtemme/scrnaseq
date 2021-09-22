#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STARSOLO }    from '../../../../subworkflows/local/starsolo'

workflow test_starsolo
{
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "10XV2"
    star_index          = null
    barcode_whitelist   = null

    STARSOLO (
        fastq,
        genome_fasta,
        gtf,
        star_index,
        protocol,
        barcode_whitelist
    )
}