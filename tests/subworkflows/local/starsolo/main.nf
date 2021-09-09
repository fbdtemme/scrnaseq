#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules                      = params.modules.clone()

def star_genomegenerate_options  = modules['star_genomegenerate']

include { STARSOLO }                   from '../../../../subworkflows/local/starsolo'
include { STAR_GENOMEGENERATE }        from '../../../../modules/nf-core/modules/star/genomegenerate/main'  addParams( options: star_genomegenerate_options )


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


workflow test_starsolo_index
{
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "10XV2"
    star_index          = STAR_GENOMEGENERATE ( genome_fasta, gtf ).index
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


workflow test_starsolo_whitelist
{
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "10XV2"
    star_index          = null
    barcode_whitelist   = file(params.test_data_scrnaseq["reference"]["tenx_V2_barcode_whitelist"], checkIfExists: true)

    STARSOLO (
        fastq,
        genome_fasta,
        gtf,
        star_index,
        protocol,
        barcode_whitelist
    )
}