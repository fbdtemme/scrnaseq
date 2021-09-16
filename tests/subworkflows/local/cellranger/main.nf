#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules          = params.modules.clone()

include { CELLRANGER }              from '../../../../subworkflows/local/cellranger'


workflow test_cellranger
{
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    fastq               = [[id:"test", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "10XV2"

    ch_fastq = Channel.from([fastq])

    CELLRANGER(
        ch_fastq,             
        genome_fasta,      
        gtf,
        protocol
    )
}
