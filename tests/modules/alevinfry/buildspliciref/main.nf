#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUILDSPLICIREF } from '../../../../modules/local/alevinfry/buildspliciref/main' addParams( options: [:] )

workflow test_buildspliciref {
    genome_fasta = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf          = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    read_length  = 100
    fastq        =  [[id:"S10_L001", single_end:false], [
                        file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                        file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]

    BUILDSPLICIREF ( genome_fasta, gtf, read_length )
}
