#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUILD_SPLICI_REF } from '../../../../modules/local/alevinfry/build_splici_ref/main' addParams( options: [:] )

workflow test_build_splici_ref {
    genome_fasta = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf          = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    read_length  = 100
    fastq        =  [[id:"S10_L001", single_end:false], [
                        file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                        file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]

    BUILD_SPLICI_REF ( genome_fasta, gtf, read_length )
}
