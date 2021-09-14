#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLRANGER_MKREF } from '../../../../modules/local/cellranger/mkref/main.nf' addParams( options: [:] )
include { CELLRANGER_COUNT } from '../../../../modules/local/cellranger/count/main.nf' addParams( options: [:] )

workflow test_cellranger_count {
    fastq         = [[id:"test", single_end:false], [
                    file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                    file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]

    genome_fasta = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf          = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)

    CELLRANGER_MKREF ( genome_fasta, gtf )

    CELLRANGER_COUNT ( fastq, CELLRANGER_MKREF.out.reference )
}
