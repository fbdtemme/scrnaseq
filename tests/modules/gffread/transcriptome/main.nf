#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFFREAD_TRANSCRIPTOME } from '../../../../modules/local/gffread/transcriptome/main.nf' addParams( options: [:] )

workflow test_gffread_transcriptome {
    genome_fasta = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf          = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)

    GFFREAD_TRANSCRIPTOME ( genome_fasta, gtf )
}
