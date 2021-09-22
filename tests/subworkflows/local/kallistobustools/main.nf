#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTO_BUSTOOLS }  from '../../../../subworkflows/local/bustools'

workflow test_kallistobustools
{
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "10XV2"
    kallisto_gene_map   = null
    kallisto_index      = null

    KALLISTO_BUSTOOLS (
        fastq,
        genome_fasta,
        gtf,
        kallisto_gene_map,
        kallisto_index,
        protocol,
    )
}
