#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules              = params.modules.clone()

def gffread_kallisto_genemap_options = modules['gffread_kallisto_genemap']
def kallistobustools_ref_options     = modules['kallistobustools_ref']

include { KALLISTO_BUSTOOLS }                   from '../../../../subworkflows/local/bustools'
include { GFFREAD as GFFREAD_KALLISTO_GENEMAP } from '../../../../modules/nf-core/modules/gffread/main'              addParams(
        options: gffread_kallisto_genemap_options )
include { KALLISTOBUSTOOLS_REF }                from '../../../../modules/nf-core/modules/kallistobustools/ref/main' addParams( 
        options: kallistobustools_ref_options )



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

workflow test_kallistobustools_genemap
{
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "10XV2"
    kallisto_index      = null
    kallisto_gene_map   = GFFREAD_KALLISTO_GENEMAP ( gtf ).gtf

    KALLISTO_BUSTOOLS (
        fastq,
        genome_fasta,
        gtf,
        kallisto_gene_map,
        kallisto_index,
        protocol,
    )
}

workflow test_kallistobustools_index
{
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "10XV2"
    kallisto_gene_map   = null
    kb_workflow         = "standard"
       
    kallisto_index = KALLISTOBUSTOOLS_REF ( genome_fasta, gtf, kb_workflow ).index

    KALLISTO_BUSTOOLS (
        fastq,
        genome_fasta,
        gtf,
        kallisto_gene_map,
        kallisto_index,
        protocol,
    )
}