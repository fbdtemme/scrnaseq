#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP }          from '../../../../modules/nf-core/modules/gunzip/main.nf'       addParams( options: [:] )
include { SALMON_INDEX }    from '../../../../modules/nf-core/modules/salmon/index/main.nf' addParams( options: [:] )
include { SALMON_ALEVIN }   from '../../../../modules/local/salmon/alevin/main.nf'          addParams( options: [:] )

workflow test_salmon_alevin {
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    transcriptome_fasta = file(params.test_data_scrnaseq["reference"]["mouse_transcriptome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    tx2gene             = file(params.test_data_scrnaseq["reference"]["mouse_tx2gene"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "chromium"
    chemistry           = "V2"
    barcode_whitelist   = file(params.test_data_scrnaseq["reference"]["tenx_V2_barcode_whitelist"], checkIfExists: true)

    SALMON_INDEX ( genome_fasta, transcriptome_fasta )
    salmon_index = SALMON_INDEX.out.index

    SALMON_ALEVIN ( 
        fastq, 
        salmon_index,
        tx2gene, 
        protocol, 
        barcode_whitelist
    )
}