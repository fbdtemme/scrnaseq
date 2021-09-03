#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules = params.modules.clone()

def star_genomegenerate_options     = modules['star_genomegenerate']
def star_align_options              = modules['star_align']

include { STAR_GENOMEGENERATE } from '../../../../modules/nf-core/modules/star/genomegenerate/main'  addParams( options: star_genomegenerate_options )
include { STAR_ALIGN }          from '../../../../modules/local/star/alignsolo/main'                 addParams( options: star_align_options )

workflow test_star_alignsolo {
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    transcriptome_fasta = file(params.test_data_scrnaseq["reference"]["mouse_transcriptome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    tx2gene             = file(params.test_data_scrnaseq["reference"]["mouse_tx2gene"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]
                          ]
    protocol            = "CB_UMI_Simple"
    chemistry           = "V2"

    barcode_whitelist   = file(params.test_data_scrnaseq["reference"]["tenx_V2_barcode_whitelist"], checkIfExists: true)

    STAR_GENOMEGENERATE ( 
        genome_fasta, 
        gtf 
    )
    star_index = STAR_GENOMEGENERATE.out.index
  
    STAR_ALIGN ( 
        fastq,
        star_index,
        gtf,
        barcode_whitelist,
        protocol
    )
}