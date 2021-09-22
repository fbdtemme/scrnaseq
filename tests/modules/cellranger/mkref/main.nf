#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules                  = params.modules.clone()
def cellranger_mkref_options = modules['cellranger_mkref']

include { CELLRANGER_MKREF } from '../../../../modules/local/cellranger/mkref/main.nf' addParams( options: cellranger_mkref_options )

workflow test_cellranger_mkref {
    genome_fasta = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf          = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)

    CELLRANGER_MKREF ( genome_fasta, gtf )
}
