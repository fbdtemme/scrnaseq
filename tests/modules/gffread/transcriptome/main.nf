#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules                        = params.modules.clone()
def gffread_transcriptome_options  = modules['gffread_transcriptome']
gffread_transcriptome_options.remove("publish_files")

include { GFFREAD_TRANSCRIPTOME }  from '../../../../modules/local/gffread/transcriptome/main' addParams( options: gffread_transcriptome_options )

workflow test_gffread_transcriptome {
    genome_fasta = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf          = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)

    GFFREAD_TRANSCRIPTOME ( genome_fasta, gtf )
}
