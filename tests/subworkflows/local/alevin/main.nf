#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules                = params.modules.clone()
def salmon_index_options   = modules["salmon_index"]

include { ALEVIN }         from '../../../../subworkflows/local/alevin'
include { SALMON_INDEX }   from '../../../../modules/nf-core/modules/salmon/index/main' addParams( options: salmon_index_options )


workflow test_alevin
{
    genome_fasta        = file(params.test_data_scrnaseq["reference"]["mouse_genome"], checkIfExists: true)
    gtf                 = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    fastq               = [[id:"S10_L001", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true),
                                file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]]
    protocol            = "10XV2"

    transcript_fasta    = null
    txp2gene            = null
    salmon_index        = null

    ALEVIN (
        fastq,
        genome_fasta,
        transcript_fasta,
        gtf,
        txp2gene,
        salmon_index,
        protocol
    )
}