#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENE_MAP }              from '../../../modules/local/genemap/main.nf'  addParams( options: [:] )

workflow test_gene_map {
    gtf = file(params.test_data_scrnaseq["reference"]["mouse_gtf"], checkIfExists: true)
    
    GENE_MAP ( gtf )
}