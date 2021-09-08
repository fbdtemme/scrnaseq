////////////////////////////////////////////////////
/* --       KALLISTO BUSTOOLS SUBWORKFLOW      -- */
////////////////////////////////////////////////////

def modules = params.modules.clone()

def kallistobustools_ref_options    = modules['kallistobustools_ref']
def kallistobustools_count_options  = modules['kallistobustools_count']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { INPUT_CHECK }                 from '../../subworkflows/local/input_check'                      addParams( options: [:] )
include { GENE_MAP }                    from '../../modules/local/genemap/main'                          addParams( options: [:] )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                      from '../../modules/nf-core/modules/gunzip/main'                 addParams( options: [:] )
include { KALLISTOBUSTOOLS_COUNT }      from '../../modules/nf-core/modules/kallistobustools/count/main' addParams( options: kallistobustools_count_options )
include { KALLISTOBUSTOOLS_REF }        from '../../modules/nf-core/modules/kallistobustools/ref/main'   addParams( options: kallistobustools_ref_options )


workflow KALLISTO_BUSTOOLS {
    take:
    reads
    genome_fasta
    gtf
    kallisto_gene_map
    kallisto_index
    protocol 

    main:

    ch_software_versions = Channel.empty()
    
    // Get the protocol parameter
    (kb_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "kallisto")
    kb_workflow = "standard"

    // Generate Kallisto Gene Map if not supplied and index is given
    // If index is given, the gene map will be generated in the 'kb ref' step 
    if (!kallisto_gene_map && kallisto_index) {
        GENE_MAP ( gtf )
        ch_kallisto_gene_map = GENE_MAP.out.gene_map
    }

    // Generate kallisto index
    if (!kallisto_index) { 
        KALLISTOBUSTOOLS_REF ( 
            genome_fasta, gtf, 
            kb_workflow 
        )
        ch_kallisto_gene_map = KALLISTOBUSTOOLS_REF.out.t2g
        ch_kallisto_index    = KALLISTOBUSTOOLS_REF.out.index
    }

    // Quantification with kallistobustools count
    KALLISTOBUSTOOLS_COUNT (
        reads,
        ch_kallisto_index,
        ch_kallisto_gene_map,
        [],
        [],
        kb_workflow,
        kb_protocol
    )

    ch_multiqc_files = Channel.empty()
    // Collect software versions
    ch_software_versions = ch_software_versions.mix(KALLISTOBUSTOOLS_COUNT.out.version.first().ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(KALLISTOBUSTOOLS_COUNT.out.count.map{it[1]}.ifEmpty([]))

    emit: 
    software_versions        = ch_software_versions
    multiqc_files            = ch_multiqc_files
}
