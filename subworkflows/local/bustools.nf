////////////////////////////////////////////////////
/* --       KALLISTO BUSTOOLS SUBWORKFLOW      -- */
////////////////////////////////////////////////////

def modules = params.modules.clone()

def kallistobustools_ref_options     = modules['kallistobustools_ref']
def kallistobustools_count_options   = modules['kallistobustools_count']
def gffread_kallisto_genemap_options = modules['gffread_kallisto_genemap']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
include { GENE_MAP }                from '../../modules/local/genemap/main'                          addParams( options: [:] )
include { POSTPROCESS }                 from '../../modules/local/postprocess/main'                  addParams( options: [:] )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////
include { GUNZIP }                              from '../../modules/nf-core/modules/gunzip/main'                 addParams( options: [:] )
include { KALLISTOBUSTOOLS_COUNT }              from '../../modules/nf-core/modules/kallistobustools/count/main' addParams( options: kallistobustools_count_options )
include { KALLISTOBUSTOOLS_REF }                from '../../modules/nf-core/modules/kallistobustools/ref/main'   addParams( options: kallistobustools_ref_options )
include { GFFREAD as GFFREAD_KALLISTO_GENEMAP } from '../../modules/nf-core/modules/gffread/main'                addParams( options: gffread_kallisto_genemap_options )


workflow KALLISTO_BUSTOOLS {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    genome_fasta       // channel: /path/to/genome.fasta
    gtf                // channel: /path/to/annotation.gtf
    kallisto_gene_map  // channel: /path/to/genemap.txt
    kallisto_index     // channel: /path/to/kallisto/index
    protocol           // channel: protocol

    main:
    ch_software_versions = Channel.empty()
    
    // Get the protocol parameter
    (kb_protocol, chemistry) = WorkflowScrnaseq.formatProtocol(protocol, "kallisto")

    // TODO make this a parameter
    kb_workflow = "standard"

    // Generate Kallisto Gene Map if not supplied and index is given
    // If index is given, the gene map will be generated in the 'kb ref' step 
    if (!kallisto_gene_map && kallisto_index) {
        GFFREAD_KALLISTO_GENEMAP ( gtf )
        ch_kallisto_gene_map = GFFREAD_KALLISTO_GENEMAP.out.gtf
        ch_software_versions = ch_software_versions.mix(GFFREAD_KALLISTO_GENEMAP.out.version.ifEmpty(null))
    } else if (kallisto_gene_map) {
        ch_kallisto_gene_map = kallisto_gene_map
    }

    // Generate kallisto index
    if (!kallisto_index) { 
        KALLISTOBUSTOOLS_REF ( 
            genome_fasta, gtf, 
            kb_workflow 
        )
        ch_kallisto_gene_map = KALLISTOBUSTOOLS_REF.out.t2g
        ch_kallisto_index    = KALLISTOBUSTOOLS_REF.out.index
    } else {
        ch_kallisto_index = kallisto_index
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

    // Reformat output
    ch_kallisto_results_files = KALLISTOBUSTOOLS_COUNT.out.count.map{it[1]}.ifEmpty([])
    // TODO there may be a cleaner way of doing this
    Channel.fromPath("$ch_kallisto_results_files/counts_filtered/cells_x_genes.mtx").set{ matrix_file }
    Channel.fromPath("$ch_kallisto_results_files/counts_filtered/cells_x_genes.genes.txt").set{ features_file }
    Channel.fromPath("$ch_kallisto_results_files/counts_filtered/cells_x_genes.barcodes.txt").set{ barcodes_file }
    POSTPROCESS ( matrix_file, features_file, barcodes_file, "Kallisto" )

    // Collect software versions
    ch_software_versions = ch_software_versions.mix(KALLISTOBUSTOOLS_COUNT.out.version.ifEmpty(null))

    emit: 
    software_versions    = ch_software_versions
    multiqc_files        = ch_kallisto_results_files
}
