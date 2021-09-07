////////////////////////////////////////////////////
/* --         KALLISTO BUSTOOLS WORKFLOW       -- */
////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
}

// Check if GTF is supplied properly
if (params.gtf) {
    Channel
    .fromPath(params.gtf)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    .set { gtf }
}

// Setup FastA channels
if (params.genome_fasta) {
    Channel
    .fromPath(params.genome_fasta)
    .ifEmpty { exit 1, "Fasta file not found: ${params.genome_fasta}" }
    .set { genome_fasta }
} 


// Check if files for index building are given if no index is specified
if (!params.kallisto_index && (!params.genome_fasta || !params.gtf)) {
    exit 1, "Must provide a genome fasta file ('--genome_fasta') and a gtf file ('--gtf') if no index is given!"
}

// Setup channel for kallisto index if specified
if (params.kallisto_index) {
    Channel
    .fromPath(params.kallisto_index)
    .ifEmpty { exit 1, "Kallisto index not found: ${params.kallisto_index}" }
    .set { ch_kallisto_index }
}

// Kallist gene map
// Check if txp2gene file has been provided
if (params.kallisto_gene_map) {
    Channel
    .fromPath(params.kallisto_gene_map)
    .set{ ch_kallisto_gene_map } 
}

if (!params.gtf && !params.kallisto_gene_map) {
    exit 1, "Must provide either a GTF file ('--gtf') or kallisto gene map ('--kallisto_gene_map') to align with kallisto bustools!"
}

// Get the protocol parameter
(protocol, chemistry) = Workflow.formatProtocol(params.protocol, "kallisto")
kb_workflow = "standard"

// Create a channel for input read files
if (params.input) { 
    ch_input = file(params.input)
} else { 
    exit 1, 'Input samplesheet file not specified!'
}

// Check AWS batch settings
// TODO use the Checks.awsBatch() function instead


////////////////////////////////////////////////////
/* --    Define command line options           -- */
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

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
def multiqc_report    = []

workflow KALLISTO_BUSTOOLS {
    take:
    reads

    main:
    ch_software_versions = Channel.empty()

    // Generate Kallisto Gene Map if not supplied and index is given
    // If index is given, the gene map will be generated in the 'kb ref' step 
    if (!params.kallisto_gene_map && params.kallisto_index) {
        GENE_MAP ( gtf )
        ch_kallisto_gene_map = GENE_MAP.out.gene_map
    }

    // Generate kallisto index
    if (!params.kallisto_index) { 
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
        protocol
    )

    ch_multiqc_files = Channel.empty()
    // Collect software versions
    ch_software_versions = ch_software_versions.mix(KALLISTOBUSTOOLS_COUNT.out.version.first().ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(KALLISTOBUSTOOLS_COUNT.out.count.map{it[1]}.ifEmpty([]))

    emit: 
    software_versions        = ch_software_versions
    multiqc_files            = ch_multiqc_files
}