
////////////////////////////////////////////////////
/* --    VALIDATE INPUTS                       -- */
////////////////////////////////////////////////////

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowScrnaseq.initialise(params, log)

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, 
    params.multiqc_config,
    params.genome_fasta,
    params.transcript_fasta, 
    params.gtf,
    params.gff,
    params.salmon_index,
    params.star_index,
    params.kallisto_index,
    params.txp2gene,
    params.kallisto_gene_map
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// Create a channel for input read files
if (params.input) { 
    ch_input = file(params.input)
}

// Parse tools to make them a list
def tools = params.tools ? params.tools.split(',').collect{ it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '') } : []

////////////////////////////////////////////////////
/* --    CONFIG FILES                          -- */
////////////////////////////////////////////////////

// Stage config files
ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs           = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images    = file("$projectDir/docs/images/", checkIfExists: true)


// Don't overwrite global params.modules, create a copy instead and use that within the main script
def modules                    = params.modules.clone()

def multiqc_options            = modules['multiqc']
def fastqc_options             = modules['fastqc'] 
def gffread_gff3togtf_options  = modules['gffread_gff3togtf']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { GET_SOFTWARE_VERSIONS }        from '../modules/local/get_software_versions'          addParams( options: [publish_files: ['tsv':'']]       )
include { INPUT_CHECK }                  from '../subworkflows/local/input_check'               addParams( options: [:] )
include { GFFREAD as GFFREAD_GFF3TOGTF } from '../modules/nf-core/modules/gffread/main'         addParams( options: gffread_gff3togtf_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

include { FASTQC  }                from '../modules/nf-core/modules/fastqc/main'  addParams( options: fastqc_options)
include { MULTIQC }                from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options )

////////////////////////////////////////////////////
/*    IMPORT LOCAL MODULES/SUBWORKFLOWS           */
////////////////////////////////////////////////////


if ("alevin" in tools) {
    include { ALEVIN }              from '../subworkflows/local/alevin'
}

if ("alevinfry" in tools) {
    include { ALEVINFRY }           from '../subworkflows/local/alevinfry'
}

if ("star" in tools) {
    include { STARSOLO }            from '../subworkflows/local/starsolo'
}

if ("kallisto" in tools) {
    include { KALLISTO_BUSTOOLS }   from '../subworkflows/local/bustools'
}

if ("cellranger" in tools) {
    include { CELLRANGER }          from '../subworkflows/local/cellranger'
}

////////////////////////////////////////////////////
/*    SCRNASEQ WORKFLOW        */
////////////////////////////////////////////////////

workflow SCRNASEQ {

    ch_software_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    // Stage input files
    INPUT_CHECK ( ch_input )
    .map {
        meta, reads -> meta.id = meta.id.split('_')[0..-2].join('_')
        [ meta, reads ]
    }
    .groupTuple(by: [0])
    .map { it -> [ it[0], it[1].flatten() ] }
    .set { ch_fastq }

    // Run FastQC
    if (!params.skip_fastqc) {
        FASTQC ( ch_fastq )

        ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { it -> it[1] }.collect())
    } 

    // Convert GFF to GTF if annotation is given as GFF

    if (params.gff) {
        GFFREAD_GFF3TOGTF ( params.gff )
        ch_gtf = GFFREAD_GFF3TOGTF.out.gtf
    } else {
        ch_gtf = Channel.fromPath(params.gtf)
    }

    // Dispatch to specified tool

    if ("alevin" in tools) {
        ALEVIN ( 
            ch_fastq,
            params.genome_fasta,
            params.transcript_fasta,
            ch_gtf,
            params.txp2gene,
            params.salmon_index,
            params.protocol,
            params.barcode_whitelist
        )

        ch_software_versions = ch_software_versions.mix(ALEVIN.out.software_versions.collect().ifEmpty([]))
        ch_multiqc_files     = ch_multiqc_files.mix(ALEVIN.out.multiqc_files.collect().ifEmpty([]))  
    }

    if ("alevinfry" in tools) {
        //ALEVINFRY( ch_fastq )

        //ch_software_versions = ch_software_versions.mix(ALEVINFRY.out.software_versions.collect().ifEmpty([]))
        //ch_multiqc_files     = ch_multiqc_files.mix(ALEVINFRY.out.multiqc_files.collect().ifEmpty([]))  
    }

    // Run STARSolo pipeline
    if ("star" in tools) {
        STARSOLO( 
            ch_fastq,
            params.genome_fasta,
            ch_gtf,
            params.star_index,
            params.protocol,
            params.barcode_whitelist
        )

        ch_software_versions = ch_software_versions.mix(STARSOLO.out.software_versions.collect().ifEmpty([]))
        ch_multiqc_files     = ch_multiqc_files.mix(STARSOLO.out.multiqc_files.collect().ifEmpty([]))    
    }

    // Run kallisto bustools pipeline
    if ("kallisto" in tools) {
        KALLISTO_BUSTOOLS( 
            ch_fastq,
            params.genome_fasta,
            ch_gtf,
            params.kallisto_gene_map,
            params.kallisto_index,
            params.protocol
        )

        ch_software_versions = ch_software_versions.mix(KALLISTO_BUSTOOLS.out.software_versions.collect().ifEmpty([]))
        ch_multiqc_files     = ch_multiqc_files.mix(KALLISTO_BUSTOOLS.out.multiqc_files.collect().ifEmpty([]))    
    }

    // Run kallisto bustools pipeline
    if ("cellranger" in tools) {
        //CELLRANGER( ch_fastq )

        //ch_software_versions = ch_software_versions.mix(KALLISTO_BUSTOOLS.out.software_versions.collect().ifEmpty([]))
        //ch_multiqc_files     = ch_multiqc_files.mix(CELLRANGER.out.multiqc_files.collect().ifEmpty()([]))    
    }
   
    // Get software versions
    GET_SOFTWARE_VERSIONS ( ch_software_versions.map{it}.collect() )

     // MultiQC
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowScrnaseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files    = ch_multiqc_files.mix(ch_multiqc_config)
        ch_multiqc_files    = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files    = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files    = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())

        MULTIQC ( ch_multiqc_files.collect() )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    }
}
