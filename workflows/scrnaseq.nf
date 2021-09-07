////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

params.summary_params = [:]

// Create a channel for input read files
if (params.input) { 
    ch_input = file(params.input)
} else { 
    exit 1, 'Input samplesheet file not specified!'
}


// Stage config files
ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs           = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images    = file("$projectDir/docs/images/", checkIfExists: true)


// Don't overwrite global params.modules, create a copy instead and use that within the main script
def modules = params.modules.clone()

def multiqc_options                     = modules['multiqc']
def fastqc_options                      = modules['fastqc'] 

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { GET_SOFTWARE_VERSIONS }    from '../modules/local/get_software_versions'          addParams( options: [publish_files: ['csv':'']]       )
include { INPUT_CHECK }              from '../subworkflows/local/input_check'               addParams( options: [:] )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

include { FASTQC  }                     from '../modules/nf-core/modules/fastqc/main'       addParams( options: fastqc_options)
include { MULTIQC }                     from '../modules/nf-core/modules/multiqc/main'      addParams( options: multiqc_options )

def tools = params.aligner ? params.aligner.split(',').collect{ it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '') } : []

////////////////////////////////////////////////////
/* CONDITIONALY IMPORT LOCAL MODULES/SUBWORKFLOWS */
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
    fastqc_zip = Channel.empty()
    if (!params.skip_fastqc) {
        FASTQC ( ch_fastq )

        ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(fastqc_zip.map { it -> it[1] }.collect())
    } 

    // Dispatch to specified tool

    if ("alevin" in tools) {
        ALEVIN ( ch_fastq )

        ch_software_versions = ch_software_versions.mix(ALEVIN.out.software_versions.collect().ifEmpty([]))
        ch_multiqc_files     = ch_multiqc_files.mix(ALEVIN.out.multiqc_files.collect().ifEmpty([]))  
    }

    if ("alevinfry" in tools) {
        ALEVINFRY( ch_fastq )

        ch_software_versions = ch_software_versions.mix(ALEVINFRY.out.software_versions.collect().ifEmpty([]))
        ch_multiqc_files     = ch_multiqc_files.mix(ALEVINFRY.out.multiqc_files.collect().ifEmpty([]))  
    }

    // Run STARSolo pipeline
    if ("star" in tools) {
        STARSOLO( ch_fastq )

        ch_software_versions = ch_software_versions.mix(STARSOLO.out.software_versions.collect().ifEmpty([]))
        ch_multiqc_files     = ch_multiqc_files.mix(STARSOLO.out.multiqc_files.collect().ifEmpty([]))    
    }

    // Run kallisto bustools pipeline
    if ("kallisto" in tools) {
        KALLISTO_BUSTOOLS( ch_fastq )

        ch_software_versions = ch_software_versions.mix(KALLISTO_BUSTOOLS.out.software_versions.collect().ifEmpty([]))
        ch_multiqc_files     = ch_multiqc_files.mix(KALLISTO_BUSTOOLS.out.multiqc_files.collect().ifEmpty([]))    
    }

  // Run kallisto bustools pipeline
    if ("cellranger" in tools) {
        CELLRANGER( ch_fastq )

        ch_software_versions = ch_software_versions.mix(KALLISTO_BUSTOOLS.out.software_versions.collect().ifEmpty([]))
        ch_multiqc_files     = ch_multiqc_files.mix(CELLRANGER.out.multiqc_files.collect().ifEmpty()([]))    
    }
   
    // Get software versions
    GET_SOFTWARE_VERSIONS ( ch_software_versions.map{it}.collect() )

     // MultiQC
    if (!params.skip_multiqc) {
        workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_config)
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())

        MULTIQC ( ch_multiqc_files.collect() )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    }
}
