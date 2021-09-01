#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/scrnaseq
========================================================================================
 nf-core/scrnaseq Analysis Pipeline.
 ----------------------------------------------------------------------------------------
    nf-core/scrnaseq:
        An open-source pipeline to compute abundances for single cell datasets. 
----------------------------------------------------------------------------------------
    @Website
    https://nf-co.re/scrnaseq
----------------------------------------------------------------------------------------
    @Documentation
    https://nf-co.re/scrnaseq/usage
----------------------------------------------------------------------------------------
    @Github
    https://github.com/nf-core/scrnaseq
----------------------------------------------------------------------------------------
    @Slack
    https://nfcore.slack.com/channels/scrnaseq
----------------------------------------------------------------------------------------
 
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/scrnaseq --input '*_R{1,2}.fastq.gz' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)


// Check the hostnames against configured profiles
//Checks.hostname(workflow, params, log)



////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    // Run salmon alevin pipeline
    if (params.aligner == "alevin") {
        include { ALEVIN } from '../workflows/alevin'
        ALEVIN()
    }

    // Run salmon alevin-fry pipeline
    if (params.aligner == "alevinfry") {
        include { ALEVINFRY } from '../workflows/alevinfry'
        ALEVINFRY()
    }

    // Run STARSolo pipeline
    if (params.aligner == "star") {
        include { STARSOLO } from '../workflows/starsolo'
        STARSOLO()
    }

    // Run kallisto bustools pipeline
    if (params.aligner == "kallisto") {
        include { KALLISTO_BUSTOOLS } from '../workflows/bustools'
        BUSTOOLS()
    }

    // Run cell ranger pipeline
    if (params.aligner == "cellranger") {
        include { CELLRANGER } from '../workflows/cellranger'
        CELLRANGER()
    }
    
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
