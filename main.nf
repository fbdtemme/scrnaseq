#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/scrnaseq
========================================================================================
    Github : https://github.com/nf-core/scrnaseq
    Website: https://nf-co.re/scrnaseq
    Slack  : https://nfcore.slack.com/channels/scrnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/


//TODO

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
WorkflowMain.initialise(workflow, params, log)


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    // Run salmon alevin pipeline
    if (params.aligner == "alevin") {
        include { ALEVIN } from './workflows/alevin'
        ALEVIN()
    }

    // Run salmon alevin-fry pipeline
    if (params.aligner == "alevinfry") {
        include { ALEVINFRY } from './workflows/alevinfry'
        ALEVINFRY()
    }

    // Run STARSolo pipeline
    if (params.aligner == "star") {
        include { STARSOLO } from './workflows/starsolo'
        STARSOLO()
    }

    // Run kallisto bustools pipeline
    if (params.aligner == "kallisto") {
        include { KALLISTO_BUSTOOLS } from './workflows/bustools'
        KALLISTO_BUSTOOLS()
    }

    // Run cellranger pipeline
    if (params.aligner == "cellranger") {
        include { CELLRANGER } from './workflows/cellranger'
        CELLRANGER()
    }
    
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
