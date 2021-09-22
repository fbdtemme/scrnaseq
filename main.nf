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


//TODO Parse genome parameters here

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

include { SCRNASEQ } from './workflows/scrnaseq'

workflow NFCORE_SCRNASEQ {
    SCRNASEQ ()
}

workflow {
    NFCORE_SCRNASEQ ()
}


////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
