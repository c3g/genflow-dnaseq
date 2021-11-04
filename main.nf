#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/dnaseq
========================================================================================
    Github : https://github.com/nf-core/dnaseq
    Website: https://nf-co.re/dnaseq
    Slack  : https://nfcore.slack.com/channels/dnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.bwa   = WorkflowMain.getGenomeAttribute(params, 'bwa')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { DNASEQ } from './workflows/dnaseq'

//
// WORKFLOW: Run main nf-core/dnaseq analysis pipeline
//
workflow GF_DNASEQ {
    DNASEQ ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    GF_DNASEQ ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
