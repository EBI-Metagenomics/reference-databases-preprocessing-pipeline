#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { REFERENCE_DATABASES_PREPROCESSING } from './workflows/pipeline.nf'

workflow {
    REFERENCE_DATABASES_PREPROCESSING()
}