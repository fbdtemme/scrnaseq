#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MEANREADLENGTH } from '../../../modules/local/meanreadlength/main' addParams( options: [:] )

workflow test_meanreadlength {
    fastq        =  [
                        [[id:"S10_L001", single_end:false], [
                            file(params.test_data_scrnaseq["testdata"]["R1"], checkIfExists: true), 
                            file(params.test_data_scrnaseq["testdata"]["R2"], checkIfExists: true)]
                        ],
                        [
                            [id:"kallisto", single_end:false], [
                                file(params.test_data_scrnaseq["testdata"]["kallisto"]["R1"], checkIfExists: true), 
                                file(params.test_data_scrnaseq["testdata"]["kallisto"]["R2"], checkIfExists: true)]
                        ]
                    ]

    reads = Channel.from(fastq)

    reads
        .map { meta, reads -> [ reads[0] ] }
        .collect()
        .map { reads -> [ ["id": "FW"], reads ] }
        .set { ch_all_fw }
    reads
        .map { meta, reads -> [ reads[1] ] }
        .collect()
        .map { reads -> [ ["id": "RV"], reads ] }
        .set { ch_all_rv }

    ch_all_reads = ch_all_fw.mix ( ch_all_rv )

    MEANREADLENGTH ( ch_all_reads )
}
