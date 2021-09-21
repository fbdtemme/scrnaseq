#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// include { READ_LENGTH } from '../../../subworkflows/local/read_length.nf' addParams( options: [:] )
include { MEAN_READ_LENGTH } from '../../../modules/local/mean_read_length/main' addParams( options: [:] )


workflow test_mean_read_length {
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

    MEAN_READ_LENGTH ( ch_all_reads )
}
