
nextflow.enable.dsl = 2

include { MEAN_READ_LENGTH } from '../../../modules/local/mean_read_length/main' addParams( options: [:] )

workflow READ_LENGTH {
    take:   
    reads           // [ Channel: [meta, reads] ]

    main:

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

    ch_all_reads
        .mix ( ch_all_fw )
        .set { ch_all_rv }

    MEAN_READ_LENGTH ( ch_all_reads )
}