// Subworkflows that checks that the input samplesheet is correct
// and then parses it to obtain channels (reads and meta)

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheetcheck/main' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    
    main:
    SAMPLESHEET_CHECK ( samplesheet )
    SAMPLESHEET_CHECK
            .out
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channels(it) }
            .set { sample_info }

    emit:
    sample_info // channel: [ val(meta), [ reads ] ]
}


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    // Check that both fastq files exist in the samplesheet
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    def array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    return array    
}