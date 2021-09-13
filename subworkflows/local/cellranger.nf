


Set available_references = [ "GRCh38", "mm10" ]
////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow CELLRANGER {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    genome              // channel: genome [ GRCh38 | mm10 ] 
    genome_fasta        // channel: /path/to/genome.fasta
    gtf                 // channel: /path/to/annotation.gtf
    prebuild_reference  // channel: /path/to/reference


    // Fetch reference if possible
    if (prebuild_reference in available_references) {
        CELLRANGER_GETREFERENCES()

        ch_reference = CELLRANGER_GETREFERENCES.out.reference
        ch_reference_version = Channel.empty()
    }  else if (!params.prebuilt_gex_reference && !params.genome) {
        CELLRANGER_MKREF( genome_fasta, gtf, ch_reference_name )
        
        ch_reference = CELLRANGER_MKREF.out.reference
        ch_reference_version = CELLRANGER_MKREF.out.version.first().ifEmpty(null)
    } else {
        ch_reference_version = Channel.empty()
    }

    ch_software_versions = ch_software_versions.mix(ch_reference_version.ifEmpty(null))

    ch_cellranger_count = ch_reads.dump(tag: 'before merge')
                                    .map{ it -> [ it[0].gem, it[0].sample, it[1] ] }
                                    .groupTuple()
                                    .dump(tag: 'gem merge')
                                    .map{ get_meta_tabs(it) }
                                    .dump(tag: 'rearr merge')


    //
    // MODULE: Cellranger count
    //
    CELLRANGER_COUNT(
        ch_cellranger_count,
        ch_reference
    )
    ch_software_versions = ch_software_versions.mix(CELLRANGER_COUNT.out.version.ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )
