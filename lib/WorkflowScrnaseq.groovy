/*
 * This file holds several functions specific to the pipeline.
 */

class WorkflowScrnaseq {
    
    static void initialise(params, log)
    {
        genomeExists(params, log)

        if (params.gtf && params.gff) {
            log.error "Must provide either a GTF('--gtf') or a GFF3 file ('--gff'), not both"
            System.exit(1)
        }

        def tools = params.tools.split(',').collect{ it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '') }

        if ('alevin' in tools) {
            
            // Check if gtf or TXP2Gene is provided for Alevin
            if (!(params.gtf || params.gff) && !params.txp2gene) {
                log.error "Must provide either a GTF/GFF3 file ('--gtf/--gff') or transcript to gene mapping ('--txp2gene') to quantify with Alevin"
                System.exit(1)
            }

            // Check if files for index building are given if no index is specified
            if (!params.salmon_index && !params.genome_fasta) {
                log.error "Must provide a genome fasta file ('--genome_fasta') or a transcript fasta ('--transcript_fasta') if no index is given!"
                System.exit(1)
            }
        }

        if ('kallisto' in tools) {
            // Check if files for index building are given if no index is specified
            if (!params.kallisto_index && (!params.genome_fasta || !(params.gtf || params.gff))) {
                log.error "Must provide a genome fasta file ('--genome_fasta') and a gtf file ('--gtf') or gff file ('--gff') if no index is given!"
                System.exit(1)
            }

            if (!(params.gtf || params.gff) && !params.kallisto_gene_map) {
                log.error "Must provide either a GTF/GFF3 file ('--gtf/--gff') or kallisto gene map ('--kallisto_gene_map') to align with kallisto bustools!"
                System.exit(1)
            }
        }

        if ('star' in tools) {
            if (!params.star_index && (!(params.gtf || params.gff) || !params.genome_fasta)) {
                log.error "STAR needs either a GTF + FASTA or a precomputed index supplied."
                System.exit(1)
            }
        }

        if ('alevinfry' in tools) {
            if (!params.alevinfry_index && (!(params.gtf || params.gff) || !params.genome_fasta)) {
                log.error "Alevinfry needs either a GTF + FASTA or a precomputed index supplied."
                System.exit(1)
            }

              // Check if gtf or TXP2Gene is provided for Alevin
            if (params.alevinfry_index  && !params.txp2gene) {
                log.error "Must provide a txp2gene file (--alevinfry_gene_map) when using a precumputed index."
                System.exit(1)
            }
        }
    }

    // Exit pipeline if incorrect --genome key provided
    static void genomeExists(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                      "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                      "  Currently, the available genome keys are:\n" +
                      "  ${params.genomes.keySet().join(", ")}\n" +
                      "==================================================================================="
            System.exit(1)
        }
    }

        /*
     * Get workflow summary for MultiQC
     */
    static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    /*
    * Format the protocol
    * Given the protocol paramter (params.protocol) and the tool (params.tools),
    * this function formats the protocol such that it is fit for the respective
    * subworkflow
    */
    static formatProtocol(protocol, tool) {
        String new_protocol = protocol
        String chemistry = ""
        
        // alevin
        if (tool == "alevin") {
            switch(protocol) {
                case "10XV1":
                    new_protocol = "chromium"
                    chemistry = "V1"
                    break
                case "10XV2":
                    new_protocol = "chromium"
                    chemistry = "V2"
                    break
                case "10XV3":
                    new_protocol = "chromiumV3"
                    chemistry = "V3"
                    break
                case "dropseq":
                    new_protocol = "dropseq"
            }
        }

        // star
        else if (tool == "star") {
            switch(protocol) {
                case "10XV1":
                    new_protocol = "CB_UMI_Simple"
                    chemistry = "V1"
                    break
                case "10XV2":
                    new_protocol = "CB_UMI_Simple"
                    chemistry = "V2"
                    break
                case "10XV3":
                    new_protocol = "CB_UMI_Simple"
                    chemistry = "V3"
                    break
                case "dropseq":
                    new_protocol = "CB_UMI_Simple"
                    break
                case "smartseq":
                    new_protocol = "SmartSeq"
            }
        }

        // kallisto bustools
        else if (tool = "kallisto" ) {
            switch(protocol) {
                case "10XV1":
                    new_protocol = "10XV1"
                    chemistry = "V1"
                    break
                case "10XV2":
                    new_protocol = "10XV2"
                    chemistry = "V2"
                    break
                case "10XV3":
                    new_protocol = "10XV3"
                    chemistry = "V3"
                    break
                case "dropseq":
                    new_protocol = "DROPSEQ"
                    break
                case "smartseq":
                    new_protocol = "SMARTSEQ"
            }
        }
        else {
        exit 1, "Tool not recognized."
        }

        return [new_protocol, chemistry]
    }
}