nextflow.enable.dsl=2

// Merge the metadata from the CFAR and Dean medicine cohorts and generate 
// sample information, patient information, viral load information, CD4 count
// information, and HLA information.
//
// PROCESS INPUTS:
// - Path to CFAR metdata file
//  (Source: cfar_meta_path)
// - Path to Dean metdata file
//  (Source: dean_meta_path)
//
// PROCESS OUTPUTS:
// - RDA file with all the combined metadata tables
// EMITTED CHANNELS:
// - meta_rda_path :- Channel metadata RDA file
//
// NOTE: 
//
// TODO: 
process process_metadata {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/rdata/", mode : "copy"
  input:
    path meta_path
  output:
    path "Study_meta_frame.rda", emit: meta_rda_path
  script:
    """
    Rscript $moduleDir/scripts/cfar/formatMetadata.R -m ${meta_path} 
    """
}
// Format the HLA information and prepare the amino acid data for external tools 
// like GLIPH2 and DeepTCR.
//
// PROCESS INPUTS:
// - Path to sampled AIRR data aggregated by amino acid sequences
//  (Source: sample_data.out.amino_table)
// - Path to CFAR metdata RDA file
//  (Source: process_metadata.out.meta_rda_path)
// - Minimum number of sequences per repertoire
//  (Source: min_count)
// - Study ID
//  (Source: study_id)
//
// PROCESS OUTPUTS:
// - RDA file with AIRR data formatted for GLIPH runs 
// - RDA file with AIRR data formatted for DeepTCR runs
// - Table with HLA data formatted for GLIPH
// - Table with HLA data formatted for DeepTCR
// EMITTED CHANNELS:
// - gliph_one_amino_table :- Channel with AIRR data formatted for GLIPH (HIV vs Normal)
// - gliph_two_amino_table :- Channel with AIRR data formatted for GLIPH (Pre-ART, Post-ART, and Normal)
// - gliph_hla_table :- Channel with HLA data formatted for GLIPH
// - deeptcr_one :- Channel with AIRR data formatted for DeepTCR (HIV vs Normal)
// - deeptcr_two :- Channel with AIRR data formatted for DeepTCR (Pre-ART, Post-ART, and Normal)
//
// NOTE: 
//
// TODO: 
process prepare_external_runs {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/rdata", mode : "copy", 
    pattern : "${study_id}_gliph_amino_table.rda"
  publishDir "$params.output.data/tables", mode : "copy", 
    pattern : "${study_id}_gliph_hla.tsv"
  publishDir "$params.output.data/deeptcr", mode : "copy",
    pattern : "deeptcr_one"
  publishDir "$params.output.data/deeptcr", mode : "copy",
    pattern : "deeptcr_two"  
  input:
    path cfar_amino_table
    path cfar_meta_path
    val min_count
    val study_id
  output:
    path "${study_id}_gliph_one_amino_table.rds", emit: gliph_one_amino_table
    path "${study_id}_gliph_two_amino_table.rds", emit: gliph_two_amino_table
    path "${study_id}_gliph_hla.tsv", emit: gliph_hla_table
    path "deeptcr_one", emit: deeptcr_one
    path "deeptcr_two", emit: deeptcr_two
  script:
    """
    Rscript $moduleDir/scripts/cfar/prepExtRuns.R -a ${cfar_amino_table} \
      -m ${cfar_meta_path} -n ${min_count} -s ${study_id}
    """
}

// Plot alluvial plots for all samples
//
// PROCESS INPUTS:
// - Study amino acid table
//  (Source: merge_immuno_tables.out.amino_table)
// - Study metadata RDA file
//  (Source: process_metadata.out.meta_rda_path)
//
// PROCESS OUTPUTS:
// - PDF files containing the alluvial plots
//
// EMITTED CHANNELS:
// - alluvial_plot : PDF files containing the alluvial plot
//
// NOTE: 
//
// TODO: 
process get_alluvials {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/figures/alluvial", mode : "copy", pattern : "alluvial.pdf"
  input:
    path amino_table
    path meta_path
    path gliph_table 
    path gliph_hla_table
    path mcpas_table
    path vdjdb_table
  output:
    path "*_alluvial.pdf", emit: alluvial_plots
  script:
    """
    Rscript $moduleDir/scripts/cfar/getAlluvials.R -a ${amino_table} \
    -m ${meta_path} -g ${gliph_table} -h ${gliph_hla_table} \
    -c ${mcpas_table} -v ${vdjdb_table}
    """
}

// Plot plots for all samples
//
// PROCESS INPUTS:
// - Study amino acid table
//  (Source: merge_immuno_tables.out.amino_table)
// - Study metadata RDA file
//  (Source: process_metadata.out.meta_rda_path)
//
// PROCESS OUTPUTS:
// - PDF files containing the alluvial plots
//
// EMITTED CHANNELS:
// - alluvial_plot : PDF files containing the alluvial plot
//
// NOTE: 
//
// TODO: 
process get_plots {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/figures/", mode : "copy", pattern: "*.pdf"
  input:
    path amino_table
    path meta_path
    path summary_tables
    path repertoire_summary
    path translation_efficacy_tables
    path public_tables
    path network_tables
    path gliph_table 
    path gliph_hla_table
    path mcpas_table
    path vdjdb_table
  output:
    path "*.pdf", emit: study_plots
  script:
    """
    Rscript $moduleDir/scripts/cfar/plotCohortFigure.R -a ${amino_table} \
    -m ${meta_path} -r ${repertoire_summary}
    Rscript $moduleDir/scripts/cfar/plotCohortComparison.R -a ${amino_table} \
    -m ${meta_path} -r ${repertoire_summary}
    Rscript $moduleDir/scripts/cfar/plotPublics.R -a ${amino_table} \
    -m ${meta_path} -r ${repertoire_summary}
    Rscript $moduleDir/scripts/cfar/plotNetworks.R -a ${amino_table} \
    -m ${meta_path} -r ${repertoire_summary}
    Rscript $moduleDir/scripts/cfar/plotLengths.R -a ${amino_table} \
    -m ${meta_path} -r ${repertoire_summary}
    Rscript $moduleDir/scripts/cfar/plotSimilarity.R -a ${amino_table} \
    -m ${meta_path} -r ${repertoire_summary}
    """
}


