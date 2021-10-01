nextflow.enable.dsl=2

// Read the ImmunoSEQ and generate different views for each sample in the 
// cohort.
//
// PROCESS INPUTS:
// - Path to ImmunoSEQ file
//  (Source: data path channel)
//
// PROCESS OUTPUTS:
// - Tuple with sample name, study table, nucleotide table, amino acid table,
//   and summary table
// EMITTED CHANNELS:
// - immuno_data_path :- Channel immune data file
//
// NOTE: 
//
// TODO: 
process read_immuno_files {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  input:
    path immuno_path
  output:
    path "${sample}*.rda", emit: sample_data
  script:
    sample = immuno_path.getSimpleName()
    """
    Rscript $moduleDir/scripts/dataReader.R -d ${immuno_path} > stdout.txt 2> stderr.txt
    """
}

// Merge LymphoSeq tables from all samples
//
// PROCESS INPUTS:
// - Path to tuple with sample names, study table, nucleotide table, amino acid,
//   and summary table.
//  (Source: read_immuno_files.toSortedList())
//
// PROCESS OUTPUTS:
// - Path to study table
// - Path to nucleotide table
// - Path to amino acid table
// - Path to summary table
//
// EMITTED CHANNELS:
// - study_table :- Channel with merged study table
// - nucleotide_table :- Channel with merged nucleotide table
// - amino_table :- Channel with merged amino acid table
// - summary_table :- Channel with merged summary table
//
// NOTE: 
//
// TODO: 
process merge_immuno_tables {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/rdata/", mode : "copy"
  input:
    path sample_data
  output:
    path "study_table.rda", emit: study_table
    path "nucleotide_table.rda", emit: nucleotide_table
    path "amino_table.rda", emit: amino_table
    path "summary_table.rda", emit: summary_table
  script:
    """
    Rscript $moduleDir/scripts/dataMerger.R > stdout.txt 2> stderr.txt
    """
}

// Sample the metadata such that each sample contains a minimum number of 
// sequences as specified by the parameter. Additionally, the normal datasets
// are sampled to match the desired sample count
//
// PROCESS INPUTS:
// - Path to RDA file with all the merged raw study table
//  (Source: merge_immuno_tables.study_table)
// - Path to RDA file with all the merged nucleotide level aggregation table
//  (Source: merge_immuno_tables.nucleotide_table)
// - Path to RDA file with all the merged amino acid level aggregation table
//  (Source: merge_immuno_tables.amino_table)
// - Path to RDA file with all the merged summary table
//  (Source: merge_immuno_tables.summary_table)
// - Count of desired normal samples
//  (Source: study_normal_count)
// - Desired number of sequences per sample
//  (Source: study_minumum_sequences)
// - Path to metadata RDA file
//  (Source: process_metadata.out.meta_rda_path)
// - Unique study ID for the run
//  (Source: study_id)
//
// PROCESS OUTPUTS:
// - RDA file with a list of samples that will be used dowstream
// EMITTED CHANNELS:
// - study_data_list :- Channel metadata RDA file
//
// NOTE: 
//
// TODO: 
process sample_data {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/rdata/", mode : "copy", pattern : "*.rda"
  publishDir "$params.output.data/tables/", mode : "copy", pattern : "*.tsv"
  input:
    path study_table
    path nucleotide_table
    path amino_table
    path summary_table
    val study_sample_count
    val study_minimum_sequences
    val study_id
  output:
    path "${study_id}_study.rda", emit: study_table 
    path "${study_id}_nucleotide.rda", emit: nucleotide_table
    path "${study_id}_amino.rda", emit: amino_table
    path "${study_id}_data_summary.tsv", emit: study_summary_path
  script:
    """
    Rscript $moduleDir/scripts/dataSelector.R -r ${study_table} \
      -n ${nucleotide_table} -a ${amino_table} -s ${summary_table} \
      -c ${study_sample_count} -t ${study_minimum_sequences} \
      -i ${study_id} > stdout.txt 2> stderr.txt
    """
}

// Split the merged data into individual sample tables that can be processed 
// by the downstream modules.
//
// PROCESS INPUTS:
// - Path to RDA file with all the merged raw study table
//  (Source: merge_immuno_tables.study_table)
// - Path to RDA file with all the merged nucleotide level aggregation table
//  (Source: merge_immuno_tables.nucleotide_table)
// - Path to RDA file with all the merged amino acid level aggregation table
//  (Source: merge_immuno_tables.amino_table)
//
// PROCESS OUTPUTS:
// - TSV files for raw data, nucleotide level aggregation, and amino acid level
//   aggregation for each sample
// EMITTED CHANNELS:
// - study_data_list :- Channel with path to all raw AIRR seq table from 
//   each sample
// - nucleotide_data_list :- Channel with path to all nucleotide level
//   aggregation tables from each sample
// - amino_data_list :- Channel with path to all amino acid level
//   aggregation tables from each sample
//
// NOTE: 
//
// TODO: 
process split_sampled_data {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/tables/by_sample/study", mode : "copy", pattern : "*_study.tsv"
  publishDir "$params.output.data/tables/by_sample/nucleotide", mode : "copy", pattern : "*_nucleotide.tsv"
  publishDir "$params.output.data/tables/by_sample/amino", mode : "copy", pattern : "*_amino.tsv"
  input:
    path study_table
    path nucleotide_table
    path amino_table
  output:
    path "*_study.tsv", emit: study_tables
    path "*_nucleotide.tsv", emit: nucleotide_tables
    path "*_amino.tsv", emit: amino_tables

  script:
    """
    Rscript $moduleDir/scripts/dataSplitter.R -r ${study_table} \
      -n ${nucleotide_table} -a ${amino_table} > stdout.txt 2> stderr.txt
    """
}

// Iterate each sample table and iteratively calculate the summary statistics
// for each sample.
//
// PROCESS INPUTS:
// - Tuple with sampled per sample study data tables
//  (Source: split_sampled_data.out.study_table)
// - Number of iterations to perform
//  (Source: study_iterations)
// - Minimum number of sequences per sample
//  (Source: study_minimum_sequences)
//
// PROCESS OUTPUTS:
// - RDA file with the summary statistics for each sample in the study
//
// EMITTED CHANNELS:
// - summary_table_list :- Channel with path to all raw summary tables
//
// NOTE: 
//
// TODO: 
process get_summary_stats {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/rdata/by_sample/summary", mode : "copy", pattern : "*_avg_summary.rda"
  input:
    path study_table
    val study_iterations
    val study_minimum_sequences
  output:
    path "*_avg_summary.rda", emit: summary_tables

  script:
    """
    Rscript $moduleDir/scripts/iterativeSummary.R -r ${study_table} \
      -i ${study_iterations} -c ${study_minimum_sequences} > stdout.txt 2> stderr.txt
    """
}

