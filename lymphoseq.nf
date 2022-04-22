nextflow.enable.dsl=2

// Filter datasets that have less than the minimum number of sequences
//
// PROCESS INPUTS:
// - Path to ImmunoSEQ file
//  (Source: data path channel)
// - Minimum number of sequences
//  (Source: min sequences channel)
// PROCESS OUTPUTS:
// - AIRR Seq files that have a minimum number of sequences
//
// EMITTED CHANNELS:
// - input_data :- Channel immune data files
//
// NOTE: 
//
// TODO: 
process filter_immuno_files {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/tables/", mode : "copy", pattern : "*.csv"
  input:
    path immuno_path
    val min_sequences
  output:
    path "filtered/*.tsv", emit: input_data
    path "Sequence_count_summary.csv"
  script:
    """
    mkdir filtered
    Rscript $moduleDir/scripts/filterDatasets.R -m ${min_sequences} > stdout.txt 2> stderr.txt
    """
}


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
// - sample_data :- Channel immune data file
//
// NOTE: 
//
// TODO: 
process read_immuno_files {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  input:
    each immuno_path
  output:
    path "${sample}*.rda", emit: sample_data
  script:
    sample = immuno_path.getSimpleName()
    """
    Rscript $moduleDir/scripts/dataReader.R -d ${immuno_path}  > stdout.txt 2> stderr.txt
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
    val study_id
  output:
    path "study_table.rda", emit: study_table
    path "nucleotide_table.rda", emit: nucleotide_table
    path "amino_table.rda", emit: amino_table
    path "summary_table.rda", emit: summary_table
    path "${study_id}_LS.rda", emit: shiny_object
  script:
    """
    Rscript $moduleDir/scripts/dataMerger.R -s ${study_id} > stdout.txt 2> stderr.txt
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
  errorStrategy 'ignore'
  cache 'true'
  publishDir "$params.output.data/rdata/by_sample/summary", mode : "copy", pattern : "*_avg_summary.rda"
  publishDir "$params.output.data/rdata/by_sample/translation_efficacy", mode : "copy", pattern : "*_translation_efficacy.rda"
  input:
    each path(study_table)
    val study_iterations
    val study_minimum_sequences
  output:
    path "*_avg_summary.rda", emit: summary_tables
    path "*_translation_efficacy.rda", emit: translation_efficacy_tables
  script:
    """
    Rscript $moduleDir/scripts/iterativeSummary.R -r ${study_table} \
      -i ${study_iterations} -c ${study_minimum_sequences} > stdout.txt 2> stderr.txt
    """
}

// Process the input AIRR-Seq tables and summarize the prevalence of public
// TRB sequences from VDJdb.
//
// PROCESS INPUTS:
// - Raw sample RDA files
//  (Source: read_data.out.sample_data)
// - Path to VDJDb slim database file
//  (Source: vdjdb_path)
//
// PROCESS OUTPUTS:
// - TSV file containing count of public sequences and cumulative frequency
//
// EMITTED CHANNELS:
// - public_amino_table : RDA file with annotated amino table
//
// NOTE: 
//
// TODO: 
process get_public_tables {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/tables/by_sample/publics", mode : "copy", 
    pattern : "*_public_summary.tsv"
  input:
    each path(sample_data)
    path vdjdb_path
  output:
    path "*_public_summary.tsv", emit: public_summaries
    path "*_public_annotated.rda", emit: public_amino_table
  script:
    """
    Rscript $moduleDir/scripts/findPublicSeq.R -d ${vdjdb_path} -m lenient \
    > stdout.txt 2> stderr.txt
    """
}

// Run GLIPH2 given amino acid and hla information
//
// PROCESS INPUTS:
// - Sampled AIRR data aggregated by amino acid
//  (Source: gliph_amino_table)
// - HLA information when available
//  (Source: gliph_hla_table)
//
// PROCESS OUTPUTS:
// - TARBALL containing GLIPH results
//
// EMITTED CHANNELS:
// - gliph_results : TARBALL file containing GLIPH results
//
// NOTE: 
//
// TODO: 
process run_gliph {
  module 'R/4.1.0-foss-2020b'
  module 'GLIPH2'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/gliph/${run_name}", mode : "copy", pattern : "*.tar.gz"
  input:
    path amino_table
    path hla_table
    path gref
    path lref 
    path vref 
    val study_id
    val run_name
  output:
    path "${study_id}_gliph.tar.gz", emit: gliph_results
    path "${study_id}_cluster.csv", emit: gliph_table 
    path "${study_id}_HLA.csv", emit: gliph_hla_table
  script:
    """
    Rscript $moduleDir/scripts/prepGliph2.R -a ${amino_table} -m ${hla_table} \
    -r ${gref} -l ${lref} -v ${vref} \
    -s ${study_id} > stdout.txt 2> stderr.txt
    irtools.centos -c ${study_id}_gliph.txt > gstdout.txt 2> gstderr.txt
    tar -czf ${study_id}_gliph.tar.gz ./${study_id}*
    """
}

// Calculate the edit distance between sequences and generate network plot
//
// PROCESS INPUTS:
// - Amino acid table for each sample in the analysis
//  (Source: split_sampled_data.out.amino_tables)
// - Study metadata RDA file
//  (Source: process_metadata.out.meta_rda_path)
// - The minimum number of sequences per sample
//  (Source: min_count)
//
// PROCESS OUTPUTS:
// - RDA file containing the graph data structure
// - PDF file containing the network plot
//
// EMITTED CHANNELS:
// - graph_structure : RDA file containing the graph data structure
// - network_plot : PDF file containing the network plot
//
// NOTE: 
//
// TODO: 
process get_networks {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/rdata", mode : "copy", pattern : "*.rda"
  publishDir "$params.output.data/figures/network", mode : "copy", 
    pattern : "*.pdf"
  input:
    each path(amino_table)
    path meta_path
    val min_count
  output:
    path "${sample_id}_graph.rda", emit: graph_structure
    path "${sample_id}_network.pdf", emit: network_plot
  script:
    sample_id = amino_table.getSimpleName().replaceAll("_amino", "")
    """
    Rscript $moduleDir/scripts/getEditDistance.R -a ${amino_table} \
    -m ${meta_path} -n ${min_count}
    """
}

// Prepare control database
//
// PROCESS INPUTS:
// - AIRR data path for each sample
//  (Source: data_path)
//
// PROCESS OUTPUTS:
// - RDA file containing the minimum field required for the database
//
// EMITTED CHANNELS:
// - database_path : RDA file containing the minimum field required for the database
//
// NOTE: 
//
// TODO: 
process prep_controls {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data", mode : "copy", pattern : "*.rda"
  input:
    each path(data_path)
  output:
    path "${sample_id}.rda", emit: graph_structure
  script:
    sample_id = data_path.getSimpleName()
    """
    Rscript $moduleDir/scripts/prepControls.R -d ${data_path}
    """
}



