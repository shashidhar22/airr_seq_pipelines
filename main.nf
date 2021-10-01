nextflow.enable.dsl=2
include {read_immuno_files ; merge_immuno_tables ; sample_data ; split_sampled_data ; get_summary_stats } from './lymphoseq'
include { process_metadata } from './cfar'

// Workflow cfarExperimentOne will sample 188 repertoire from the CFAR cohort
// and 188 repertoires from the Dean cohort. For each repertoire we will 
// itratively calculate the summary statistics over 100 iterations and average
// the results. We will also select the top 5000 sequences and track the 
// the prevelance of public sequences, understand the network structure of 
// the repertoire and run DeepTCR and GLIPH2 to identify the antigenic 
// specificities of the sequen
workflow cfarExperimentOne {
  cfar_meta_path = Channel.fromPath(params.meta.cfar_meta)
  norm_meta_path = Channel.fromPath(params.meta.dean_meta)
  cfar_data_path = Channel.fromPath(params.inputs.cfar_data)
  norm_data_path = Channel.fromPath(params.inputs.dean_data)
  normal_count = Channel.value(params.run_one.normal_count)
  min_count = Channel.value(params.meta.min_sequences)
  iterations = Channel.value(params.meta.iterations)
  study_id = Channel.value(params.run_one.study_id)
  main:
    // Process the metadata to get the HLA data in the same format
    // and process the clinical and patient level information
    process_metadata(cfar_meta_path, norm_meta_path)
    // Read all the AIRR sequence data
    read_immuno_files(cfar_data_path.concat(norm_data_path))
    // Merge all the AIRR sequence data into one table
    merge_immuno_tables(read_immuno_files.out.sample_data.flatten().toSortedList())
    // Filter the merged table to only include the samples that have
    // enough sequences to be considered and the desired number of normal
    // samples
    sample_data(merge_immuno_tables.out.study_table, 
      merge_immuno_tables.out.nucleotide_table, 
      merge_immuno_tables.out.amino_table,
      merge_immuno_tables.out.summary_table, 
      normal_count, min_count, study_id)
    // Split the sampled merge data table unto individual sample table for 
    // further downstream analysis
    split_sampled_data(sample_data.out.study_table, 
      sample_data.out.nucleotide_table, sample_data.out.amino_table)
    // Calculate average summary statistics for each repertoire
    get_summary_stats(split_sampled_data.out.study_tables, iterations, 
      min_count)  

}

workflow cfarExperimentTwo {
  cfar_meta_path = Channel.fromPath(params.meta.cfar_meta)
  norm_meta_path = Channel.fromPath(params.meta.dean_meta)
  cfar_data_path = Channel.fromPath(params.inputs.cfar_data)
  norm_data_path = Channel.fromPath(params.inputs.dean_data)
  normal_count = Channel.value(params.run_two.normal_count)
  min_count = Channel.value(params.meta.min_sequences)
  study_id = Channel.value(params.run_two.study_id)
  main:
    process_metadata(cfar_meta_path, norm_meta_path)
    sample_data(normal_count, min_count, study_id, 
      process_metadata.out.meta_rda_path, cfar_data_path.toSortedList(), 
      norm_data_path.toSortedList())
}
