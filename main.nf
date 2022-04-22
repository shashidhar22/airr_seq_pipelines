nextflow.enable.dsl=2
include {filter_immuno_files ; read_immuno_files ; merge_immuno_tables ; 
  sample_data ; split_sampled_data ; get_summary_stats ; get_public_tables ;
  get_networks ; run_gliph as run_gliph_one ; 
  run_gliph as run_gliph_two ; prep_controls} from './lymphoseq'
include { process_metadata ; prepare_external_runs ; get_alluvials ; get_plots } from './cfar'

// Workflow cfarExperimentOne will sample 188 repertoire from the CFAR cohort
// and 188 repertoires from the Dean cohort. For each repertoire we will 
// itratively calculate the summary statistics over 100 iterations and average
// the results. We will also select the top 5000 sequences and track the 
// the prevelance of public sequences, understand the network structure of 
// the repertoire and run DeepTCR and GLIPH2 to identify the antigenic 
// specificities of the sequen
workflow cfarExperimentOne {
  meta_path = Channel.fromPath(params.meta.meta_path)
  vdjdb_table = Channel.fromPath(params.meta.vdjdb_path)
  mcpas_table = Channel.fromPath(params.meta.mcpas_path)
  cfar_data_path = Channel.fromPath(params.inputs.cfar_data)
  norm_data_path = Channel.fromPath(params.inputs.dean_data)
  normal_count = Channel.value(params.run_one.normal_count)
  min_count = Channel.value(params.meta.min_sequences)
  iterations = Channel.value(params.meta.iterations)
  study_id = Channel.value(params.run_one.study_id)
  gliph_gref = Channel.fromPath(params.gliph.gref)
  gliph_vref = Channel.fromPath(params.gliph.vref)
  gliph_lref = Channel.fromPath(params.gliph.lref)
  main:
    // Process the metadata to get the HLA data in the same format
    // and process the clinical and patient level information
    process_metadata(meta_path)
    // Read all the AIRR sequence data
    read_immuno_files(cfar_data_path.concat(norm_data_path))
    // Merge all the AIRR sequence data into one table
    merge_immuno_tables(read_immuno_files.out.sample_data.flatten().toSortedList(), study_id)
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
    // Find public sequences from the samples
    get_public_tables(read_immuno_files.out.sample_data, vdjdb_table)
    // Prepare amino acid table and HLA information for GLIPH and DeepTCR
    prepare_external_runs(sample_data.out.amino_table, 
      process_metadata.out.meta_rda_path, min_count, study_id)
    // Run GLIPH2 to identify antigenic specificities between HIV and normal
    run_gliph_one(prepare_external_runs.out.gliph_one_amino_table, 
      prepare_external_runs.out.gliph_hla_table, gliph_gref, gliph_vref, 
      gliph_lref, study_id, "gliph_one")
    // Run GLIPH2 to identify antigenic specificities between Pre-ART, Post-ART
    // and HIV
    run_gliph_two(prepare_external_runs.out.gliph_two_amino_table, 
      prepare_external_runs.out.gliph_hla_table, gliph_gref, gliph_vref, 
      gliph_lref, study_id, "gliph_two")
    // Generate network plots for each sample
    get_networks(split_sampled_data.out.amino_tables, 
      process_metadata.out.meta_rda_path, min_count)
    // Generate alluvial plots for each sample
    get_alluvials(sample_data.out.amino_table, 
      process_metadata.out.meta_rda_path, 
      run_gliph_one.out.gliph_table,
      run_gliph_one.out.gliph_hla_table,
      mcpas_table,
      vdjdb_table)
    get_plots(sample_data.out.amino_table, 
      process_metadata.out.meta_rda_path, 
      get_summary_stats.out.summary_tables.toSortedList(),
      sample_data.out.study_summary_path,
      get_summary_stats.out.translation_efficacy_tables.toSortedList(),
      get_public_tables.out.public_summaries.toSortedList(),
      get_networks.out.graph_structure.toSortedList(),
      run_gliph_one.out.gliph_table,
      run_gliph_one.out.gliph_hla_table,
      mcpas_table,
      vdjdb_table)
    
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


workflow uFiftyFour {
  uff_data_path = Channel.fromPath(params.inputs.uff_data)
  min_count = Channel.value(params.meta.min_sequences)
  iterations = Channel.value(params.meta.iterations)
  study_id = Channel.value(params.run_one.study_id)
  main:
    // Parse all input AIRR data and filter samples that do not meet the 
    // minimum sequence requirement
    filter_immuno_files(uff_data_path.toSortedList(), min_count)
    // Read all the AIRR sequence data
    read_immuno_files(filter_immuno_files.out.input_data)
    // Merge all the AIRR sequence data into one table
    merge_immuno_tables(read_immuno_files.out.sample_data.flatten().toSortedList(), study_id)
    // Split the sampled merge data table unto individual sample table for 
    // further downstream analysis
    split_sampled_data(merge_immuno_tables.out.study_table, 
      merge_immuno_tables.out.nucleotide_table, merge_immuno_tables.out.amino_table)
    // Calculate average summary statistics for each repertoire
    get_summary_stats(split_sampled_data.out.study_tables, iterations, 
      min_count)  
    
}

workflow kaposiSarcoma {
  ks_data_path = Channel.fromPath(params.inputs.ks_data)
  min_count = Channel.value(params.meta.min_sequences)
  iterations = Channel.value(params.meta.iterations)
  study_id = Channel.value(params.run_one.study_id)
  main:
    // Parse all input AIRR data and filter samples that do not meet the 
    // minimum sequence requirement
    filter_immuno_files(ks_data_path.toSortedList(), min_count)
    // Read all the AIRR sequence data
    read_immuno_files(filter_immuno_files.out.input_data)
    // Merge all the AIRR sequence data into one table
    merge_immuno_tables(read_immuno_files.out.sample_data.flatten().toSortedList(), study_id)
    // Split the sampled merge data table unto individual sample table for 
    // further downstream analysis
    split_sampled_data(merge_immuno_tables.out.study_table, 
      merge_immuno_tables.out.nucleotide_table, merge_immuno_tables.out.amino_table)
    // Calculate average summary statistics for each repertoire
    get_summary_stats(split_sampled_data.out.study_tables, iterations, 
      min_count)  
    
}


workflow lungVA {
  lung_data_path = Channel.fromPath(params.inputs.lung_data)
  min_count = Channel.value(params.meta.min_sequences)
  iterations = Channel.value(params.meta.iterations)
  study_id = Channel.value(params.run_one.study_id)
  main:
    // Read all the AIRR sequence data
    read_immuno_files(lung_data_path, min_count)
    // Merge all the AIRR sequence data into one table
    merge_immuno_tables(read_immuno_files.out.sample_data.flatten().toSortedList(), study_id)
    // Split the sampled merge data table unto individual sample table for 
    // further downstream analysis
    split_sampled_data(merge_immuno_tables.out.study_table, 
      merge_immuno_tables.out.nucleotide_table, merge_immuno_tables.out.amino_table)
    // Calculate average summary statistics for each repertoire
    get_summary_stats(split_sampled_data.out.study_tables, iterations, 
      min_count)  
    
}

workflow emersonControl {
  emerson_data_path = Channel.fromPath(params.inputs.emerson_data)
  study_id = Channel.value(params.run_one.study_id)
  main:
    // Prepare all the AIRR dataset to create a reference DB
    prep_controls(emerson_data_path)
}
