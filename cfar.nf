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
    path cfar_meta_path
    path norm_meta_path
  output:
    path "Study_meta_frame.rda", emit: meta_rda_path
  script:
    """
    Rscript $moduleDir/scripts/cfar/formatMetadata.R -c ${cfar_meta_path} \
      -n ${norm_meta_path} 
    """
}

