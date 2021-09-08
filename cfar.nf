nextflow.enable.dsl=2

process select_datasets {
  module 'R/4.1.0-foss-2020b'
  label 'low_mem'
  errorStrategy 'retry'
  publishDir "$params.output.data/hiv_positive", mode : "copy", pattern : "*_CFAR.tsv"
  publishDir "$params.output.data/normal", mode : "copy", pattern : "HIP*.tsv"
  input:
    path cfar_path
    path norm_path
  output:
    path "*_CFAR.tsv", emit: cfar_data
    path "HIP*.tsv", emit: norm_data
  script:
    """
    Rscript $moduleDir/scripts/cfar/dataSelector.R -c ${cfar_path} -n ${norm_path}
    """
}


