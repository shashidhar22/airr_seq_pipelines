nextflow.enable.dsl=2

process summarize_inputs {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.rdata/summary", mode : "copy"
  input:
    each tcr_path
  output:
    path "${sample}_summary.rda", emit: summary_files
  script:
    sample = tcr_path.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(LymphoSeq2)

    tcr_path <-  "${tcr_path}"
    tcr_table <- readImmunoSeq(tcr_path)
    tcr_summary <- clonality(tcr_table)
    summary_path <- "${sample}_summary.rda"
    save(tcr_summary, file = summary_path)
    """
}

process study_summary {
  module 'R/4.1.0-foss-2020b'
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.tables/summary", mode : "copy"
  input:
    path summary_inputs
  output:
    path "study_summary.tsv", emit: summarized_study
  script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(LymphoSeq2)

    summary_paths <- list.files(path = "./", 
      pattern = "*_summary.rda", 
      full.names = TRUE)
    
    loadSummary <- function(summary_path) {
      load(summary_path)
      return(tcr_summary)
    }

    tcr_summaries <- summary_paths %>% 
      map(loadSummary) %>%
      bind_rows()
    
    summary_path <- "study_summary.tsv"
    write_tsv(tcr_summaries, summary_path)
    """
}

process gliph_runner {
  echo false
  publishDir "$params.output.tables/Gliph2/", mode : "copy", pattern : "*_cluster.csv"
  publishDir "$params.output.tables/Gliph2/", mode : "copy", pattern : "*_HLA.csv"
  label 'mid_mem'
  module 'R/4.1.0-foss-2020b'
  module 'GLIPH2'

  input: 
  path atable
  path htable
  val study_id
  
  output: 
  path "${study_id}_cluster.csv", emit: gliph_cluster
  path "${study_id}_HLA.csv", emit: gliph_hla
  path "${study_id}_gliph.txt", emit: gliph_t0ext

  script:
  """
  Rscript $moduleDir/scripts/prepGliph2.R -a ${atable} -m ${htable} \
  -s ${study_id} -r $params.meta.gref -l $params.meta.lref \
  -v $params.meta.vref 
  irtools.centos -c ${study_id}_gliph.txt
  """
}

process deeptcr_runner {
  echo false
  publishDir "$params.output.tables/Gliph2/", mode : "copy", pattern : "*_cluster.csv"
  publishDir "$params.output.tables/Gliph2/", mode : "copy", pattern : "*_HLA.csv"
  executor 'slurm'
  errorStrategy 'retry'
  maxRetries 3
  time '4d'
  cpus 36
  memory '720 GB'
  queue 'campus-new'
  conda 'DeepTCR'

  input: 
  path dpath
  val study_id

  output: 
  path "${study_id}_Results", emit: deep_results

  script:
  """
  python $moduleDir/scripts/deepRunner.py -d ${dpath} -s ${study_id}
  """
}
