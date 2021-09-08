nextflow.enable.dsl=2
include { summarize_inputs; study_summary; gliph_runner } from './lymphoseq'
include { select_datasets } from './cfar'
workflow summarizeData {
  tcr_path = Channel.fromPath(params.tcr_file + '/*.tsv')
  main:
    summarize_inputs(tcr_path)
    study_summary(summarize_inputs.out.summary_files.collect{ "$it"})
}

workflow cfarReanalysis {
  cfar_path = Channel.fromPath(params.data.cfar_data_path)
  norm_path = Channel.fromPath(params.data.norm_data_path)
  main:
    select_datasets(cfar_path, norm_path)
}

workflow gliphRunner {
  airr_path = Channel.fromPath(params.data.gliph_data_path)
  hla_path = Channel.fromPath(params.data.hla_data_path)
  study_id = Channel.from(params.meta.study_id)
  main:
    gliph_runner(airr_path, hla_path, study_id)
}

workflow runClassificationTools {
  gliph_path = Channel.fromPath(params.data.gliph_data_path)
  dtcr_twoclass_path = Channel.fromPath(params.data.dtcr_twoclass_data_path)
  dtcr_twentyclass_path = Channel.fromPath(params.data.dtcr_twentyclass_data_path)
  hla_path = Channel.fromPath(params.data.hla_data_path)
  study_id = Channel.from(params.meta.study_id)
  dtcr_mode = Channel.from(params.meta.dtcr_mode)

  main:
    gliph_runner(gliph_path, hla_path, study_id)
    

}
