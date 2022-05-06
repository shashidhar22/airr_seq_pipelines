# AIRR Seq pipelines

A group of Nextflow workflows and processes for the analysis of AIRR-Seq data from the Adaptive ImmunoSEQ platform and the 10X VDJ platform. The processes mirror the utilities from LymphoSeq2. 
AIRR Seq pipelines also contain wokrflows for publication specific code. The workflows are designed to 
run using the module framework on slurm clusters. Future updates will try to make these pipelines more
compatible cloud services like AWS or Azure. This repository contains the following workflows:

1. [Serial Analysis of the T-Cell Receptor β-Chain Repertoire in People Living With HIV Reveals Incomplete Recovery After Long-Term Antiretroviral Therapy](https://doi.org/10.3389/fimmu.2022.879190)


## General requirements

All Nextflow workflows in this repository reuquire at least Nextflow version `20.10.0`. The workflows use Nextflow DSL2 syntax. It is configured for slurm clusters and relies on the modules package manager 
framework for the time being. The following modules are required:
1. Anaconda3
2. Nextflow
3. Singularity
4. R/4.1.0
5. Gliph2

The following R packages are required:
1. [tidyverse](https://www.tidyverse.org/packages/)
2. [lubridate](https://lubridate.tidyverse.org/)
3. [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html)
4. [readxl](https://readxl.tidyverse.org/)
5. [optparse](https://cran.r-project.org/web/packages/optparse/index.html)
6. [LymphoSeq2](https://github.com/shashidhar22/LymphoSeq2)
7. [igraph](https://igraph.org/r/)
8. [tidygraph](https://tidygraph.data-imaginist.com/)
9. [ggraph](https://www.data-imaginist.com/2017/ggraph-introduction-layouts/)
10. [ggalluvial](https://corybrunson.github.io/ggalluvial/)
11. [patchwork](https://patchwork.data-imaginist.com/)
12. [scales](https://scales.r-lib.org/)
13. [DT](https://rstudio.github.io/DT/)
14. [ComplexHeatmap](https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
15. [ggpubr](https://rpkgs.datanovia.com/ggpubr/)
16. [circlize](https://jokergoo.github.io/circlize_book/book/)
17. [viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)

## Serial Analysis of the T-Cell Receptor β-Chain Repertoire in People Living With HIV Reveals Incomplete Recovery After Long-Term Antiretroviral Therapy

The workflow for the analysis of T-cell Receptor β-chain repertoires are defined in `main.nf` under the entry point `cfarExperimentOne`. A template parameter files is provided here [`cfar_template.yaml`](https://github.com/shashidhar22/airr_seq_pipelines/blob/main/params/cfar_template.yaml). The required metadata file is provided in the Github repo under [`metadata/cfar/CFAR_Dean_metadata.xlsx`](https://github.com/shashidhar22/airr_seq_pipelines/blob/main/metadata/cfar/CFAR_Dean_metadata.xlsx). TRB repertoire sequencing data from the 192 peripheral blood samples is publicly available for access on the Adaptive immunoSEQ Analyzer portal (https://clients.adaptivebiotech.com/pub/towlerton-2022-hiv). The control dataset was sampled from the 786 TRB repertoire sequencing datasets from bone-marrow donors available on the Adaptive immunoSEQ Analyzer portal (Adaptive Biotechnologies, Seattle, WA; https://clients.adaptivebiotech.com/pub/Dean-2015-GenomeMed).

Before running the workflow, please edit the `run.sh` to specify the following variable:
```
SCRATCH=/path/to/scratch/directory # Nextflow caches a lot, not specifying a scratch directory will quickly fill up your run directory.

PARAMS=/path/to/params/file.yaml # Please make sure to edit the cfar_template.yaml with the correct paths for all the inputs.
```

You will also need to edit the `nextflow.config` with the correct `queue` name for your slurm configuration.

To run the workflow, run the following command:
```
./run.sh
```


