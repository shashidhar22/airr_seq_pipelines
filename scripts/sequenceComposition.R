library(tidyverse)
library(patchwork)
library(ggalluvial)
library(lubridate)
library(ggplot2)
library(viridis)
library(readxl)
library(scales)
library(lubridate)
library(ComplexHeatmap)
library(tidyHeatmap)
library(ggpubr)
devtools::load_all("~/projects/LymphoSeq2")
option_list <- list(
    optparse::make_option(c("-s", "--stable_data"), 
                type="character",  
                help="Path to study data RDA file",
                dest="sdata"),
    optparse::make_option(c("-m", "--meta_data"), 
                type="character",  
                help="Path to meta data RDA file",
                dest="mdata"),
    optparse::make_option(c("-p", "--patient_id"),
                          type="numeric",
                          help="Patient ID",
                          dest="pid"),
    optparse::make_option(c("-r", "--seed"),
                          type="numeric",
                          help="Seed value for random sampling of the repertoire",
                          dest="seed")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                               print_help_and_exit = TRUE)

load(parser$sdata)
load(parser$mdata)

getRarefiedCounts <- function(study_table, patient_table) {
    #Get the smallest repertoire size from the cohort
    rep_size <- study_table %>% 
                dplyr::group_by(repertoire_id) %>% 
                dplyr::summarise(rep_size = sum(duplicate_count))  %>% 
                ungroup()
    #Uncount all junction sequences and randomly sample them to smallest rep size, finally group by unique rearrangement and add count and frequency of rearrangement
    set.seed(parser$seed)
    rarefied_table <- study_table %>% 
            tidyr::uncount(duplicate_count) %>% 
            dplyr::group_by(repertoire_id) %>% 
            dplyr::sample_n(size = min(rep_size$rep_size)) %>% 
            dplyr::group_by_all() %>% 
            dplyr::summarise(duplicate_count = n()) %>% 
            dplyr::ungroup() %>% 
            dplyr::group_by(repertoire_id) %>% 
            dplyr::mutate(duplicate_frequency = duplicate_count/ sum(duplicate_count)) %>% 
            dplyr::ungroup()
    #Get productive sequences
    rarefied_amino_table <- rarefied_table %>% 
                            LymphoSeq2::productiveSeq() %>%
                            left_join(sample_table, by = "repertoire_id")
    return(rarefied_amino_table)
}

#Get rarefied count through random sampling
stable <- stable %>% 
          filter(patientID == parser$pid)
rarefied_amino_table <- getRarefiedCounts(stable, sample_table)

sorensen_matrix <- scoringMatrix(rarefied_amino_table, mode = "PSI")
samp_table <- sample_table %>% filter(patientID == parser$pid)
samp_table <- left_join(samp_table, patient_table, by = "patientID") %>% 
                group_by(patientID) %>% 
                mutate(relTime = (dateART %--% dateTCR) %/% days(1)) %>% 
                ungroup()
row_dates <- samp_table %>% pull(relTime)
names(row_dates) <- samp_table %>%  pull(repertoire_id)
col_fun = colorRamp2(c(0, 50, 100), c("#f7fbff", "#6baed6", "#08306b"))
sorensen_plot <- Heatmap(sorensen_matrix, 
                         col = col_fun, 
                         column_split = samp_table$time_point, 
                         row_split = samp_table$time_point, 
                         row_labels = row_dates[rownames(sorensen_matrix)], 
                         cluster_columns = FALSE, 
                         cluster_rows = FALSE, 
                         show_column_names = FALSE, 
                         cell_fun = function(j, i, x, y, width, height, fill) {
                                    grid.text(sprintf("%.0f", sorensen_matrix[i, j]), x, y, gp = gpar(fontsize = 10))}, 
                         name = "PSI", 
                         column_title = paste("Sample ID ", parser$pid, sep = ":"))
out_file <- paste(parser$pid, "simPlot.pdf", sep = "_")
save_pdf(sorensen_plot, out_file, width = 4.5, height = 3, units = "in")

