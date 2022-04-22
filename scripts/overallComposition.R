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
    optparse::make_option(c("-a", "--atable_data"), 
                type="character",  
                help="Path to amino data RDA file",
                dest="adata"),
    optparse::make_option(c("-m", "--meta_data"), 
                type="character",  
                help="Path to meta data RDA file",
                dest="mdata"),
    optparse::make_option(c("-r", "--seed"),
                          type="numeric",
                          help="Seed value for random sampling of the repertoire",
                          dest="seed")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                               print_help_and_exit = TRUE)

load(parser$sdata)
load(parser$mdata)
load(parser$adata)
getRarefiedCounts <- function(study_table, sample_table) {

    #Get the smallest repertoire size from the cohort
    rep_size <- study_table %>% 
                dplyr::group_by(repertoire_id) %>% 
                dplyr::summarise(rep_size = sum(duplicate_count)) 
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
rarefied_atable <- stable %>% 
          getRarefiedCounts(sample_table)

#Get add normal labels and get labels
rarefied_atable <- rarefied_atable %>%  mutate(patientID = as.character(patientID), 
                                               patientID = if_else(is.na(patientID), "Normal", patientID), 
                                               patientID = as_factor(patientID),
                                               time_point = as.character(time_point), 
                                               time_point = if_else(is.na(time_point), "Normal", time_point), 
                                               time_point = as_factor(time_point))

#Get PSI matrix
psi_matrix <- scoringMatrix(rarefied_atable, mode = "PSI")

#Prepare labels
sample_table <- atable %>% 
                select(repertoire_id, patientID, dateTCR, time_point) %>% 
                distinct() %>% 
                left_join(patient_table, by = "patientID") %>% 
                mutate(relTime =  (dateART %--% dateTCR) %/% days(1))
row_dates <- sample_table %>% pull(relTime)
names(row_dates) <- sample_table %>% pull(repertoire_id)
col_fun = colorRamp2(c(0, 50, 100), c("#f7fbff", "#6baed6", "#08306b"))

#Plot rarefied plots
psi_plot <- Heatmap(psi_matrix, 
                         col = col_fun, 
                         column_split = sample_table$patientID, 
                         row_split = sample_table$patientID, 
                         row_labels = row_dates[rownames(psi_matrix)], 
                         cluster_columns = FALSE, 
                         cluster_rows = FALSE, 
                         row_title_rot = 0,
                         row_title_gp = gpar(fontsize = 20),
                         column_title_gp = gpar(fontsize = 24),
                         show_column_names = FALSE,
                         show_row_names = FALSE,
                         name = "PSI", 
                         column_title = "CFAR cohort")
out_file <- "CFAR_simPlot.pdf"
save_pdf(psi_plot, out_file, width = 14, height = 14, units = "in")


#Get add normal labels and get labels
atable <- atable %>%  mutate(patientID = as.character(patientID), 
                             patientID = if_else(is.na(patientID), "Normal", patientID), 
                             patientID = as_factor(patientID),
                             time_point = as.character(time_point), 
                             time_point = if_else(is.na(time_point), "Normal", time_point), 
                             time_point = as_factor(time_point))

#Get PSI matrix
psi_matrix_raw <- scoringMatrix(atable, mode = "PSI")

psi_plot_raw <- Heatmap(psi_matrix_raw, 
                         col = col_fun, 
                         column_split = sample_table$patientID, 
                         row_split = sample_table$patientID, 
                         row_labels = row_dates[rownames(psi_matrix_raw)], 
                         cluster_columns = FALSE, 
                         cluster_rows = FALSE, 
                         row_title_rot = 0,
                         row_title_gp = gpar(fontsize = 20),
                         column_title_gp = gpar(fontsize = 24),
                         show_column_names = FALSE,
                         show_row_names = FALSE,
                         name = "PSI", 
                         column_title = "CFAR cohort")
out_file_raw <- "CFAR_simPlot_raw.pdf"
save_pdf(psi_plot_raw, out_file_raw, width = 14, height = 14, units = "in")

