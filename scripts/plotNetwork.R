library(tidyverse)     
library(lubridate)
library(readxl)
library(optparse)
library(TCellPack)
library(patchwork)
library(ggforce)
# Accept arguments
option_list <- list(
    optparse::make_option(c("-a", "--amino_data"), 
                          type="character",  
                          help="Path to amino acid RDA file", 
                          dest="adata"),
    optparse::make_option(c("-s", "--study_data"), 
                          type="character",  
                          help="Path to study RDA file", 
                          dest="sdata"),
    optparse::make_option(c("-p", "--patient_id"),
                          type="numeric",
                          help="Patient ID",
                          dest="pid")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                               print_help_and_exit = TRUE)
load(parser$adata)
load(parser$sdata)
specTable <- function(atable) {
  atable <- atable %>% 
            dplyr::mutate(yearsRel = (dateART %--% dateTCR) %/% lubridate::days(1)) 

  repertoire_ids<- atable %>%
                   pull(repertoire_id) %>%
                   unique()
  patient_ids <- atable %>%
                pull(patientID) %>%
                unique()
  year_arts <- atable %>% 
              pull(yearsRel) %>%
              unique()
  gpath <- paste(repertoire_ids, "tsv", sep = ".")
  gtable <- read_tsv(gpath, col_names = c("spec_count", "sgroup", "sequences")) %>%
            separate_rows(sequences, sep = " ") %>%
            left_join(atable, by = c("sequences" = "junction_aa"))
            
  return(gtable)
}

getGini <- function(gtable) {
  sgini <- gtable %>% 
           dplyr::filter(spec_count > 10) %>%   
           arrange(desc(spec_count))  %>% 
           select(sgroup, sequences, duplicate_count) %>%
           distinct() %>%
           group_by(sgroup) %>%
           mutate(spec_freq = duplicate_count/sum(duplicate_count)) %>% 
           summarize(gini_index = ineq::Gini(spec_freq),
                     count = n()) %>%
           ungroup() %>%
           arrange(desc(count)) %>% 
           pull(gini_index) %>%
           mean()
  cgini <- gtable %>%
           select(sequences, duplicate_frequency) %>%
           distinct() %>%
           pull(duplicate_frequency) %>%
           ineq::Gini()
  patient_id <- gtable %>%
                pull(patientID) %>%
                unique()
  repertoire_id <- gtable %>%
                   pull(repertoire_id) %>%
                   unique()
  time_point <- gtable %>%
                pull(timePoint) %>%
                unique()
  date_tct <- gtable %>%
              pull(dateTCR) %>% 
              unique()
  years_rel <- gtable %>% 
               pull(yearsRel) %>% 
               unique()
  gini_table <- tibble(repertoire_id = repertoire_id,
                       patientID = patient_id,
                       timePoint = time_point,
                       cluster_gini = sgini,
                       vertex_gini = cgini,
                       yearsRel = years_rel)
  return(gini_table)
}

plotPack <- function(atable) {
  repertoire_id <- atable %>%
                   dplyr::pull(repertoire_id) %>%
                   unique()
  patient_id <- atable %>%
                dplyr::mutate(patientID = as.character(patientID)) %>%
                dplyr::pull(patientID) %>%
                unique()
  year_art <- atable %>% 
              dplyr::mutate(yearsRel = (dateART %--% dateTCR) %/% lubridate::days(1)) %>%
              dplyr::pull(yearsRel) %>%
              unique()
  atable <- atable %>%
            dplyr::select(junction_aa, duplicate_frequency) %>%
            dplyr::rename(clonotype = junction_aa, 
                          frequency = duplicate_frequency)
  
  gpath <- paste(repertoire_id, "tsv", sep = ".")
  pplot <- TCellPack::PlotTCellPack(gliph = gpath,
                                    clonotype.data = atable,
                                    specificity.color =  if(year_art < 0) "#bcbddc" else "#fdae6b",
                                    clonotype.color = if(year_art < 0) "#7570b3" else "#d95f02",
                                    line.color = if(year_art < 0) "#bcbddc" else "#fdae6b") 
  out_path <- paste(patient_id, paste(year_art, "network.pdf", sep = "_"), sep = "_")
  ggsave(out_path, plot = pplot, width = 3, height = 3, units = "in", device = "pdf", limitsize = FALSE)
}

gplot <- atable %>%
         dplyr::filter(patientID == parser$pid) %>%
         dplyr::group_by(repertoire_id) %>%
         dplyr::group_split() %>%
         purrr::map(plotPack) 

gtable <- stable %>% 
          dplyr::filter(patientID == parser$pid) %>%
          group_by(repertoire_id) %>% 
          group_split() %>% 
          purrr::map(specTable) %>% 
          bind_rows()


top_ten <- gtable %>% 
           group_by(sgroup) %>% 
           summarize(max_count = max(spec_count)) %>% 
           arrange(desc(max_count)) %>% 
           slice_head(n = 10) %>% 
           pull(sgroup)

names_ten <- seq(length(top_ten)) 
names_ten <- paste("SG", names_ten, sep = "")
names(names_ten) <- top_ten


gfiltered <- gtable %>% 
             filter(sgroup %in% top_ten) %>%
             group_by(sgroup, repertoire_id, yearsRel) %>% 
             mutate(spec_freq = duplicate_count/sum(duplicate_count)) %>% 
             summarize(spec_count = first(spec_count),
                       spec_gini = ineq::Gini(spec_freq)) %>%
             arrange(spec_count) %>% 
             mutate(time_point = if_else(yearsRel < 0, "Pre-ART", "Post-ART"),
                    yearsRel = as_factor(yearsRel)) %>% 
             ungroup() %>% 
             mutate(spec_rank = names_ten[sgroup]) %>%
             mutate(spec_rank = as_factor(spec_rank),
                    spec_rank = fct_relevel(spec_rank, c("SG1", "SG2", "SG3", "SG4", "SG5", "SG6", "SG7", "SG8", "SG9", "SG10")))
ggini <- gtable %>% 
         mutate(timePoint = if_else(yearsRel < 0, "Pre-ART", "Post-ART")) %>% 
         group_by(repertoire_id) %>% 
         group_split() %>% 
         purrr::map(getGini) %>% 
         bind_rows() 
top_plot <- gfiltered %>% 
            mutate(spec_rank = fct_rev(spec_rank)) %>% 
            ggplot(aes(x= yearsRel, y = spec_rank)) + 
            ggforce::geom_link2(aes(size = spec_gini, color = spec_gini, group = spec_rank), lineend = 'round', n = 500) +
            ggplot2::theme_classic() +
            ggplot2::xlab("Days relative to ART initiation") + 
            ggplot2::ylab("Specificity group") +
            ggplot2::labs(fill = "Gini index") + 
            ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
top_path <- paste(parser$pid, "topPlot.pdf", sep = "_")
ggsave(top_path, plot = top_plot, width = 4.5, height = 3, units = "in", device = "pdf", limitsize = FALSE)


# Scatter plot colored by groups ("Species")
ggini <- ggini %>% 
         mutate(yearsRel = as_factor(yearsRel))
gini_plot <- ggplot(ggini, aes(x = cluster_gini, y = vertex_gini, label = yearsRel, color = timePoint)) +
             geom_point() +
             geom_text(hjust = 0, nudge_x = 0.05) + 
             ggplot2::theme_classic() +
             ggplot2::xlim(0, 1) + 
             ggplot2::ylim(0, 1) + 
             ggplot2::xlab("Mean Specificity Group\nGini Co-efficient") + 
             ggplot2::ylab("Global Repertoire\nGini Co-efficient") +
             ggplot2::labs(color = "Treatment status") + 
             ggplot2::scale_color_manual(values = c("Pre-ART" = "#7570b3", "Post-ART" = "#d95f02")) + 
             ggplot2::theme(legend.position = "none")

gini_path <- paste(parser$pid, "giniPlot.pdf", sep = "_")
ggsave(gini_path, plot = gini_plot, width = 3, height = 3, units = "in", device = "pdf", limitsize = FALSE)

