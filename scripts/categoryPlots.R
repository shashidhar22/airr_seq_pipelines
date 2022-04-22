library(tidyverse)
library(ggpubr)
library(lubridate)
library(patchwork)
library(coin)

devtools::load_all("~/projects/LymphoSeq2")
option_list <- list(
    optparse::make_option(c("-m", "--meta_data"), 
                type="character",  
                help="Path to meta data RDA file",
                dest="mdata"),
    optparse::make_option(c("-a", "--amino_data"), 
                type="character",  
                help="Path to amino data RDA file",
                dest="adata"),
    optparse::make_option(c("-s", "--study_data"), 
                type="character",  
                help="Path to study data RDA file",
                dest="sdata"),
    optparse::make_option(c("-i", "--immuno_data"), 
                type="character",  
                help="Path to ImmunoSeq summary data RDA file",
                dest="idata")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                               print_help_and_exit = TRUE)
# Load data
load(parser$mdata)
load(parser$adata)
load(parser$idata)
load(parser$sdata)

# Define functions
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
           arrange(desc(spec_count))  %>% 
           dplyr::filter(spec_count > 1) %>% 
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


# Create rarefied counts table
rep_size <- stable %>% 
            dplyr::group_by(repertoire_id) %>% 
            dplyr::summarise(rep_size = sum(duplicate_count)) 
rstable <- stable %>% 
           tidyr::uncount(duplicate_count) %>% 
           dplyr::group_by(repertoire_id) %>% 
           dplyr::sample_n(size = min(rep_size$rep_size)) %>% 
           dplyr::group_by(repertoire_id, junction, junction_aa, v_call, j_call, d_call, v_family, d_family, j_family) %>% 
           dplyr::summarise(duplicate_count = n(),
                            reading_frame = first(reading_frame),
                            patientID = first(patientID),
                            dateTCR = first(dateTCR),
                            time_point = first(time_point),
                            sex = first(sex),
                            birthYear = first(birthYear),
                            deathYear = first(deathYear),
                            dateART = first(dateART),
                            HIV_status = first(HIV_status),
                            group = first(group),
                            response = first(response),
                            HIVsubtype = first(HIVsubtype),
                            dateMinEvidenceHIV = first(dateMinEvidenceHIV),
                            age_at_art = first(age_at_art)) %>% 
           dplyr::ungroup() %>% 
           dplyr::group_by(repertoire_id) %>% 
           dplyr::mutate(duplicate_frequency = duplicate_count/ sum(duplicate_count)) %>% 
           dplyr::ungroup()
              
ratable <- LymphoSeq2::productiveSeq(rstable)
ritable <- LymphoSeq2::clonality(rstable) 
ritable <- dplyr::left_join(ritable, sample_table, by = "repertoire_id")

supression_table <- vl_table %>% 
                    mutate(vl_suppresion = if_else(result < 200, TRUE, FALSE)) %>%
                    filter(time_point == "Post-ART") %>% 
                    group_by(patientID, vl_suppresion) %>% 
                    summarise(count = n()) %>% 
                    pivot_wider(names_from = vl_suppresion, values_from = count, values_fill = 0) %>% 
                    mutate(supression_percent = `TRUE` * 100/(`FALSE` + `TRUE`),
                           n = `FALSE` + `TRUE`,
                           supression = case_when(supression_percent > 50 & supression_percent <= 100 ~ "High",
                                                         supression_percent <= 50 ~ "Low"))
cc_table <- cc_table %>% 
            mutate(aids = if_else(result < 200, TRUE, FALSE))

ritable <- ritable %>% 
           select(repertoire_id, total_sequences) %>% 
           rename(rarefied_count = total_sequences)

gtable <- atable %>% 
          filter(!is.na(dateART)) %>% 
          group_by(repertoire_id) %>% 
          group_split() %>% 
          purrr::map(specTable) %>% 
          bind_rows()
  
gini_table  <-  gtable %>% 
                mutate(timePoint = if_else(yearsRel < 0, "Pre-ART", "Post-ART")) %>% 
                group_by(repertoire_id) %>% 
                group_split() %>% 
                purrr::map(getGini) %>% 
                bind_rows() 

itable <- itable %>% 
          filter(!is.na(time_point)) %>% 
          left_join(supression_table, by = "patientID") %>% 
          left_join(ritable, by = "repertoire_id")
itable <- itable %>% 
          left_join(gini_table, by = c("repertoire_id", "patientID")) 

cc_table <- cc_table %>% 
            left_join(supression_table, by = "patientID")

cc_2class_plot <- ggboxplot(cc_table, 
                            x = "time_point", 
                            y = "result",
                            color = "time_point", 
                            palette = c("Pre-ART" = "#7570b3", "Post-ART" = "#d95f02"),
                            add = "jitter", 
                            notch = TRUE, 
                            xlab = "Treatment status", 
                            ylab = "CD4 counts (cells/uL)") +
                  stat_compare_means(comparisons = list(c("Pre-ART", "Post-ART")), 
                                     label = "p.adj", 
                                     label.y = 1400) +
                  labs(color = "Treatment status") 

cc_2class_plot <- ggpar(cc_2class_plot, ylim = c(0, 2000)) 

ggplot2::ggsave(filename = "cc2class.pdf", 
                plot = cc_2class_plot, 
                width = 4, 
                height = 4, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)


normalized_count <- ggboxplot(itable, 
                             x = "time_point", 
                             y = "rarefied_count",
                             color = "time_point", 
                             palette = c("Pre-ART" = "#7570b3", "Post-ART" = "#d95f02"),
                             add = "jitter",
                             notch = TRUE,
                             xlab = "Treatment status",
                             ylab = "Rarefied counts") +
                    stat_compare_means(comparisons = list(c("Pre-ART", "Post-ART")),
                                       label = "p.adj",
                                       label.y = 1400) +
                    labs(color = "Treatment status")

normalized_count <- ggpar(normalized_count, ylim = c(0, 1500)) 

ggplot2::ggsave(filename = "nc2class.pdf", 
                plot = normalized_count, 
                width = 4, 
                height = 4, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)


raw_gini <- ggboxplot(itable, 
                      x = "time_point", 
                      y = "gini_coefficient",
                      color = "time_point", 
                      palette = c("Pre-ART" = "#7570b3", "Post-ART" = "#d95f02"),
                      add = "jitter", 
                      notch = TRUE, 
                      xlab = "Treatment status", 
                      ylab = "Global repertoire\nGini coefficient") +
            stat_compare_means(comparisons = list(c("Pre-ART", "Post-ART")), 
                               label = "p.adj", 
                               label.y = 0.9) +
            labs(color = "Treatment status")

raw_gini <- ggpar(raw_gini, ylim = c(0, 1)) 


ggplot2::ggsave(filename = "rg2class.pdf", 
                plot = raw_gini, 
                width = 4, 
                height = 4, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)

cluster_gini <- ggboxplot(itable, 
                          x = "time_point",
                          y = "cluster_gini",
                          color = "time_point",
                          palette = c("Pre-ART" = "#7570b3", "Post-ART" = "#d95f02"),
                          add = "jitter",
                          notch = TRUE,
                          xlab = "Treatment status",
                          ylab = "Mean specificity cluster\nGini coefficient") +
                stat_compare_means(comparisons = list(c("Pre-ART", "Post-ART")), 
                                   label = "p.adj") +
                labs(color = "Treatment status")

ggplot2::ggsave(filename = "cg2class.pdf", 
                plot = cluster_gini, 
                width = 4, 
                height = 4, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)

cc_4class_plot <- cc_table %>% 
                  mutate(time_point = fct_relevel(time_point, c("Pre-ART", "Post-ART"))) %>%
                  ggboxplot(x = "time_point",
                            y = "result", 
                            color = "time_point",
                            palette = "jco",
                            add = "jitter", 
                            shape = "supression", 
                            facet.by = "supression", 
                            order = c("Pre-ART", "Post-ART"),
                            xlab = "Viral load supression rate", 
                            ylab = "CD4 counts (cells/uL)", 
                            label.y = 0.8) +
                  stat_compare_means(comparisons = list(c("Pre-ART", "Post-ART")), 
                                     label = "p.adj") +
                  labs(color = "Viral load supression rate")

cc_4class_plot <- ggpar(cc_4class_plot, ylim = c(0, 2000))

ggplot2::ggsave(filename = "cc4class.pdf", 
                plot = cc_4class_plot, 
                width = 4, 
                height = 4, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)

supression_cluster <- ggboxplot(itable, 
                                x = "time_point", 
                                y = "cluster_gini", 
                                color = "time_point", 
                                palette = "jco",
                                add = "jitter", 
                                shape = "supression", 
                                facet.by = "supression", 
                                order = c("Pre-ART", "Post-ART"),
                                xlab = "Viral load supression rate", 
                                ylab = "Mean specificity cluster\nGini coefficient", 
                                label.y = 0.8) +
                      stat_compare_means(comparisons = list(c("Pre-ART", "Post-ART")), 
                                         label = "p.adj") +
                      labs(color = "Viral load supression rate")

ggplot2::ggsave(filename = "cg4class.pdf", 
                plot = supression_cluster, 
                width = 4, 
                height = 4, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)


supression_vertex <- ggboxplot(itable, 
                               x = "time_point", 
                               y = "vertex_gini", 
                               color = "time_point", 
                               palette = "jco",
                               add = "jitter", 
                               shape = "supression", 
                               facet.by = "supression", 
                               order = c("Pre-ART", "Post-ART"),
                               xlab = "Viral load supression rate", 
                               ylab = "Global repertoire\nGini coefficient", 
                               label.y = 0.9) +
                     stat_compare_means(comparisons = list(c("Pre-ART", "Post-ART")), 
                                        label = "p.adj") +                                                        
                     labs(color = "Viral load supression rate")

supression_vertex <- ggpar(supression_vertex, ylim = c(0, 1)) 

ggplot2::ggsave(filename = "vg4class.pdf", 
                plot = supression_vertex, 
                width = 4, 
                height = 4, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)

supression_normalized <- ggboxplot(itable, 
                                   x = "time_point", 
                                   y = "rarefied_count", 
                                   color = "time_point", 
                                   palette = "jco",
                                   add = "jitter", 
                                   shape = "supression", 
                                   facet.by = "supression", 
                                   order = c("Pre-ART", "Post-ART"),
                                   xlab = "Viral load supression rate", 
                                   ylab = "Rarefied counts") +
                         stat_compare_means(comparisons = list(c("Pre-ART", "Post-ART")), 
                                            label = "p.adj", 
                                            label.y = 1400)+
                         labs(color = "Viral load supression rate")  

supression_normalized <- ggpar(supression_normalized, ylim = c(0, 1500))

ggplot2::ggsave(filename = "nc4class.pdf", 
                plot = supression_normalized, 
                width = 4, 
                height = 4, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)

final_design <- "AAAABBBBCCCCDDDD
                 AAAABBBBCCCCDDDD
                 AAAABBBBCCCCDDDD
                 AAAABBBBCCCCDDDD
                 EEEEFFFFGGGGHHHH
                 EEEEFFFFGGGGHHHH
                 EEEEFFFFGGGGHHHH
                 EEEEFFFFGGGGHHHH"
final_summary <- cc_2class_plot + normalized_count + raw_gini + cluster_gini + cc_4class_plot + supression_normalized + supression_vertex + supression_cluster +
                 plot_layout(design = final_design, guides = "collect") & theme(legend.position = "bottom")

ggplot2::ggsave(filename = "classification_plot.pdf", 
                plot = final_summary, 
                width = 16, 
                height = 8, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)

getIndTest <- function(itable, cc_table) {
  pre_cctable <- cc_table %>% 
                 filter(time_point == "Pre-ART") %>% 
                 select(supression, result) %>% 
                 mutate(supp = if_else(supression == "High", 1, 0)) %>% 
                 select(supp, result)
  post_cctable <- cc_table %>% 
                 filter(time_point == "Post-ART") %>% 
                 select(supression, result) %>% 
                 mutate(supp = if_else(supression == "High", 1, 0)) %>% 
                 select(supp, result)
  pre_itable <- itable %>% 
                filter(time_point == "Pre-ART") %>% 
                select(supression, rarefied_count, cluster_gini, vertex_gini) %>% 
                mutate(supp = if_else(supression == "High", 1, 0)) %>% 
                select(supp, rarefied_count, cluster_gini, vertex_gini)
  post_itable <- itable %>% 
                filter(time_point == "Post-ART") %>% 
                select(supression, rarefied_count, cluster_gini, vertex_gini) %>% 
                mutate(supp = if_else(supression == "High", 1, 0)) %>% 
                select(supp, rarefied_count, cluster_gini, vertex_gini)
  pre_src <- pvalue(independence_test(supp ~ rarefied_count, data = pre_itable))
  pre_svg <- pvalue(independence_test(supp ~ vertex_gini, data = pre_itable))
  pre_scg <- pvalue(independence_test(supp ~ cluster_gini, data = pre_itable))
  pre_scc <- pvalue(independence_test(supp ~ result, data = pre_cctable))
  post_src <- pvalue(independence_test(supp ~ rarefied_count, data = post_itable))
  post_svg <- pvalue(independence_test(supp ~ vertex_gini, data = post_itable))
  post_scg <- pvalue(independence_test(supp ~ cluster_gini, data = post_itable))
  post_scc <- pvalue(independence_test(supp ~ result, data = post_cctable))

  ind_itable <- itable %>% 
                mutate(tmp = if_else(time_point == "Pre-ART", 0, 1)) %>% 
                select(tmp, rarefied_count, cluster_gini, vertex_gini)
  ind_cctable <- cc_table %>% 
                 mutate(tmp = if_else(time_point == "Pre-ART", 0, 1)) %>% 
                 select(tmp, result)
  ind_rc <- pvalue(independence_test(tmp ~ rarefied_count, data = ind_itable))
  ind_vg <- pvalue(independence_test(tmp ~ vertex_gini, data = ind_itable))
  ind_cg <- pvalue(independence_test(tmp ~ cluster_gini, data = ind_itable))
  ind_cc <- pvalue(independence_test(tmp ~ result, data = ind_cctable))

  test_tibble <- tibble(group = c("Pre-ART: High Vs Low", "Pre-ART: High Vs Low", "Pre-ART: High Vs Low", "Pre-ART: High Vs Low", "Post-ART: High Vs Low", "Post-ART: High Vs Low", "Post-ART: High Vs Low", "Post-ART: High Vs Low", "Pre-ART Vs Post-ART", "Pre-ART Vs Post-ART", "Pre-ART Vs Post-ART", "Pre-ART Vs Post-ART"),
                        measure = c("Rarefied count", "Vertex gini", "Cluster_gini", "CD4 count", "Rarefied count", "Vertex gini", "Cluster_gini", "CD4 count", "Rarefied count", "Vertex gini", "Cluster_gini", "CD4 count"),
                        pvalue = c(pre_src, pre_svg, pre_scg, pre_scc, post_src, post_svg, post_scg, post_scc, ind_rc, ind_vg, ind_cg, ind_cc)) 
  
}

test_table <- getIndTest(itable, cc_table)
write_csv(test_table, "Summary_test_table.tsv")