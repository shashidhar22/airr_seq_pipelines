library(tidyverse)
library(lubridate)
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

# Transform data

## Update infection table and create an infection dictionary to standardize the nomenclature
intable <- infection_table %>%
           tidyr::pivot_longer(-patientID, 
                        names_to = "infection", 
                        values_to = "dateInfection")
inpal <- c("EBV" = "EBV",
           "HBV" = "HBV",
           "HCV" = "HCV",
           "HSV1" = "HSV1",
           "HSV2" = "HSV2",
           "Gonnorrhea" = "Gonnorrhea",
           "Syphillis" = "Syphillis",
           "Kaposi_sarcoma" = "Kaposi's Sarcoma",
           "Anal_cancer" = "Anal Cancer", 
           "Colorectal_cancer"  = "Colorectal Cancer",
           "Cervical_cancer" =  "Cervical Cancer",
           "Skin_cancer_squamous_cell_carcinoma" = "Squamous Cell Carcinoma",
           "Esophagus_cancer" = "Esophagus Cancer")

## Transform CD4 counts to be relative to time point nearest to date of treament
    min_cd4 <- cc_table %>% 
            group_by(patientID) %>% 
            dplyr::filter(dateCollection %within% interval(dateART - months(6), dateART)) %>%
            dplyr::arrange(patientID, desc(dateCollection)) %>% 
            dplyr::summarize(result = first(result)) %>%
            dplyr::ungroup()

min_vec <- min_cd4 %>% 
           dplyr::pull(result) 
names(min_vec) <- min_cd4 %>% 
                  dplyr::pull(patientID)

cc_table <- cc_table %>% 
            dplyr::mutate(resultsRel = result - min_vec[as.character(patientID)], 
                          patientID = forcats::fct_rev(patientID),
                          yearsRel = (dateART %--% dateCollection) %/% lubridate::days(1),
                          time_point = dplyr::if_else(dateCollection <= dateART, "Pre-ART", "Post-ART"))

## Calculate days relative ART treatment for ImmunoSeq and Viral load table
vl_table <- vl_table %>%
            dplyr::mutate(patientID = forcats::fct_rev(patientID),
                          yearsRel = (dateART %--% dateCollection) %/% lubridate::days(1),
                          time_point = dplyr::if_else(dateCollection <= dateART, "Pre-ART", "Post-ART"))
itable <- itable %>%
          dplyr::filter(!is.na(dateTCR)) %>%
          dplyr::mutate(patientID = forcats::fct_rev(patientID),
                        yearsRel = (dateART %--% dateTCR) %/% lubridate::days(1), 
                        time_point = dplyr::if_else(dateTCR <= dateART, "Pre-ART", "Post-ART"))

#Calculate ranking

rank_table <- vl_table %>% 
              filter(time_point == "Post-ART") %>%
               mutate(vl_supression = if_else(result < 200, TRUE, FALSE)) %>% 
              group_by(patientID, vl_supression) %>% 
              summarise(count = n()) %>% 
              pivot_wider(names_from = vl_supression, values_from = count, values_fill = 0) %>% 
              mutate(supression_rate = `TRUE`/(`FALSE` + `TRUE`)) %>% 
              select(supression_rate)
# Generate min and max date range
date_range <- unique(c(vl_table$yearsRel, 
                       cc_table$yearsRel, 
                       itable$yearsRel))
min_date <- min(date_range, na.rm=TRUE) - 200
max_date <- max(date_range, na.rm=TRUE) + 200
date_lim <- c(min_date, max_date)      

# Plot sex data
plot_gender <- patient_table %>%
               left_join(rank_table, by = "patientID") %>%
               dplyr::mutate(patientID = forcats::fct_reorder(patientID, supression_rate),
                             gender =  1) %>%
               ggplot2::ggplot(aes(y= patientID, x = gender, fill = sex)) + 
               ggplot2::geom_tile(color = "white") + 
               ggplot2::coord_equal() +
               ggplot2::theme_classic(base_size = 75) +
               ggplot2::xlab("Sex") + 
               ggplot2::ylab("Subject ID") + 
               ggplot2::labs(fill = "Sex") +
               ggplot2::scale_x_discrete(position = "top") +
               ggplot2::scale_fill_manual(values = c("F" = "#e41a1c", "M" = "#377eb8")) +
               ggplot2::theme(legend.position = "bottom",
                              axis.ticks = ggplot2::element_blank(),
                              axis.title.x = ggplot2::element_text(angle = 90, vjust = 0, hjust=1),
                              legend.key.size = unit(4, "cm"),
                              axis.line = ggplot2::element_blank(),
                              axis.text.y = ggplot2::element_text(size = 75, face = "bold"),
                              axis.title.y = ggplot2::element_text(size = 75, face = "bold"))
ggplot2::ggsave(filename = "gender.pdf", 
                plot = plot_gender, 
                width = 15, 
                height = 70, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)

# Plot age data
plot_age <- patient_table %>%
            left_join(rank_table, by = "patientID") %>%
            dplyr::mutate(patientID = forcats::fct_reorder(patientID, supression_rate),
                          age = 1) %>%
            ggplot2::ggplot(aes( y = patientID, x = age, fill = age_at_art)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::coord_equal() +
            ggplot2::theme_classic(base_size = 75) +
            ggplot2::xlab("Age at ART") +
            ggplot2::ylab("Subject ID") +
            ggplot2::labs(fill = "Age") +
            ggplot2::scale_x_discrete(position = "top") +
            ggplot2::scale_fill_viridis_b() + 
            ggplot2::theme(legend.position = "bottom",
                           axis.ticks = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.line = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_text(angle = 90, vjust = 0, hjust=1),
                           axis.ticks.y = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank(),
                           legend.key.size = unit(4, "cm"))
ggplot2::ggsave(filename = "age.pdf", 
       plot = plot_age, 
       width = 15, 
       height = 70, 
       device = "pdf", 
       units = "in", 
       limitsize = FALSE)

# Plot viral load
plot_vl <- vl_table %>% 
           left_join(rank_table, by = "patientID") %>%
           dplyr::mutate(patientID = forcats::fct_reorder(patientID, supression_rate)) %>%
           ggplot2::ggplot(aes(x = yearsRel, y = patientID)) +
           ggplot2::geom_point(aes(color = time_point, size = log10(result))) +
           ggplot2::geom_segment(x=-6500, 
                                 xend=5500, 
                                 aes(y=patientID, 
                                 yend=patientID)) +
           ggplot2::theme_classic(base_size = 75) +
           ggplot2::xlab("Viral load (cells/mL)\nDays relative to ART initiation") + 
           ggplot2::ylab("Patient ID") +
           ggplot2::labs(color = "Treatment status",
                         size = "Log10(VL)") +
           ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size=24))) +
           ggplot2::scale_x_continuous(position = "top", 
                                       limits = date_lim,
                                       breaks = seq(-6500, 5500, 1000)) + 
           ggplot2::scale_size_continuous(range = c(5, 35)) +
           ggplot2::scale_color_manual(values = c("Pre-ART" = "#7570b3", "Post-ART" = "#d95f02")) +
           ggplot2::theme(legend.position = "bottom",
                          axis.ticks = ggplot2::element_blank(),
                          axis.title.y = ggplot2::element_blank(),
                          axis.line = ggplot2::element_blank(),
                          axis.ticks.y = ggplot2::element_blank(),
                          axis.text.y = ggplot2::element_blank(),
                          legend.key.size = unit(4, "cm"),
                          axis.text.x.top = ggplot2::element_text(size = 64, angle = 90, vjust = 0, hjust=0))
ggplot2::ggsave(filename = "viralload.pdf", 
                plot = plot_vl, 
                width = 40, 
                height = 70, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)

# Plot CD4 count variation
plot_cc <- cc_table %>% 
           left_join(rank_table, by = "patientID") %>%
           dplyr::mutate(patientID = forcats::fct_reorder(patientID, supression_rate)) %>%
           ggplot2::ggplot(aes(x = yearsRel, y = patientID)) +
           ggplot2::geom_point(aes(color = time_point, size = resultsRel)) +
           ggplot2::geom_segment(x=-6500, 
                                 xend=5500, 
                                 aes(y=patientID, 
                                 yend=patientID)) +
           ggplot2::theme_classic(base_size = 75) +
           ggplot2::xlab("CD4 count (cell/uL)\nDays relative to ART initiation") + 
           ggplot2::ylab("Patient ID") +
           ggplot2::labs(color = "Treatment status",
                         size = "CD4 Count") +
           ggplot2::guides(color = FALSE) +
           ggplot2::scale_x_continuous(position = "top", 
                                       limits = date_lim,
                                       breaks = seq(-6500, 5500, 1000)) + 
           ggplot2::scale_size_continuous(range = c(5, 45)) +
           ggplot2::scale_color_manual(values = c("Pre-ART" = "#7570b3", "Post-ART" = "#d95f02")) +
           ggplot2::theme(legend.position = "bottom",
                          axis.ticks = ggplot2::element_blank(),
                          axis.title.y = ggplot2::element_blank(),
                          axis.line = ggplot2::element_blank(),
                          axis.ticks.y = ggplot2::element_blank(),
                          axis.text.y = ggplot2::element_blank(),
                          legend.key.size = unit(4, "cm"),
                          axis.text.x.top = ggplot2::element_text(size = 64, angle = 90, vjust = 0, hjust=0))
ggplot2::ggsave(filename = "cdfcount.pdf", 
                plot = plot_cc, 
                width = 40, 
                height = 70, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)

# Plot clonality variation 
plot_cl<-  itable %>%
           left_join(rank_table, by = "patientID") %>%
           dplyr::mutate(patientID = forcats::fct_reorder(patientID, supression_rate)) %>%
           #dplyr::mutate(patientID = fct_rev(patientID)) %>%
           ggplot2::ggplot(aes(x = yearsRel, y = patientID)) +
           ggplot2::geom_point(aes(color = time_point, size = clonality)) +
           ggplot2::geom_segment(x=-6500, 
                                 xend=5500, 
                                 aes(y=patientID, 
                                     yend=patientID)) +
           ggplot2::theme_classic(base_size = 75) +
           ggplot2::xlab("T-cell repertoire sequencing\nDays relative to ART initiations") + 
           ggplot2::ylab("Patient ID") +
           ggplot2::labs(color = "Treatment status",
                         size = "Clonality of \n TCR repertoire") +
           ggplot2::guides(color = FALSE) +
           ggplot2::scale_x_continuous(position = "top", 
                                       limits = date_lim,
                                       breaks = seq(-6500, 5500, 1000)) + 
           ggplot2::scale_size_continuous(range = c(5, 45)) +
           ggplot2::scale_color_manual(values = c("Pre-ART" = "#7570b3", "Post-ART" = "#d95f02")) +
           ggplot2::theme(legend.position = "bottom",
                          axis.ticks = ggplot2::element_blank(),
                          axis.title.y = ggplot2::element_blank(),
                          axis.line = ggplot2::element_blank(),
                          axis.ticks.y = ggplot2::element_blank(),
                          axis.text.y = ggplot2::element_blank(),
                          legend.key.size = unit(4, "cm"),
                          axis.text.x.top = ggplot2::element_text(size = 64, angle = 90, vjust = 0, hjust=0))
ggplot2::ggsave(filename = "clonality.pdf", 
                plot = plot_cl, 
                width = 40, 
                height = 70, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)

# Plot complication and co-infection grid plot
plot_in <-  intable %>%
            left_join(rank_table, by = "patientID") %>%
            dplyr::mutate(patientID = forcats::fct_reorder(patientID, supression_rate),
                          infected = if_else(is.na(dateInfection), "Absent", "Present")) %>%
            ggplot2::ggplot(aes( y = patientID, x = infection, fill = infected)) +
            ggplot2::geom_tile(color = "black", size =2 ) +
            ggplot2::theme_classic(base_size = 75) +
            ggplot2::xlab("Co-infections / Complications") +
            ggplot2::ylab("Patient ID") +
            ggplot2::guides(fill = FALSE) +
            ggplot2::scale_fill_manual(values = c("Absent" = "#ffffff",
                                                  "Present" = "#969696"))+
            ggplot2::scale_x_discrete(labels = inpal,
                                      position = "top") +
            ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.line = ggplot2::element_blank(),
                           axis.ticks.y = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank(),
                           legend.key.size = unit(4, "cm"),
                           axis.text.x.top = ggplot2::element_text(size = 64, angle = 90, vjust = 0, hjust=0))
ggplot2::ggsave(filename = "infections.pdf", 
                plot = plot_in, 
                width = 30, 
                height = 70, 
                device = "pdf", 
                units = "in", 
                limitsize = FALSE)        

# Compose plots
comp_plot <- patchwork::wrap_plots(plot_gender, 
                                   plot_age,  
                                   plot_in, 
                                   plot_vl, 
                                   plot_cc, 
                                   plot_cl, 
                                   nrow = 1, 
                                   widths = c(10,10,15,45,45,45)) + 
             patchwork::plot_layout()
ggplot2::ggsave("CFAR_summary_raw.pdf", 
                plot = comp_plot, 
                width = 170, 
                height = 90, 
                units = "in", 
                device = "pdf", 
                limitsize = FALSE)