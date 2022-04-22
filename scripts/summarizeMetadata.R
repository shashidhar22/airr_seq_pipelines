library(LymphoSeq2)
library(tidyverse)     
library(lubridate)
library(readxl)
library(patchwork)
library(lubridate)
library(readxl)
library(viridis)
library(fossil)
library(optparse)
# Collect arguments
option_list <- list(
    make_option(c("-m", "--meta_data"), 
                type="character",  
                help="Path to RDA file containing metadata information",
                dest="meta")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
# Load metadata
load(parser$meta)
plot_gender <- patient_table %>%
                dplyr::mutate(patientID = fct_rev(patientID),
                                gender =  1) %>%
                ggplot2::ggplot(aes(y= patientID, x = gender, fill = sex)) + 
                ggplot2::geom_tile(color = "white") + 
                ggplot2::coord_equal() +
                ggplot2::theme_classic(base_size = 46) +
                ggplot2::xlab("Sex") + 
                ggplot2::ylab("Patient ID") + 
                ggplot2::labs(fill = "Sex") +
                ggplot2::scale_x_discrete(position = "top") +
                ggplot2::scale_fill_manual(values = c("#d95f02", "#1b9e77")) +
                ggplot2::theme(axis.ticks = element_blank(),
                                legend.key.size = unit(2, "cm"))
plot_age <- patient_table %>%
            dplyr::mutate(patientID = fct_rev(patientID),
                        age = 1) %>%
            ggplot2::ggplot(aes( y = patientID, x = age, fill = age_at_art)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::coord_equal() +
            ggplot2::theme_classic(base_size = 46) +
            ggplot2::xlab("Age at \n ART") +
            ggplot2::ylab("Patient ID") +
            ggplot2::labs(fill = "Age") +
            ggplot2::scale_x_discrete(position = "top") +
            ggplot2::scale_fill_viridis_b() + 
            ggplot2::theme(axis.ticks = element_blank(),
                        axis.title.y = element_blank(),
                        axis.line.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.text.y = element_blank(),
                        legend.key.size = unit(2, "cm"))
# Create plot for HLA alleles
hla_plot <-  hla_table %>%
            dplyr::mutate(patientID = fct_rev(patientID)) %>%
            ggplot2::ggplot(aes(x = allele, y = patientID)) +
            ggplot2::facet_grid(cols = vars(HLA), as.table =TRUE, margins = FALSE) +
            ggplot2::geom_tile(color = "white", aes(fill = HLA)) +
            ggplot2::geom_text(aes(label = value), size = 12, fontface = "bold") + 
            ggplot2::theme_classic(base_size = 46) +
            ggplot2::xlab("HLA") +
            ggplot2::ylab("Patient ID") +
            ggplot2::scale_x_discrete(position = "top") +
            ggplot2::scale_fill_manual(values = c("#fbb4ae", "#b3cde3", "#ccebc5",
                                                "#decbe4", "#fed9a6", "#ffffcc",
                                                "#e5d8bd", "#fddaec", "#f2f2f2")) +
            ggplot2::theme(axis.ticks = element_blank(),
                        axis.title.y = element_blank(),
                        axis.line.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.title.x = element_blank(),
                        legend.position = "none", 
                        strip.background = element_rect(colour=NA, fill=NA),
                        strip.placement = "outside",
                        strip.text = element_text(size = 44),
                        panel.grid = element_blank())
# Create ImmunoSeq plot
samp_table  <- sample_table %>% 
            dplyr::filter(timePoint != "Normal") %>% 
            dplyr::group_by(patientID, timePoint) %>% 
            dplyr::summarize(count = n()) %>%
            dplyr::mutate(patientID = as.integer(patientID), 
                            patientID = as.factor(patientID))

sample_plot <- samp_table %>%
            dplyr::mutate(patientID = fct_rev(patientID)) %>%
            ggplot2::ggplot(aes(y = patientID, x= timePoint, fill= timePoint)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::coord_equal() +
            ggplot2::labs(fill = "Time point") +
            ggplot2::xlab("ImmunoSeq") +
            ggplot2::scale_x_discrete(position = "top") +
            ggplot2::geom_text(aes(label = count), size = 12, fontface = "bold") +
            ggplot2::theme_classic(base_size = 46) +
            ggplot2::theme(axis.ticks = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.x = element_blank(),
                            axis.line.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.text.y = element_blank(),
                            legend.key.size = unit(2, "cm"))
# Combine plot
meta_summary <- plot_gender + 
                plot_age + 
                hla_plot  + 
                plot_layout(guides = "collect")
ggplot2::ggsave(meta_summary, 
                filename = "CFAR_study_summary.pdf", 
                height = 30, 
                width = 65, 
                units = "in", 
                device = "pdf",
                limitsize = FALSE)
seq_summary <- plot_gender + 
                plot_age + 
                sample_plot  + 
                plot_layout(guides = "collect")
ggplot2::ggsave(seq_summary, 
                filename = "CFAR_immunoSeq_summary.pdf", 
                height = 30, 
                width = 15, 
                units = "in", 
                device = "pdf", 
                limitsize = FALSE)