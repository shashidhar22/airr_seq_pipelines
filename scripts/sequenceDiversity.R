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
devtools::load_all("~/projects/LymphoSeq2")
# Accept arguments
option_list <- list(
    optparse::make_option(c("-a", "--amino_data"), 
                          type="character",  
                          help="Path to productive sequence RDA file",
                          dest="adata"),
    optparse::make_option(c("-i", "--immuno_data"), 
                          type="character",  
                          help="Path to ImmunoSeq summary RDA file",
                          dest="idata"),
    optparse::make_option(c("-m", "--meta_data"), 
                          type="character",  
                          help="Path to metadta RDA file",
                          dest="mdata"),
    optparse::make_option(c("-p", "--patient_id"),
                          type="numeric",
                          help="Patient ID",
                          dest="pid")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                               print_help_and_exit = TRUE)
load(parser$adata)
load(parser$idata)
load(parser$mdata)
# Calculate rank of sequence
getRank <- function(atable) {
    atable <- atable %>% 
              dplyr::group_by(repertoire_id) %>%
              dplyr::arrange(desc(duplicate_count)) %>%
              dplyr::mutate(rank = dplyr::row_number(desc(duplicate_count))) %>%
              dplyr::ungroup()
    return(atable)
}
# Generate rank abundance plots
plotRank <- function(rtable, itable) {
    # Prepare rank data for plotting
    patient <- rtable %>%
               dplyr::pull(patientID) %>%
               unique()
    rtable  <- rtable %>%
               dplyr::mutate(yearsRel = (dateART %--% dateTCR) %/% lubridate::days(1),
                      yearsRel = forcats::as_factor(yearsRel),
                      yearsRel = forcats::fct_rev(yearsRel),
                      time_point = dplyr::if_else(as.integer(yearsRel) < 0, "Pre-ART", "Post-ART"))
    itable <- itable %>% 
              dplyr::mutate(yearsRel = (dateART %--% dateTCR) %/% lubridate::days(1),
                      yearsRel = forcats::as_factor(yearsRel),
                      yearsRel = forcats::fct_rev(yearsRel),
                      time_point = dplyr::if_else(as.integer(yearsRel) < 0, "Pre-ART", "Post-ART"))
    # Generate rank abundance line graphs
    art_palette <- colorRampPalette(c("#d95f02", "#7570b3"))
    times <- rtable %>%
             dplyr::select(patientID, yearsRel) %>%
             dplyr::mutate(yearsRel = as.integer(yearsRel)) %>%
             dplyr::arrange(yearsRel) %>% 
             dplyr::pull(yearsRel) %>%
             unique()
    treatment_timepoint <- rtable %>% 
                           dplyr::mutate(year_ints = as.integer(yearsRel)) %>%
                           dplyr::arrange(year_ints) %>%
                           dplyr::mutate(yearsRel = forcats::fct_reorder(yearsRel, year_ints)) %>%
                           dplyr::pull(yearsRel) %>%
                           unique()
    treatment_palette <- art_palette(length(treatment_timepoint))
    names(treatment_palette) <- treatment_timepoint
    rank_plot <- rtable %>%
                 filter(rank <= 1000) %>%
                 ggplot2::ggplot(aes(x=rank, y=duplicate_frequency, color = yearsRel)) +
                 ggplot2::geom_point(size = 0.5) +
                 ggplot2::geom_line(size = 0.5) +
                 ggplot2::theme_classic(base_size = 12) + 
                 ggplot2::ylab("CDR3 frequency") + 
                 ggplot2::xlab("Rank") +
                 ggplot2::scale_color_manual(values = treatment_palette) +
                 ggplot2::scale_x_log10() +
                 ggplot2::scale_y_log10() +
                 ggplot2::labs(color = "Days from ART") +
                 ggplot2::theme(legend.position = "none")
    # Generate rarefaction curves
    rep_year <- rtable %>%
                dplyr::select(repertoire_id, yearsRel)
    rarefaction_tables <- rtable %>% 
                          dplyr::group_by(yearsRel) %>% 
                          dplyr::group_split() %>% 
                          purrr::map(runINext) %>% 
                          dplyr::bind_rows()
    rarefaction_tables <- dplyr::left_join(rarefaction_tables, rep_year, by = "repertoire_id")
    rarefaction_tables <- rarefaction_tables %>%
                          dplyr::mutate(method = dplyr::recode(method, 
                                                               observed = "Interpolated", 
                                                               interpolated = "Interpolated", 
                                                               extrapolated = "Extrapolated"),
                                        yearsRel = forcats::as_factor(yearsRel),
                                        yearsRel = forcats::fct_rev(yearsRel),
                                        time_point = if_else(as.integer(yearsRel) < 0, "Pre-ART", "Post-ART")) 
    points_table  <- rarefaction_tables %>% 
                     group_by(repertoire_id, method) %>% 
                     arrange(m) %>% 
                     slice_tail(n = 1) %>% 
                     ungroup()
    rarefaction_curves <- ggplot2::ggplot(rarefaction_tables, 
                                          aes(x = m, y = qD, color = yearsRel, fill= yearsRel)) + 
                          ggplot2::geom_line(size = 0.5) + 
                          ggplot2::geom_point(data = points_table, aes(x = m, y = qD, shape = method, color = yearsRel, fill= yearsRel), size = 1) +
                          ggplot2::geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha=0.2) +
                          ggplot2::theme_classic(base_size = 12) +
                          ggplot2::xlab("Total number of sequences") +
                          ggplot2::ylab("TCR diversity") +
                          ggplot2::scale_color_manual(values = treatment_palette) +
                          ggplot2::scale_fill_manual(values = treatment_palette) +
                          ggplot2::scale_y_continuous(labels=scales::label_number(), breaks = scales::pretty_breaks(4)) +
                          ggplot2::scale_x_continuous(labels=scales::label_number(), breaks = scales::pretty_breaks(4)) +
                          ggplot2::labs(fill = "Treatment status", 
                                        color="Treatment status", 
                                        shape="Count status") +
                          ggplot2::guides(color = FALSE,
                                          fill = FALSE,
                                          shape = ggplot2::guide_legend(order=2)) +
                          ggplot2::theme(legend.position = c(0.25, 0.85)) 
    # Generate diveristy bar plots
    prod_graph <- itable %>% 
                  mutate(yearsRel = forcats::fct_rev(yearsRel)) %>% 
                  ggplot2::ggplot(aes(x = yearsRel, y= as.numeric(unique_productive_sequences), color= yearsRel, fill = yearsRel, group = patientID)) +
                  ggplot2::geom_line(size = 1) +
                  ggplot2::geom_point(size = 2) +
                  ggplot2::scale_color_manual(values = treatment_palette) + 
                  ggplot2::theme_classic(base_size = 12) +
                  geom_vline(xintercept = 0) +
                  ggplot2::ylab("Unique productive sequences") +
                  ggplot2::xlab("Days from ART initiation") +
                  ggplot2::scale_y_continuous(labels=scales::label_number(), breaks = scales::pretty_breaks(4)) +
                  ggplot2::labs(fill = "Pre or Post ART", color="Pre or Post ART") +
                  ggplot2::guides(color = guide_legend(nrow = 1), fill = FALSE) +
                  ggplot2::theme(legend.position = "bottom")
    #Write component plots
    prod_out <- paste(patient, "prod.pdf", sep = "_")
    rare_out <- paste(patient, "rarefaction.pdf", sep = "_")
    rank_out <- paste(patient, "rankplot.pdf", sep = "_")
    ggsave(plot = prod_graph, 
           filename = prod_out,
           width = 4,
           height = 4,
           device = "pdf",
           units = "in")
    ggsave(plot = rarefaction_curves, 
           filename = rare_out,
           width = 9,
           height = 8,
           device = "pdf",
           units = "in")
    ggsave(plot = rank_plot, 
           filename = rank_out,
           width = 9,
           height = 8,
           device = "pdf",
           units = "in")
    # Assemble and save the plot
    design <- "111112222233333
               111112222233333
               111112222233333"
    
    rarefaction_curves <- rarefaction_curves + 
                          prod_graph + 
                          rank_plot + 
                          plot_layout(design=design, widths=c(5, 5, 5), heights = 3) +
                          plot_annotation(title = patient,
                                          theme = theme(plot.title = element_text(face="bold")))
    out_file <- paste(patient, "diversity.pdf", sep = "_")
    ggsave(rarefaction_curves,
            filename = out_file, 
            height = 3.5,
            width = 15,
            device = "pdf",
            units = "in")
}

# Generate rank abundance compound plots
atable <- atable %>%
          dplyr::filter(patientID == parser$pid)
itable <- itable %>% 
          dplyr::filter(patientID == parser$pid)
rank_table <- getRank(atable) 
plotRank(rank_table, itable)
