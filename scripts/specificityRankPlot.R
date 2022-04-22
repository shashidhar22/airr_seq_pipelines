library(tidyverse)
library(lubridate)
devtools::load_all("~/projects/LymphoSeq2")
option_list <- list(
    optparse::make_option(c("-a", "--amino_data"), 
                type="character",  
                help="Path to amino data RDA file",
                dest="adata"),
    optparse::make_option(c("-p", "--patient_id"),
                          type="numeric",
                          help="Patient ID",
                          dest="pid")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                               print_help_and_exit = TRUE)


load(parser$adata)

readGliph <- function(atable) {
    repertoire_id <- atable %>%
                   dplyr::pull(repertoire_id) %>%
                   unique() 
    gpath <- paste(repertoire_id, "tsv", sep = ".")
    gtable <- readr::read_tsv(gpath, col_names = c("count", "sgroup", "sequences")) %>%
              tidyr::separate_rows(sequences, sep =" ") %>%
              dplyr::select(count, sgroup, sequences)
    atable <- left_join(atable, gtable, by = c("junction_aa" = "sequences"))
    atable <- atable %>% 
              dplyr::group_by(sgroup) %>% 
              dplyr::arrange(desc(duplicate_count)) %>%
              dplyr::mutate(rank = dplyr::row_number(desc(duplicate_count))) %>%
              dplyr::ungroup()
}

plotRank <- function(rtable) {
    # Prepare rank data for plotting
    patient <- rtable %>%
               dplyr::pull(patientID) %>%
               unique()
    spec <- rtable %>% 
            dplyr::pull(sgroup) %>% 
            unique()
    rtable  <- rtable %>%
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
                 filter(rank <= 100) %>%
                 ggplot2::ggplot(aes(x=rank, y=duplicate_frequency, color = yearsRel)) +
                 ggplot2::geom_point() +
                 ggplot2::geom_line() +
                 ggplot2::theme_classic() + 
                 ggplot2::ylab("CDR3 frequency") + 
                 ggplot2::xlab("Rank") +
                 ggplot2::scale_color_manual(values = treatment_palette) +
                 ggplot2::scale_x_log10() +
                 ggplot2::scale_y_log10() +
                 ggplot2::labs(color = "Days from ART") +
                 ggplot2::ggtitle(paste("Patient ID:", patient, "; Specificity group:", spec)) +
                 ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=12, face="bold"))
    return(rank_plot)
}

atable <- atable %>% 
          dplyr::filter(patientID == parser$pid) %>% 
          dplyr::group_by(repertoire_id) %>% 
          dplyr::group_split() %>% 
          purrr::map(readGliph) %>% 
          dplyr::bind_rows() 

top_spec <- atable %>% 
            dplyr::group_by(sgroup) %>% 
            dplyr::summarize(scount = max(count)) %>% 
            dplyr::ungroup() %>% 
            dplyr::arrange(desc(scount)) %>% 
            dplyr::top_n(n = 5) %>% 
            dplyr::pull(sgroup) 

rplots <- atable %>% 
          dplyr::filter(sgroup %in% top_spec) %>% 
          dplyr::group_by(sgroup) %>% 
          dplyr::group_split() %>%  
          purrr::map(plotRank) %>% 
          patchwork::wrap_plots() +
          patchwork::plot_layout(guides = "collect")

out_path <- paste(parser$pid, "specrank.pdf", sep = "_")
ggsave(out_path, plot = rplots, width = 18, height = 12, units = "in", device = "pdf", limitsize = FALSE)
