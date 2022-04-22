library(tidyverse)
library(patchwork)
library(lubridate)
library(readxl)
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

atable <- atable %>%
          dplyr::filter(patientID == parser$pid) %>%
          dplyr::mutate(yearsRel = (dateART %--% dateTCR) %/% lubridate::days(1),
                        yearsRel = forcats::as_factor(yearsRel),
                        time_point = dplyr::if_else(as.integer(yearsRel) < 0, "Pre-ART", "Post-ART"))
ttens <- atable %>%
         topSeqs(top = 10) %>% 
         dplyr::group_by(yearsRel) %>% 
         dplyr::summarize(cummulative_freq = sum(duplicate_frequency), top = 10, time_point = first(time_point)) %>% 
         dplyr::ungroup()
thund <- atable %>%
         topSeqs(top = 100) %>% 
         dplyr::group_by(yearsRel) %>% 
         dplyr::summarize(cummulative_freq = sum(duplicate_frequency), top = 100, time_point = first(time_point)) %>% 
         dplyr::ungroup()
tthou <- atable %>%
         topSeqs(top =  1000) %>% 
         dplyr::group_by(yearsRel) %>% 
         dplyr::summarize(cummulative_freq = sum(duplicate_frequency), top = 1000, time_point = first(time_point)) %>% 
         dplyr::ungroup()

top_table <- dplyr::bind_rows(ttens, thund) %>%
             dplyr::bind_rows(tthou) 

top_plot <- top_table %>%
            dplyr::mutate(top = forcats::as_factor(top)) %>% 
            ggplot2::ggplot(aes(x = yearsRel, y = cummulative_freq, color = top, group = top)) +
            ggplot2::geom_point() +
            ggplot2::geom_line() + 
            ggplot2::xlab("Days relative to ART initiation") + 
            ggplot2::ylab("Cumulative frequency") +
            ggplot2::labs(color = "Treatment status",
                         group = "Top n sequences") +
            ggplot2::theme_classic(base_size = 14) +
            ggplot2::ggtitle(paste("Subject ID", parser$pid, sep = " : ")) 

out_path <- paste(parser$pid, "topSeqsPercent.pdf", sep = "_")
ggplot2::ggsave(out_path, 
                top_plot, 
                width = 8,
                height = 8,
                unit = "in",
                device = "pdf")


