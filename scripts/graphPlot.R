library(LymphoSeq2)
library(tidyverse)     
library(lubridate)
library(readxl)
library(patchwork)
library(viridis)
library(fossil)
library(tidygraph)
library(ggraph)
library(optparse)

# Accept arguments
option_list <- list(
        make_option(c("-t", "--plot_type"), 
                    type="character",  
                    help="Type of plot",
                    dest="ptype")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
# Load data
graph_files <- list.files(path = ".", pattern = "*_NT.tsv")
# Summarize graph files 
summarizeGraph <- function(rpath) {
    rgraph <- read_tsv(rpath)
    rgraph <- rgraph %>%
              dplyr::select(repertoire_id, patientID, dateTCR, occupency) %>%
              dplyr::group_by(repertoire_id) %>%
              dplyr::summarize(patientID = first(patientID),
                               dateTCR = ymd(first(dateTCR)),
                               total_nodes = n(),
                               max_node = max(occupency)) %>%
              dplyr::ungroup()
    return(rgraph)
}
# Plot sequence expansion
plotExpansion <- function(rgraph) {
    patient <- rgraph %>%
               pull(patientID) %>%
               unique()
    rgraph <- rgraph %>%
              mutate(node_ratio = max_node / total_nodes)
    rplot <- ggplot(rgraph, aes(x = as.factor(dateTCR), y = node_ratio, fill = as.character(dateTCR))) +
             geom_bar(stat = "identity") +
             theme_classic(base_size = 18) +
             scale_fill_viridis_d() +
             ggtitle(patient) +
             theme(plot.title = element_text(face="bold", hjust = 0.5))
}
rgraph <- graph_files %>%
          purrr::map(summarizeGraph) %>%
          dplyr::bind_rows()

rplot <- rgraph %>%
         dplyr::group_by(patientID) %>%
         dplyr::group_split() %>%
         purrr::map(plotExpansion) %>%
         patchwork::wrap_plots()

out_file <- paste("Sequence", parser$ptype, "expansion.pdf", sep = "_")
ggsave(rplot, filename = out_file, width = 60, height = 50, units = "in", device = "pdf", limitsize = FALSE)
