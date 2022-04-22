#' Generate the cohort network figure. 
#' @description
#' Generate the nextwork plots comparing the control and PLHIV cohort
#' @return Cohort summary figure
library(tidyverse)
library(LymphoSeq2)
library(patchwork)
library(tidygraph)
library(lubridate)
library(ggraph)
library(ggpubr)
library(ggalluvial)
library(readxl)
library(knitr)
library(scales)
library(optparse)
library(DT)
library(ComplexHeatmap)
set.seed(12357) # Set the seed for reproducibility
option_list <- list(
    make_option(c("-a", "--amino_table"),
        type = "character",
        help = "Path to the sampled amino acid table",
        dest = "amino_table"),
    make_option(c("-m", "--meta_table"),
        type = "character",
        help = "Path to the study meta data in RDA format",
        dest = "meta_table"),
    make_option(c("-r", "--repertoire_summary_table"),
        type = "character",
        help = "Path to the repertoire summary tables",
        dest = "rep_table")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list, add_help_option=FALSE), 
    print_help_and_exit = TRUE) # Parse the arguments

load(parser$amino_table)
load(parser$meta_table)
data_summary <- read_tsv(parser$rep_table)

comp_palette <- c("Post-ART" = "#f1a340", "Pre-ART" = "#998ec3", "Normal" = "#000000")
comp_one <- list(c("Pre-ART", "Post-ART"), c("Pre-ART", "Normal"), c("Post-ART", "Normal"))
comp_two <- list(c("Pre-ART", "Post-ART"))

annotate_table <- amino_table %>% 
    left_join(study_sample_table, by = "repertoire_id") %>% 
    left_join(study_patient_table, by = "patientID") %>%
    mutate(lengths = nchar(junction_aa)) %>% 
    select(repertoire_id, time_group, time_point, junction_aa, lengths) 

readGraphs <- function(graph_table) {
    load(graph_table)
    repertoire_id <- graph_table %>% 
        basename() %>% 
        str_remove("_graph.rda")
    number_of_cluster <- graph %>% 
        as_tibble() %>% 
        pull(community) %>% 
        unique() %>% 
        length()
    max_cluster_size <- graph %>% 
        as_tibble() %>% 
        pull(occupency) %>% 
        max()
    graph_table <- tibble(repertoire_id = repertoire_id,
        number_of_clusters = number_of_cluster,
        max_cluster_size = max_cluster_size)
    return(graph_table)
}

getDegreeFrequency <- function(graph_table) {
  load(graph_table) 
  repertoire_id <- graph_table %>% 
    basename() %>% 
    str_remove("_graph.rda") 
  degree_table <- graph %>% 
    activate(nodes) %>% 
    mutate(degree = centrality_degree()) %>% 
    arrange(desc(degree)) %>% 
    as_tibble() %>% 
    mutate(repertoire_id = repertoire_id)
  return(degree_table)
}

getNetworks <- function(rep_id, study_sample_table) {
  graph_table <- paste("./", rep_id, "_graph.rda", sep = "")
  time_group <- study_sample_table %>% 
    filter(repertoire_id == rep_id) %>% 
    pull(years_rel_to_art)
  time_point <- study_sample_table %>% 
    filter(repertoire_id == rep_id) %>% 
    pull(time_point)
  load(graph_table) 
  set.seed(12357)
  graph <- graph %>% 
    activate(nodes) %>% 
    mutate(degree = centrality_degree(),
           stong_community = as_factor(group_components(type="strong")))
  comm_number <- graph %>% 
    activate(nodes) %>% 
    pull(stong_community) %>% 
    unique() %>% 
    length()
  if (time_point == "Pre-ART") {
    gradient_scale <- colorRampPalette(c("#cab2d6", "#998ec3"))
  } else if (time_point == "Post-ART") {
    gradient_scale <- colorRampPalette(c("#fdbf6f", "#f1a340"))
  } else if (time_point == "Normal") {
    gradient_scale <- colorRampPalette(c("#525252", "#d9d9d9"))
  }
  gradient_pals <- gradient_scale(comm_number)
  
  network <- ggraph(graph, layout = "stress") +
    geom_edge_link() +
    geom_node_point(aes(colour = stong_community)) +
    scale_colour_manual(values = gradient_pals)  +
    coord_fixed() +
    theme_void(base_size = 20) +
    ggtitle(as.character(time_group)) +
    theme(legend.position = "none")
  return(network)
}
graph_table <- list.files(path = "./", pattern = "_graph.rda", full.names = TRUE) %>% 
    purrr::map(readGraphs) %>%
    bind_rows()

degree_table <- list.files(path = "./", pattern = "_graph.rda", full.names = TRUE) %>% 
    purrr::map(getDegreeFrequency) %>%
    bind_rows() 

graph_table <- annotate_table %>% 
    select(repertoire_id, time_group, time_point) %>% 
    distinct() %>% 
    right_join(graph_table, by = "repertoire_id")

degree_table <- annotate_table %>% 
  select(repertoire_id, time_group, time_point) %>% 
  distinct() %>% 
  right_join(degree_table, by = "repertoire_id")

network_structure <- degree_table %>% 
  group_by(repertoire_id, degree) %>% 
  summarize(count = n(), time_point = first(time_point)) %>% 
  ungroup() %>% 
  group_by(repertoire_id) %>% 
  mutate(frequency = (count/sum(count)) * 100) %>% 
  ggplot(aes(x = degree, y = frequency,color = time_point )) +
  geom_point(alpha = 0.05) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs") ) +
  scale_color_manual(values = comp_palette) + 
  scale_x_log10() +
  scale_y_log10() + 
  theme_classic(base_size = 20) +
  labs(y = "Cummulative\nFrequency", x = "Degree (k)", color = "Time point") +
  theme(strip.background = element_blank(), 
        strip.placement = "outside",
        legend.key = element_rect(colour = NA, fill = NA),
        legend.position = c(0.8, 0.8)) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 22) ) )



network_repertoires <- c("015V09001205_CFAR", "015V09003735_CFAR", "015V17000496_CFAR", "HIP09062")
annotated_sample_table <- study_sample_table %>% 
  left_join(study_patient_table, by = "patientID") %>% 
  mutate(time_from_art =  years_rel_to_art)
preplot <- getNetworks("015V09003724_CFAR", annotated_sample_table)
postplot <- getNetworks("015V10002397_CFAR", annotated_sample_table)
lastplot <- getNetworks("015V16001483_CFAR", annotated_sample_table)
normplot <- getNetworks("HIP09062", annotated_sample_table)
design <- "
  11
  23
  45
"
network_plot <- network_structure + preplot + postplot + lastplot + normplot +  plot_layout(design = design , heights = c(6, 6,6,6,6))
ggsave("Figure7_CFAR_networks.pdf", network_plot, device = "pdf", width = 14,
    height = 12, units = "in", limitsize = FALSE)