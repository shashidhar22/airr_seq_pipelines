#' Calculate edit distance and build network of all CDR3 sequences
#' @description
#' Generate a matrix of all the sample amino acid sequences and compute the 
#' edit pairwise edit distance between all sequences. Filter sequences that 
#' are within an edit distance of one. And generate a network of CDR3 sequences.
#' @param amino_table Path to the sampled amino acid table for each repertoire 
#' @return Repertoire network of related sequences
library(tidyverse)
library(LymphoSeq2)
library(readxl)
library(optparse)
library(igraph)
library(tidygraph)
library(ggraph)
set.seed(12357) # Set the seed for reproducibility
option_list <- list(
    optparse::make_option(c("-a", "--amino_table"),
        type = "character",
        help = "Path to the sampled amino acid table",
        dest = "amino_table"),
    optparse::make_option(c("-m", "--meta_table"),
        type = "character",
        help = "Path to the study meta data in RDA format",
        dest = "meta_table"),
    optparse::make_option(c("-n", "--min_sequences"),
        type="integer",
        help="Minimum number of sequences to be considered",
        dest="min_sequences")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE) # Parse the arguments
# Load meta data
load(parser$meta_table)
# Read the amino acid table
amino_table <- read_tsv(parser$amino_table)
# Sample the top 5000 sequences from the dataset
amino_table <- amino_table %>% 
    arrange(desc(duplicate_frequency)) %>%
    slice_head(n = 5000)
# Get meta information
repertoire_name <- amino_table %>% 
    pull(repertoire_id) %>%
    unique()
time_point <- study_sample_table %>% 
    filter(repertoire_id == repertoire_name) %>%
    pull(time_group)
patient_id <- study_sample_table %>% 
    filter(repertoire_id == repertoire_name) %>%
    pull(patientID) 
plot_title <- paste("Patient ID : ", patient_id, " | Time Point : ", 
    time_point, sep="")
# Compute the edit distance between all sequences. Convert the distance matrix
# to a tibble, and filter out sequences that are within an edit distance not
# equal to one.
edge_table <- stringdist::stringdistmatrix(amino_table$junction_aa, 
    amino_table$junction_aa, method = "lv")
colnames(edge_table) <- amino_table$junction_aa
edge_table <- edge_table %>% 
    as_tibble() %>% 
    mutate(junction_one = amino_table$junction_aa) %>% 
    select(junction_one, everything()) %>%
    pivot_longer(-junction_one, names_to = "junction_two", 
        values_to = "distance") %>% 
    filter(distance == 1) %>% 
    select(junction_one, junction_two) %>% 
    rename(from = junction_one, to = junction_two)
node_table <- unique(c(edge_table$from, edge_table$to)) %>% 
    as_tibble() %>%
    rename(junction_aa = value) %>% 
    left_join(amino_table, by = "junction_aa") %>%
    select(junction_aa, duplicate_count) %>% 
    rename(name = junction_aa, size = duplicate_count)
# Generate the network using tidygraph
graph <- tidygraph::tbl_graph(nodes = node_table,
    edges = edge_table)  %>% 
    mutate(relatedness = centrality_degree(mode = 'in'), 
        community = as_factor(group_components(type="strong")))  %>%
    group_by(community) %>% 
    mutate(occupency = n()) %>%
    ungroup() 
# Plot the network
graph_plot <- graph %>% 
    filter(occupency >= 4) %>%
    ggraph(layout = 'fr') +
    geom_edge_link(aes(alpha = ..index..),
        show.legend = FALSE) +
    geom_node_point(aes(size = size, color = community)) + 
    theme_void() +
    theme(legend.position = "none") +
    labs(title = plot_title)

# Generate output files
graph_file <- paste(repertoire_name, "graph.rda", sep = "_")
plot_file <- paste(repertoire_name, "network.pdf", sep = "_")
save(graph, file = graph_file)
ggsave(plot_file, graph_plot, width =  16 , height = 16, units = "in", 
    device = "pdf")
