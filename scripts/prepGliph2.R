library(tidyverse)     
library(lubridate)
library(readxl)
library(optparse)
# Accept arguments
option_list <- list(
    optparse::make_option(c("-a", "--amino_data"), 
                type="character",  
                help="Path to annotated amino acid collapsed RDA file",
                dest="apath"),
    optparse::make_option(c("-m", "--hla_data"), 
            type="character",  
            help="Path to HLA table, the first column should correspond to sample name",
            dest="hpath"),
    optparse::make_option(c("-r", "--ref_file"), 
                type="character",  
                help="Path to GLIPH2 reference",
                dest="gref"),
    optparse::make_option(c("-l", "--lref_file"), 
                type="character",  
                help="Path to GLIPH2 length distribution file",
                dest="lref"),
    optparse::make_option(c("-v", "--vref_file"), 
                type="character",  
                help="Path to GLIPH2 vgene distribution file",
                dest="vref"),
    optparse::make_option(c("-s", "--study_id"), 
                type="character",  
                help="Study ID",
                dest="sid")            
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE)
gliph_amino_table <- readRDS(parser$apath)

atable <- gliph_amino_table %>% 
    mutate(cdrA = NA) %>%
    select(junction_aa, v_call, j_call, cdrA, gliph_id, duplicate_frequency) 

htable <- read_tsv(parser$hpath) # HLA table
sid <- parser$sid
gliph_in <- paste(sid, "GLIPH2", "in.tsv", sep = "_")
hla_in <- paste(sid, "HLA", "in.tsv", sep = "_")
gref_in <- paste(sid, "Gref", "in.tsv", sep = "_")
vref_in <- paste(sid, "Vref", "in.tsv", sep = "_")
lref_in <- paste(sid, "Lref", "in.tsv", sep = "_")

write_tsv(atable, gliph_in, col_names = FALSE)
write_tsv(htable, hla_in, col_names = FALSE)

gliph_config <- paste(sid, "gliph.txt", sep = "_")

#            paste("hla_file", hla_in, sep = "="),

config <- c(paste("out_prefix", sid, sep = "="),
            paste("cdr3_file", gliph_in, sep = "="),
            paste("hla_file", hla_in, sep = "="),
            paste("refer_file", parser$gref, sep = "="),
            paste("v_usage_freq_file", parser$vref, sep = "="),
            paste("cdr3_length_freq_file", parser$lref , sep = "="),
            "local_min_pvalue=0.001",
            "p_depth = 1000",
            "global_convergence_cutoff = 1",
            "simulation_depth=1000",
            "kmer_min_depth=3",
            "local_min_OVE=10",
            "algorithm=GLIPH2",
            "all_aa_interchangeable=1")
writeLines(config, gliph_config)
