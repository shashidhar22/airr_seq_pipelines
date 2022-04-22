library(tidyverse)
library(lubridate)
devtools::load_all("~/projects/LymphoSeq2")
option_list <- list(
    optparse::make_option(c("-a", "--amino_data"), 
                type="character",  
                help="Path to amino data RDA file",
                dest="adata"),
    optparse::make_option(c("-g", "--gliph"),
                          type="numeric",
                          help="Gliph path",
                          dest="gpath")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                               print_help_and_exit = TRUE)

writeGliph <- function(gtable) {
    repertoire_id <- gtable %>% 
                     pull(repertoire_id) %>% 
                     unique()
    out_path <- paste(repertoire_id, "tsv", sep = ".")
    gtable <- gtable %>% 
              group_by(pattern) %>% 
              summarize(scount = n(), 
                        junction_aa = paste(TcRb, collapse = " ")) %>% 
              select(scount, pattern, junction_aa)
    write_tsv(gtable, out_path, col_names = FALSE)
}


gtable <- read_csv(parser$gpath) %>% 
          select(number_unique_cdr3, pattern, TcRb, Sample) %>% 
          separate(Sample, c("subject_id", "days_rel"), sep = ":") %>% 
          filter(subject_id != "NA" & days_rel != "NA") %>%
          mutate(days_rel = as.double(days_rel)) %>%  
          select(number_unique_cdr3, pattern, TcRb, subject_id, days_rel)

load(parser$adata)

atable <- atable %>% 
          mutate(days_rel = (lubridate::as_date(dateART) %--% dateTCR) %/% lubridate::days(1)) %>% 
          select(patientID, days_rel, repertoire_id, junction_aa) %>%
          distinct() 

gtable <- left_join(gtable, atable, by = c("subject_id" = "patientID", "days_rel", "TcRb" = "junction_aa")) 

gtable %>% 
group_by(repertoire_id) %>% 
group_split() %>% 
purrr::map(writeGliph)