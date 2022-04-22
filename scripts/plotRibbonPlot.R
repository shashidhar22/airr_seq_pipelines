library(tidyverse)
library(patchwork)
library(ggalluvial)
library(lubridate)
library(ggplot2)
library(viridis)
library(readxl)
library(scales)
library(lubridate)
devtools::load_all("~/projects/LymphoSeq2")
option_list <- list(
    optparse::make_option(c("-u", "--public_path"), 
                type="character",  
                help="Path to public TRB sequences",
                dest="ppath"),
    optparse::make_option(c("-a", "--amino_data"), 
                type="character",  
                help="Path to amino data RDA file",
                dest="adata"),
    optparse::make_option(c("-m", "--meta_data"), 
                type="character",  
                help="Path to meta data RDA file",
                dest="mdata"),
    optparse::make_option(c("-i", "--immuno_data"), 
                type="character",  
                help="Path to ImmunoSeq summary data RDA file",
                dest="idata"),
    optparse::make_option(c("-p", "--patient_id"),
                          type="numeric",
                          help="Patient ID",
                          dest="pid")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                               print_help_and_exit = TRUE)
load(parser$adata)
load(parser$mdata)
load(parser$idata)

# Function definitions
readPublic <- function(file_path) {
       ptable <- read_csv(file_path)
       sample_name <- basename(file_path) %>% str_extract("^\\w+")
       if ("aminoAcid.beta" %in% colnames(ptable)) {
              ptable <- ptable %>%
                        pivot_longer(cols= c(aminoAcid.beta),
                                     names_to = "Chain",
                                     values_to = "aminoAcid") %>%
                        filter(Chain == "aminoAcid.beta") %>%
                        mutate(antigen = sample_name) %>%
                        select(antigen, aminoAcid)
       } else if ("CDR3" %in% colnames(ptable)){
              ptable <- ptable %>%
                        mutate(antigen = sample_name) %>%
                        rename(aminoAcid = CDR3) %>%
                        select(antigen, aminoAcid)
       } else {
              ptable <- ptable %>%
                        mutate(antigen = sample_name) %>%
                        select(antigen, aminoAcid)              
       }
       ptable <- ptable %>% drop_na()
       return(ptable)
}

plotSankey <- function(sample_table, vtable, cdtable, dtable, cltable, pbfactor, ptid) {
    set.seed(12345)
    # Filter samples
    sample_table <- sample_table %>%
          dplyr::filter(patientID == ptid)
    vtable <- vtable %>%
                dplyr::filter(patientID == ptid)
    cdtable <- cdtable %>%
                dplyr::filter(patientID == ptid)
    dtable <- dtable %>%
                       dplyr::filter(patientID == ptid)
    cltable <- cltable %>%
              dplyr::filter(patientID == ptid)
    # Get top sequences
    sample_table <- sample_table %>%
                    LymphoSeq2::cloneTrack() 
    # Collect repertoire_id information and response variable
    sample_name <- sample_table %>% 
                   dplyr::pull(repertoire_id) %>% 
                   unique()
    patient_name <- sample_table %>%
                    dplyr::pull(patientID) %>%
                    unique()
    date_art <- sample_table %>%
                dplyr::pull(dateART) %>%
                unique()
    # Get viral load 
    vtable <- vtable %>% 
              dplyr::mutate(result = dplyr::case_when(result == 0 ~ 1,
                                    TRUE ~ as.numeric(result)),
                     value_type = "Viral Load",
                     days_rel = ( lubridate::as_date(date_art) %--% dateCollection) %/% lubridate::days(1))
    # Get CD4 count
    cdtable <- cdtable %>%
               mutate(value_type = "CD4 count",
                      days_rel = ( lubridate::as_date(date_art) %--% dateCollection) %/% lubridate::days(1))
    # Get infection table
    cval <- c("ART", "EBV", "HBV", "HCV", "HSV1", "HSV2", "Chlamydia",
              "Gonnorrhea", "Syphilis", "Kaposi's sarcoma", "Anal cancer", "Colorectal cancer",
              "Skin cancer: squamous cell carcinoma", "Esophagus cancer", "Cervical cancer")
    cpal <- c("#e31a1c", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", 
              "#33a02c", "#ff7f00", "#6a3d9a", "#b15928", "#140063", "#7a6fa6", 
              "#20154d", "#037b80", "#33edf5")
    cdict <- c("ART", "EBV", "HBV", "HCV", "HSV1", "HSV2", "Chlamydia", 
               "Gonnorrhea", "Syphillis", "Kaposi_sarcoma", "Anal_cancer", "Colorectal_cancer",
               "Skin_cancer_squamous_cell_carcinoma", "Esophagus_cancer", "Cervical_cancer")
    names(cpal) <- cval
    names(cval) <- cdict
    dtable <- dtable %>%
              tidyr::pivot_longer(-patientID, names_to = "Condition", values_to = "Day") %>%
              tibble::add_row(patientID = as.factor(patient_name), 
                              Condition = "ART", 
                              Day = lubridate::as_date(date_art)) %>%
              dplyr::mutate(Day = lubridate::as_date(Day),
                            days_rel = (lubridate::as_date(date_art) %--% Day) %/% lubridate::days(1))
    # Gather all dates and create a date range for the samples
    date_range <- unique(c(sample_table$yearsRel, 
                           cdtable$days_rel, 
                           vtable$days_rel, 
                           dtable$days_rel))
    min_date <- min(date_range, na.rm=TRUE) - 200
    max_date <- max(date_range, na.rm=TRUE) + 200
    date_lim <- c(min_date, max_date)
    # Color sequences by occurence
    fcommon <- sample_table %>% 
               dplyr::filter(seen > 1 & duplicate_count > 1) %>%
               dplyr::group_by(junction_aa) %>%
               dplyr::summarize(dateEarliest = min(dateTCR)) %>%
               dplyr::ungroup() %>%
               dplyr::mutate(foundStatus = dplyr::if_else(dateEarliest <= lubridate::as_date(date_art), "Pre-ART", "Post-ART"),
                             foundStatus = dplyr::if_else(junction_aa %in% names(pbfactor), pbfactor[junction_aa], foundStatus))
    # Now that we have identified the status for each sequence and the earliest date of occurrence. Let pull up the to 50 sequences
    atable <- sample_table
    sample_table <- sample_table %>%
                    topSeqs(top = 50)
    fcommon <- fcommon %>%
               filter(junction_aa %in% sample_table$junction_aa)
    all_seq <- unique(sample_table$junction_aa)
    pre_seq <- fcommon %>% filter(foundStatus == "Pre-ART") %>% dplyr::pull(junction_aa)
    pre_pal <- if (length(pre_seq) > 0) p1(length(pre_seq)) else NULL
    pos_seq <- fcommon %>% filter(foundStatus == "Post-ART") %>% dplyr::pull(junction_aa)
    pos_pal <- if (length(pos_seq) > 0) p2(length(pos_seq)) else NULL
    pub_seq <- fcommon %>% filter(foundStatus %in% pbfactor) %>% dplyr::pull(junction_aa)
    pub_pal <- if (length(pub_seq) > 0) viridis::viridis(length(pub_seq)) else NULL
    fnd_seq <- c(pre_seq, pos_seq, pub_seq)
    mis_seq <- setdiff(all_seq, fnd_seq)
    mis_pal <- p3(length(mis_seq))
    tot_seq <- c(pre_seq, pos_seq, pub_seq, mis_seq)
    tot_pal <- c(pre_pal, pos_pal, pub_pal, mis_pal)
    names(tot_pal) <- tot_seq    

    pub_lab <- paste(pbfactor[pub_seq], pub_seq, sep=" : ")
    # Generate clonaltiy table
    clonality_table <- itable %>% 
                       filter(patientID == patient_name) %>%
                       mutate(days_rel = (as_date(date_art) %--% dateTCR) %/% days(1)) %>%
                       select(days_rel, clonality)

    

    # Generate the alluvial plot to trace the frequency variation of top n Amino acid sequences
    sample_table <- dplyr::left_join(sample_table, 
                                     clonality_table, 
                                     by=c("yearsRel" = "days_rel"))
    sample_limits <- sample_table %>% 
                     dplyr::mutate(yearsRel = forcats::as_factor(yearsRel)) %>% 
                     dplyr::pull(yearsRel) %>%
                     unique()
    sample_names <- sample_table %>%
                    dplyr::pull(repertoire_id) %>%
                    unique()
    names(sample_limits) <- sample_names
    sample_table <- sample_table %>% 
                    dplyr::mutate(repertoire_id = yearsRel)
    sankey <- plotTrack(sample_table, alist = tot_seq, apal = tot_pal) +
              ggalluvial::geom_stratum(aes(y=duplicate_frequency), width=120) + 
              ggplot2::geom_vline(xintercept = 0, color="red", size=1) + 
              ggplot2::scale_y_continuous(name = "Frequency",
                                          sec.axis = sec_axis( ~ ., name = "Clonality"),
                                          limits = c(0, 0.5)) +
              ggplot2::theme_classic(base_size=12) + 
              ggplot2::scale_x_continuous(limits = date_lim,
                                          breaks = seq(date_lim[1], date_lim[2], 1000)) +
              ggplot2::guides(fill = ggplot2::guide_legend(nrow=length(pub_lab), byrow=TRUE)) +
              ggplot2::scale_fill_manual(values = tot_pal, 
                                         name="Public TCR", 
                                         drop=FALSE, 
                                         breaks=pub_seq, 
                                         labels=pub_lab) +
              ggplot2::xlab("Date") +
              ggplot2::geom_point(aes(y=clonality), size=2, color="purple", fill="purple") +
              ggplot2::ylab("Frequency") +
              ggplot2::theme(legend.position="none",
                             axis.line.y = ggplot2::element_line(colour = "black", size = 1),
                             axis.text.x = ggplot2::element_text(size=10),
                             axis.title.x = ggplot2::element_text(size=12),
                             axis.title.y = ggplot2::element_text(size=12),
                             axis.text.y = ggplot2::element_text(size=10),
                             legend.text = ggplot2::element_text(size=10, face="bold"),
                             legend.title = ggplot2::element_text(size=12, face="bold"),
                             legend.key.size = unit(.2,"cm")) 

    # Generate connected line graph to show Viral load calculations
    vmax <- max(vtable$result, na.rm=TRUE) + 10
    vline <- ggplot2::ggplot(vtable, aes(x= days_rel)) +
             ggplot2::geom_line(aes(y=result), size=1, color="black") + 
             ggplot2::geom_point(aes(y=result),size=2, color="black") + 
             ggplot2::geom_vline(xintercept = 0, color="red", size=1) + 
             ggplot2::theme_classic(base_size = 12) + 
             ggplot2::scale_y_continuous(name = "Viral Load\n(copies/mL)",
                                         labels = number_format(),
                                         trans = "log10",
                                         limits = c(1, vmax),
                                         sec.axis = sec_axis(~.)) + 
             ggplot2::scale_x_continuous(limits = date_lim,
                                         breaks = seq(date_lim[1], date_lim[2], 1000)) +
             ggplot2::ylab("Viral Load") +
             ggplot2::theme(legend.position = "none",
                            axis.line.x= ggplot2::element_blank(),
                            axis.text.y.right = ggplot2::element_blank(),
                            axis.ticks.y.right = ggplot2::element_blank(),
                            axis.line.y = ggplot2::element_line(colour = "black", size = 1),
                            axis.ticks.x = ggplot2::element_blank(), 
                            axis.text.x = ggplot2::element_blank(), 
                            axis.title.x = ggplot2::element_blank(),
                            axis.title.y = ggplot2::element_text(size=12),
                            axis.text.y = ggplot2::element_text(size=10))

    # Generate connected line graph to show CD4 count
    cmax <- max(cdtable$result, na.rm=TRUE) + 10
    cline <- ggplot2::ggplot(cdtable, aes(x= days_rel))+
             ggplot2::geom_line( aes(y=result), size=1, color="grey") +
             ggplot2::geom_point(aes(y=result),size=2, color="grey") +
             ggplot2::geom_vline(xintercept = 0, color="red", size=1) + 
             ggplot2::theme_classic(base_size = 12) + 
             ggplot2::scale_y_continuous(name = "CD4 count\n(count/uL)",
                                         labels = number_format(),
                                         limits = c(0, 1200),
                                         sec.axis = sec_axis(~.)) + 
             ggplot2::scale_x_continuous(limits = date_lim,
                                         breaks = seq(date_lim[1], date_lim[2], 1000)) +
             ggplot2::ylab("CD4 Count") +
             ggplot2::theme(legend.position = "none",
                            axis.line.x = ggplot2::element_blank(),
                            axis.text.y.right = ggplot2::element_blank(),
                            axis.ticks.y.right = ggplot2::element_blank(),
                            axis.line.y = ggplot2::element_line(colour = "black", size = 1),
                            axis.ticks.x = ggplot2::element_blank(), 
                            axis.text.x = ggplot2::element_blank(), 
                            axis.title.x = ggplot2::element_blank(),
                            axis.title.y = ggplot2::element_text(size=12),
                            axis.text.y = ggplot2::element_text(size=10))
   
    # Generate point plot to display clinical events
    dlist <- dtable %>% 
             dplyr::filter(!is.na(days_rel)) %>%
             dplyr::mutate(disease = cval[Condition]) %>%
             dplyr::pull(disease) %>% 
             unique()
    dline <- dtable %>%
             dplyr::mutate(disease = cval[Condition]) %>%
             ggplot2::ggplot( aes(y= 0, x= days_rel)) +
             ggplot2::geom_point(aes( color=disease), size = 4) +
             ggplot2::scale_x_continuous(limits = date_lim,
                                         breaks = seq(date_lim[1], date_lim[2], 1000)) +
             ggplot2::scale_y_continuous(limits = c(0,0), 
                                         breaks = 0,
                                         sec.axis = sec_axis(~.)) +
             ggplot2::ylab("Clinical\nEvents") + 
             ggplot2::labs(color = "Condition") +
             ggplot2::theme_classic(base_size = 12) + 
             ggplot2::scale_color_manual(values=cpal,
                                         breaks = dlist) +
             ggplot2::theme(legend.position="top", 
                            axis.ticks.y = ggplot2::element_blank(), 
                            axis.line.x = ggplot2::element_blank(),
                            axis.text.y.right = ggplot2::element_blank(),
                            axis.ticks.y.right = ggplot2::element_blank(),
                            axis.line.y = ggplot2::element_line(colour = "black", size = 1),
                            axis.text.y = ggplot2::element_blank(),
                            axis.title.x = ggplot2::element_blank(),
                            axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank(),
                            legend.text = ggplot2::element_text(size=10))
    # Plot the top sequences
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
    tmax <- max(top_table$cummulative_freq, na.rm=TRUE) + .1
    top_plot <- top_table %>%
                dplyr::mutate(top = forcats::as_factor(top)) %>% 
                ggplot2::ggplot(aes(x = yearsRel, y = cummulative_freq, color = top, group = top)) +
                ggplot2::geom_point(size = 2) +
                ggplot2::geom_line(size = 1) + 
                ggplot2::labs(color = "Top n sequences",
                              group = "Top n sequences") +
                ggplot2::geom_vline(xintercept = 0, color="red", size=1) + 
                ggplot2::scale_color_grey() + 
                ggplot2::scale_x_continuous(limits = date_lim,
                                         breaks = seq(date_lim[1], date_lim[2], 1000)) +
                ggplot2::scale_y_continuous(name = "Cumulative\nfrequency",
                                         limits = c(0, 1),
                                         sec.axis = sec_axis(~.)) + 
                ggplot2::theme_classic(base_size = 12) +
                ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                            axis.text.y.right = ggplot2::element_blank(),
                            axis.ticks.y.right = ggplot2::element_blank(),
                            axis.line.y = ggplot2::element_line(colour = "black", size = 1),
                            axis.ticks.x = ggplot2::element_blank(), 
                            axis.text.x = ggplot2::element_blank(), 
                            axis.title.x = ggplot2::element_blank(),
                            axis.title.y = ggplot2::element_text(size=12),
                            axis.text.y = ggplot2::element_text(size=10),
                            legend.position = c(0.2, .8))
    # Generate final plot layout using Patchwork
    layout <- "
    AAAAAAAAA
    BBBBBBBBB
    CCCCCCCCC
    EEEEEEEEE
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    "
    final_plot <- dline + vline + cline + sankey + top_plot + 
                  plot_layout(design = layout) + 
                  plot_annotation(title = patient_name)
}

# Load public TRB data
public_list <- list.files(parser$ppath, full.names = TRUE, all.files = FALSE, 
                          recursive = FALSE, pattern = ".csv", 
                          include.dirs = FALSE)
pbtable <- public_list %>% 
           purrr::map(readPublic) %>% 
           dplyr::bind_rows() %>% 
           dplyr::rename(junction_aa = aminoAcid)
mcppath <- paste(parser$ppath, "publicDB/McPAS-TCR.csv", sep = "/") 
mptable <- readr::read_csv(mcppath) %>% 
           dplyr::rename(junction_aa = CDR3.beta.aa, antigen = Pathology) %>% 
           dplyr::select(antigen, junction_aa)

antigen_table <- c(EBV = "Epstein Barr virus (EBV)", 
                   CMV = "Cytomegalovirus (CMV)",
                   flu = "Influenza",
                   HSV = "Herpes simplex virus 1 (HSV1)",
                   `Human herpes virus 1` = "Herpes simplex virus 1 (HSV1)",
                   `Hepatitis C virus` = "Hepatitis C virus (HCV)",
                   HIV = "Human immunodeficiency virus (HIV)",
                   KS = "Kaposi Sarcoma",
                   `Psoriatic arthritis` = "Psoriatic Arthritis",
                   `HIV-1` = "Human immunodeficiency virus (HIV)",
                   InfluenzaA = "Influenza",
                   HPV = "Human papilloma virus",
                   HCV = "Hepatitis C virus (HCV)",
                   `HSV-2` = "Herpes simplex virus 2 (HSV2)",
                   LCMV = "Lymphocytic choriomeningitis virus (LCMV)",
                   YFV = "Yellow fever virus")

filter_species <- c("MusMusculus", "HomoSapiens", "GallusGallus")
vdjpath <- paste(parser$ppath, "publicDB/vdjdb_full.txt", sep = "/") 
vdjtable <- readr::read_delim(vdjpath, , delim="\t" ) %>%
            dplyr::select(cdr3.beta, antigen.species) %>%
            dplyr::rename(junction_aa = cdr3.beta, antigen = antigen.species) %>%
            dplyr::filter(!antigen %in% filter_species)

pbtable <- dplyr::bind_rows(pbtable, mptable) %>%
           dplyr::bind_rows(vdjtable) %>%
           dplyr::mutate(antigen = recode_factor(antigen, !!!antigen_table)) %>%
           dplyr::distinct_all() %>%
           dplyr::group_by(junction_aa) %>%
           dplyr::summarize(antigen = paste(antigen, collapse = ",")) %>%
           dplyr::ungroup()

pbnames <- pbtable$junction_aa
pbfactor <- pbtable$antigen
names(pbfactor) <- pbnames

atable <- atable %>%
          dplyr::mutate(yearsRel = (dateART %--% dateTCR) %/% lubridate::days(1),
                        time_point = dplyr::if_else(as.integer(yearsRel) < 0, "Pre-ART", "Post-ART"))
p1 <- colorRampPalette(c("#7570b3", "#cab2d6"))
p2 <- colorRampPalette(c("#d95f02", "#fdbf6f"))
p3 <- colorRampPalette(c("#737373", "#525252"))

sankey_plot <- plotSankey(atable, vl_table, cc_table, infection_table, itable, pbfactor, parser$pid)
out_path <- paste(parser$pid, "topFreq.pdf", sep = "_")
ggsave(out_path, sankey_plot, width = 8, height = 12, device = "pdf", units = "in")