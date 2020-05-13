#### LOAD LIBRARIES #### 
pacman::p_load(shiny, shinythemes, dplyr, plotly, conflicted, tidyr, ggplot2)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT)

# set conflict preference
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")

setwd("/Users/biplabendudas/Documents/GitHub/Shiny_apps/TC5_zplots/app")

## Load results database -------------
# load the database for TC5
db <- dbConnect(RSQLite::SQLite(), "/Users/biplabendudas/Documents/GitHub/R-scripts_zombie_ant_lab/RSQLite/sql_dbs/TC5_data.db")
# src_dbi(db)

## Load blast database ---------------
blast_db <- dbConnect(RSQLite::SQLite(),"/Users/biplabendudas/Documents/GitHub/R-scripts_zombie_ant_lab/RSQLite/sql_dbs/blast_data.db")
# src_dbi(blast_db)



# Your function here ------------------------------------------------------
### Data
## load core dataset
cflo.annots.exp <- tbl(db, "annot_fpkm") %>% collect()

## load zscore datasets
cflo.zscores.for <- tbl(db, "zscore_for") %>% collect()
cflo.zscores.nur <- tbl(db, "zscore_nur") %>% collect()

## load eJTK results
for.ejtk <- tbl(db, "ejtk_all") %>% filter(caste=="for") %>% collect()
nur.ejtk <- tbl(db, "ejtk_all") %>% filter(caste=="nur") %>% collect()

## load DEG results
cflo_Ian_cuffdiff <- tbl(db, "ophio_all_DEG") %>% collect()

## sig DEGs for biting vs control
sig.degs.biting.control <- tbl(db, "ophio_biting_control") %>% collect() %>% pull(gene)


### Functions
# source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/plot_zscores.R")
# source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/vennDia.R")
# source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/enrichment_analysis.R")
source("./data/Functions/theme_publication.R")


# List of expressed genes -------------------------------------------------

dat.f <- cflo.annots.exp %>% 
  dplyr::select(gene_name, X2F:X24F) %>% 
  # omit the genes that have NAs
  na.omit() %>% 
  # remove genes that have at least 0/0.5/1 fpkm for at least one time point
  filter_at(vars(starts_with("X")), any_vars(. >= 1)) 


dat.n <- cflo.annots.exp %>% 
  dplyr::select(gene_name, X2N:X24N) %>% 
  # omit the genes that have NAs
  na.omit() %>% 
  # remove genes that have at least 0.5 fpkm for at least one time point
  filter_at(vars(starts_with("X")), any_vars(. >= 1))  

# cflo.annots.exp %>% 
#   select(gene_name, X2N:X24N) %>% 
#   # omit the genes that have NAs
#   na.omit() %>% 
#   # remove genes that have at least 0.5 fpkm for at least one time point
#   filter_at(vars(starts_with("X")), any_vars(. > 1)) %>% 
#   nrow()
# 
# venndiagram(x= cflo.annots.exp %>% 
#               select(gene_name, X2F:X24F) %>% 
#               # omit the genes that have NAs
#               na.omit() %>% 
#               # remove genes that have at least 0.5 fpkm for at least one time point
#               filter_at(vars(starts_with("X")), any_vars(. > 0)) %>% pull(gene_name),
#             y= cflo.annots.exp %>% 
#               select(gene_name, X2N:X24N) %>% 
#               # omit the genes that have NAs
#               na.omit() %>% 
#               # remove genes that have at least 0.5 fpkm for at least one time point
#               filter_at(vars(starts_with("X")), any_vars(. > 0)) %>% pull(gene_name),  
#             unique=T, title="", 
#             labels=c("Foragers", "Nurses"), 
#             lines=1:2, lcol=1:2, tcol=c(3,1,2), diacol=1, plot=T, type="2")



cntNonZero.f<- apply(dat.f[-1], 1, function(x) sum(x >= 1))
cntNonZero.n<- apply(dat.n[-1], 1, function(x) sum(x >= 1))

# subset the data and only keep the filtered genes
gs1 <- dat.f[which(cntNonZero.f >=6),] %>% pull(gene_name)
gs2 <- dat.n[which(cntNonZero.n >=6),] %>% pull(gene_name)


# set background genes
bg.genes = unique(c(dat.f %>% pull(gene_name),
                    dat.n %>% pull(gene_name)))


# Function to plot z-scores
z.plot <- function(gene_names, caste = "both", lwd=1.5, alpha=0.9) {
  # load the required libraries
  # library(tidyverse)
  # library(gridExtra)
  # library(ggplot2)
  # load the core datasets
  # load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
  # load(file = "./data/TC5_core_datasets.RData")
  
  #save the current directory
  current.dir <- getwd()
  #change the directory to the folder where we have our zscores
  #setwd("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data")
  # save the list of gene names in a character vector
  g <- gene_names
  # Read the zscores for caste
  # Read the zscore file
  #foragers_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_foragers.csv", header = T, stringsAsFactors = F)
  #nurses_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_nurses.csv", header = T, stringsAsFactors = F)
  foragers_zscores <- cflo.zscores.for
  nurses_zscores <- cflo.zscores.nur

  # Read the annotation file for description of the genes
  #all_genes <- read.csv("~/R-scripts_zombie_ant_lab/Enrichment_analysis/Ants/cflo_annotations_expression_v2.csv", header = T, stringsAsFactors = FALSE)
  all_genes <- cflo.annots.exp
  all_genes_annots <- all_genes %>%
    dplyr::select(gene_name, annot = old_annotation)

  # Transforming the data to be able to use ggplot
  # Let's try the logic on a dummy subset
  dummy.for <- foragers_zscores %>%
    filter(gene_name %in% g) %>%
    gather(ZT, zscore, -1) %>%
    arrange(gene_name) %>%
    mutate(caste = "for") %>%
    mutate(ZT = readr::parse_number(ZT)) %>%
    left_join(all_genes_annots, by = "gene_name")

  dummy.nur <- nurses_zscores %>%
    filter(gene_name %in% g) %>%
    gather(ZT, zscore, -1) %>%
    arrange(gene_name) %>%
    mutate(caste = "nur") %>%
    mutate(ZT = readr::parse_number(ZT)) %>%
    left_join(all_genes_annots, by = "gene_name")

  # make an if else statement to modify the dataset as per the "caste" parameter
  if (caste == "for" | caste == "forager" | caste == "foragers" | caste == "f"){
    dummy <- dummy.for
    col.scheme <- c("#F23030")
  }
  else if (caste == "nur" | caste == "nurse" | caste == "nurses" | caste == "n") {
    dummy <- dummy.nur
    col.scheme <- c("#1A80D9")
  } else if (caste == "all" | caste == "both" | caste == "a" | caste == "b") {
    dummy <- rbind(dummy.for, dummy.nur)
    col.scheme <- c("#F23030","#1A80D9")
  } else {
    print("Invalid value for caste. Use one of the following options: for, nur, all.")
    # stop the function and give the user an error.
    stop();
  }

  # make the gene_name and caste columns as factors
  dummy[[1]] <- as.factor(dummy[[1]])
  dummy[[4]] <- as.factor(dummy[[4]])
  dummy[[5]] <- as.factor(dummy[[5]])

  # # Initialize a list to save the plots
  # l <- list()
  
  ZT <- as.numeric(as.character(dummy[dummy$gene_name==g,]$ZT)) 
  Z.score <- round(dummy[dummy$gene_name==g,]$zscore,2)
  
  # Let's plot
  library(ggplot2)
  pd <- position_dodge(0.2)
  l <- 
    # lapply(sort(unique(dummy[[1]])), function(i) {
    ggplot(dummy[dummy$gene_name==g,], aes(x=ZT, y=Z.score)) +
      #geom_errorbar(aes(ymin=Corrected_Exp-se, ymax=Corrected_Exp+se), width=.1, position=pd) +
      # geom_hline(yintercept=0, color = "red", size = 2) +
      ## if you need highlighting parts of the graph (dark phase in my case)
      # geom_rect(aes(xmin = 1, xmax = 23.5, ymin = -Inf, ymax = Inf),
      #           fill = "lightgrey", alpha = 0.8, color=NA) +
      #facet_grid(~gene_name, scales = "free") +
      facet_wrap(~ gene_name) +
      #ggtitle("TC5 - Z-scores") +
      xlab("") +
      ylab("") +
      theme_Publication() +
      scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
      scale_y_continuous(limits = c(min(dummy[[3]])-0.5,max(dummy[[3]])+0.5)) +
      geom_line(data = dummy[dummy$gene_name==g,], position=pd,
                aes(col=caste), size=lwd, alpha=alpha) +
    
      # # let's add points:
      # geom_point(size=2.5, color="black", alpha=0.6) +
      # geom_point(position=pd,
      #            aes(fill=as.factor(caste), shape = as.factor(caste)), size=2.5, show.legend = F,
      #            color="black", pch=21) +
      # 
      scale_fill_manual(values =col.scheme) +
      scale_color_manual(values=col.scheme) +
      theme(text = element_text(size = 20, colour = 'black'),
            legend.position = "none"
            # axis.title.x=element_blank(),
            # axis.text.x=element_blank()
            ) +
      # set transparency
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA)) +
      theme(axis.text.y=element_text(color=c("transparent","black","transparent","transparent",
                                             "black","transparent","transparent","transparent",
                                             "transparent","black")))

  # # Let's name the plots with their resp. gene names
  # names(l) <- sort(g);
  # # let's save the data frame with annotation of the genes into the list as well
  # annotations <- all_genes_annots %>%
  #   filter(gene_name %in% g)
  # 
  # # save the list of plots (l) and the annotation table (annotations) into a list
  # l2 <- list()
  # l2$plots <- l
  # l2$annots <- annotations

  # set the working directory to the pre-existing one
  setwd(current.dir)

  # return the list with all the plots and the data frame containing
  return(l);

}



# log.plot ----------------------------------------------------------------

log.plot <- function(gene_names, caste = "both", log=T, lwd=2, alpha=0.9) {
  # # load the required libraries
  # library(tidyverse)
  # library(gridExtra)
  # library(ggplot2)
  # load the core datasets
  # load(file = "./data/TC5_core_datasets.RData")
  #save the current directory 
  current.dir <- getwd()
  #change the directory to the folder where we have our zscores
  #setwd("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data")
  # save the list of gene names in a character vector
  g <- gene_names
  # Read the zscores for caste
  # Read the zscore file
  #foragers_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_foragers.csv", header = T, stringsAsFactors = F)
  #nurses_zscores <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/Zscore_data/TC5_zscores_nurses.csv", header = T, stringsAsFactors = F)
  foragers <- cflo.annots.exp %>% dplyr::select(gene_name, X2F:X24F)
  nurses <- cflo.annots.exp %>% dplyr::select(gene_name, X2N:X24N)
  
  # Read the annotation file for description of the genes
  #all_genes <- read.csv("~/R-scripts_zombie_ant_lab/Enrichment_analysis/Ants/cflo_annotations_expression_v2.csv", header = T, stringsAsFactors = FALSE)
  all_genes <- cflo.annots.exp
  all_genes_annots <- all_genes %>% 
    dplyr::select(gene_name, annot = old_annotation)
  
  # Transforming the data to be able to use ggplot
  # Let's try the logic on a dummy subset
  
  dummy.for <- foragers %>%
    filter(gene_name %in% g) %>% 
    gather(ZT, exp, -1) %>% 
    # mutate(log.exp = log10(exp+1)) %>% 
    # dplyr::select(-exp) %>% 
    arrange(gene_name) %>% 
    mutate(caste = "for") %>% 
    mutate(ZT = readr::parse_number(ZT)) %>% 
    left_join(all_genes_annots, by = "gene_name")
  
  dummy.nur <- nurses %>%
    filter(gene_name %in% g) %>% 
    gather(ZT, exp, -1) %>% 
    # mutate(log.exp = log10(exp+1)) %>% 
    # dplyr::select(-exp) %>% 
    arrange(gene_name) %>% 
    mutate(caste = "nur") %>% 
    mutate(ZT = readr::parse_number(ZT)) %>% 
    left_join(all_genes_annots, by = "gene_name")
  
  # make an if else statement to modify the dataset as per the "caste" parameter
  if (caste == "for" | caste == "forager" | caste == "foragers" | caste == "f"){
    dummy <- dummy.for
    col.scheme <- c("#F23030")
  } 
  else if (caste == "nur" | caste == "nurse" | caste == "nurses" | caste == "n") {
    dummy <- dummy.nur
    col.scheme <- c("#1A80D9")
  } else if (caste == "all" | caste == "both" | caste == "a" | caste == "b") {
    dummy <- rbind(dummy.for, dummy.nur)
    col.scheme <- c("#F23030","#1A80D9")
  } else {
    print("Invalid value for caste. Use one of the following options: for, nur, all.")
    # stop the function and give the user an error.
    stop();
  }
  
  # make the gene_name and caste columns as factors
  dummy[[1]] <- as.factor(dummy[[1]])
  dummy[[4]] <- as.factor(dummy[[4]])
  dummy[[5]] <- as.factor(dummy[[5]])
  
  # Initialize a list to save the plots
  # l <- list()
  # Let's plot
  # library(ggplot2)
  pd <- position_dodge(0.2)
  
  zt <- as.numeric(as.character(dummy[dummy$gene_name==g,]$ZT)) 
  
  if(log==T) {
    
    log2.exp <- round(log2(dummy[dummy$gene_name==g,]$exp+1),2)
    
    l <- 
      # lapply(sort(unique(dummy[[1]])), function(i) {
      ggplot(dummy[dummy$gene_name==g,], aes(x=zt, y=log2.exp)) + 
        #geom_errorbar(aes(ymin=Corrected_Exp-se, ymax=Corrected_Exp+se), width=.1, position=pd) +
        # geom_hline(yintercept=0, color = "red", size = 2) +
        ## if you need highlighting parts of the graph (dark phase in my case)
        geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),    
                  fill = "lightgrey", alpha = 0.02, color=NA) +
        #facet_grid(~gene_name, scales = "free") +
        facet_wrap(~ gene_name, scales = "free_y") +
        #ggtitle("TC5 - Z-scores") +
        xlab("") +
        ylab("") +
        theme_Publication() +
        scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
        # scale_y_continuous(limits = c(min(dummy[[3]]),max(dummy[[3]]))) +  
        geom_line(data = dummy[dummy$gene_name==g,], position=pd, 
                  aes(col=caste), size=lwd, alpha=alpha) +
        # geom_point(position=pd, 
        #            aes(fill=as.factor(caste), shape = as.factor(caste)), size=5, show.legend = F,
        #            color="black", pch=21) +
        #scale_fill_manual(values = c("#F2CB05","#0FBF67")) +
      scale_fill_manual(values =col.scheme) +
      scale_color_manual(values=col.scheme) +
      theme(text = element_text(size = 20, colour = 'black'),
            legend.position = "none"
            # axis.title.x=element_blank(),
            # axis.text.x=element_blank()
      ) +
      # set transparency
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA)) +
      theme(axis.text.y=element_text(color=c("transparent","black","transparent","transparent",
                                             "black","transparent","transparent","transparent",
                                             "transparent","black")))
    
    }
  
  else if (log==F) {
    
    expn = round(dummy[dummy$gene_name==g,]$exp,2)
    
    l <- 
      # lapply(sort(unique(dummy[[1]])), function(i) {
      ggplot(dummy[dummy$gene_name==g,], aes(x=zt, y=expn)) + 
        #geom_errorbar(aes(ymin=Corrected_Exp-se, ymax=Corrected_Exp+se), width=.1, position=pd) +
        # geom_hline(yintercept=0, color = "red", size = 2) +
        ## if you need highlighting parts of the graph (dark phase in my case)
        geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),    
                  fill = "lightgrey", alpha = 0.02, color=NA) +
        #facet_grid(~gene_name, scales = "free") +
        facet_wrap(~ gene_name, scales = "free_y") +
        #ggtitle("TC5 - Z-scores") +
        xlab("") +
        ylab("") +
        theme_Publication() +
        scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
        # scale_y_continuous(limits = c(min(dummy[[3]]),max(dummy[[3]]))) +  
        geom_line(data = dummy[dummy$gene_name==g,], position=pd, 
                  aes(col=caste), size=lwd, alpha=alpha) +
        # geom_point(position=pd, 
        #            aes(fill=as.factor(caste), shape = as.factor(caste)), size=5, show.legend = F,
        #            color="black", pch=21) +
        #scale_fill_manual(values = c("#F2CB05","#0FBF67")) +
      scale_fill_manual(values =col.scheme) +
      scale_color_manual(values=col.scheme) +
      theme(text = element_text(size = 20, colour = 'black'),
            legend.position = "none"
            # axis.title.x=element_blank(),
            # axis.text.x=element_blank()
      ) +
      # set transparency
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA)) +
      theme(axis.text.y=element_text(color=c("transparent","black","transparent","transparent",
                                             "black","transparent","transparent","transparent",
                                             "transparent","black")))
    
  }
  
  # # Let's name the plots with their resp. gene names
  # names(l) <- sort(g);
  # # let's save the data frame with annotation of the genes into the list as well
  # annotations <- all_genes_annots %>% 
  #   filter(gene_name %in% g)
  # 
  # # save the list of plots (l) and the annotation table (annotations) into a list
  # l2 <- list()
  # l2$plots <- l
  # l2$annots <- annotations
  
  # set the working directory to the pre-existing one
  setwd(current.dir)
  
  # return the list with all the plots and the data frame containing   
  return(l);
  
}


# Annotations to display --------------------------------------------------

annot <- function(gene) {
  # load the core datasets
  # load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
  # load(file = "./data/TC5_core_datasets.RData")
  
  df <- cflo.annots.exp %>% 
    dplyr::select(gene_name, Blast = old_annotation, GO = GOs, PFAM = pfams, signalP, TMHMM) %>% 
    filter(gene_name %in% gene)
  
  df
  
}



# rhythmicity ------------------------------

rhy.results <- function(gene) {
  
  ultra.8h.df <- tbl(db, "ejtk_8h_all") %>% 
    filter(gene_name %in% gene) %>% 
    select(gene_name, caste, GammaP_8h = GammaP)
  
  ultra.12h.df <- tbl(db, "ejtk_12h_all") %>% 
    filter(gene_name %in% gene) %>% 
    select(gene_name, caste, GammaP_12h = GammaP)
  
  df <- tbl(db, "ejtk_all") %>% 
    filter(gene_name %in% gene) %>% 
    select(gene_name, caste:GammaP) %>% 
    left_join(ultra.8h.df, by=c("gene_name","caste")) %>% 
    left_join(ultra.12h.df, by=c("gene_name","caste")) %>% 
    select(-1)
    
  df
}



# Display DGE information -------------------------------------------------

cflo.cuffdiff <- 
  cflo_Ian_cuffdiff %>% 
  mutate(treatment_1 = ifelse(sample_1 == "Alive", "Biting", sample_1),
         treatment_2 = ifelse(sample_2 == "Alive", "Biting", sample_2)) %>% 
  select(gene_name = gene, 
         # locus,
         annotation,
         treatment_1, treatment_2,
         significant,
         value_1, 
         value_2,
         logFC = log2.fold_change.,
         test_stat, 
         p_value,
         q_value) %>% 
  arrange(gene_name, treatment_1, treatment_2) %>% 
  mutate_at(c(7), as.numeric) %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% 
  mutate(DEG_biting = ifelse(gene_name %in% sig.degs.biting.control, "yes", "no"))


cuffdiff.results <- function(gene) {
  df <- cflo.cuffdiff %>% 
    filter(gene_name %in% gene) %>% 
    select(-1)
  df
}


