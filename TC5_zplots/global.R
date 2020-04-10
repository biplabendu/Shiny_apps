#### LOAD LIBRARIES #### 
pacman::p_load(shiny, shinythemes, deSolve, dplyr, plotly)


#### DEFINE CUSTOM FUNCTIONS ####
#' @title ODE System
#' 
#' @description This function represents the ODE system to simulate the 
#'  collective selection of a food source by ant colonies.  
#' 
#' @param t The time step to compute the ODE at.
#' @param y The values of the state variables of the ODE at the time step: 
#'  number of ants at the nest (N) and at each food source (S1...Sn), quantity 
#'  of pheromone on the routes to each food source (Q1...Qn).
#' @param parms Parameters of the ODE system:
#'  \itemize{
#'   \item l: Distance to each food source in seconds.
#'   \item q: Quantity of pheromone deposited by the ants on the way back from 
#'    each food source.
#'   \item qe: Quantity of pheromone deposited by the ants on the way to the 
#'    food sources.
#'   \item alpha, beta, gamma, eta: Regulates recruitment at the nest as a 
#'    function of the quantity of pheromone at the nest entrance.
#'   \item psi: Probability for an ant to leave the food source.
#'   \item rho: Rate of evaporation of the pheromone.
#'   \item k: Intrinsic attractiveness of each route.
#'   \item n: Degree of nonlinearity.
#'  }
#' 
#' @return A list of the changes in \code{y} suitable for \code{\link{dede}} 
#'  function. 
#' 
#' @author Simon Garnier: \email{garnier@@njit.edu}, 
#'  \link[https://twitter.com/sjmgarnier]{@@sjmgarnier}
#' 
# ode_sys <- function(t, y, parms) {
#   nS <- length(parms$l)
#   
#   Q <- sum(y[(nS + 2):(2 * nS + 1)])
#   phi <- (parms$alpha * (parms$beta + Q) ^ parms$eta) / ((parms$beta + Q) ^ parms$eta + parms$gamma)
#   
#   dN <- -phi*y[1]
#   dS <- rep(0,nS)
#   dQ <- rep(0,nS)
#   
#   ylag <- {}
#   
#   for (i in 1:nS) {
#     if (t < parms$l[i]) {
#       ylag <- c(parms$iniN, parms$iniS, parms$iniQ)
#       Qlag <- 0
#       philag <- 0
#     } else {
#       ylag <- lagvalue(t - parms$l[i])
#       Qlag <- sum(ylag[(nS + 2):(2 * nS + 1)])
#       philag <- (parms$alpha * (parms$beta + Qlag) ^ parms$eta) / ((parms$beta + Qlag) ^ parms$eta + parms$gamma)
#     }
#     
#     P <- ((parms$k[i] + y[i + nS + 1]) ^ parms$n) / sum((parms$k + y[(nS + 2):(2 * nS + 1)]) ^ parms$n)
#     Plag <- ((parms$k[i] + ylag[i + nS + 1]) ^ parms$n) / sum((parms$k + ylag[(nS + 2):(2 * nS + 1)]) ^ parms$n)
#     
#     dN <- dN + parms$psi*ylag[i + 1]
#     
#     dS[i] <- -parms$psi*y[i + 1] + philag * ylag[1] * Plag
#     
#     dQ[i] <- parms$qe * P * phi * y[1] + parms$q[i] * parms$psi * ylag[i + 1] - parms$rho * y[i + nS + 1]
#     
#   }
#   
#   return(list(c(dN, dS, dQ)))
# }


# Your function here ------------------------------------------------------

load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
# my functions
# source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/plot_zscores.R")
source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/vennDia.R")
source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/enrichment_analysis.R")
source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/theme_publication.R")


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
  library(tidyverse)
  library(gridExtra)
  library(ggplot2)
  # load the core datasets
  load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
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


# Annotations to display --------------------------------------------------

annot <- function(gene) {
  # load the core datasets
  load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
  
  df <- cflo.annots.exp %>% 
    dplyr::select(gene_name, Blast = old_annotation, GO = GOs, PFAM = pfams, signalP, TMHMM) %>% 
    filter(gene_name %in% gene)
  
  df
  
}



# Results from eJTK -------------------------------------------------------

# eJTK - Run 1 (all genes with non-zero expression) -------------------------------------------------

for.run1 <- read.table("~/OneDrive - University of Central Florida/BD-TC5/1212LD_24h/eJTK/Results/eJTK_Run1/output/for/tc5_for_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt", header = T)
nur.run1 <- read.table("~/OneDrive - University of Central Florida/BD-TC5/1212LD_24h/eJTK/Results/eJTK_Run1/output/nur/tc5_nur_cos24_ph00-22_by2_a02-22_by2_OTHERTEXT_jtkout_GammaP.txt", header = T)

# identify rhythmicity ------------------------------
pbh.fdr = 0.1

for.pbh <- for.run1 %>% 
  ## keep the "expressed" genes
  filter(ID %in% gs1) %>% 
  ## Do p.adjust
  mutate(empP_BH = p.adjust(empP, method = "BH")) %>% 
  mutate(caste = "for") %>% 
  dplyr::select(ID, caste, empP, empP_BH, GammaBH, Phase, Mean, Std_Dev, Max_Amp)

nur.pbh <- nur.run1 %>% 
  ## keep the "expressed" genes
  filter(ID %in% gs2) %>% 
  ## Do p.adjust
  mutate(empP_BH = p.adjust(empP, method = "BH")) %>% 
  mutate(caste = "nur") %>% 
  dplyr::select(ID, caste, empP, empP_BH, GammaBH, Phase, Mean, Std_Dev, Max_Amp)

for.nur.pbh <- rbind(for.pbh, nur.pbh)

rhy.results <- function(gene) {
  df <- for.nur.pbh %>% 
    filter(ID %in% gene) %>% 
    dplyr::select(-1) %>% 
    mutate("FDR_10%" = ifelse(empP_BH < 0.1, "sig. rhy", "ns"))
  df
}



# # COVID_19 code -----------------------------------------------------------
# 
# library(dplyr)
# library(tidyr)
# f1 = list(family="Courier New, monospace", size=12, color="rgb(30,30,30)")
# 
# # Read the curated dataset from my GitHub account
# 
# urlfile <- "https://github.com/biplabendu/homepage/raw/master/covid19_data_india/curated_data/curated_data_BD.csv"
# 
# 
# daysSinceLastUpdate = function(fileName) {
#   (as.numeric(as.POSIXlt(Sys.time())) - as.numeric(file.info(fileName)$ctime))/60/60/24
# }
# 
# loadData = function(fileName) {
#   if(!file.exists(fileName) || daysSinceLastUpdate(fileName) > 0.5) {
#     data = read.table(url(urlfile), sep = ",",
#                       header = T, stringsAsFactors=FALSE) %>%
#       mutate(date=as.Date(date))
#     # select(-Lat, -Long) %>%
#     # pivot_longer(-(1:2), names_to="date", values_to=columnName) %>%
#     # mutate(
#     #   date=as.Date(date, format="%m/%d/%y"),
#     #   `Country/Region`=if_else(`Country/Region` == "", "?", `Country/Region`),
#     #   `Province/State`=if_else(`Province/State` == "", "<all>", `Province/State`))
# 
#     save(data, file=fileName)
#   } else {
#     load(file=fileName)
#   }
#   return(data)
# }
# 
# # allData =
# #   loadData("time_series_covid19_confirmed_global.csv", "CumConfirmed") %>%
# #   inner_join(loadData("time_series_covid19_deaths_global.csv", "CumDeaths")) %>%
# #   inner_join(loadData("time_series_covid19_recovered_global.csv", "CumRecovered"))
# #
# # df.curated.global <- allData
# #
# 
# df <- loadData("curated_data_BD.csv")
# 
# datetime <- paste0(file.info("curated_data_BD.csv")$ctime %>% as.character()," EDT")
# # df <- read.csv(url(urlfile),
# #                   header = T, stringsAsFactors = F)
# 
# df <- df[-1]
# names(df)[1] <- "State"
# names(df)[2] <- "Country"
# 
# allData <- df
# 
