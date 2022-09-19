#FINAL PROJECT
# written by Eden Meidan and Racheli Cohen

######################################
install.packages("ggplot2")
install.packages("airr")
install.packages("ape")
install.packages("graphics")
install.packages("grid")
install.packages("ape")
install.packages("igraph")
install.packages("Matrix")
install.packages("methods")
install.packages("progress")
install.packages("Rcpp")
install.packages("readr")
install.packages("rlang")
install.packages("scales")
install.packages("seqinr")
install.packages("stats")
install.packages("ape")
install.packages("stringi")
install.packages("tibble")
install.packages("tidyr")
install.packages("utils")
install.packages("IRanges")
install.packages("data.table")
install.packages("ggpubr")
install.packages("dplyr")
install.packages("stringr")
install.packages("foreach")
install.packages("BiocGenerics")
install.packages("GenomicAlignments")
install.packages("GenomicAlignments")
install.packages("ggplot2")
install.packages("alakazam")
install.packages("Biostrings")
install.packages("diptest")
install.packages("KernSmooth")
install.packages("doParallel")
install.packages("gtools")
install.packages("lazyeval")
#install.packages("shazam")
######################################
suppressMessages(library(seqinr))
suppressMessages(library(tools))
suppressMessages(library(tigger))
suppressMessages(library(optparse))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(alakazam))
suppressMessages(library(reshape2))
suppressMessages(library(shazam))
suppressMessages(library(rabhit))
suppressMessages(library(utils))
suppressMessages(library(splitstackshape))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(alakazam))
suppressMessages(library(parallel))
suppressMessages(library(foreach))    
suppressMessages(library(stringi))
library(ggplot2)
library(airr)
library(ape)
library(graphics)
library(grid)
library(ape)
library(igraph)
library(Matrix)
library(methods)
library(progress)
library(Rcpp)
library(readr)
library(rlang)
library(scales)
library(seqinr)
library(stats)
library(ape)
library(stringi)
library(tibble)
library(tidyr)
library(utils)
library(IRanges)
library(data.table)
library(ggpubr)
library(dplyr)
library(stringr)
library(foreach)
library(BiocGenerics)
library(GenomicAlignments)
library(alakazam)
library(Biostrings)
library(diptest)
library(KernSmooth)
library(doParallel)
library(gtools)
library(lazyeval)
library(shazam)
library(tigger)
library(wesanderson)
library(ComplexHeatmap)
setwd("C:/Users/user/OneDrive - Bar-Ilan University - Students/Documents/ביואינפורמטיקה/project gur yaari/heatmap")
source("functions_tigger_1.R")

### refrence #######################################################################################################################
germline <- readIgFasta("IGHV_gap_full.fasta")
germline_a <- names(germline)

### genotype seq data ##############################################################################################################
# get rid of multiple assignment, starts without ., no N, consensus_count = 2
data_geno_seq <- readRDS("data_tigger_p1_p11.rds")
print("loaded data_tigger_p1_p11.rds")
data_geno_seq  <- data_geno_seq[!grepl(",", data_geno_seq$v_call),]
data_geno_seq  <- data_geno_seq %>% filter(!grepl("^[.]", sequence_alignment), consensus_count>=2) %>% rowwise() %>% mutate(v_seq = substr(sequence_alignment, 1, 318)) %>% filter(!grepl("N", v_seq), !grepl("-", v_seq))
# add v_alleles and v_gene column
data_geno_seq <- data_geno_seq %>% mutate(v_alleles = getAllele(v_call, strip_d = F, first = T)) %>% mutate(v_gene = getGene(v_call, strip_d = F)) 
data_geno_seq <- filter(data_geno_seq, grepl("P1_", sample)) # filter out P11_
data_geno_seq <- filter(data_geno_seq, v_alleles %in% germline_a)
print("after clean up - data_tigger_p1_p11.rds")
#data_geno_seq<-read.csv("data_filter.csv")

### Upload genotype data   ##########################################################################################################
geno_tigger <- readRDS("genos_tigger_relative.rda")
geno_tigger <- rbindlist(geno_tigger, idcol = "subject")
geno_tigger$genotyped_alleles <- geno_tigger$alleles
geno_tigger$genotyped_imgt_alleles <- sapply(1:nrow(geno_tigger),function(i){
  alleles <- geno_tigger$genotyped_alleles[i]
  alleles <- strsplit(alleles,",")[[1]]
  paste0(geno_tigger$gene[i],"*",alleles, collapse = ",")
})
# filter out P11_
geno_tigger <- filter(geno_tigger, grepl("P1_", subject))
print("after clean up - genos_tigger_relative.rda")

## split data to 9 parts ##############################################################################
data_geno_seq1 <- filter(data_geno_seq, sample %in% c("P1_I1_S1","P1_I10_S1","P1_I100_S1","P1_I11_S1","P1_I12_S1","P1_I13_S1","P1_I14_S1","P1_I15_S1","P1_I16_S1","P1_I17_S1","P1_I18_S1"))
seq <- data_geno_seq1

data_geno_seq2 <- filter(data_geno_seq, sample %in% c("P1_I19_S1","P1_I2_S1","P1_I20_S1","P1_I21_S1","P1_I22_S1","P1_I23_S1", "P1_I24_S1","P1_I25_S1","P1_I26_S1","P1_I27_S1","P1_I28_S1"))
seq <- data_geno_seq2
#
data_geno_seq3 <- filter(data_geno_seq, sample %in% c("P1_I29_S1","P1_I3_S1","P1_I30_S1","P1_I31_S1","P1_I32_S1","P1_I33_S1","P1_I34_S1","P1_I35_S1","P1_I37_S1","P1_I38_S1","P1_I39_S1"))
seq <- data_geno_seq3

data_geno_seq4 <- filter(data_geno_seq, sample %in% c("P1_I4_S1","P1_I40_S1","P1_I41_S1","P1_I42_S1","P1_I43_S1","P1_I44_S1","P1_I45_S1","P1_I46_S1","P1_I47_S1","P1_I48_S1","P1_I49_S1"))
seq <- data_geno_seq4

data_geno_seq5 <- filter(data_geno_seq, sample %in% c("P1_I5_S1","P1_I50_S1","P1_I51_S1","P1_I52_S1","P1_I53_S1","P1_I54_S1","P1_I55_S1","P1_I56_S1","P1_I57_S1","P1_I58_S1","P1_I59_S1"))
seq <- data_geno_seq5

data_geno_seq6 <- filter(data_geno_seq, sample %in% c("P1_I6_S1","P1_I60_S1","P1_I61_S1","P1_I62_S1","P1_I63_S1","P1_I64_S1", "P1_I65_S1","P1_I66_S1","P1_I67_S1","P1_I68_S1","P1_I69_S1"))
seq <- data_geno_seq6

data_geno_seq7 <- filter(data_geno_seq, sample %in% c("P1_I7_S1","P1_I70_S1","P1_I71_S1","P1_I72_S1","P1_I73_S1","P1_I74_S1","P1_I75_S1","P1_I76_S1","P1_I77_S1","P1_I78_S1","P1_I79_S1" ))
seq <- data_geno_seq7

data_geno_seq8 <- filter(data_geno_seq, sample %in% c("P1_I8_S1","P1_I80_S1","P1_I81_S1","P1_I82_S1","P1_I83_S1","P1_I84_S1","P1_I85_S1","P1_I86_S1","P1_I87_S1","P1_I88_S1","P1_I89_S1"))
seq <- data_geno_seq8

data_geno_seq9 <- filter(data_geno_seq, sample %in% c("P1_I9_S1","P1_I90_S1","P1_I91_S1","P1_I92_S1","P1_I93_S1","P1_I94_S1","P1_I95_S1","P1_I96_S1","P1_I97_S1","P1_I98_S1","P1_I99_S1" ))
seq <- data_geno_seq9
#######################################################################################################


# numCores ######################################
numCores <- parallel::detectCores() # Requires library(parallel)
numCores

##  data swap #########################################################################################
geno_novel_list_df_all_def <- data.frame()
geno_novel_list_df_all_new <- data.frame()

people<-unique(seq$sample)
registerDoParallel(2)
foreach (p = people) %dopar% {
  genes_p <- head(unique(seq$v_gene[seq$sample == p]))
  geno_final_df <- c()
  novel_def <- c()
  novel_new <- c()
  geno_novel_list_def <- list()
  geno_novel_list_new <- list()
  for(g in genes_p){
    geno_alleles <- unique(seq$v_alleles[seq$v_gene == g])
    geno_novel_list_def[[g]] <- c()
    geno_novel_list_new[[g]] <- c()
    for(a in geno_alleles) {
      print(paste0(p, " ", g, " ", a))
      geno_sub <- seq[seq$v_gene == g, ] #data of g
      geno_sub$v_call <- a # swap all the alleles of g to a
      
      novel_def <- findNovelAlleles(geno_sub, germline[a] ,v_call="v_call", j_call="j_call",
                                    junction="junction", junction_length="junction_length",
                                    germline_min=1, seq="sequence_alignment", min_seqs = 1, 
                                    pos_range = 1:318,mut_range = 1:10, alpha=0.05,j_max = 0.15)
      novel_new <- findNovelAlleles(geno_sub, germline[a] ,v_call="v_call", j_call="j_call",
                                    junction="junction", junction_length="junction_length",
                                    germline_min=1, seq="sequence_alignment", min_seqs = 1, 
                                    pos_range = 1:318, mut_range = 1:9, alpha= 0.1, j_max = 0.3)
    # # #
    geno_novel_list_def[[g]] <- dplyr::bind_rows(geno_novel_list_def[[g]], novel_def)
    geno_novel_list_new[[g]] <- dplyr::bind_rows(geno_novel_list_new[[g]], novel_new)
    }
  }
  geno_novel_list_df_def <- data.table::rbindlist(geno_novel_list_def)
  geno_novel_list_df_new <- data.table::rbindlist(geno_novel_list_new)
  geno_novel_list_df_def <- geno_novel_list_df_def %>% mutate(subject = p)
  geno_novel_list_df_new <- geno_novel_list_df_new %>% mutate(subject = p)
  geno_novel_list_df_all_def <- dplyr::bind_rows(geno_novel_list_df_all_def,geno_novel_list_df_def)
  geno_novel_list_df_all_new <- dplyr::bind_rows(geno_novel_list_df_all_new,geno_novel_list_df_new)
  # keep only found alleles
  geno_novel_list_df_all_def <- geno_novel_list_df_all_def %>% filter(note == "Novel allele found!")
  geno_novel_list_df_all_new <- geno_novel_list_df_all_new %>% filter(note == "Novel allele found!")
  tmp_d <- paste0("geno_novel_list_df_",p,"_def.csv")
  tmp_n <- paste0("geno_novel_list_df_",p,"_new.csv")
  write.csv(geno_novel_list_df_all_def ,tmp_d, row.names = FALSE)
  write.csv(geno_novel_list_df_all_new ,tmp_n, row.names = FALSE)
  print(paste0("DONE P ", p, " !!!!!!!!!!!!!!!!!!!!!!!!!!!"))
}
stopCluster(cl)

# # # # # # #

## csv merge ###########################################################################
library(readr)
## run on all people = full data = can run only if we have 99*2 csv######################
geno_novel_list_df_all_def <- data.frame()
geno_novel_list_df_all_new <- data.frame()
# for(p in pepole) {}
# peeps <- c("P1_I19_S1","P1_I20_S1")
peeps <- c("P1_I9_S1")
for(p in peeps) {
  tmp_d <- paste0("geno_novel_list_df_",p,"_def.csv")
  tmp_n <- paste0("geno_novel_list_df_",p,"_new.csv")
  geno_novel_list_df_tmp_def <-  data.frame()
  geno_novel_list_df_tmp_new <-  data.frame()
  geno_novel_list_df_tmp_def <- read_csv(tmp_d, show_col_types = FALSE)
  geno_novel_list_df_tmp_new <- read_csv(tmp_n, show_col_types = FALSE)
  
  geno_novel_list_df_all_def <- dplyr::bind_rows(geno_novel_list_df_all_def,geno_novel_list_df_tmp_def)
  geno_novel_list_df_all_new <- dplyr::bind_rows(geno_novel_list_df_all_new,geno_novel_list_df_tmp_new)
}
#View(geno_novel_list_df_all_def)
#View(geno_novel_list_df_all_new)


## run on part data - run through csv files  ##############################################
## def.csv - list of all def file names
geno_novel_list_df_all_def <- data.frame()
def_list_tmp <- list.files(path=".", pattern="_def.csv", all.files=TRUE, full.names=TRUE)
def_list <- character(0)
for(d in def_list_tmp) {
  dd <- substr(d,3,nchar(d))
  def_list <- append(def_list, dd)
}
for(d in def_list) {
  geno_novel_list_df_tmp_def <-  data.frame()
  geno_novel_list_df_tmp_def <- read_csv(d, show_col_types = FALSE)
  geno_novel_list_df_all_def <- dplyr::bind_rows(geno_novel_list_df_all_def,geno_novel_list_df_tmp_def)
}
View(geno_novel_list_df_all_def)
## new_csv - list of all new file names
new_list_tmp <- list.files(path=".", pattern="S1_new.csv", all.files=TRUE, full.names=TRUE) ## CHANGE get rid of S1
new_list <- character(0)
for(n in new_list_tmp) {
  nn <- substr(n,3,nchar(n))
  new_list <- append(new_list, nn)
}
geno_novel_list_df_all_new <- data.frame()
for(n in new_list) {
  geno_novel_list_df_tmp_new <-  data.frame()
  geno_novel_list_df_tmp_new <- read_csv(n, show_col_types = FALSE)
  geno_novel_list_df_all_new <- dplyr::bind_rows(geno_novel_list_df_all_new,geno_novel_list_df_tmp_new)
}
View(geno_novel_list_df_all_new)


# allele list found novel #################################################
alleles_list_default <- unique(geno_novel_list_df_all_def$germline_call)
alleles_list_new <- unique(geno_novel_list_df_all_new$germline_call)

# in data_geno_seq - genotype data  #######################################
seq<-data_geno_seq
alleles <- unique(seq$v_alleles)
n_alleles <- length(alleles)
samples <- unique(geno_novel_list_df_all_new$subject) 
n_samples <- length(samples)

# matrix creation for heatmap  ############################################
allele_geno_mat <- matrix(0, nrow = n_samples, n_alleles ,dimnames = list(samples, alleles))
legend_p <- c(none = 0, defalut = 1, new = 2, both = 3)
for(pp in samples) {
  alleles_list_default <- unique(geno_novel_list_df_all_def$germline_call[geno_novel_list_df_all_def$subject == pp])
  alleles_list_new <- geno_novel_list_df_all_new$germline_call[geno_novel_list_df_all_new$subject == pp]
  g_aa <- as.double(as.numeric(alleles %in% alleles_list_default))
  t_aa <- as.double(as.numeric(alleles %in% alleles_list_new) * 2)
  aa_in <- g_aa + t_aa
  allele_geno_mat[pp,] <- aa_in
}
allele_geno_df <- as.data.frame(allele_geno_mat)
#View(allele_geno_df)
#t_df<-t(allele_geno_df) # transpose

## AYELET HEATMAP ###################################################################################
colors_pal1 <- setNames(c("#FFFACD","#FFFFF0","#EE3B3B", "#FFA07A"),(c('none','default','new','both')))
plot_heatmap_geno <- function(allele_geno_mat, colors_pal1, nrow_legend = NULL, ncol_legend = 1, legend_loc = "right", to_pdf = F, filename, draw_p = F){
  legend_p <- c(none = 0, defalut = 1, new = 2, both = 3)
  allele_geno_df <- as.data.frame(allele_geno_mat)
  allele_geno_df$subject <- rownames(allele_geno_df)
  allele_geno_df <- reshape2::melt(allele_geno_df, id.var = "subject")
  allele_geno_df$value <- sapply(allele_geno_df$value, function(x) names(legend_p)[legend_p==x])
  set.seed(123)
  allele_geno_df_bar <- allele_geno_df %>% dplyr::group_by(variable, value) %>% dplyr::count()
  allele_geno_df_bar$variable <- as.character(allele_geno_df_bar$variable)
  vals_bar <- list(both = c(), default = c(), new = c())
  for(a in unique(alleles)){
    for(cat in unique(allele_geno_df_bar$value)){
      v <- allele_geno_df_bar$n[allele_geno_df_bar$variable==a & allele_geno_df_bar$value==cat]
      if(length(v)==0) v <- 0
      vals_bar[[cat]] <- c(vals_bar[[cat]], v)
    }
  }
  
  ############################################
  #m_bar = matrix(0, nrow = 3, ncol = n_alleles)
  #m_bar[2,] <- vals_bar$default
  #m_bar[3,] <- vals_bar$new
  #m_bar[1,] <- vals_bar$both
  
  #anno = anno_barplot(t(m_bar), gp = gpar(fill = colors_pal1[c("both","default","new")]), bar_width = 1, height = unit(6, "cm"), axis_param = list(
  #gp=gpar(fontsize=18)))
  # getFuntionalGroup <- function(x, pat = "IGHV[0-9]F[0-9]-") strsplit(gsub(pat,"",x),"[*]")[[1]][1]
  #col_fun = setNames(sapply(colnames(allele_geno_mat), function(x) group_color_vector[getFuntionalGroup(x,"F[0-9]-G")], USE.NAMES = F),colnames(allele_geno_mat))
  #column_ha = HeatmapAnnotation(bar1 = anno, show_annotation_name = F,
  #                              foo = colnames(allele_geno_mat), 
  #                              col = list(foo = col_fun), show_legend = F,
  #                              annotation_label = gpar(fontsize = 20))
  ############################################
  
  #order matrix
  if(to_pdf) pdf(filename, width = 37, height = 25)
  
  groups_heatmap <- getGene(colnames(allele_geno_mat), strip_d = F, omit_nl = F, sep = "[*]")
  
  ha = Heatmap(allele_geno_mat, name = "mat", 
               row_names_gp = gpar(fontsize = 16), 
               column_names_gp = gpar(fontsize = 16), 
               #top_annotation = column_ha, 
               cluster_rows = F, 
               cluster_columns = F, 
               col = unname(colors_pal1), 
               heatmap_legend_param = list(title = "genotype",
                                           labels = c('none','default','new','both'),
                                           title_gp = gpar(fontsize = 18, fontface = "bold"),
                                           labels_gp = gpar(fontsize = 18),
                                           nrow = nrow_legend,
                                           ncol = ncol_legend,
                                           border = "black"
               )
  )
  
  if(draw_p) draw(ha, heatmap_legend_side = legend_loc) 
  else{
    ha = Heatmap(allele_geno_mat, name = "mat", 
                 row_names_gp = gpar(fontsize = 16), 
                 column_names_gp = gpar(fontsize = 5), #16
                 column_title = "genotype",
                 column_title_side="bottom",
                 row_title = "subject",
                 row_title_side = "right",
                 #top_annotation = column_ha, 
                 cluster_rows = F, 
                 cluster_columns = F, 
                 col = unname(colors_pal1), 
                 #show_heatmap_legend =F,
                 heatmap_legend_param = list(title = "Novel Allele Found",
                                             labels =rev(c('none','default','new','both')),
                                             title_gp = gpar(fontsize = 18, fontface = "bold"),
                                             labels_gp = gpar(fontsize = 12),
                                             nrow = nrow_legend,
                                             ncol = ncol_legend,
                                             border = "black",
                                             title_position = ifelse(!is.null(nrow_legend),"leftcenter",NULL),
                                             heatmap_legend_side = legend_loc
                 )
    )
    lgd = Legend(labels = names(colors_pal1), legend_gp = gpar(fill = unname(colors_pal1)), 
                 title = "genotype",
                 title_gp = gpar(fontsize = 18, fontface = "bold"),
                 labels_gp = gpar(fontsize = 18),
                 nrow = nrow_legend,
                 ncol = ncol_legend,
                 border = "black",
                 title_position = ifelse(!is.null(nrow_legend),"leftcenter",NULL)
    )
    #return(list(heatmap = ha, legend = lgd))
    return(draw(ha, heatmap_legend_side = legend_loc))
  }
  
  if(to_pdf) dev.off()
}

### heatmap  ####################################################################################################
ha = plot_heatmap_geno(allele_geno_mat = allele_geno_mat, colors_pal1, to_pdf = F, nrow_legend = 1, ncol_legend = NULL, legend_loc = "top", draw_p = F)
#################################################################################################################