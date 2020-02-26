Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

library(sleuth, lib.loc="./scripts/R")
library(RColorBrewer, lib.loc="./scripts/R")
library(ggplot2, lib.loc="./scripts/R")
library(reshape, lib.loc="./scripts/R")
library(gplots, lib.loc="./scripts/R")

#design <- "experimental_design_test.txt"
#transcript_map <- "transcript_to_gene_map.txt"
#project <- "Sinergia"

args<-commandArgs(TRUE)

DEanalysis<-function(design, transcript_map, project){
  ##Read experimental design file
  expDesign <- read.table(design, header=TRUE, sep="\t", colClasses = c("character", "character"))
  
  ##Read transcripts to genes map
  ##Must have two columns: "target_id" and "gene_name"; tabulation separated fields
  Ids_mapping <- read.table(transcript_map, header = T)

  ###------------VVV----------###
  ### Count based aggregation ###
  ##Prepare model
  s2c <- read.table(file.path(design), header = TRUE, stringsAsFactors=FALSE)
  paths <- paste0(project, "/kallisto/", s2c$sample)
  s2c <- dplyr::select(s2c, sample = 1, condition)
  s2c <- dplyr::mutate(s2c, path = paths)
  #print(s2c)

  ##Fit model and calculate differential expression
  #exp_filter <- function(row, min_reads = 5, min_prop = 0.125){
  #                     mean(row >= min_reads) >= min_prop
  #            }
  # so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = Ids_mapping, aggregation_column = 'gene_name', gene_mode=TRUE, num_cores = 1, filter_fun = exp_filter)
  so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = Ids_mapping, aggregation_column = 'gene_name', gene_mode=TRUE)
  so <- sleuth_fit(so, ~condition, 'full')
  so <- sleuth_fit(so, ~1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE)
  
  #Write table into file
  counts_file=paste0(getwd(), "/", project, '/sleuth/est_counts.tsv')
  write.table(sleuth_table, file=counts_file, sep=',', col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  ##DEGs
  sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
  
  #Write table of significant genes
  counts_file_significant=paste0(getwd(),  project, '/sleuth/est_counts_significant.tsv')
  write.table(sleuth_significant, file=counts_file_significant, sep=',', col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  ##Scaled reads per base (expression) Matrix
  ExpressionMatrix = sleuth_to_matrix(obj = so,which_df = "obs_norm",which_units = "tpm")
  expMat=paste0(getwd(),  project, '/sleuth/expression_matrix.csv')
  write.table(ExpressionMatrix, file=expMat, sep=',', col.names = NA, row.names = TRUE, quote = FALSE)

  ##Plot pca
  p <- plot_pca(so, color_by = "condition")
  pdf(project+"/sleuth/PC.pdf")
  print(p, ntop=500)
  dev.off()

}

#args[1]=Experimental_design; args[2]=transcript_map; args[3]=project_name
DEanalysis(args[1], args[2], args[3])



