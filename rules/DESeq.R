## Perfoms DE analysis using kallisto output (est_counts) and plots:
# PC analysis
# heatmap of 100 most differentially expressed genes
# MA plot

args<-commandArgs(TRUE)

CreateDDS<-function(design, name.dds,
                    name.rld, name.heatmap, name.pca, name.maplot){

  library(DESeq2)
  library(biomaRt)
  library(RColorBrewer)
  library(gplots)
  
  #design<-"/home/carlos/remoteDir/data7/projects/p283_rna_and_disease/projects/sinergia/sinergia/sinergia_main/experimental_design_test.txt"
  #name.dds<-"/home/carlos/remoteDir/data7/projects/p283_rna_and_disease/projects/sinergia/sinergia/sinergia_main/DESeq/dds.out"
  #name.rld<-"/home/carlos/remoteDir/data7/projects/p283_rna_and_disease/projects/sinergia/sinergia/sinergia_main/DESeq/rld.out"
  #name.heatmap<-"/home/carlos/remoteDir/data7/projects/p283_rna_and_disease/projects/sinergia/sinergia/sinergia_main/DESeq/heatmap.pdf"
  #name.pca<-"/home/carlos/remoteDir/data7/projects/p283_rna_and_disease/projects/sinergia/sinergia/sinergia_main/DESeq/pca.pdf"
  #name.maplot<-"/home/carlos/remoteDir/data7/projects/p283_rna_and_disease/projects/sinergia/sinergia/sinergia_main/DESeq/maplot.pdf"
  
  # Get path to each sample quantification done by kallisto and create counts table -------------------------------------------------------------
  # col #4 in kallisto output
  setwd(getwd())
  samples_path <- c()
  expDesign<-read.table(design, header=FALSE, sep="\t", colClasses = c("character", "character"))
  samples <- expDesign[,1]
  samples_path <- paste0(getwd(),"/kallisto/", samples, "_kallisto_quantification/abundance.tsv")

  #sample1_path <- paste0("/kallisto/", samples[1], "kallisto_quantification/abundance.tsv")
  s1 <- read.table(samples_path[1], sep = '\t', header = TRUE)
  #tpms <- cbind(s1[,5])  #for tpm
  est_counts <- cbind(s1[,4])  #for est_counts
  colnames(est_counts)<-samples[1]
  rownames(est_counts)<-s1[,"target_id"]

  for (s in samples_path){
    if (s == samples_path[1]) {
      next
    }
    tsv <- read.table(file = s, sep = '\t', header = TRUE)
    est_counts <- cbind(est_counts, tsv[,4])
  }
  colnames(est_counts)<-samples
  
  # write est_counts into files -------------------------------------------------------------
  est_counts<-as.data.frame(est_counts)
  counts_file=paste0(getwd(), '/DESeq/est_counts.tsv')
  write.table(est_counts, file=counts_file, sep='\t', col.names = NA, row.names = TRUE)

  # Prepare data object and experimental design -------------------------------------------------------------
  data <- read.table(counts_file, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  d <- as.data.frame(apply(data, 2, as.integer))
  rownames(d) <- rownames(data)
  data<-d

  expDesign<-read.table(design, header=FALSE, sep="\t", colClasses = c("character", "character"))
  colnames(expDesign)<-c("Sample", "Group")
  expDesign$Group <- as.factor(expDesign$Group)
  
  expDesign<-expDesign[match(colnames(data), as.character(expDesign$Sample)),]
  diff<-DESeqDataSetFromMatrix(data, expDesign, ~Group)  
  
  # DE analysis -------------------------------------------------------------
  dds<-DESeq(diff, betaPrior = TRUE)
  #keep <- rowSums(counts(dds)) >= 4

  cntNonZero <- apply(counts(dds), 1, function(x) sum(x != 0))
  n <- round(ncol(counts(dds))/2)-1 #half number of columns minus 1
  cntNonZero <- cntNonZero >= n
  
  dds <- dds[cntNonZero,]
  saveRDS(dds, name.dds)
  
  rld<-rlog(dds, blind=TRUE) # apply a regularized log transformation, ignoring information about experimental groups
  saveRDS(rld, name.rld)  
  
  # PCA -------------------------------------------------------------------
  p<-plotPCA(rld, intgroup=c("Group"))
  pdf(name.pca)
  print(p, ntop=500)
  dev.off()
  
  # Heatmap of most highly expressed genes ----------------------------------
  #select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:500]
  res <- results(dds)
  n<- nrow(res)/4
  if (n > 100) { n<- 100 }
  select <- order(res[,6])[1:n]
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  pdf(name.heatmap)
  heatmap.2(assay(rld)[select,], col = hmcol, trace="none", margin=c(10, 6),
            labCol=colnames(dds), cexRow = 0.4)
  dev.off()
  
  # MA plot -----------------------------------------------------------------
  pdf(name.maplot)
  plotMA(dds)
  dev.off()
  
}


CreateDDS(args[1], args[2], args[3], args[4], args[5], args[6])
  
  
  
  
  
  
  
  
  
  
    