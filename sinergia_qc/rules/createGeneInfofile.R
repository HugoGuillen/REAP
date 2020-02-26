## The script expects two arguments:
## 1) ensembl dataset to be used, e.g. hsapiens_gene_ensembl. The script does not check if the provided value is valid. Biomart will give an informative error if it isn't
## 2) Table of counts from featureCounts. The script checks that there is a column with ensembl IDs


args<-commandArgs(TRUE)



#########################################################################################################
# FUNCTIONS
#########################################################################################################

createGeneTable<-function(ensembl_dataset, count_table){
  ## Assumes that count_table contains a column with ensemblIDs with header "Geneid" (as produced by featureCounts)
  
  library(biomaRt)
  # take ensembl IDs for all genes in the experiment from the count table
  # get gene length, gc content and biotype from biomart
  # Note that this will use the LATEST genome version, e.g. hg20
  ensembl<-useMart("ENSEMBL_MART_ENSEMBL")
  ensembl<-useDataset(ensembl_dataset, mart=ensembl)
  attributes<-c("ensembl_gene_id", "gene_biotype", "transcript_length", "percentage_gc_content")
  # some of these may not exist depending on the species
  missingColumns <- attributes[!is.element(attributes, listAttributes(ensembl)$name)]
  attributes <- attributes[is.element(attributes, listAttributes(ensembl)$name)]
  
  geneinfo<-getBM(attributes, filters="ensembl_gene_id", values=data$Geneid, mart=ensembl)
  
  # Retain only the longest value for transcript length for each gene
  maxRows<-by(geneinfo, as.factor(geneinfo$ensembl_gene_id), function(x) x[which.max(x$transcript_length),])
  geneinfo<-do.call("rbind", maxRows)
  
  # if there were missing attributes, add column containing specific value (qualimap counts cannot handle NA)
  # Neither transcript length nor GC-content will be used in the final report
  for (i in missingColumns) {
    switch( i,
            "ensembl_gene_id" = { stop("Table with gene info for qualimap counts cannot be created because ensembl gene id does not exist") },
            "gene_biotype" = { 
              geneinfo$gene_biotype <- "artifact"
              print("Gene biotype is not known. All genes will be put in category 'other' for qualimap counts report")
              }, # will be converted to biotype "other"
            "transcript_length" = { 
              geneinfo$transcript_length <- sample(1:1000, nrow(geneinfo), replace=TRUE)
              print("Transcript length not known. Genes will be assigned a random value between 1 and 1000 for qualimap counts report")
              },
            "percentage_gc_content" = { 
              geneinfo$percentage_gc_content <- sample(seq(1,100,by=0.1), nrow(geneinfo), replace=TRUE)
              print("GC-content not known. Genes will be assigned a random value between 1 and 100 for qualimap counts report")
              }
    )
    
  }
  
  temp<-merge(x=as.data.frame(data$Geneid), y=geneinfo, by.x="data$Geneid", by.y="ensembl_gene_id", all.x=TRUE, all.y=FALSE)
  invisible(temp)
}
  
  

simplifyBiotype<-function(df){
  # set up a named vector translationTable based on
  # http://www.ensembl.org/Help/Faq?id=468 and
  # https://www.gencodegenes.org/gencode_biotypes.html
  translationTable<-c("protein_coding", "protein_coding", "protein_coding", "protein_coding", "protein_coding", "protein_coding", "protein_coding", "protein_coding", "protein_coding", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "rRNA", "tRNA", "short_noncoding_RNA", "short_noncoding_RNA", "rRNA", "short_noncoding_RNA", "short_noncoding_RNA", "short_noncoding_RNA", "short_noncoding_RNA", "short_noncoding_RNA", "short_noncoding_RNA", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "protein_coding", "protein_coding", "protein_coding", "protein_coding", "protein_coding", "long_noncoding_RNA", "long_noncoding_RNA", "long_noncoding_RNA", "long_noncoding_RNA", "long_noncoding_RNA", "long_noncoding_RNA", "long_noncoding_RNA", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "pseudogene", "other", "long_noncoding_RNA", "long_noncoding_RNA", "long_noncoding_RNA", "pseudogene", "short_noncoding_RNA", "long_noncoding_RNA") 
  names(translationTable)<-c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_D_gene", "IG_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "TR_V_pseudogene", "TR_J_pseudogene", "Mt_rRNA", "Mt_tRNA", "miRNA", "misc_RNA", "rRNA", "scRNA", "snRNA", "snoRNA", "ribozyme", "sRNA", "scaRNA", "Mt_tRNA_pseudogene", "tRNA_pseudogene", "snoRNA_pseudogene", "snRNA_pseudogene", "scRNA_pseudogene", "rRNA_pseudogene", "misc_RNA_pseudogene", "miRNA_pseudogene", "TEC", "nonsense_mediated_decay", "non_stop_decay", "retained_intron", "protein_coding", "processed_transcript", "non_coding", "ambiguous_orf", "sense_intronic", "sense_overlapping", "antisense", "known_ncrna", "pseudogene", "processed_pseudogene", "polymorphic_pseudogene", "retrotransposed", "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "artifact", "lincRNA", "macro_lncRNA", "3prime_overlapping_ncRNA", "disrupted_domain", "vaultRNA", "bidirectional_promoter_lncRNA")
  
  df$reduced_biotype<-translationTable[df$gene_biotype]
  # if a gene_biotype does not exist in the lookup table, reduced_biotype will be NA
  invisible(df)
}


reformatTable<-function(df){
  temp<-data.frame(biotypes=df$reduced_biotype, 
                   length=df$transcript_length,
                   gc=df$percentage_gc_content)
  rownames(temp)<-df$`data$Geneid`
  invisible(temp)
}
  
#########################################################################################################


data<-read.table(args[2], header=TRUE, check.names=FALSE)

# Check if the counts table contains a column called "Geneid"
if (!("Geneid" %in% colnames(data))){
    # checks if there is a column containing only entries that start with "ENS"
    ens.ids<-apply(data, 2, function(x) {length(x)==length(grep("ENS", x))})
    # if yes, this column will be given the header "Geneid"
    if(sum(ens.ids)==1){
      colnames(data)[ens.ids]<-"Geneid"
    } else {
      stop("The count table does not contain a column with ensembl IDs")
    }
  }

temp<-createGeneTable(args[1], data) 
temp<-simplifyBiotype(temp)
temp<-reformatTable(temp)
write.table(temp, "geneinfofile.txt", quote=FALSE, sep="\t", row.names=TRUE)
