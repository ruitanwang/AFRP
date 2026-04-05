# bulkseq_rawdata_process

## GSE124439
library(GEOquery)
library(org.Hs.eg.db)
library(tidyverse)
library(Biobase)

rawcounts <- read.csv("./analysiswork/rawdata/bulkRNAseq/GSE124439/GSE124439_raw_counts_GRCh38.p13_NCBI.tsv",sep = "")

gene_ids <- rawcounts$GeneID %>% as.character()

symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

rawcounts$GeneID <- symbols

rawcounts <- rawcounts[!is.na(rawcounts$GeneID), ] %>% as.data.frame()

rawcounts <- rawcounts %>%
  group_by(GeneID) %>%
  summarise(across(where(is.numeric),  mean)) %>% as.data.frame()

rawcounts <- rawcounts %>% 
  column_to_rownames(var = "GeneID")

softdir <- "./analysiswork/rawdata/bulkRNAseq/GSE124439/GSE124439_family.soft.gz"
metasoft <- getGEO(filename = softdir,GSEMatrix = TRUE)
gsmlist <- GSMList(metasoft)
metalist <- list()
for (i in 1:length(gsmlist)) {
  metavec <- gsmlist[[i]]@header[[2]]
  metalist[[i]] <- metavec
  names(metalist)[i] <- gsmlist[[i]]@header[["geo_accession"]]
}

metadf <- as.data.frame(do.call(rbind,  metalist))
rawmetadf <- tibble::rownames_to_column(metadf,var = "Sample")
cleanmetadf <- rawmetadf[,c(1,3,5)]
colnames(cleanmetadf) <- c("sample","tissue","disease")

cleanmetadf$sample <- sub("^.*?:", "", cleanmetadf$sample) %>% trimws()
cleanmetadf$tissue <- sub("^.*?:", "", cleanmetadf$tissue) %>% trimws()
cleanmetadf$disease <- sub("^.*?:", "", cleanmetadf$disease) %>% trimws()

rownames(cleanmetadf) <- cleanmetadf$sample

GSE124439_eset <- ExpressionSet(
  assayData = rawcounts %>% as.matrix(),
  phenoData = AnnotatedDataFrame(cleanmetadf)
  )

save(GSE124439_eset,file = "./analysiswork/cleandata/GSE124439_eset.RData")






######### GSE153960
library(GEOquery)
library(org.Hs.eg.db)
library(tidyverse)
library(Biobase)

rawcounts <- read.csv("./analysiswork/rawdata/bulkRNAseq/GSE153960/GSE153960_raw_counts_GRCh38.p13_NCBI.tsv",sep = "")
gene_ids <- rawcounts$GeneID %>% as.character()

symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

rawcounts$GeneID <- symbols

rawcounts <- rawcounts[!is.na(rawcounts$GeneID), ] %>% as.data.frame()

rawcounts <- rawcounts %>%
  group_by(GeneID) %>%
  summarise(across(where(is.numeric),  mean)) %>% as.data.frame()

rawcounts <- rawcounts %>% 
  column_to_rownames(var = "GeneID")

softdir <- "./analysiswork/rawdata/bulkRNAseq/GSE153960/GSE153960_family.soft.gz"
metasoft <- getGEO(filename = softdir,GSEMatrix = TRUE)
gsmlist <- GSMList(metasoft)
metalist <- list()
for (i in 1:length(gsmlist)) {
  metavec <- gsmlist[[i]]@header[[2]]
  metalist[[i]] <- metavec
  names(metalist)[i] <- gsmlist[[i]]@header[["geo_accession"]]
}


metadf <- as.data.frame(do.call(rbind,  metalist))
cleanmetadf <- metadf[,c(2,4)]
colnames(cleanmetadf) <- c("tissue","disease")


cleanmetadf$tissue <- sub("^.*?:", "", cleanmetadf$tissue) %>% trimws()
cleanmetadf$disease <- sub("^.*?:", "", cleanmetadf$disease) %>% trimws()
cleanmetadf$sample <- rownames(cleanmetadf)
shared_sample <- intersect(rownames(cleanmetadf),colnames(rawcounts))
cleanmetadf <- cleanmetadf[shared_sample,]


GSE153960_eset <- ExpressionSet(
  assayData = rawcounts %>% as.matrix(),
  phenoData = AnnotatedDataFrame(cleanmetadf)
)

save(GSE153960_eset,file = "./analysiswork/cleandata/GSE153960_eset.RData")






## GSE272624
library(GEOquery)
library(org.Hs.eg.db)
library(tidyverse)
library(Biobase)
library(archive)

# process and obtain the read counts and metadata of GSE272624
location_dir <- "./analysiswork/rawdata/bulkRNAseq/GSE272624/counts"

file_list <- list.files(location_dir,  full.names  = TRUE, pattern = "\\.txt|\\.csv|\\.tsv")
reads_list <- list()

for (i in 1:length(file_list)) {
  sample_id <- sub(".*/([^_]+)_.*", "\\1", file_list[i])
  
  df <- read_delim(
    file = gzfile(file_list[i]), 
    delim = "\t", 
    locale = locale(encoding = "GB18030"),
    guess_max = 100000
  )
  
  colnames(df) <- c("Gene",sample_id)
  
  reads_list[[i]] <- df
}

merged_df <- purrr::reduce(
  reads_list, 
  function(x, y) {full_join(x, y, by = "Gene")}
)
merged_df <- merged_df %>% column_to_rownames(var = "Gene")

softdir <- "./analysiswork/rawdata/bulkRNAseq/GSE272624/metadata/GSE272624_family.soft.gz"
metasoft <- getGEO(filename = softdir,GSEMatrix = TRUE)
gsmlist <- GSMList(metasoft)
metalist <- list()

for (i in 1:length(gsmlist)) {
  metavec <- gsmlist[[i]]@header[[2]]
  metalist[[i]] <- metavec
  names(metalist)[i] <- gsmlist[[i]]@header[["geo_accession"]]
}

metadf <- as.data.frame(do.call(rbind,  metalist))
rawmetadf <- tibble::rownames_to_column(metadf,var = "Sample")
rownames(rawmetadf) <- rawmetadf$Sample
cleanmetadf <- rawmetadf[,c(1,2,3,4)]
colnames(cleanmetadf) <- c("sample","tissue","disease","gender")


cleanmetadf$tissue <- sub("^.*?:", "", cleanmetadf$tissue) %>% trimws()
cleanmetadf$disease <- sub("^.*?:", "", cleanmetadf$disease) %>% trimws()
cleanmetadf$gender <- sub("^.*?:", "", cleanmetadf$gender) %>% trimws()

GSE272624_eset <- ExpressionSet(
  assayData = merged_df %>% as.matrix(),
  phenoData = AnnotatedDataFrame(cleanmetadf)
)

save(GSE272624_eset,file = "./analysiswork/cleandata/GSE272624_eset.RData")





## GSE272626
library(GEOquery)
library(org.Hs.eg.db)
library(tidyverse)
library(Biobase)
library(archive)

# process and obtain the read counts and metadata of GSE272626
location_dir <- "./analysiswork/rawdata/bulkRNAseq/GSE272626/counts"

file_list <- list.files(location_dir,  full.names  = TRUE, pattern = "\\.txt|\\.csv|\\.tsv")
reads_list <- list()

for (i in 1:length(file_list)) {
  sample_id <- sub(".*/([^_]+)_.*", "\\1", file_list[i])
  
  df <- read_delim(
    file = gzfile(file_list[i]), 
    delim = "\t", 
    locale = locale(encoding = "GB18030"),
    guess_max = 100000
  )
  
  colnames(df) <- c("Gene",sample_id)
  
  reads_list[[i]] <- df
}

merged_df <- purrr::reduce(
  reads_list, 
  function(x, y) {full_join(x, y, by = "Gene")}
)
merged_df <- merged_df %>% column_to_rownames(var = "Gene")

softdir <- "./analysiswork/rawdata/bulkRNAseq/GSE272626/meta/GSE272626_family.soft.gz"
metasoft <- getGEO(filename = softdir,GSEMatrix = TRUE)
gsmlist <- GSMList(metasoft)
metalist <- list()

for (i in 1:length(gsmlist)) {
  metavec <- gsmlist[[i]]@header[[2]]
  metalist[[i]] <- metavec
  names(metalist)[i] <- gsmlist[[i]]@header[["geo_accession"]]
}

metadf <- as.data.frame(do.call(rbind,  metalist))
rawmetadf <- tibble::rownames_to_column(metadf,var = "Sample")
rownames(rawmetadf) <- rawmetadf$Sample
cleanmetadf <- rawmetadf[,c(1,2,3)]
colnames(cleanmetadf) <- c("sample","tissue","disease")

cleanmetadf$sample <- sub("^.*?:", "", cleanmetadf$sample) %>% trimws()
cleanmetadf$tissue <- sub("^.*?:", "", cleanmetadf$tissue) %>% trimws()
cleanmetadf$disease <- sub("^.*?:", "", cleanmetadf$disease) %>% trimws()

GSE272626_eset <- ExpressionSet(
  assayData = merged_df %>% as.matrix(),
  phenoData = AnnotatedDataFrame(cleanmetadf)
)

save(GSE272626_eset,file = "./analysiswork/cleandata/GSE272626_eset.RData")



# 
# ## GSE116622
# rm(list = ls())
# library(GEOquery)
# library(org.Hs.eg.db)
# library(tidyverse)
# library(Biobase)
# library(GenomicFeatures)
# rawcounts <- read.csv("./analysiswork/rawdata/bulkRNAseq/GSE116622/GSE116622_norm_counts_TPM_GRCh38.p13_NCBI.tsv",sep = "")
# txdb <- makeTxDbFromGFF("./analysiswork/helpdata/gencode.v49.chr_patch_hapl_scaff.annotation (1).gtf.gz", format = "gtf")
# exons_by_gene <- exonsBy(txdb, by = "gene")
# gene_exon_lengths <- reduce(exons_by_gene)
# gene_lengths <- sum(width(gene_exon_lengths))
# gene_ids <- rawcounts$GeneID %>% as.character()
# 
# gene_length_df <- data.frame(
#   gene_id = names(gene_lengths),
#   exon_length = as.numeric(gene_lengths)
# )
# 
# gene_length_df$gene_symbol <- mapIds(
#   org.Hs.eg.db,
#   keys = sub("\\..*", "", gene_length_df$gene_id),
#   keytype = "ENSEMBL",
#   column = "SYMBOL"
# )
# 
# ncbi_ids <- gene_ids
# ensembl_ids <- mapIds(org.Hs.eg.db, keys = ncbi_ids, keytype = "ENTREZID", column = "ENSEMBL")
# 
# ssa <- getGEO(GEO = "GSE116622")
