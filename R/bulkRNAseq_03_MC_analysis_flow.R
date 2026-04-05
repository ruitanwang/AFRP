### MC flow
library(base)
library(GEOquery)
library(tidyverse)
library(ggplot2)
library(qs2)
library(Biobase)
library(RColorBrewer)
library(ggpie)
library(stringr)
library(Seurat)
library(monocle3)
library(PCAtools)
library(RColorBrewer)
library(ggnewscale)
library(cowplot)
library(DESeq2)
library(edgeR)
library(limma)
library(ggvenn)
library(ComplexUpset)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GseaVis)

load("./analysiswork/cleandata/GSE124439_eset.RData")
load("./analysiswork/cleandata/GSE153960_eset.RData")
load("./analysiswork/cleandata/GSE272624_eset.RData")
load("./analysiswork/cleandata/GSE272626_eset.RData")

seuobj <- qs_read("./analysiswork/cleandata/snRNA_MCA_seurat.qs2")
pseudobulk <- AggregateExpression(
  object = seuobj,
  assays = "RNA",
  group.by = "sample_ID"
)


peseudobulk_df <- pseudobulk[["RNA"]] %>% as.data.frame()
sample_meta_df <- qs_read("./analysiswork/cleandata/snRNAseq_splited_cohort_Seurat_sample_meta.qs2")
select_meta_df <- sample_meta_df[,c("sample_ID","Region","sub_disease_type","Sex","Batch")]
colnames(select_meta_df) <- c("sample","tissue","disease","gender","batch")
select_meta_df_with_rownames <-column_to_rownames(select_meta_df,var = "sample")
select_meta_df_with_rownames$sample <- rownames(select_meta_df_with_rownames)
select_meta_df_sorted <- select_meta_df_with_rownames[colnames(peseudobulk_df),]
select_meta_df_sorted <- select_meta_df_sorted[,c(5,1,2,3,4)]
peseudobulk_eset <- ExpressionSet(
  assayData = peseudobulk_df %>% as.matrix(),
  phenoData = AnnotatedDataFrame(select_meta_df_sorted)
)

peseudobulk_eset@phenoData@dimLabels <- c("rowNames","columnNames")

for (esetname in ls(pattern = "^GSE")) {
  esetob <- get(esetname)
  pData <- Biobase::pData(esetob)
  pData$batch <- str_replace(esetname,"_.*","")
  phenoData(esetob) <- AnnotatedDataFrame(pData)
  assign(esetname,esetob)
}

common_genes <- Reduce(
  intersect,
  list(
    GSE124439_genes <- rownames(GSE124439_eset),
    GSE153960_genes <- rownames(GSE153960_eset),
    GSE272624_genes <- rownames(GSE272624_eset),
    GSE272626_genes <- rownames(GSE272626_eset),
    peseudobulk_genes <- rownames(peseudobulk_eset)
  )
)


filter_GSE124439_eset <- GSE124439_eset[common_genes,]
filter_GSE153960_eset <- GSE153960_eset[common_genes,]
filter_GSE272624_eset <- GSE272624_eset[common_genes,]
filter_GSE272626_eset <- GSE272626_eset[common_genes,]
filter_peseudobulk_eset <- peseudobulk_eset[common_genes,]

combined_eset <- Biobase::combine(
  filter_GSE124439_eset, 
  filter_GSE153960_eset, 
  filter_GSE272624_eset, 
  filter_GSE272626_eset,
  filter_peseudobulk_eset
)

length(unique(rownames(combined_eset)))

eset <- combined_eset
zero_prop <- rowMeans(Biobase::exprs(eset) == 0)
keep_genes <- zero_prop < 0.5
filtered_eset <- eset[keep_genes, ]
length(unique(rownames(filtered_eset)))
final_bulkseq_eset <- filtered_eset

eset <- final_bulkseq_eset
phenedf <- Biobase::pData(eset)
unique(phenedf$disease)
table(phenedf$disease)
phenedf <- phenedf %>%
  mutate(
    group = case_when(
      disease == "ALS Spectrum MND" ~ "ALS",
      disease == "C9ALS" ~ "ALS",
      disease == "SALS" ~ "ALS",
      disease == "C9FTLD" ~ "FTLD",
      disease == "SFTLD" ~ "FTLD",
      disease == "Non-Neurological Control" ~ "CON",
      disease == "Other Neurological Disorders" ~ "OND",
      disease == "ALS Spectrum MND, Other Neurological Disorders" ~ "ALS+OND",
      disease == "Pre-fALS, Other Neurological Disorders" ~ "Pre-fALS",
      disease == "Unknown" ~ "Unknown",
      disease == "Other MND" ~ "Other MND",
      disease == "Pre-fALS" ~ "Pre-fALS",
      disease == "Other Neurological DIsorders" ~ "OND",
      disease == "ALS" ~ "ALS",
      disease == "FTD" ~ "FTLD",
      disease == "CTL" ~ "CON",
      disease == "CON" ~ "CON",
      disease == "NEURCTL" ~ "CON",
      disease == "ALS+FTD" ~ "FTLD",
      TRUE ~ "General"
    )
  )
unique(phenedf$group)

unique(phenedf$tissue)
table(phenedf$tissue)
phenedf <- phenedf %>%
  mutate(
    region = case_when(
      tissue == "MC" ~ "MC",
      tissue == "FC" ~ "FC",
      tissue == "Frontal Cortex" ~ "FC",
      tissue == "Motor Cortex (Medial)" ~ "MC",
      tissue == "Motor Cortex (Lateral)" ~ "MC",
      tissue == "Motor Cortex" ~ "MC",
      tissue == "Spinal Cord Lumbar" ~ "SC",
      tissue == "Spinal Cord Cervical" ~ "SC",
      tissue == "Cerebellum" ~ "Cere",
      tissue == "Cortex Motor Lateral" ~ "MC",
      tissue == "Cortex Frontal" ~ "FC",
      tissue == "Cortex Motor Unspecified" ~ "MC",
      tissue == "Cortex Motor Medial" ~ "MC",
      tissue == "Cortex Sensory" ~ "SenC",
      tissue == "Cortex Occipital" ~ "OC",
      tissue == "Spinal Cord Thoracic" ~ "SC",
      tissue == "Cortex Temporal" ~ "TC",
      tissue == "Hippocampus" ~ "Hipp",
      tissue == "motor cortex" ~ "MC",
      tissue == "frontal cortex" ~ "FC",
      tissue == "medial motor cortex" ~ "MC",
      tissue == "lateral motor cortex" ~ "MC",
      tissue == "lumbar spinal cord" ~ "SC",
      tissue == "thoracic spinal cord" ~ "SC",
      tissue == "cervical spinal cord" ~ "SC",
      tissue == "spinal cord" ~ "SC",
      TRUE ~ "General"
    )
  )
unique(phenedf$region)

phenoData(eset) <- AnnotatedDataFrame(phenedf)
final_bulkseq_eset <- eset
eset <- final_bulkseq_eset
mc_samples <- rownames(Biobase::pData(eset))[Biobase::pData(eset)$region == "MC"]
eset_mc <- eset[,mc_samples]
Biobase::pData(eset_mc)$mergegroup <- case_when(
  Biobase::pData(eset_mc)$group %in% c("ALS","CON") ~ Biobase::pData(eset_mc)$group,
  TRUE ~ "Other"
)
qs_save(eset_mc,"./analysiswork/cleandata/bulkRNAseq_mc_eset.qs2")


eset <- qs_read("./analysiswork/cleandata/bulkRNAseq_mc_eset.qs2")
phenodf <- Biobase::pData(eset)
select_disease <- c("ALS","FTLD","CON")
p1 <- group_bar(
  plotdf = phenodf[phenodf$group %in% select_disease,],
  groupcol = group,
  plot_title = "Motor Cortex Samples",
  legend_ncol = 1,
  disease_color = c("#f09ba0","#b7b7eb","#9bbbee")
)
p1


p2 <- ggpie3D(
  data = phenodf[phenodf$group %in% select_disease,],
  group_key = "group",
  count_type = "full",
  fill_color = c("#f09ba0","#b7b7eb","#9bbbee"),
  start_degrees = 0,
  tilt_degrees = -20,
  height = 0.2,
  darken = 0.15,
  camera_eye = c(0, 3, 3),
  camera_look_at = c(0, 0, 0),
  show_label = TRUE,
  label_info = "ratio",
  label_split = "[[:space:]]+",
  label_len = 40,
  label_size = 5
)

p2

combined_plot <- cowplot::plot_grid(
  p1,
  p2,
  nrow = 1,
  ncol = 2,
  labels = c("A","B"),
  rel_widths = c(5,5)
)


ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_001.tiff",
  plot = combined_plot,
  dpi = 1000,
  width = 9,
  height = 5,
  units = "in"
)


eset <- qs_read("./analysiswork/cleandata/bulkRNAseq_mc_eset.qs2")
select_disease <- c("ALS","FTLD","CON")
select_samples <- rownames(Biobase::pData(eset))[Biobase::pData(eset)$group %in% select_disease]
eset <- eset[,select_samples]
count_matrix <- Biobase::exprs(eset)
count_matrix[is.na(count_matrix)] <- 0
logcpm <- cpm(count_matrix, log = TRUE, prior.count = 1)

pca_res <- pca(
  mat = logcpm, 
  metadata = Biobase::pData(eset),
  center = TRUE,
  scale = TRUE,
  rank = 10
)


p1 <- pca_plot(input_batch = Biobase::pData(eset)$batch,plot_title = "Before remove batch") +
  theme(
    legend.position = "right"
  )
p1


corrected_counts <- sva::ComBat_seq(Biobase::exprs(eset), batch = Biobase::pData(eset)$batch,group = Biobase::pData(eset)$group)

rm_batch_mc_bulkseq_eset <- ExpressionSet(
  assayData = corrected_counts %>% as.matrix(),
  phenoData = Biobase::AnnotatedDataFrame(Biobase::pData(eset))
)

qs_save(rm_batch_mc_bulkseq_eset,"./analysiswork/cleandata/bulkRNAseq_MC_rm_batch_eset.qs2")
sum(is.na(rownames(corrected_counts)))

eset <- rm_batch_mc_bulkseq_eset
count_matrix <- Biobase::exprs(eset)
count_matrix[is.na(count_matrix)] <- 0
logcpm <- cpm(count_matrix, log = TRUE, prior.count = 1)

pca_res <- PCAtools::pca(
  mat = logcpm, 
  metadata = Biobase::pData(eset),
  center = TRUE,
  scale = TRUE,
  rank = 10
)


p2 <- pca_plot(input_batch = Biobase::pData(eset)$batch,plot_title = "After remove batch") +
  theme(
    legend.position = "right"
  )
p2

common_legend <- get_legend(
  p2 + theme(
    legend.box = "vertical",
    legend.position = "top",
    legend.box.just = "center",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,0),
    legend.spacing.x = unit(1,"cm")
  ) +
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        barwidth = 15,
        barheight = 0.8,
        frame.color = "black",
        frame.linewidth = 0.5
      ),
      size = guide_legend(
        direction = "horizontal",
        nrow = 1
      )
    )
)


combined_plot <- cowplot::plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  nrow = 1,
  ncol = 2,
  labels = c("A","B"),
  rel_widths = c(5,5)
)

combined_plot <- cowplot::plot_grid(
  combined_plot,
  common_legend,
  nrow = 2,
  ncol = 1,
  rel_heights = c(10,1)
)


ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_002.tiff",
  plot = combined_plot,
  dpi = 1000,
  width = 9,
  height = 5,
  units = "in"
)


### 进行差异分析
eset_mc <- qs_read("./analysiswork/cleandata/bulkRNAseq_MC_rm_batch_eset.qs2")
counts <- Biobase::exprs(eset_mc)
metadata <- Biobase::pData(eset_mc)
group <- factor(metadata$group,levels = c("CON","ALS","FTLD"))
batch <- metadata$batch

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(group = group),
  design = ~ group 
)

dds <- DESeq(dds)
resultsNames(dds)
deseq2_res_ALS_vs_CON <- results(dds, contrast = c("group", "ALS", "CON"))
deseq2_res_FTLD_vs_CON <- results(dds, contrast = c("group", "FTLD", "CON"))


dds_result <- deseq2_res_ALS_vs_CON
deseq2_result_table <- dds_result@listData %>% as.data.frame()
rownames(deseq2_result_table) <- dds_result@rownames
deseq2_result_table <- tibble::rownames_to_column(deseq2_result_table,var = "gene")
deseq2_result_table$compare <- "ALS_vs_CON"
deseq2_res_ALS_vs_CON <- deseq2_result_table


dds_result <- deseq2_res_FTLD_vs_CON
deseq2_result_table <- dds_result@listData %>% as.data.frame()
rownames(deseq2_result_table) <- dds_result@rownames
deseq2_result_table <- tibble::rownames_to_column(deseq2_result_table,var = "gene")
deseq2_result_table$compare <- "FTLD_vs_CON"
deseq2_res_FTLD_vs_CON <- deseq2_result_table



# edgeR
dge <- DGEList(counts = counts, group = group)
dge <- calcNormFactors(dge)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
dge <- estimateDisp(dge, design)
fit_edgeR <- glmQLFit(dge, design)
contrast_ALS_vs_CON <- makeContrasts(
  ALS - CON,
  levels = design
)
qlf <- glmQLFTest(fit_edgeR, contrast=contrast_ALS_vs_CON)
edgeR_res_ALS_vs_CON <- qlf
edgeR_res_ALS_vs_CON <- topTags(edgeR_res_ALS_vs_CON,n=Inf,adjust.method = "BH")$table %>% 
  as.data.frame()

## FTLD
contrast_FTLD_vs_CON <- makeContrasts(
  FTLD - CON,
  levels = design
)
qlf <- glmQLFTest(fit_edgeR, contrast=contrast_FTLD_vs_CON)
edgeR_res_FTLD_vs_CON <- qlf
edgeR_res_FTLD_vs_CON <- topTags(edgeR_res_FTLD_vs_CON,n=Inf,adjust.method = "BH")$table %>% 
  as.data.frame()


# limma
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(factor(group))
v <- voom(dge, design = design, plot = TRUE)

fit_limma <- lmFit(v, design)
contrast <- makeContrasts(ALS - CON,levels = design)
fit_limma <- contrasts.fit(fit_limma,contrast)
fit_limma <- eBayes(fit_limma)
limma_res_ALS_vs_CON <- fit_limma
limma_res_ALS_vs_CON <- topTable(limma_res_ALS_vs_CON,coef = 1,number = Inf,adjust.method = "BH")


fit_limma <- lmFit(v, design)
contrast <- makeContrasts(FTLD - CON,levels = design)
fit_limma <- contrasts.fit(fit_limma,contrast)
fit_limma <- eBayes(fit_limma)
limma_res_FTLD_vs_CON <- fit_limma
limma_res_FTLD_vs_CON <- topTable(limma_res_FTLD_vs_CON,coef = 1,number = Inf,adjust.method = "BH")


### 合并结果
select_deseq2_res_ALS_vs_CON <- deseq2_res_ALS_vs_CON[,c(1,3,7,8)]
select_deseq2_res_FTLD_vs_CON <- deseq2_res_FTLD_vs_CON[,c(1,3,7,8)]
colnames(select_deseq2_res_ALS_vs_CON) <- c("gene","log2FC","padj","compare")
colnames(select_deseq2_res_FTLD_vs_CON) <- c("gene","log2FC","padj","compare")

edgeR_res_ALS_vs_CON <- rownames_to_column(edgeR_res_ALS_vs_CON,var = "gene")
edgeR_res_FTLD_vs_CON <- rownames_to_column(edgeR_res_FTLD_vs_CON,var = "gene")
edgeR_res_ALS_vs_CON$compare <- "ALS_vs_CON"
edgeR_res_FTLD_vs_CON$compare <- "FTLD_vs_CON"
select_edgeR_res_ALS_vs_CON <- edgeR_res_ALS_vs_CON[,c(1,2,6,7)]
select_edgeR_res_FTLD_vs_CON <- edgeR_res_FTLD_vs_CON[,c(1,2,6,7)]
colnames(select_edgeR_res_ALS_vs_CON) <- c("gene","log2FC","padj","compare")
colnames(select_edgeR_res_FTLD_vs_CON) <- c("gene","log2FC","padj","compare")

limma_res_ALS_vs_CON <- rownames_to_column(limma_res_ALS_vs_CON,var = "gene")
limma_res_FTLD_vs_CON <- rownames_to_column(limma_res_FTLD_vs_CON,var = "gene")
limma_res_ALS_vs_CON$compare <- "ALS_vs_CON"
limma_res_FTLD_vs_CON$compare <- "FTLD_vs_CON"
select_limma_res_ALS_vs_CON <- limma_res_ALS_vs_CON[,c(1,2,6,8)]
select_limma_res_FTLD_vs_CON <- limma_res_FTLD_vs_CON[,c(1,2,6,8)]
colnames(select_limma_res_ALS_vs_CON) <- c("gene","log2FC","padj","compare")
colnames(select_limma_res_FTLD_vs_CON) <- c("gene","log2FC","padj","compare")

select_deseq2_res_ALS_vs_CON$diffmethod <- "DESeq2"
select_deseq2_res_FTLD_vs_CON$diffmethod <- "DESeq2"
select_edgeR_res_ALS_vs_CON$diffmethod <- "edgeR"
select_edgeR_res_FTLD_vs_CON$diffmethod <- "edgeR"
select_limma_res_ALS_vs_CON$diffmethod <- "limma"
select_limma_res_FTLD_vs_CON$diffmethod <- "limma"

diff_result_list <- list(
  deseq2_res_ALS_vs_CON = select_deseq2_res_ALS_vs_CON,
  deseq2_res_FTLD_vs_CON = select_deseq2_res_FTLD_vs_CON,
  edgeR_res_ALS_vs_CON = select_edgeR_res_ALS_vs_CON,
  edgeR_res_FTLD_vs_CON = select_edgeR_res_FTLD_vs_CON,
  limma_res_ALS_vs_CON = select_limma_res_ALS_vs_CON,
  limma_res_FTLD_vs_CON = select_limma_res_FTLD_vs_CON
)

diff_result_df <- do.call(rbind,diff_result_list)
qs_save(diff_result_df,"./analysiswork/cleandata/bulkseq_MC_diff_result_df.qs2")


## Venn交集分析
diff_result_df <- qs_read("./analysiswork/cleandata/bulkseq_MC_diff_result_df.qs2")

res_ALS_vs_CON_df <- diff_result_df[diff_result_df$compare == "ALS_vs_CON",]


p1 <- multi_volcano(
  plotdata = res_ALS_vs_CON_df,
  group_order = c("DESeq2","edgeR","limma"),
  group_color = c("#6ca3d4","#e8b86c","#e1807e"),
  log2FC_col = log2FC,
  padj_col = padj,
  group_col = diffmethod,
  upcolor = "#b42b22",
  downcolor = "#315a89",
  rect_width = 0.45,
  rect_height = 0.545,
  jitter_size = 0.8,
  jitter_alpha = 0.8,
  plot_title = "ALS_vs_CON",
  plot_title_size = 16,
  trunc_line = 2.5,
  phold = 0.05,
  log2fchold = 0.585
)
p1

res_FTLD_vs_CON_df <- diff_result_df[diff_result_df$compare == "FTLD_vs_CON",]

p2 <- multi_volcano(
  plotdata = res_FTLD_vs_CON_df,
  group_order = c("DESeq2","edgeR","limma"),
  group_color = c("#6ca3d4","#e8b86c","#e1807e"),
  log2FC_col = log2FC,
  padj_col = padj,
  group_col = diffmethod,
  upcolor = "#b42b22",
  downcolor = "#315a89",
  rect_width = 0.45,
  rect_height = 0.545,
  jitter_size = 0.8,
  jitter_alpha = 0.8,
  plot_title = "FTLD_vs_CON",
  plot_title_size = 16,
  trunc_line = 2.5,
  phold = 0.05,
  log2fchold = 0.585
)
p2

combined_plot <- cowplot::plot_grid(
  p1,p2,
  nrow = 2,
  labels = c("A","B"),
  label_fontfamily = "Times New Roman",
  label_size = 14
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_003.tiff",
  plot = combined_plot,
  dpi = 1000,
  width = 6,
  height = 8,
  units = "in"
)


# ALS_vs_CON up_venn_plot
ALS_vs_CON_deseq2_sig_up_genes <- res_ALS_vs_CON_df[res_ALS_vs_CON_df$diffmethod == "DESeq2" & res_ALS_vs_CON_df$log2FC > 0.585 & res_ALS_vs_CON_df$padj < 0.05,"gene"]
ALS_vs_CON_edger_sig_up_genes <- res_ALS_vs_CON_df[res_ALS_vs_CON_df$diffmethod == "edgeR" & res_ALS_vs_CON_df$log2FC > 0.585 & res_ALS_vs_CON_df$padj < 0.05,"gene"]
ALS_vs_CON_limma_sig_up_genes <- res_ALS_vs_CON_df[res_ALS_vs_CON_df$diffmethod == "limma" & res_ALS_vs_CON_df$log2FC > 0.585 & res_ALS_vs_CON_df$padj < 0.05,"gene"]


ALS_vs_CON_up_venn_list <- list(
  DESeq2 = ALS_vs_CON_deseq2_sig_up_genes,
  edgeR = ALS_vs_CON_edger_sig_up_genes,
  limma = ALS_vs_CON_limma_sig_up_genes
)
ALS_vs_CON_up_common_genes <- Reduce(intersect,ALS_vs_CON_up_venn_list)


# ALS_vs_CON down_venn_plot
ALS_vs_CON_deseq2_sig_down_genes <- res_ALS_vs_CON_df[res_ALS_vs_CON_df$diffmethod == "DESeq2" & res_ALS_vs_CON_df$log2FC < -0.585 & res_ALS_vs_CON_df$padj < 0.05,"gene"]
ALS_vs_CON_edger_sig_down_genes <- res_ALS_vs_CON_df[res_ALS_vs_CON_df$diffmethod == "edgeR" & res_ALS_vs_CON_df$log2FC < -0.585 & res_ALS_vs_CON_df$padj < 0.05,"gene"]
ALS_vs_CON_limma_sig_down_genes <- res_ALS_vs_CON_df[res_ALS_vs_CON_df$diffmethod == "limma" & res_ALS_vs_CON_df$log2FC < -0.585 & res_ALS_vs_CON_df$padj < 0.05,"gene"]

# ALS_vs_CON down_venn_plot
ALS_vs_CON_down_venn_list <- list(
  DESeq2 = ALS_vs_CON_deseq2_sig_down_genes,
  edgeR = ALS_vs_CON_edger_sig_down_genes,
  limma = ALS_vs_CON_limma_sig_down_genes
)
ALS_vs_CON_down_common_genes <- Reduce(intersect,ALS_vs_CON_down_venn_list)


# FTLD_vs_CON up_venn_plot
FTLD_vs_CON_deseq2_sig_up_genes <- res_FTLD_vs_CON_df[res_FTLD_vs_CON_df$diffmethod == "DESeq2" & res_FTLD_vs_CON_df$log2FC > 0.585 & res_FTLD_vs_CON_df$padj < 0.05,"gene"]
FTLD_vs_CON_edger_sig_up_genes <- res_FTLD_vs_CON_df[res_FTLD_vs_CON_df$diffmethod == "edgeR" & res_FTLD_vs_CON_df$log2FC > 0.585 & res_FTLD_vs_CON_df$padj < 0.05,"gene"]
FTLD_vs_CON_limma_sig_up_genes <- res_FTLD_vs_CON_df[res_FTLD_vs_CON_df$diffmethod == "limma" & res_FTLD_vs_CON_df$log2FC > 0.585 & res_FTLD_vs_CON_df$padj < 0.05,"gene"]


FTLD_vs_CON_up_venn_list <- list(
  DESeq2 = FTLD_vs_CON_deseq2_sig_up_genes,
  edgeR = FTLD_vs_CON_edger_sig_up_genes,
  limma = FTLD_vs_CON_limma_sig_up_genes
)
FTLD_vs_CON_up_common_genes <- Reduce(intersect,FTLD_vs_CON_up_venn_list)


# FTLD_down_venn_plot
FTLD_vs_CON_deseq2_sig_down_genes <- res_FTLD_vs_CON_df[res_FTLD_vs_CON_df$diffmethod == "DESeq2" & res_FTLD_vs_CON_df$log2FC < -0.585 & res_FTLD_vs_CON_df$padj < 0.05,"gene"]
FTLD_vs_CON_edger_sig_down_genes <- res_FTLD_vs_CON_df[res_FTLD_vs_CON_df$diffmethod == "edgeR" & res_FTLD_vs_CON_df$log2FC < -0.585 & res_FTLD_vs_CON_df$padj < 0.05,"gene"]
FTLD_vs_CON_limma_sig_down_genes <- res_FTLD_vs_CON_df[res_FTLD_vs_CON_df$diffmethod == "limma" & res_FTLD_vs_CON_df$log2FC < -0.585 & res_FTLD_vs_CON_df$padj < 0.05,"gene"]

# FTLD_vs_CON down_venn_plot
FTLD_vs_CON_down_venn_list <- list(
  DESeq2 = FTLD_vs_CON_deseq2_sig_down_genes,
  edgeR = FTLD_vs_CON_edger_sig_down_genes,
  limma = FTLD_vs_CON_limma_sig_down_genes
)
FTLD_vs_CON_down_common_genes <- Reduce(intersect,FTLD_vs_CON_down_venn_list)

sig_venn_list <- list(
  ALS_vs_CON_up = ALS_vs_CON_up_venn_list,
  ALS_vs_CON_down = ALS_vs_CON_down_venn_list,
  FTLD_vs_CON_up = FTLD_vs_CON_up_venn_list,
  FTLD_vs_CON_down = FTLD_vs_CON_down_venn_list
)

qs_save(sig_venn_list,"./analysiswork/cleandata/bulkseq_MC_sig_venn_list.qs2")

sig_venn_list <- qs_read("./analysiswork/cleandata/bulkseq_MC_sig_venn_list.qs2")
## 依次画图
p1 <- ggvenn(
  data = sig_venn_list$ALS_vs_CON_up,
  columns = c("DESeq2","edgeR","limma"),
  show_elements = FALSE,
  show_set_totals = "none",
  show_stats = "c",
  show_counts = TRUE,
  show_percentage = FALSE,
  digits = 1,
  label_sep = ",",
  count_column = NULL,
  show_outside = "auto",
  auto_scale = FALSE,
  fill_color = c("#6ca3d4","#e8b86c","#e1807e"),
  fill_alpha = 0.6,
  stroke_color = "gray95",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = c("#6ca3d4","#e8b86c","#e1807e"),
  set_name_size = 6,
  text_color = "black",
  text_size = 4,
  comma_sep = FALSE,
  padding = 0.2,
  max_elements = 6,
  text_truncate = TRUE
)
p1

p2 <- ggvenn(
  data = sig_venn_list$ALS_vs_CON_down,
  columns = c("DESeq2","edgeR","limma"),
  show_elements = FALSE,
  show_set_totals = "none",
  show_stats = "c",
  show_counts = TRUE,
  show_percentage = FALSE,
  digits = 1,
  label_sep = ",",
  count_column = NULL,
  show_outside = "auto",
  auto_scale = FALSE,
  fill_color = c("#6ca3d4","#e8b86c","#e1807e"),
  fill_alpha = 0.6,
  stroke_color = "gray95",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = c("#6ca3d4","#e8b86c","#e1807e"),
  set_name_size = 6,
  text_color = "black",
  text_size = 4,
  comma_sep = FALSE,
  padding = 0.2,
  max_elements = 6,
  text_truncate = TRUE
)
p2

p3 <- ggvenn(
  data = sig_venn_list$FTLD_vs_CON_up,
  columns = c("DESeq2","edgeR","limma"),
  show_elements = FALSE,
  show_set_totals = "none",
  show_stats = "c",
  show_counts = TRUE,
  show_percentage = FALSE,
  digits = 1,
  label_sep = ",",
  count_column = NULL,
  show_outside = "auto",
  auto_scale = FALSE,
  fill_color = c("#6ca3d4","#e8b86c","#e1807e"),
  fill_alpha = 0.6,
  stroke_color = "gray95",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = c("#6ca3d4","#e8b86c","#e1807e"),
  set_name_size = 6,
  text_color = "black",
  text_size = 4,
  comma_sep = FALSE,
  padding = 0.2,
  max_elements = 6,
  text_truncate = TRUE
)
p3


p4 <- ggvenn(
  data = sig_venn_list$FTLD_vs_CON_down,
  columns = c("DESeq2","edgeR","limma"),
  show_elements = FALSE,
  show_set_totals = "none",
  show_stats = "c",
  show_counts = TRUE,
  show_percentage = FALSE,
  digits = 1,
  label_sep = ",",
  count_column = NULL,
  show_outside = "auto",
  auto_scale = FALSE,
  fill_color = c("#6ca3d4","#e8b86c","#e1807e"),
  fill_alpha = 0.6,
  stroke_color = "gray95",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = c("#6ca3d4","#e8b86c","#e1807e"),
  set_name_size = 6,
  text_color = "black",
  text_size = 4,
  comma_sep = FALSE,
  padding = 0.2,
  max_elements = 6,
  text_truncate = TRUE
)
p4

combined_plot <- cowplot::plot_grid(
  p1,p2,p3,p4,
  nrow = 1,
  byrow = TRUE,
  labels = c("A","B","C","D"),
  label_fontfamily = "Times New Roman",
  label_size = 14
)


ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_004.tiff",
  plot = combined_plot,
  dpi = 500,
  width = 15,
  height = 4,
  units = "in"
)

shared_list <- list(
  ALS_shared_up = Reduce(intersect,sig_venn_list$ALS_vs_CON_up),
  ALS_shared_down = Reduce(intersect,sig_venn_list$ALS_vs_CON_down),
  FTLD_shared_up = Reduce(intersect,sig_venn_list$FTLD_vs_CON_up),
  FTLD_shared_down = Reduce(intersect,sig_venn_list$FTLD_vs_CON_down)
)

p <- ggvenn(
  data = shared_list,
  columns = c("ALS_shared_up","ALS_shared_down","FTLD_shared_up","FTLD_shared_down"),
  show_elements = FALSE,
  show_set_totals = "none",
  show_stats = "c",
  show_counts = TRUE,
  show_percentage = FALSE,
  digits = 1,
  label_sep = ",",
  count_column = NULL,
  show_outside = "auto",
  auto_scale = FALSE,
  fill_color = c("#72b063","#badab5","#4a5f7e","#719aac"),
  fill_alpha = 0.6,
  stroke_color = "gray95",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = c("#72b063","#badab5","#4a5f7e","#719aac"),
  set_name_size = 6,
  text_color = "black",
  text_size = 4,
  comma_sep = FALSE,
  padding = 0.2,
  max_elements = 6,
  text_truncate = TRUE
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_005.tiff",
  plot = p,
  dpi = 500,
  width = 7,
  height = 7,
  units = "in"
)

shared_up_list <- list(
  ALS_shared_up = Reduce(intersect,sig_venn_list$ALS_vs_CON_up),
  FTLD_shared_up = Reduce(intersect,sig_venn_list$FTLD_vs_CON_up)
)

shared_down_list <- list(
  ALS_shared_down = Reduce(intersect,sig_venn_list$ALS_vs_CON_down),
  FTLD_shared_down = Reduce(intersect,sig_venn_list$FTLD_vs_CON_down)
)





# alsod_up_upset_plot
load("./analysiswork/helpdata/ALSOD_Genes.RData")

sig_venn_list <- qs_read("./analysiswork/cleandata/bulkseq_MC_sig_venn_list.qs2")
als_alsod_up_venn_list <- list(
  DESeq2 = sig_venn_list$ALS_vs_CON_up$DESeq2,
  edgeR = sig_venn_list$ALS_vs_CON_up$edgeR,
  limma = sig_venn_list$ALS_vs_CON_up$limma,
  ALSOD = ALSOD_Genes %>% unlist() %>% as.vector()
)

als_alsod_down_venn_list <- list(
  DESeq2 = sig_venn_list$ALS_vs_CON_down$DESeq2,
  edgeR = sig_venn_list$ALS_vs_CON_down$edgeR,
  limma = sig_venn_list$ALS_vs_CON_down$limma,
  ALSOD = ALSOD_Genes %>% unlist() %>% as.vector()
)

fltd_alsod_up_venn_list <- list(
  DESeq2 = sig_venn_list$FTLD_vs_CON_up$DESeq2,
  edgeR = sig_venn_list$FTLD_vs_CON_up$edgeR,
  limma = sig_venn_list$FTLD_vs_CON_up$limma,
  ALSOD = ALSOD_Genes %>% unlist() %>% as.vector()
)

fltd_alsod_down_venn_list <- list(
  DESeq2 = sig_venn_list$FTLD_vs_CON_down$DESeq2,
  edgeR = sig_venn_list$FTLD_vs_CON_down$edgeR,
  limma = sig_venn_list$FTLD_vs_CON_down$limma,
  ALSOD = ALSOD_Genes %>% unlist() %>% as.vector()
)

base::Reduce(intersect,als_alsod_up_venn_list)
base::Reduce(intersect,als_alsod_down_venn_list)
base::Reduce(intersect,fltd_alsod_up_venn_list)
base::Reduce(intersect,fltd_alsod_down_venn_list)

p1 <- ggvenn(
  data = als_alsod_up_venn_list,
  columns = c("DESeq2","edgeR","limma","ALSOD"),
  show_elements = FALSE,
  show_set_totals = "none",
  show_stats = "c",
  show_counts = TRUE,
  show_percentage = FALSE,
  digits = 1,
  label_sep = ",",
  count_column = NULL,
  show_outside = "auto",
  auto_scale = FALSE,
  fill_color = c("#6ca3d4","#e8b86c","#e1807e","#b7b7eb"),
  fill_alpha = 0.6,
  stroke_color = "gray95",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = c("#6ca3d4","#e8b86c","#e1807e","#b7b7eb"),
  set_name_size = 5,
  text_color = "black",
  text_size = 4,
  comma_sep = FALSE,
  padding = 0.2,
  max_elements = 6,
  text_truncate = TRUE
)


p2 <- ggvenn(
  data = als_alsod_down_venn_list,
  columns = c("DESeq2","edgeR","limma","ALSOD"),
  show_elements = FALSE,
  show_set_totals = "none",
  show_stats = "c",
  show_counts = TRUE,
  show_percentage = FALSE,
  digits = 1,
  label_sep = ",",
  count_column = NULL,
  show_outside = "auto",
  auto_scale = FALSE,
  fill_color = c("#6ca3d4","#e8b86c","#e1807e","#b7b7eb"),
  fill_alpha = 0.6,
  stroke_color = "gray95",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = c("#6ca3d4","#e8b86c","#e1807e","#b7b7eb"),
  set_name_size = 5,
  text_color = "black",
  text_size = 4,
  comma_sep = FALSE,
  padding = 0.2,
  max_elements = 6,
  text_truncate = TRUE
)

p3 <- ggvenn(
  data = fltd_alsod_up_venn_list,
  columns = c("DESeq2","edgeR","limma","ALSOD"),
  show_elements = FALSE,
  show_set_totals = "none",
  show_stats = "c",
  show_counts = TRUE,
  show_percentage = FALSE,
  digits = 1,
  label_sep = ",",
  count_column = NULL,
  show_outside = "auto",
  auto_scale = FALSE,
  fill_color = c("#6ca3d4","#e8b86c","#e1807e","#b7b7eb"),
  fill_alpha = 0.6,
  stroke_color = "gray95",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = c("#6ca3d4","#e8b86c","#e1807e","#b7b7eb"),
  set_name_size = 5,
  text_color = "black",
  text_size = 4,
  comma_sep = FALSE,
  padding = 0.2,
  max_elements = 6,
  text_truncate = TRUE
)

p4 <- ggvenn(
  data = fltd_alsod_down_venn_list,
  columns = c("DESeq2","edgeR","limma","ALSOD"),
  show_elements = FALSE,
  show_set_totals = "none",
  show_stats = "c",
  show_counts = TRUE,
  show_percentage = FALSE,
  digits = 1,
  label_sep = ",",
  count_column = NULL,
  show_outside = "auto",
  auto_scale = FALSE,
  fill_color = c("#6ca3d4","#e8b86c","#e1807e","#b7b7eb"),
  fill_alpha = 0.6,
  stroke_color = "gray95",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = c("#6ca3d4","#e8b86c","#e1807e","#b7b7eb"),
  set_name_size = 5,
  text_color = "black",
  text_size = 4,
  comma_sep = FALSE,
  padding = 0.2,
  max_elements = 6,
  text_truncate = TRUE
)


combined_plot <- cowplot::plot_grid(
  p1,p2,p3,p4,
  nrow = 1,
  byrow = TRUE,
  labels = c("A","B","C","D"),
  label_fontfamily = "Times New Roman",
  label_size = 14
)


ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_006.tiff",
  plot = combined_plot,
  dpi = 500,
  width = 15,
  height = 4,
  units = "in"
)

shared_list_with_alsod <- list(
  ALS_up = shared_list$ALS_shared_up,
  ALS_down = shared_list$ALS_shared_down,
  FTLD_up = shared_list$FTLD_shared_up,
  FTLD_down = shared_list$FTLD_shared_down,
  ALSOD = ALSOD_Genes %>% unlist() %>% as.vector()
)

Reduce(intersect,gene_list)

library(UpSetR)
tiff("./analysiswork/result/plot/bulkRNA_MC_combined_007.tiff",width = 8, height = 6, units = "in",res = 1000)
upset(
  fromList(shared_list_with_alsod),
  keep.order = F,
  order.by = "freq",
  main.bar.color = "#e3625d",
  sets.bar.color = c(ALS_up = "#e2735f",ALS_down = "#f5b99d",FTLD_up = "#79b4d7",FTLD_down = "#c6e2ed",ALSOD = "#a87ec0")
)
dev.off()




###
sig_venn_list <- qs_read("./analysiswork/cleandata/bulkseq_MC_sig_venn_list.qs2")
als_up_shared_genes <- Reduce(intersect,sig_venn_list$ALS_vs_CON_up)
als_down_shared_genes <- Reduce(intersect,sig_venn_list$ALS_vs_CON_down)
fltd_up_shared_genes <- Reduce(intersect,sig_venn_list$FTLD_vs_CON_up)
fltd_down_shared_genes <- Reduce(intersect,sig_venn_list$FTLD_vs_CON_down)

als_up_kegg_result <- KEGGflow(als_up_shared_genes,qvaluehold = 1)
p1 <- KEGGbar(als_up_kegg_result,plot_title = "ALS Up Shared Genes KEGG Result",numsize = 2)

als_down_kegg_result <- KEGGflow(als_down_shared_genes,qvaluehold = 1)
p2 <- KEGGbar(als_down_kegg_result,plot_title = "ALS Down Shared Genes KEGG Result",numsize = 2)

fltd_up_kegg_result <- KEGGflow(fltd_up_shared_genes,qvaluehold = 1)
p3 <- KEGGbar(fltd_up_kegg_result,plot_title = "FTLD Up Shared Genes KEGG Result",numsize = 2)

fltd_down_kegg_result <- KEGGflow(fltd_down_shared_genes,qvaluehold = 1)
p4 <- KEGGbar(fltd_down_kegg_result,plot_title = "FTLD Down Shared Genes KEGG Result",numsize = 2)

p1 <- p1 + theme(legend.position = "none",aspect.ratio = 1)
p2 <- p2 + theme(legend.position = "none",aspect.ratio = 1)
p3 <- p3 + theme(legend.position = "none",aspect.ratio = 1)
p4 <- p4 + theme(legend.position = "none",aspect.ratio = 1)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_008.tiff",
  plot = p1,
  dpi = 500,
  width = 10,
  height = 6,
  units = "in"
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_009.tiff",
  plot = p2,
  dpi = 500,
  width = 10,
  height = 6,
  units = "in"
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0010.tiff",
  plot = p3,
  dpi = 500,
  width = 10,
  height = 6,
  units = "in"
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0011.tiff",
  plot = p4,
  dpi = 500,
  width = 10,
  height = 6,
  units = "in"
)



## go
als_up_go_result <- GOflow(als_up_shared_genes,qvaluehold = 1)
p1 <- GObar(als_up_go_result,plot_title = "ALS Up Shared Genes GO Result")

als_down_go_result <- GOflow(als_down_shared_genes,qvaluehold = 1)
p2 <- GObar(als_down_go_result,plot_title = "ALS Down Shared Genes GO Result")

fltd_up_go_result <- GOflow(fltd_up_shared_genes,qvaluehold = 1)
p3 <- GObar(fltd_up_go_result,plot_title = "FTLD Up Shared Genes GO Result")

fltd_down_go_result <- GOflow(fltd_down_shared_genes,qvaluehold = 1)
p4 <- GObar(fltd_down_go_result,plot_title = "FTLD Down Shared Genes GO Result")

p1 <- p1 + theme(legend.position = "none",aspect.ratio = 1)
p2 <- p2 + theme(legend.position = "none",aspect.ratio = 1) ## 没有富集到
p3 <- p3 + theme(legend.position = "none",aspect.ratio = 1)
p4 <- p4 + theme(legend.position = "none",aspect.ratio = 1)


ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0012.tiff",
  plot = p1,
  dpi = 500,
  width = 10,
  height = 5,
  units = "in"
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0013.tiff",
  plot = p2,
  dpi = 500,
  width = 10,
  height = 5,
  units = "in"
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0014.tiff",
  plot = p3,
  dpi = 500,
  width = 10,
  height = 5,
  units = "in"
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0015.tiff",
  plot = p4,
  dpi = 500,
  width = 10,
  height = 5,
  units = "in"
)







### GSEA
diff_result_df <- qs_read("./analysiswork/cleandata/bulkseq_MC_diff_result_df.qs2")
select_diff <- diff_result_df[diff_result_df$compare == "ALS_vs_CON" & diff_result_df$diffmethod == "DESeq2",]
data_sort <- select_diff %>%
  arrange(desc(log2FC))

gene_list <- data_sort$log2FC
names(gene_list) <- data_sort$gene

### gsea_KEGG
id_map <- bitr(
  geneID = names(gene_list), 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)

geneList_entrez <- gene_list[id_map$SYMBOL]
names(geneList_entrez) <- id_map$ENTREZID

gsea_kegg <- gseKEGG(
  geneList = geneList_entrez,
  organism = "hsa",
  keyType = "kegg",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea"
)

gsea_kegg_up <- gsea_kegg[gsea_kegg@result$NES > 0,]
gsea_kegg_down <- gsea_kegg[gsea_kegg@result$NES < 0,]

gsea_kegg_readable <- setReadable(
  x = gsea_kegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)



gsea_kegg_sig_hsanum <- gsea_kegg_up[["ID"]]
all_mc_als_gsea_upsig_kegg <- gsea_kegg_sig_hsanum
qs_save(all_mc_als_gsea_upsig_kegg,"./analysiswork/cleandata/all_mc_als_gsea_upsig_kegg.qs2")
gsea_kegg_sig_hsanum <- gsea_kegg_down[["ID"]]
all_mc_als_gsea_downsig_kegg <- gsea_kegg_sig_hsanum
qs_save(all_mc_als_gsea_downsig_kegg,"./analysiswork/cleandata/all_mc_als_gsea_downsig_kegg.qs2")


select_diff <- diff_result_df[diff_result_df$compare == "FTLD_vs_CON" & diff_result_df$diffmethod == "DESeq2",]
data_sort <- select_diff %>%
  arrange(desc(log2FC))

gene_list <- data_sort$log2FC
names(gene_list) <- data_sort$gene

### gsea_KEGG
id_map <- bitr(
  geneID = names(gene_list), 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)

geneList_entrez <- gene_list[id_map$SYMBOL]
names(geneList_entrez) <- id_map$ENTREZID

gsea_kegg <- gseKEGG(
  geneList = geneList_entrez,
  organism = "hsa",
  keyType = "kegg",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea"
)


gsea_kegg_up <- gsea_kegg[gsea_kegg@result$NES > 0,]
gsea_kegg_down <- gsea_kegg[gsea_kegg@result$NES < 0,]

gsea_kegg_readable <- setReadable(
  x = gsea_kegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)



gsea_kegg_sig_hsanum <- gsea_kegg_up[["ID"]]
all_mc_ftld_gsea_upsig_kegg <- gsea_kegg_sig_hsanum
qs_save(all_mc_ftld_gsea_upsig_kegg,"./analysiswork/cleandata/all_mc_ftld_gsea_upsig_kegg.qs2")
gsea_kegg_sig_hsanum <- gsea_kegg_down[["ID"]]
all_mc_ftld_gsea_downsig_kegg <- gsea_kegg_sig_hsanum
qs_save(all_mc_ftld_gsea_downsig_kegg,"./analysiswork/cleandata/all_mc_ftld_gsea_downsig_kegg.qs2")




als_up_sig_kegg_hsanum <- rownames(als_up_kegg_result[als_up_kegg_result$p.adjust < 0.05,])
als_down_sig_kegg_hsanum <- rownames(als_down_kegg_result[als_down_kegg_result$p.adjust < 0.05,])
ftld_up_sig_kegg_hsanum <- rownames(fltd_up_kegg_result[fltd_up_kegg_result$p.adjust < 0.05,])
ftld_down_sig_kegg_hsanum <- rownames(fltd_down_kegg_result[fltd_down_kegg_result$p.adjust < 0.05,])

shared_mc_als_kegg_hsanum <- intersect(gsea_kegg_sig_hsanum,unique(c(als_up_sig_kegg_hsanum,als_down_sig_kegg_hsanum)))
shared_mc_ftld_kegg_hsanum <- intersect(gsea_kegg_sig_hsanum,unique(c(ftld_up_sig_kegg_hsanum),ftld_down_sig_kegg_hsanum))


gsea_plist <- list()
for (hsaname in shared_kegg_hsanum) {
  p <- gseaNb(
    object = gsea_kegg_readable,
    geneSetID = hsaname,
    newGsea = T,
    rmHt = F,
    addPoint = T,
    addPval = T,
    pvalX = 0.8,
    pvalY = 0.8,
    pCol = 'black',
    newHtCol = c("#a244a1","white", "#668099"),
    pHjust = 0,
    subPlot = 2,
    termWidth = 50,
    legend.position = c(0.8,0.8)
  )
  gsea_plist[[hsaname]] <- p
}


combined_plot <- plot_grid(
  plotlist = gsea_plist,
  ncol = 1,
  byrow = TRUE
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0016.tiff",
  plot = combined_plot,
  dpi = 500,
  width = 10,
  height = 10,
  units = "in"
)



### gsea_GO
gsea_go <- gseGO(
  geneList = gene_list,
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)


gsea_go_readable <- setReadable(
  x = gsea_go,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL"
)


gsea_go_sig_hsanum <- gsea_go@result[["ID"]]
als_up_sig_go_hsanum <- rownames(als_up_go_result[als_up_go_result$p.adjust < 0.05,])
als_down_sig_go_hsanum <- rownames(als_down_go_result[als_down_go_result$p.adjust < 0.05,])
ftld_up_sig_go_hsanum <- rownames(fltd_up_go_result[fltd_up_go_result$p.adjust < 0.05,])
ftld_down_sig_go_hsanum <- rownames(fltd_down_go_result[fltd_down_go_result$p.adjust < 0.05,])


shared_go_als_up_num <- intersect(gsea_go_sig_hsanum,als_up_sig_go_hsanum)
shared_go_als_down_num <- intersect(gsea_go_sig_hsanum,als_down_sig_go_hsanum)
shared_go_ftld_up_num <- intersect(gsea_go_sig_hsanum,ftld_up_sig_go_hsanum)
shared_go_ftld_down_num <- intersect(gsea_go_sig_hsanum,ftld_down_sig_go_hsanum)

shared_go_num_list <- list(
  als_up = shared_go_als_up_num,
  als_down = shared_go_als_down_num,
  ftld_up = shared_go_ftld_up_num,
  ftld_down = shared_go_ftld_down_num
)

qs_save(shared_go_num_list,"./analysiswork/cleandata/shared_mc_go_num.qs2")

gsea_plist <- list()
for (hsaname in shared_go_num_list$als_up) {
  p <- gseaNb(
    object = gsea_go_readable,
    geneSetID = hsaname,
    newGsea = F,
    lineSize = 1,
    rmHt = T,
    addPoint = T,
    addPval = T,
    htCol = c("#5861ac","#f28080"),
    pvalX = 0.6,
    pvalY = 0.8,
    pCol = 'black',
    pHjust = 0,
    subPlot = 3,
    termWidth = 50,
    legend.position = "top",
    curveCol = RColorBrewer::brewer.pal(3,"Set1")
  )
  gsea_plist[[hsaname]] <- p
}

combined_plot1 <- plot_grid(
  plotlist = gsea_plist,
  ncol = 2,
  byrow = TRUE
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0017.tiff",
  plot = combined_plot1,
  dpi = 500,
  width = 12,
  height = 5,
  units = "in"
)




p2 <- gseaNb(
  object = gsea_go_readable,
  geneSetID = shared_go_num_list$als_down,
  newGsea = F,
  lineSize = 1,
  rmHt = T,
  addPoint = T,
  addPval = F,
  htCol = c("#5861ac","#f28080"),
  pvalX = 0.1,
  pvalY = 0.9,
  pCol = 'black',
  pHjust = 0,
  subPlot = 3,
  termWidth = 50,
  legend.position = "top",
  curveCol = RColorBrewer::brewer.pal(6,"Set1")
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0018.tiff",
  plot = p2,
  dpi = 500,
  width = 12,
  height = 10,
  units = "in"
)



p3 <- gseaNb(
  object = gsea_go_readable,
  geneSetID = shared_go_num_list$ftld_up,
  newGsea = F,
  lineSize = 1,
  rmHt = T,
  addPoint = T,
  addPval = F,
  htCol = c("#5861ac","#f28080"),
  pvalX = 0.1,
  pvalY = 0.9,
  pCol = 'black',
  pHjust = 0,
  subPlot = 3,
  termWidth = 50,
  legend.position = "top",
  curveCol = RColorBrewer::brewer.pal(11,"Set3")
)

ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0019.tiff",
  plot = p3,
  dpi = 500,
  width = 15,
  height = 10,
  units = "in"
)

p4 <- gseaNb(
  object = gsea_go_readable,
  geneSetID = shared_go_num_list$ftld_down,
  newGsea = F,
  lineSize = 1,
  rmHt = T,
  addPoint = T,
  addPval = F,
  htCol = c("#5861ac","#f28080"),
  pvalX = 0.1,
  pvalY = 0.9,
  pCol = 'black',
  pHjust = 0,
  subPlot = 3,
  termWidth = 50,
  legend.position = "top",
  curveCol = RColorBrewer::brewer.pal(12,"Set3")
)


ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_combined_0020.tiff",
  plot = p4,
  dpi = 500,
  width = 15,
  height = 10,
  units = "in"
)



### 4.GSVA
library(GSVA)
library(limma)
library(GSEABase)
library(msigdbr)
library(pheatmap)
library(ggplot2)
library(Biobase)

eset <- qs_read("./analysiswork/cleandata/bulkRNAseq_MC_rm_batch_eset.qs2")
metadata <- Biobase::pData(eset)
expr_matrix <- Biobase::exprs(eset)
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

expr_matrix <- expr_matrix[shared_gens,]

gsva_param <- gsvaParam(
  exprData = expr_matrix,
  geneSets = gene_sets,
  minSize = 5,
  maxSize = 500,
  kcdf = "Gaussian"
)

gsva_results <- gsva(gsva_param)




library(ggplot2)
library(tidyverse)
library(Biobase)
library(ggpubr)
library(gghalves)
library(ggsignif)
### 绘制火山图

### 交集基因箱线图
eset <- qs_read("./analysiswork/")
shared_genes <- Reduce(intersect,venn_genes_list$alsod_up_venn_list)
count_matrix <- exprs(eset)
shared_count_matrix <- count_matrix[shared_genes,]
shared_count_matrix[is.na(shared_count_matrix)] <- 0
logcpm <- cpm(shared_count_matrix, log = TRUE, prior.count = 1)
t_logcpm <- logcpm %>% t()
group <- pData(eset)$mergegroup %>% as.data.frame()
colnames(group) <- "group"
final_df <- cbind(t_logcpm %>% as.data.frame(),group)
wide_final_df <- final_df %>% rownames_to_column(var = "sample")
long_final_df <- wide_final_df %>% 
  pivot_longer(
    cols = -c(sample,group),
    names_to = "gene",
    values_to = "expression"
  ) %>% 
  select(
    sample,group,expression,everything()
  )

long_final_df <- long_final_df %>% 
  filter(group != "Other")
plotdata <- long_final_df
group_color <- setNames(c("#E41A1C","#984EA3"),c("ALS","Control"))

ggplot(plotdata) +
  stat_boxplot(
    mapping = aes(x = gene,y = expression, fill = group),
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  geom_boxplot(
    mapping = aes(x = gene,y = expression, fill = group), 
    errorbar.length = 0.2,
    width = 0.4,
    alpha = 0.7,
    position = position_dodge(width = 0.6),
    outlier.color = NA
  ) +
  scale_fill_manual(values = group_color) +
  scale_y_continuous(expand = expansion(0,0.1)) +
  labs(x = "", y = "Logcpm of counts", title = "") +
  theme(
    plot.background = element_blank(),
    plot.title = element_text(),
    panel.background = element_blank(),
    
    axis.line.x = element_line(colour = "black",linewidth = 1),
    axis.line.y = element_line(color = "black",linewidth = 1),
    
    axis.text.x = element_text(size = 12,face = "bold",colour = "black",angle = 45,hjust = 1,vjust = 1),
    axis.text.y = element_text(size = 12,face = "bold",colour = "black"),
    
    axis.title.x = element_text(size = 20,face = "bold",colour = "black"),
    axis.title.y = element_text(size = 20,face = "bold",colour = "black"),
    axis.ticks.length = unit(3,"mm"),
    axis.ticks = element_line(linewidth = 1),
    legend.position = c(0.5,1),
    panel.grid.major = element_line(color = "grey80"),
    legend.title = element_text(hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(
      ncol = 2
    )
  ) +
  coord_equal(ratio = 1/3)



### 交集基因热图
library(pheatmap)
