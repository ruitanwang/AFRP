#### 

FC_diff_genes <- qs_read("./analysiswork/cleandata/snRNA_VAL_FC_allcelltype_diff_genes.qs2")
C9ALS_FC_diff_genes <- FC_diff_genes[FC_diff_genes$compare == "C9ALS_vs_CON",]
C9FTLD_FC_diff_genes <- FC_diff_genes[FC_diff_genes$compare == "C9FTLD_vs_CON",]
SALS_FC_diff_genes <- FC_diff_genes[FC_diff_genes$compare == "SALS_vs_CON",]
SFTLD_FC_diff_genes <- FC_diff_genes[FC_diff_genes$compare == "SFTLD_vs_CON",]


## Oligo
C9ALS_oligo <- C9ALS_FC_diff_genes[C9ALS_FC_diff_genes$celltype == "Oligo",]
C9ALS_oligo <- C9ALS_oligo %>% 
  arrange(desc(avg_log2FC))

gene_list <- C9ALS_oligo$avg_log2FC
names(gene_list) <- C9ALS_oligo$gene

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

C9ALS_oligo_gsea_go_result <- gsea_go[,c("Description","NES","qvalue")]

SALS_oligo <- SALS_FC_diff_genes[SALS_FC_diff_genes$celltype == "Oligo",]
SALS_oligo <- SALS_oligo %>% 
  arrange(desc(avg_log2FC))

gene_list <- SALS_oligo$avg_log2FC
names(gene_list) <- SALS_oligo$gene

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

SALS_oligo_gsea_go_result <- gsea_go[,c("Description","NES","qvalue")]


C9FTLD_oligo <- C9FTLD_FC_diff_genes[C9FTLD_FC_diff_genes$celltype == "Oligo",]
C9FTLD_oligo <- C9FTLD_oligo %>% 
  arrange(desc(avg_log2FC))

gene_list <- C9FTLD_oligo$avg_log2FC
names(gene_list) <- C9FTLD_oligo$gene

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
C9FTLD_oligo_gsea_go_result <- gsea_go[,c("Description","NES","qvalue")]


SFTLD_oligo <- SFTLD_FC_diff_genes[SFTLD_FC_diff_genes$celltype == "Oligo",]
SFTLD_oligo <- SFTLD_oligo %>% 
  arrange(desc(avg_log2FC))

gene_list <- SFTLD_oligo$avg_log2FC
names(gene_list) <- SFTLD_oligo$gene

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

SFTLD_oligo_gsea_go_result <- gsea_go[,c("Description","NES","qvalue")]

oligo_gsea_go_result <- list(
  C9ALS = C9ALS_oligo_gsea_go_result,
  C9FTLD = C9FTLD_oligo_gsea_go_result,
  SALS = SALS_oligo_gsea_go_result,
  SFTLD = SFTLD_oligo_gsea_go_result
)
oligo_gsea_go_result <- do.call(rbind,oligo_gsea_go_result)
oligo_gsea_go_result$disease <- sub("\\..*","",rownames(oligo_gsea_go_result))

shared_go <- oligo_gsea_go_result %>% 
  add_count(Description) %>% 
  filter(n == 4)

plotdata <- shared_go
ggplot(data = plotdata,mapping = aes(x = disease,y = Description)) +
  geom_point(mapping = aes(size = -log10(qvalue),color = NES)) +
  scale_color_gradientn(colors = c("#f4eff7","#83559e")) +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black")
  )





seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_seurat_annofinal.qs2")
all_df <- qs_read("./analysiswork/cleandata/snRNA_VAL_MC_allcelltype_diff_genes.qs2")

all_df_sig <- all_df %>% 
  filter(p_val_adj < 0.001 & abs(avg_log2FC) > 1)

result_longer <- all_df_sig %>% 
  dplyr::count(compare,celltype,updown)

colnames(result_longer) <- c("compare","sub_cell_type","updown","num")

cellinfo <- seuobj@meta.data[,c("cell_type","sub_cell_type")] %>% unique()
colnames(cellinfo) <- c("main_cell_type","sub_cell_type")

final_df <- result_longer %>% 
  left_join(cellinfo,by = "sub_cell_type")


multigroup_bar <- function(
    plotdata = df,
    plot_title = "",
    compare_color = "",
    x_title = "",
    x_title_size = 15,
    sub_celltype_color = all_sub_type_color,
    main_celltype_color = main_cell_type_color
) {
  
  all_maintype_color <- main_celltype_color[names(main_celltype_color) %in% unique(plotdata$main_cell_type)]
  all_subtype_color <- sub_celltype_color[names(sub_celltype_color) %in% unique(plotdata$sub_cell_type)]
  plotdata$main_cell_type <- factor(plotdata$main_cell_type,levels = names(all_maintype_color))
  plotdata$sub_cell_type <- factor(plotdata$sub_cell_type,levels = names(all_subtype_color))
  
  p <- ggplot(
    data = plotdata,
    mapping = aes(x = num, y = interaction(sub_cell_type,main_cell_type))
  ) +
    geom_vline(xintercept = 0,color = "gray50",linewidth = 0.2,linetype = 2) +
    geom_col(
      mapping = aes(group = updown,fill = updown),
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    scale_x_continuous(
      labels = function(x) abs(x)
    ) +
    scale_fill_manual(
      values = c("UP" = "#f09ba0","DOWN" = "#9bbbe1")
    ) +
    facet_wrap(
      ~ compare,
      scales = "free_x",
      nrow = 1
    ) +
    labs(
      x = x_title,
      y = "",
      title = plot_title
    ) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.background = element_blank(),
      plot.title = element_text(),
      panel.grid.major = element_line(color = "gray90"),
      
      panel.background = element_blank(),
      
      axis.line.x = element_line(colour = "black",linewidth = 1),
      axis.line.y = element_line(color = "black",linewidth = 1),
      
      axis.text.x = element_text(size = 12,face = "bold",colour = "black",angle = 45,hjust = 1),
      axis.text.y = element_text(size = 12,face = "bold"),
      
      axis.title.x = element_text(size = x_title_size,face = "bold",colour = "black"),
      axis.title.y = element_text(size = 20,face = "bold",colour = "black"),
      
      strip.background = element_blank(),
      strip.text = element_text(color = "gray30",face = "bold",size = 10)
    ) +
    guides(
      y = legendry::guide_axis_nested(
        key = legendry::key_range_auto(sep = "\\."),
        levels_text = list(
          element_text(color = all_subtype_color,face = "bold"),
          element_text(color = all_maintype_color,face = "bold")
        )
      )
    ) +
    theme(
      legend.position = "none"
    )
  
  return(p)
}

multigroup_bar(final_df,compare_color = up_compare_color,x_title = "Regulate gene number") +
  theme(legend.position = "right")



all_df <- qs_read("./analysiswork/cleandata/snRNA_VAL_MC_allcelltype_diff_genes.qs2")
all_df_sig <- all_df %>% 
  filter(p_val_adj < 0.001 & abs(avg_log2FC) > 1)
load("./analysiswork/helpdata/ALSOD_Genes.RData")
ALSOD_Genes <- ALSOD_Genes %>% unlist() %>% as.vector()

all_df_sig_select_ALSOD_Genes <- all_df_sig[all_df_sig$gene %in% ALSOD_Genes,]
result_longer <- all_df_sig_select_ALSOD_Genes %>% 
  dplyr::count(compare,celltype,updown)

colnames(result_longer) <- c("compare","sub_cell_type","updown","num")

cellinfo <- seuobj@meta.data[,c("cell_type","sub_cell_type")] %>% unique()
colnames(cellinfo) <- c("main_cell_type","sub_cell_type")

final_df <- result_longer %>% 
  left_join(cellinfo,by = "sub_cell_type")


multigroup_bar(final_df,compare_color = up_compare_color,x_title = "Regulate gene number") +
  theme(legend.position = "right")

unique(all_df_sig_select_ALSOD_Genes$gene)

all_df_sig_select_ALSOD_Genes_C9ALS <- all_df_sig_select_ALSOD_Genes[all_df_sig_select_ALSOD_Genes$compare == "C9ALS_vs_CON",]
gene_count <- all_df_sig_select_ALSOD_Genes_C9ALS %>% 
  dplyr::count(gene,sort = TRUE)
gene_count$gene <- factor(gene_count$gene,levels = gene_count$gene)
p1gene <- gene_count$gene
p1 <- ggplot(gene_count,aes(x = gene,y = n,group = 1)) +
  geom_line(color = "#377eb8",linewidth = 1) +
  geom_point(size = 2,color = "#e41a1c") +
  theme_bw() +
  theme(
    text = element_text(color = "black",family = "Times New Roman",face = "bold"),
    axis.text.x = element_text(angle = 45,hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "Gene",
    y = "Count",
    title = "C9ALS_vs_CON"
  )


all_df_sig_select_ALSOD_Genes_C9FTLD <- all_df_sig_select_ALSOD_Genes[all_df_sig_select_ALSOD_Genes$compare == "C9FTLD_vs_CON",]
gene_count <- all_df_sig_select_ALSOD_Genes_C9FTLD %>% 
  count(gene,sort = TRUE)
gene_count$gene <- factor(gene_count$gene,levels = gene_count$gene)
p2gene <- gene_count$gene
p2 <- ggplot(gene_count,aes(x = gene,y = n,group = 1)) +
  geom_line(color = "#377eb8",linewidth = 1) +
  geom_point(size = 2,color = "#e41a1c") +
  theme_bw() +
  theme(
    text = element_text(color = "black",family = "Times New Roman",face = "bold"),
    axis.text.x = element_text(angle = 45,hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "Gene",
    y = "Count",
    title = "C9FTLD_vs_CON"
  )

all_df_sig_select_ALSOD_Genes_SALS <- all_df_sig_select_ALSOD_Genes[all_df_sig_select_ALSOD_Genes$compare == "SALS_vs_CON",]
gene_count <- all_df_sig_select_ALSOD_Genes_SALS %>% 
  count(gene,sort = TRUE)
gene_count$gene <- factor(gene_count$gene,levels = gene_count$gene)
p3gene <- gene_count$gene
p3 <- ggplot(gene_count,aes(x = gene,y = n,group = 1)) +
  geom_line(color = "#377eb8",linewidth = 1) +
  geom_point(size = 2,color = "#e41a1c") +
  theme_bw() +
  theme(
    text = element_text(color = "black",family = "Times New Roman",face = "bold"),
    axis.text.x = element_text(angle = 45,hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "Gene",
    y = "Count",
    title = "SALS_vs_CON"
  )

all_df_sig_select_ALSOD_Genes_SFTLD <- all_df_sig_select_ALSOD_Genes[all_df_sig_select_ALSOD_Genes$compare == "SFTLD_vs_CON",]
gene_count <- all_df_sig_select_ALSOD_Genes_SFTLD %>% 
  count(gene,sort = TRUE)
gene_count$gene <- factor(gene_count$gene,levels = gene_count$gene)
p4gene <- gene_count$gene
p4 <- ggplot(gene_count,aes(x = gene,y = n,group = 1)) +
  geom_line(color = "#377eb8",linewidth = 1) +
  geom_point(size = 2,color = "#e41a1c") +
  theme_bw() +
  theme(
    text = element_text(color = "black",family = "Times New Roman",face = "bold"),
    axis.text.x = element_text(angle = 45,hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "Gene",
    y = "Count",
    title = "SFTLD_vs_CON"
  )


combined_plot <- cowplot::plot_grid(
  p1,p2,p3,p4,
  nrow = 4,
  byrow = TRUE,
  labels = c("A","B","C","D"),
  label_fontfamily = "Times New Roman",
  label_size = 14
)


ggsave(
  filename = "./analysiswork/result/plot/cell_diff_combined_001.tiff",
  plot = combined_plot,
  dpi = 500,
  width = 15,
  height = 15,
  units = "in"
)

gene_list <- list(
  C9ALS = p1gene %>% as.character(),
  C9FTLD = p2gene %>% as.character(),
  SALS = p3gene %>% as.character(),
  SFTLD = p4gene %>% as.character()
)

Reduce(intersect,gene_list)
library(UpSetR)
tiff("./analysiswork/result/plot/cell_diff_combined_002.tiff",width = 8, height = 6, units = "in",res = 600)
upset(
  fromList(gene_list),
  keep.order = F,
  order.by = "freq",
  main.bar.color = "#a87ec0",
  sets.bar.color = c(C9ALS = "#e2735f",C9FTLD = "#f5b99d",SALS = "#79b4d7",SFTLD = "#c6e2ed")
)
dev.off()



