library(donutsk)
library(tidyverse)
library(ggpie)
library(ggplot2)
library(cowplot)
library(qs2)
library(Seurat)
library(ggh4x)

ex_cell_type_color <- setNames(c("#CC247C","#E95351","#F7A24F","#7e4909","#9bdfdf","#7c9895","gray80"),c("L2-3 IT","L5 IT","L6 CT","L6 IT","L6 IT Car3","L5-6 mixed","Excambiguous"))
in_cell_type_color <- setNames(c("#376e8e","#84badb","#c2b4d3","#839dd1","gray80"),c("PVALB","VIP","SST","LAMP5","Inhambiguous"))
noneu_cell_type_color <- setNames(c("#3c234a","#925eb0","#89b780","#eed0e0","#dbc9b3","#7ec0c3","gray80"),c("Oligo","OPC","Astro","Endo","Micro","VLMC","Neuambiguous"))
all_sub_type_color <- c(ex_cell_type_color,in_cell_type_color,noneu_cell_type_color)



#### 001 大类注释umap
seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_seurat_annofinal.qs2")
### combined plot
main_cell_type_color <- setNames(c("#c0321a","#547bb4","#eab883"),c("Excitatory neuron","Inhibitory neuron","Non neuron"))
p1 <- umap_group(
  seurat_obj = seuobj,
  groupcol = cell_type,
  groupcolor = main_cell_type_color,
  plot_title = "Main cell type",
  point_size = 0.08,
  point_alpha = 0.5,
  color_title = "",
  axis_title_size = 10,
  axis_text_size = 8,
  plot_title_size = 15,
  legend_text_size = 8,
  legend_size = 5,
  label_text_size = 4,
  legend_position = c(0.8,0.1)
) +
  theme(
    panel.border = element_rect(colour = "black",size = 1),
    axis.line = element_blank()
  )

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_001.tiff",
  plot = p1,
  dpi = 500,
  width = 6,
  height = 6,
  units = "in"
)



#### 002. 分组批次umap
disease_color <- setNames(c("#a45c72","#3c789d","#9369ab"),c("ALS","FTLD","CON"))
p1 <- umap_group(
  seurat_obj = seuobj,
  groupcol = disease_type,
  plot_title = "Disease batch removed",
  groupcolor = disease_color,
  color_title = "",
  remove_name = c("ALS","FTLD","CON"),
  point_size = 0.08,
  point_alpha = 0.5,
  axis_title_size = 10,
  axis_text_size = 8,
  plot_title_size = 13,
  legend_text_size = 8,
  legend_size = 5
)


sub_disease_color <- setNames(c("#e2735f","#f5b99d","#9369ab","#4f9dcb","#c6e2ed"),c("C9ALS","SALS","CON","C9FTLD","SFTLD"))
p2 <- umap_group(
  seurat_obj = seuobj,
  groupcol = sub_disease_type,
  plot_title = "Sub disease batch removed",
  groupcolor = sub_disease_color,
  color_title = "",
  remove_name = c("C9ALS","SALS","CON","C9FTLD","SFTLD"),
  point_size = 0.08,
  point_alpha = 0.5,
  axis_title_size = 10,
  axis_text_size = 8,
  plot_title_size = 13,
  legend_text_size = 8,
  legend_size = 5
)

### 疾病批次
sex_color <- setNames(c("#0787c3","#a51c36"),c("M","F"))
p3 <- umap_group(
  seurat_obj = seuobj,
  groupcol = Sex,
  plot_title = "Sex batch removed",
  groupcolor = sex_color,
  color_title = "",
  remove_name = c("M","F"),
  point_size = 0.08,
  point_alpha = 0.5,
  axis_title_size = 10,
  axis_text_size = 8,
  plot_title_size = 13,
  legend_text_size = 8,
  legend_size = 5
)

Batch_color <- setNames(c("#37939a","#304e7e"),c("GSE219280","SYN53421586"))
p4 <- umap_group(
  seurat_obj = seuobj,
  groupcol = Batch,
  plot_title = "Project batch removed",
  groupcolor = Batch_color,
  color_title = "",
  remove_name = c("GSE219280","SYN53421586"),
  point_size = 0.08,
  point_alpha = 0.5,
  axis_title_size = 10,
  axis_text_size = 8,
  plot_title_size = 13,
  legend_text_size = 8,
  legend_size = 5
)

p1 <- p1 + theme(legend.position = c(0.85,0.15))
p2 <- p2 + theme(legend.position = c(0.85,0.15))
p3 <- p3 + theme(legend.position = c(0.85,0.15))
p4 <- p4 + theme(legend.position = c(0.85,0.15))

pcom <- plot_grid(
  p1,p2,p3,p4,
  labels = c('A','B',"C","D"),
  nrow = 2, 
  byrow = TRUE,
  label_size = 14,
  label_fontfamily = "Times New Roman",
  rel_heights = c(1,1)
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_002.tiff",
  plot = pcom,
  dpi = 500,
  width = 10,
  height = 10,
  units = "in"
)



sub_disease_color <- setNames(c("#e2735f","#f5b99d","#9369ab","#4f9dcb","#c6e2ed"),c("C9ALS","SALS","CON","C9FTLD","SFTLD"))

meta_df <- seuobj@meta.data
umap_df <- Embeddings(
  object = seuobj,
  reduction = "umap"
)
merged_df <- merge(umap_df,meta_df,by = "row.names")
merged_df <- column_to_rownames(merged_df,var = "Row.names")
plotdata <- merged_df

p <- ggplot() +
  geom_point(
    data = plotdata,
    mapping = aes(x = umap_1,y = umap_2,color = cell_type),
    size = 0.001
  ) +
  scale_color_manual(
    values = main_cell_type_color
  ) +
  facet_grid(
    rows = vars(Region),
    cols = vars(sub_disease_type),
    scale = "fixed"
  ) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "",
    title = ""
  ) +
  theme(
    plot.title = element_text(size = 15,face = "bold",color = "black",hjust = 0.5),
    text = element_text(family = "Times New Roman"),
    panel.background = element_rect(fill = "white",color = "black",size = 1),
    panel.border = element_rect(color = "black",size = 1),
    panel.grid.minor = element_line(color = "gray80"),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 10,face = "bold",color = "black"),
    strip.text.x = element_text(size = 10,face = "bold",color = "black"),
    strip.text.y = element_text(size = 10,face = "bold",color = "black"),
    legend.key = element_rect(color = NA,fill = NA),
    legend.key.size = unit(0.3,"cm"),
    legend.key.width = unit(0.1,"cm"),
    legend.key.height = unit(0.2,"cm"),
    legend.key.spacing.y = unit(0.05,"cm"),
    legend.background = element_rect(fill = "transparent")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5,alpha = 0.8)
      )
  ) +
  coord_fixed() +
  theme(legend.position = "bottom")

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_003.tiff",
  plot = p,
  dpi = 500,
  width = 10,
  height = 5,
  units = "in"
)




select_main_marker_list <- list(
  `Excitatory neuron` = c("SLC17A7","SATB2","NRGN","CNKSR2","CAMK2A","RORB"),
  `Inhibitory neuron` = c("GAD1","GAD2","LHX6","DLX6-AS1"),
  `Non neuron` = c("SLC1A3","SOX9")
)
main_markers_vec <- select_main_marker_list %>% unlist() %>% unname()
seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_seurat_annofinal.qs2")

p <- multi_marker_umap(
  marker_vec = main_markers_vec,
  seuratobj = seuobj,
  color_vec = c("#fae5b8","#992224")
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_004.tiff",
  plot = p,
  dpi = 500,
  width = 10,
  height = 10,
  units = "in"
)


## 005.大类marker bubble图
p <- bubble_genes(
  seurat_obj = seuobj,
  groupcol = cell_type,
  grouporder = c("Excitatory neuron","Inhibitory neuron","Non neuron"),
  groupcolor = main_cell_type_color,
  plotmarker = main_markers_vec,
  x_text_size = 10,
  y_text_size = 10,
  fill_palette = c("#fffec8","#fc913c","#800025")
)

p1 <- p + theme(legend.position = "none")

legend <- get_legend(
  p + theme(
    legend.box = "horizontal",
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
  p1,
  legend,
  nrow = 2,
  rel_heights = c(9,1)
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_005.tiff",
  plot = combined_plot,
  dpi = 500,
  width = 8,
  height = 3,
  units = "in"
)


ex_seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_ex_seurat_annotated.qs2")
p1 <- umap_group(
  seurat_obj = ex_seuobj,
  groupcol = sub_cell_type,
  groupcolor = ex_cell_type_color,
  plot_title = "Exc cell subtype",
  color_title = "",
  point_size = 0.2,
  point_alpha = 0.5,
  axis_title_size = 10,
  axis_text_size = 8,
  plot_title_size = 15,
  legend_text_size = 8,
  legend_size = 5,
  remove_name = c("Neuambiguous","Excambiguous","Inhambiguous")
)

### inh
in_seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_in_seurat_annotated.qs2")
p2 <- umap_group(
  seurat_obj = in_seuobj,
  groupcol = sub_cell_type,
  groupcolor = in_cell_type_color,
  plot_title = "Inh cell subtype",
  color_title = "",
  point_size = 0.2,
  point_alpha = 0.5,
  legend_position = "none",
  axis_title_size = 10,
  axis_text_size = 8,
  plot_title_size = 15,
  legend_text_size = 8,
  legend_size = 5,
  remove_name = c("Neuambiguous","Excambiguous","Inhambiguous")
)

## neuro
noneu_seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_noneu_seurat_annotated.qs2")

p3 <- umap_group(
  seurat_obj = noneu_seuobj,
  groupcol = sub_cell_type,
  groupcolor = noneu_cell_type_color,
  plot_title = "Noneu cell subtype",
  color_title = "",
  point_size = 0.2,
  point_alpha = 0.5,
  legend_position = "none",
  axis_title_size = 10,
  axis_text_size = 8,
  plot_title_size = 15,
  legend_text_size = 8,
  legend_size = 5,
  remove_name = c("Neuambiguous","Excambiguous","Inhambiguous")
)

## 所有亚类
seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_seurat_annofinal.qs2")
p4 <- umap_group(
  seurat_obj = seuobj,
  groupcol = sub_cell_type,
  groupcolor = all_sub_type_color,
  plot_title = "All cell subtype",
  color_title = "",
  point_size = 0.1,
  point_alpha = 0.5,
  axis_title_size = 10,
  axis_text_size = 8,
  plot_title_size = 15,
  legend_text_size = 8,
  legend_size = 8,
  remove_name = c("Neuambiguous","Excambiguous","Inhambiguous")
)

common_legend <- get_legend(
  p4 + theme(
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.box.just = "center",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,0)
  ) +
    guides(
      color = guide_legend(
        nrow = 2,
        override.aes = list(size = 5,alpha = 0.8))
    )
)


p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")

pcom <- plot_grid(
  p1,p2,p3,p4,
  labels = c('A','B',"C","D"),
  nrow = 2, 
  byrow = TRUE,
  label_size = 14,
  label_fontfamily = "Times New Roman",
  rel_heights = c(1,1),
  rel_widths = c(1,1)
)

pcom_addlegend <- plot_grid(
  pcom,
  common_legend,
  nrow = 2, 
  byrow = TRUE,
  label_fontfamily = "Times New Roman",
  rel_heights = c(9,1)
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_006.tiff",
  plot = pcom_addlegend,
  dpi = 500,
  width = 10,
  height = 10,
  units = "in"
)



### 007
p <- ggplot() +
  geom_point(
    data = plotdata,
    mapping = aes(x = umap_1,y = umap_2,color = sub_cell_type),
    size = 0.005
  ) +
  scale_color_manual(values = all_sub_type_color) +
  facet_grid(
    rows = vars(Region),
    cols = vars(sub_disease_type),
    scale = "fixed"
  ) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "",
    title = ""
  ) +
  theme(
    plot.title = element_text(size = 15,face = "bold",color = "black",hjust = 0.5),
    text = element_text(family = "Times New Roman"),
    panel.background = element_rect(fill = "white",color = "black",size = 1),
    panel.border = element_rect(color = "black",size = 1),
    panel.grid.minor = element_line(color = "gray80"),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 10,face = "bold",color = "black"),
    strip.text.x = element_text(size = 10,face = "bold",color = "black"),
    strip.text.y = element_text(size = 10,face = "bold",color = "black"),
    legend.key = element_rect(color = NA,fill = NA),
    legend.key.size = unit(0.3,"cm"),
    legend.key.width = unit(0.1,"cm"),
    legend.key.height = unit(0.2,"cm"),
    legend.key.spacing.y = unit(0.05,"cm"),
    legend.background = element_rect(fill = "transparent")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5,alpha = 0.8)
    )
  ) +
  coord_fixed() +
  theme(legend.position = "none")



common_legend <- get_legend(
  p + theme(
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.box.just = "center",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,0)
  ) +
    guides(
      color = guide_legend(
        nrow = 2,
        override.aes = list(size = 5,alpha = 0.8))
    )
)

pcom_addlegend <- plot_grid(
  p,
  common_legend,
  nrow = 2, 
  byrow = TRUE,
  label_fontfamily = "Times New Roman",
  rel_heights = c(9,1)
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_007.tiff",
  plot = pcom_addlegend,
  dpi = 500,
  width = 10,
  height = 5,
  units = "in"
)



### 008
select_main_marker_list <- list(
  `Excitatory neuron` = c("SLC17A7","SATB2","NRGN","CNKSR2","CAMK2A","RORB"),
  `Inhibitory neuron` = c("GAD1","GAD2","LHX6","DLX6-AS1"),
  `Non neuron` = c("SLC1A3","SOX9")
)
main_markers_vec <- select_main_marker_list %>% unlist() %>% unname()
markers_vec <- main_markers_vec
cell_type_color <- all_sub_type_color


p <- bubble_genes(
  seurat_obj = seuobj,
  groupcol = sub_cell_type,
  grouporder = names(cell_type_color),
  groupcolor = cell_type_color,
  plotmarker = markers_vec,
  x_text_size = 10,
  y_text_size = 10,
  fill_palette = c("#fffec8","#fc913c","#800025")
)

p1 <- p + theme(legend.position = "none")

legend <- get_legend(
  p + theme(
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
  p1,
  legend,
  nrow = 2,
  rel_heights = c(9,1)
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_008.tiff",
  plot = combined_plot,
  dpi = 500,
  width = 8,
  height = 8,
  units = "in"
)


### 009
ex_marker_list <- list(
  `L2-3 IT` = c("CUX2"),
  `L5 IT `= c("RORB"),
  `L6 CT` = c("HS3ST4","SEMA5A"),
  `L6 IT` = c("THEMIS"),
  `L6 IT Car3` = c("RGS12","ATP10A")
)

in_marker_list <- list(
  PVALB = c("PVALB"),
  SST = c("SST"),
  VIP = c("VIP"),
  LAMP5 = c("LAMP5")
)


noneu_marker_list <- list(
  Oligo = c("ENPP2","ST18","CDK18","MBP","ENPP6","LRRC63","LINC01378"),
  Micro = c("DOCK8","APBB1IP","CD74","P2RY12"),
  OPC = c("VCAN","STK32A","COL20A1"),
  Astro = c("COL5A3"),
  Endo = c("SRGN","CLDN5"),
  VLMC = c("COLEC12")
)

all_list <- c(ex_marker_list,in_marker_list,noneu_marker_list)


markers_vec <- all_list %>% unlist() %>% unname()
cell_type_color <- all_sub_type_color

select_order <- c("L2-3 IT","L5 IT","L6 CT","L6 IT","L6 IT Car3","L5-6 mixed","Excambiguous","PVALB","SST","VIP","LAMP5","Inhambiguous","Oligo","OPC","Astro","Endo","Micro","VLMC","Neuambiguous")
p <- bubble_genes(
  seurat_obj = seuobj,
  groupcol = sub_cell_type,
  grouporder = select_order,
  groupcolor = cell_type_color,
  plotmarker = markers_vec,
  x_text_size = 10,
  y_text_size = 10,
  fill_palette = c("#fffec8","#fc913c","#800025")
)


ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_009.tiff",
  plot = p,
  dpi = 500,
  width = 8,
  height = 5,
  units = "in"
)



### 010 
ex_seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_ex_seurat_annotated.qs2")
in_seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_in_seurat_annotated.qs2")
noneu_seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_noneu_seurat_annotated.qs2")
seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_seurat_annofinal.qs2")

ex_seuobj@graphs <- list()
in_seuobj@graphs <- list()
noneu_seuobj@graphs <- list()
seuobj@graphs <- list()

ex_seuobj@reductions <- list()
in_seuobj@reductions <- list()
noneu_seuobj@reductions <- list()
seuobj@reductions <- list()

ex_seurat_list <- SplitObject(ex_seuobj,split.by = "Region")
in_seurat_list <- SplitObject(in_seuobj,split.by = "Region")
noneu_seurat_list <- SplitObject(noneu_seuobj,split.by = "Region")
seurat_list <- SplitObject(seuobj,split.by = "Region")

p1 <- percent_bar(
  seurat_obj = ex_seurat_list$MC,
  group1 = new_sample,
  group2 = sub_cell_type,
  groupcolor = ex_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Exc cell subtype"
)


p2 <- percent_bar(
  seurat_obj = in_seurat_list$MC,
  group1 = new_sample,
  group2 = sub_cell_type,
  groupcolor = in_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Inh cell subtype"
)



p3 <- percent_bar(
  seurat_obj = noneu_seurat_list$MC,
  group1 = new_sample,
  group2 = sub_cell_type,
  groupcolor = noneu_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Noneu cell subtype"
)



p4 <- percent_bar(
  seurat_obj = seurat_list$MC,
  group1 = new_sample,
  group2 = sub_cell_type,
  groupcolor = all_sub_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "All cell subtype"
)


p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")

common_legend <- get_legend(
  p4 + theme(
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.box.just = "center",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,0)
  ) +
    guides(
      fill = guide_legend(
        nrow = 2,
        override.aes = list(size = 5,alpha = 0.8))
    )
)

pcom <- plot_grid(
  p1,p2,p3,p4,
  labels = c('A','B',"C","D"),
  nrow = 4, 
  byrow = TRUE,
  label_size = 14,
  label_fontfamily = "Times New Roman",
  rel_heights = c(1,1,1,1),
  rel_widths = c(1,1,1,1)
)

pcom_addlegend <- plot_grid(
  pcom,
  common_legend,
  nrow = 2, 
  byrow = TRUE,
  label_fontfamily = "Times New Roman",
  rel_heights = c(9,1)
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_010.tiff",
  plot = pcom_addlegend,
  dpi = 500,
  width = 10,
  height = 15,
  units = "in"
)




## 011.亚类分组流式条形图
p1 <- percent_bar(
  seurat_obj = ex_seurat_list$MC,
  group1 = sub_disease_type,
  group2 = sub_cell_type,
  groupcolor = ex_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Exc cell subtype"
)

p2 <- percent_bar(
  seurat_obj = in_seurat_list$MC,
  group1 = sub_disease_type,
  group2 = sub_cell_type,
  groupcolor = in_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Inh cell subtype"
)

p3 <- percent_bar(
  seurat_obj = noneu_seurat_list$MC,
  group1 = sub_disease_type,
  group2 = sub_cell_type,
  groupcolor = noneu_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Noneu cell subtype"
)

p4 <- percent_bar(
  seurat_obj = seurat_list$MC,
  group1 = sub_disease_type,
  group2 = sub_cell_type,
  groupcolor = all_sub_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "All cell subtype"
)


p1 <- p1 + theme(legend.position = "none",axis.text.x = element_text(color = "black",size = 10),axis.text.y = element_text(color = "black",size = 15))
p2 <- p2 + theme(legend.position = "none",axis.text.x = element_text(color = "black",size = 10),axis.text.y = element_text(color = "black",size = 15))
p3 <- p3 + theme(legend.position = "none",axis.text.x = element_text(color = "black",size = 10),axis.text.y = element_text(color = "black",size = 15))
p4 <- p4 + theme(legend.position = "none",axis.text.x = element_text(color = "black",size = 10),axis.text.y = element_text(color = "black",size = 15))

common_legend <- get_legend(
  p4 + theme(
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.box.just = "center",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,0)
  ) +
    guides(
      fill = guide_legend(
        nrow = 2,
        override.aes = list(size = 5,alpha = 0.8))
    )
)

pcom <- plot_grid(
  p1,p2,p3,p4,
  labels = c('A','B',"C","D"),
  nrow = 2, 
  byrow = TRUE,
  label_size = 14,
  label_fontfamily = "Times New Roman",
  rel_heights = c(1,1,1,1),
  rel_widths = c(1,1,1,1)
)

pcom_addlegend <- plot_grid(
  pcom,
  common_legend,
  nrow = 2, 
  byrow = TRUE,
  label_fontfamily = "Times New Roman",
  rel_heights = c(9,1)
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_011.tiff",
  plot = pcom_addlegend,
  dpi = 500,
  width = 10,
  height = 10,
  units = "in"
)


### 012
p1 <- percent_bar(
  seurat_obj = ex_seurat_list$FC,
  group1 = new_sample,
  group2 = sub_cell_type,
  groupcolor = ex_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Exc cell subtype"
)


p2 <- percent_bar(
  seurat_obj = in_seurat_list$FC,
  group1 = new_sample,
  group2 = sub_cell_type,
  groupcolor = in_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Inh cell subtype"
)



p3 <- percent_bar(
  seurat_obj = noneu_seurat_list$FC,
  group1 = new_sample,
  group2 = sub_cell_type,
  groupcolor = noneu_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Noneu cell subtype"
)



p4 <- percent_bar(
  seurat_obj = seurat_list$FC,
  group1 = new_sample,
  group2 = sub_cell_type,
  groupcolor = all_sub_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "All cell subtype"
)


p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")

common_legend <- get_legend(
  p4 + theme(
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.box.just = "center",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,0)
  ) +
    guides(
      fill = guide_legend(
        nrow = 2,
        override.aes = list(size = 5,alpha = 0.8))
    )
)

pcom <- plot_grid(
  p1,p2,p3,p4,
  labels = c('A','B',"C","D"),
  nrow = 4, 
  byrow = TRUE,
  label_size = 14,
  label_fontfamily = "Times New Roman",
  rel_heights = c(1,1,1,1),
  rel_widths = c(1,1,1,1)
)

pcom_addlegend <- plot_grid(
  pcom,
  common_legend,
  nrow = 2, 
  byrow = TRUE,
  label_fontfamily = "Times New Roman",
  rel_heights = c(9,1)
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_012.tiff",
  plot = pcom_addlegend,
  dpi = 500,
  width = 10,
  height = 15,
  units = "in"
)




## 013
p1 <- percent_bar(
  seurat_obj = ex_seurat_list$FC,
  group1 = sub_disease_type,
  group2 = sub_cell_type,
  groupcolor = ex_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Exc cell subtype"
)

p2 <- percent_bar(
  seurat_obj = in_seurat_list$FC,
  group1 = sub_disease_type,
  group2 = sub_cell_type,
  groupcolor = in_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Inh cell subtype"
)

p3 <- percent_bar(
  seurat_obj = noneu_seurat_list$FC,
  group1 = sub_disease_type,
  group2 = sub_cell_type,
  groupcolor = noneu_cell_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "Noneu cell subtype"
)

p4 <- percent_bar(
  seurat_obj = seurat_list$FC,
  group1 = sub_disease_type,
  group2 = sub_cell_type,
  groupcolor = all_sub_type_color,
  legend_title = "Cell Type",
  ifflow = TRUE,
  axis_text_size = 6,
  plot_title = "All cell subtype"
)


p1 <- p1 + theme(legend.position = "none",axis.text.x = element_text(color = "black",size = 10),axis.text.y = element_text(color = "black",size = 15))
p2 <- p2 + theme(legend.position = "none",axis.text.x = element_text(color = "black",size = 10),axis.text.y = element_text(color = "black",size = 15))
p3 <- p3 + theme(legend.position = "none",axis.text.x = element_text(color = "black",size = 10),axis.text.y = element_text(color = "black",size = 15))
p4 <- p4 + theme(legend.position = "none",axis.text.x = element_text(color = "black",size = 10),axis.text.y = element_text(color = "black",size = 15))

common_legend <- get_legend(
  p4 + theme(
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.box.just = "center",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,0)
  ) +
    guides(
      fill = guide_legend(
        nrow = 2,
        override.aes = list(size = 5,alpha = 0.8))
    )
)

pcom <- plot_grid(
  p1,p2,p3,p4,
  labels = c('A','B',"C","D"),
  nrow = 2, 
  byrow = TRUE,
  label_size = 14,
  label_fontfamily = "Times New Roman",
  rel_heights = c(1,1,1,1),
  rel_widths = c(1,1,1,1)
)

pcom_addlegend <- plot_grid(
  pcom,
  common_legend,
  nrow = 2, 
  byrow = TRUE,
  label_fontfamily = "Times New Roman",
  rel_heights = c(9,1)
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_013.tiff",
  plot = pcom_addlegend,
  dpi = 500,
  width = 10,
  height = 10,
  units = "in"
)



### 细胞特异性pseudobulk分析
seuobj <- qs_read("./analysiswork/cleandata/snRNA_VAL_seurat_annofinal.qs2")
seuobj$region_disease_cell <- paste(seuobj$Region,seuobj$sub_disease_type,seuobj$sub_cell_type, sep = "_")

# seuobj$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
# Idents(ifnb) <- "celltype.stim"
# mono.de <- FindMarkers(ifnb, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", verbose = FALSE)
# 
# pseudo_seuobj <- AggregateExpression(
#   object = seuobj, 
#   assays = "RNA", 
#   return.seurat = T, 
#   group.by = c("Region","sub_disease_type","new_sample","sub_cell_type")
# )
# 
# pseudo_seuobj$cell_disease <- paste(pseudo_seuobj$Region,pseudo_seuobj$sub_cell_type, pseudo_seuobj$sub_disease_type, sep = "_")

Idents(seuobj) <- "region_disease_cell"
all_celltype <- unique(seuobj$sub_cell_type)
remove_celltype <- c("Excambiguous","Inhambiguous","Neuambiguous")
keep_celltype <- all_celltype[!all_celltype %in% remove_celltype]

compare_group <- c("C9ALS","C9FTLD","SALS","SFTLD")     
con_group <- "CON"

result_list <- list()
select_region <- "MC"
for (compare in compare_group) {
  con <- con_group
  for (celltype in keep_celltype) {
    ident1 <- paste0(select_region,"_",compare,"_",celltype)
    ident2 <- paste0(select_region,"_",con,"_",celltype)
    
    diffresult <- FindMarkers(
      object = seuobj, 
      ident.1 = ident1, 
      ident.2 = ident2,
      test.use = "wilcox"
    )
    
    diffresult <- rownames_to_column(diffresult,var = "gene")
    diffresult$compare <- paste0(compare,"_vs_",con)
    diffresult$celltype <- celltype
    
    result_list <- append(result_list,list(diffresult))
  }
}

all_df <- do.call(rbind,result_list)
all_df <- all_df %>% 
  mutate(
    updown = case_when(
      avg_log2FC > 0 ~ "UP",
      avg_log2FC < 0 ~ "DOWN"
    )
  )


qs_save(all_df,"./analysiswork/cleandata/snRNA_VAL_MC_allcelltype_diff_genes.qs2")




result_list <- list()
select_region <- "FC"
for (compare in compare_group) {
  con <- con_group
  for (celltype in keep_celltype) {
    ident1 <- paste0(select_region,"_",compare,"_",celltype)
    ident2 <- paste0(select_region,"_",con,"_",celltype)
    
    diffresult <- FindMarkers(
      object = seuobj, 
      ident.1 = ident1, 
      ident.2 = ident2,
      test.use = "wilcox"
    )
    
    diffresult <- rownames_to_column(diffresult,var = "gene")
    diffresult$compare <- paste0(compare,"_vs_",con)
    diffresult$celltype <- celltype
    
    result_list <- append(result_list,list(diffresult))
  }
}
all_df <- do.call(rbind,result_list)
all_df <- all_df %>% 
  mutate(
    updown = case_when(
      avg_log2FC > 0 ~ "UP",
      avg_log2FC < 0 ~ "DOWN"
    )
  )


qs_save(all_df,"./analysiswork/cleandata/snRNA_VAL_FC_allcelltype_diff_genes.qs2")


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



ex_cell_type_color <- setNames(c("#CC247C","#E95351","#F7A24F","#7e4909","#9bdfdf","#7c9895","gray80"),c("L2-3 IT","L5 IT","L6 CT","L6 IT","L6 IT Car3","L5-6 mixed","Excambiguous"))
in_cell_type_color <- setNames(c("#376e8e","#84badb","#c2b4d3","#839dd1","gray80"),c("PVALB","VIP","SST","LAMP5","Inhambiguous"))
noneu_cell_type_color <- setNames(c("#3c234a","#925eb0","#89b780","#eed0e0","#dbc9b3","#7ec0c3","gray80"),c("Oligo","OPC","Astro","Endo","Micro","VLMC","Neuambiguous"))
all_sub_type_color <- c(ex_cell_type_color,in_cell_type_color,noneu_cell_type_color)

up_compare_color <- setNames(c("#e2735f","#e2735f","#e2735f","#e2735f"),c("C9ALS_vs_CON","SALS_vs_CON","C9FTLD_vs_CON","SFTLD_vs_CON"))
down_compare_color <- setNames(c("#4f9dcb","#4f9dcb","#4f9dcb","#4f9dcb"),c("C9ALS_vs_CON","SALS_vs_CON","C9FTLD_vs_CON","SFTLD_vs_CON"))



multigroup_bar <- function(
    plotdata = df,
    plot_title = "",
    compare_color = "",
    x_title = "",
    x_title_size = 15
) {
  p <- ggplot(
    data = plotdata,
    mapping = aes(x = num, y = interaction(sub_cell_type,main_cell_type))
  ) +
    geom_vline(xintercept = 0,color = "gray50",linewidth = 0.2,linetype = 2) +
    geom_col(
      mapping = aes(group = compare,fill = compare),
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    scale_x_continuous(
      labels = function(x) abs(x)
    ) +
    scale_fill_manual(
      values = compare_color
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
          element_text(color = all_sub_type_color,face = "bold"),
          element_text(color = main_cell_type_color,face = "bold")
        )
      )
    ) +
    theme(
      legend.position = "none"
    )
  
  return(p)
}

p1 <- multigroup_bar(final_df[final_df$updown == "UP",],compare_color = up_compare_color,x_title = "UP regulate gene number")
p2 <- multigroup_bar(final_df[final_df$updown == "DOWN",],compare_color = down_compare_color,x_title = "Down regulate gene number")

pcom <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  labels = c("A","B"),
  nrow = 2
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_014.tiff",
  plot = pcom,
  width = 10,
  height = 10,
  units = "in"
)




all_df <- qs_read("./analysiswork/cleandata/snRNA_VAL_FC_allcelltype_diff_genes.qs2")

all_df_sig <- all_df %>% 
  filter(p_val_adj < 0.001 & abs(avg_log2FC) > 1)

result_longer <- all_df_sig %>% 
  dplyr::count(compare,celltype,updown)

colnames(result_longer) <- c("compare","sub_cell_type","updown","num")

cellinfo <- seuobj@meta.data[,c("cell_type","sub_cell_type")] %>% unique()
colnames(cellinfo) <- c("main_cell_type","sub_cell_type")

final_df <- result_longer %>% 
  left_join(cellinfo,by = "sub_cell_type")



ex_cell_type_color <- setNames(c("#CC247C","#E95351","#F7A24F","#7e4909","#9bdfdf","#7c9895","gray80"),c("L2-3 IT","L5 IT","L6 CT","L6 IT","L6 IT Car3","L5-6 mixed","Excambiguous"))
in_cell_type_color <- setNames(c("#376e8e","#84badb","#c2b4d3","#839dd1","gray80"),c("PVALB","VIP","SST","LAMP5","Inhambiguous"))
noneu_cell_type_color <- setNames(c("#3c234a","#925eb0","#89b780","#eed0e0","#dbc9b3","#7ec0c3","gray80"),c("Oligo","OPC","Astro","Endo","Micro","VLMC","Neuambiguous"))
all_sub_type_color <- c(ex_cell_type_color,in_cell_type_color,noneu_cell_type_color)

up_compare_color <- setNames(c("#f09ba0","#f09ba0","#f09ba0","#f09ba0"),c("C9ALS_vs_CON","SALS_vs_CON","C9FTLD_vs_CON","SFTLD_vs_CON"))
down_compare_color <- setNames(c("#9bbbe1","#9bbbe1","#9bbbe1","#9bbbe1"),c("C9ALS_vs_CON","SALS_vs_CON","C9FTLD_vs_CON","SFTLD_vs_CON"))



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
      mapping = aes(group = compare,fill = compare),
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    scale_x_continuous(
      labels = function(x) abs(x)
    ) +
    scale_fill_manual(
      values = compare_color
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

p1 <- multigroup_bar(final_df[final_df$updown == "UP",],compare_color = up_compare_color,x_title = "UP regulate gene number")
p2 <- multigroup_bar(final_df[final_df$updown == "DOWN",],compare_color = down_compare_color,x_title = "Down regulate gene number")

pcom <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  labels = c("A","B"),
  nrow = 2
)

ggsave(
  filename = "./analysiswork/result/plot/snRNA_VAL_combined_015.tiff",
  plot = pcom,
  width = 10,
  height = 10,
  units = "in"
)
