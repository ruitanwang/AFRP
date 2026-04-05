### marker cruted mannual list
human_cellmarker_genes <- read.xlsx("./analysiswork/helpdata/Cell_marker_Human.xlsx")
human_brain_normal_cellmarker_genes <- human_cellmarker_genes %>% 
  filter(tissue_class == "Brain" & cancer_type == "Normal")

human_brain_normal_cellmarker_genes <- human_brain_normal_cellmarker_genes %>% 
  filter(species == "Human") %>% 
  group_by(cell_name) %>% 
  summarise(gene_set = list(unique(marker))) %>% 
  deframe()


human_main_ref_marker <- list(
  `Excitatory neuron` = human_brain_normal_cellmarker_genes$`Excitatory neuron`,
  `Inhibitory neuron` = human_brain_normal_cellmarker_genes$`Inhibitory neuron`,
  `Non neuron` = human_brain_normal_cellmarker_genes$`Non-neuron`
)


main_marker_list <- list(
  `Excitatory neuron` = c("SLC17A7","SATB2","NRGN","CNKSR2","CAMK2A","RORB","STMN2"),
  `Inhibitory neuron` = c("GAD1","GAD2","LHX6","DLX6-AS1","CALB2","NR2F2"),
  `Non neuron` = c("SLC1A3","SOX9")
)

main_markers_vec <- main_marker_list %>% unlist() %>% unname()

## 兴奋性神经元marker
ref_sub_anno <- read.csv("./analysiswork/helpdata/ref_sub_anno.csv")
ref_sub_anno <- ref_sub_anno[,c(1,8:12)]
df_to_marker_list <- function(data){
  cell_col <- data[,1,drop = TRUE]
  marker_cols <- data[,-1,drop = FALSE]
  marker_list <- apply(marker_cols,1,function(row){
    unique(na.omit(as.character(row)))
  })
  
  names(marker_list) <- cell_col
  
  marker_list <- marker_list[sapply(marker_list,length) > 0]
  
  return(marker_list)
}


ref_sub_anno_list <- df_to_marker_list(ref_sub_anno)
ex_main_ref_marker_list <- ref_sub_anno_list[6:13]


### 抑制性神经元marker
ref_sub_anno <- read.csv("./analysiswork/helpdata/ref_sub_anno.csv")
ref_sub_anno <- ref_sub_anno[,c(1,8:12)]
df_to_marker_list <- function(data){
  cell_col <- data[,1,drop = TRUE]
  marker_cols <- data[,-1,drop = FALSE]
  marker_list <- apply(marker_cols,1,function(row){
    unique(na.omit(as.character(row)))
  })
  
  names(marker_list) <- cell_col
  
  marker_list <- marker_list[sapply(marker_list,length) > 0]
  
  return(marker_list)
}


ref_sub_anno_list <- df_to_marker_list(ref_sub_anno)


target <- c("LAMP5","PVALB","SST","VIP")
inh_marker_manual_list <- lapply(target, function(x) {
  new_samll_list <- ref_sub_anno_list[grep(paste(x,collapse = "|"),names(ref_sub_anno_list))]
  
  part1 <- unique(unlist(new_samll_list))
  part2 <- unique(sapply(strsplit(names(new_samll_list)," "),tail,1))
  part_final <- unique(c(part1,part2))
  part_final <- part_final[part_final != ""]
  
  return(part_final)
})

names(inh_marker_manual_list) <- target




### 非神经元marker
ref_sub_anno <- read.csv("./analysiswork/helpdata/ref_sub_anno.csv")
ref_sub_anno <- ref_sub_anno[,c(1,8:12)]
df_to_marker_list <- function(data){
  cell_col <- data[,1,drop = TRUE]
  marker_cols <- data[,-1,drop = FALSE]
  marker_list <- apply(marker_cols,1,function(row){
    unique(na.omit(as.character(row)))
  })
  
  names(marker_list) <- cell_col
  
  marker_list <- marker_list[sapply(marker_list,length) > 0]
  
  return(marker_list)
}


ref_sub_anno_list <- df_to_marker_list(ref_sub_anno)
target <- c("Oligo","Astro","Micro","Micro-PVM","OPC","Endo","VLMC")
noneu_marker_manual_list <- lapply(target, function(x) {
  new_samll_list <- ref_sub_anno_list[grep(paste(x,collapse = "|"),names(ref_sub_anno_list))]
  
  part1 <- unique(unlist(new_samll_list))
  part2 <- unique(sapply(strsplit(names(new_samll_list)," "),tail,1))
  part_final <- unique(c(part1,part2))
  part_final <- part_final[part_final != ""]
  
  return(part_final)
})
names(noneu_marker_manual_list) <- target



ex_marker <- qs_read("./analysiswork/cleandata/snRNA_MCA_ex_seurat_marker.qs2")
in_marker <- qs_read("./analysiswork/cleandata/snRNA_MCA_in_seurat_marker.qs2")
noneu_marker <- qs_read("./analysiswork/cleandata/snRNA_MCA_noneu_seurat_marker.qs2")




ex_marker_list <- ex_marker %>% 
  group_by(celltype) %>% 
  summarise(gene_set = list(unique(intersect_genes))) %>% 
  deframe()

ex_marker_list_unique <- lapply(ex_marker_list, function(x){
  all_genes <- unlist(strsplit(x,","))
  unique(all_genes)
})
genes_counts <- table(unlist(ex_marker_list_unique))
unique_genes <- names(genes_counts)[genes_counts == 1]

ex_marker_list_unique_final <- lapply(ex_marker_list_unique, function(genes){
  intersect(genes,unique_genes)
})


in_marker_list <- in_marker %>% 
  group_by(celltype) %>% 
  summarise(gene_set = list(unique(intersect_genes))) %>% 
  deframe()

in_marker_list_unique <- lapply(in_marker_list, function(x){
  all_genes <- unlist(strsplit(x,","))
  unique(all_genes)
})
genes_counts <- table(unlist(in_marker_list_unique))
unique_genes <- names(genes_counts)[genes_counts == 1]

in_marker_list_unique_final <- lapply(in_marker_list_unique, function(genes){
  intersect(genes,unique_genes)
})



### 手动挑选
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

