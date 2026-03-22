#' 绘制单细胞UMAP分组可视化图
#'
#' 该函数基于Seurat对象生成高质量的UMAP降维图，根据指定分组进行颜色编码，支持子集筛选和自动聚类标注功能。
#'
#' @param seurat_obj Seurat对象（必需），包含UMAP嵌入和元数据
#' @param groupcol 分组列（必需），指定用于着色的元数据列（如细胞类型）
#' @param groupcolor 颜色向量（必需），长度需与分组数匹配，格式：`c("类型1"="#FF0000", "类型2"="#00FF00")`
#' @param plot_title 图像主标题（默认：空字符串）
#' @param point_size 点大小（默认：0.008），建议范围0.001-0.1
#' @param point_alpha 点透明度（默认：0.3），范围0（全透明）~1（不透明）
#' @param color_title 图例标题（默认：空字符串），通常设置为分组变量名
#' @param subcol 子集筛选列（可选），指定用于数据筛选的元数据列
#' @param subneed 子集筛选值（可选），当提供subcol时，指定需保留的子集值
#' @param num_size 细胞计数文本大小（默认：5）
#' @param legend_position 图例位置（默认："right"），可选："none"/"left"/"right"/"bottom"/"top"
#' @param label_text_size 聚类标签文本大小（默认：5）
#' @param legend_text_size 图例文本大小（默认：6）
#' @param legend_size 图例色标大小（默认：5），单位mm
#' @param legend_title_size 图例标题大小（默认：6）
#' @param axis_title_size 坐标轴标题大小（默认：7）
#' @param axis_text_size 坐标轴刻度文本大小（默认：6）
#' @param plot_title_size 主标题文本大小（默认：10）
#' @param family_theme 字体类型（默认："Times New Roman"），推荐"Arial"
#' @param autolabel 是否自动标注聚类（默认：TRUE），在分组中心添加标签
#' @param remove_name 需排除标注的分组名（默认："Excambiguous"）
#' @param text_bg_color 标签背景色（默认："black"），增强标签可读性
#'
#' @return 返回ggplot2对象
#' @export
#' @import Seurat ggplot2 purrr dplyr tibble DescTools rlang shadowtext forcats ggalluvial lubridate purrr readr stringr tidyr
#'
umap_group <- function(
    seurat_obj,
    groupcol,
    groupcolor,
    plot_title = "",
    point_size = 0.008,
    point_alpha = 0.3,
    color_title = "",
    subcol,
    subneed,
    num_size = 5,
    legend_position = "right",
    label_text_size = 5,
    legend_text_size = 6,
    legend_size = 5,
    legend_title_size = 6,
    axis_title_size = 7,
    axis_text_size = 6,
    plot_title_size = 10,
    family_theme = "Times New Roman",
    autolabel = TRUE,
    remove_name = "Excambiguous",
    text_bg_color = "black"
) {
  groupname <- rlang::enquo(groupcol)
  subcolname <- rlang::enquo(subcol)


  umap_df <- Embeddings(
    object = seurat_obj,
    reduction = "umap"
  )

  seurat_meta <- seurat_obj@meta.data

  if (missing(subcol)) {
    groupinfo <- seurat_meta[,rlang::as_label(groupname),drop = FALSE]
    colnames(groupinfo) <- "groupinfo"

    umapob <- merge(umap_df,groupinfo,by = "row.names")
    umapob <- column_to_rownames(umapob,var = colnames(umapob)[1])
    colnames(umapob) <- c("UMAP_1","UMAP_2","groupinfo")
  } else {
    groupinfo <- seurat_meta[,c(rlang::as_label(groupname),rlang::as_label(subcolname)),drop = FALSE]

    colnames(groupinfo) <- c("groupinfo","subgroup")
    umapob <- merge(umap_df,groupinfo,by = "row.names")
    umapob <- column_to_rownames(umapob,var = colnames(umapob)[1])
    colnames(umapob) <- c("UMAP_1","UMAP_2","groupinfo","subgroup")

    umapob <- umapob[umapob$subgroup == subneed,]
  }


  # groupinfo <- seurat_obj@meta.data %>%
  #   dplyr::pull(!!groupname) %>%
  #   as.data.frame()

  plotdata <- umapob
  cell_types <- unique(plotdata$groupinfo)
  cell_color <- groupcolor


  p <- ggplot(plotdata, aes(x = UMAP_1, y = UMAP_2, color = groupinfo)) +
    geom_point(size = point_size,alpha = point_alpha,shape = 16) +
    scale_color_manual(values = cell_color) +
    theme_bw() +
    labs(
      x = "UMAP1",
      y = "UMAP2",
      title = plot_title,
      color = color_title
    ) +
    theme(
      text = element_text(color = "black",family = family_theme),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black",linewidth = 0.5),
      axis.ticks = element_line(colour = "black",linewidth = 0.5),
      axis.text = element_text(size = axis_text_size,face = "bold",family = family_theme),
      plot.title = element_text(size = plot_title_size,hjust = 0.5,face = "bold",family = family_theme),
      axis.title = element_text(size = axis_title_size,face = "bold",family = family_theme),
      legend.title = element_text(size = legend_title_size,face = "bold",family = family_theme),
      legend.text = element_text(size = legend_text_size,face = "bold",family = family_theme),
      legend.key.size = unit(0.3,"cm"),
      legend.key.width = unit(0.1,"cm"),
      legend.key.height = unit(0.2,"cm"),
      legend.key.spacing.y = unit(0.05,"cm"),
      legend.background = element_rect(fill = "transparent"),
      aspect.ratio = 1:1
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1,0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0.13,0.1))) +
    annotate(
      geom = "text",
      x = min(plotdata[,1]) + 0.1*(DescTools::Range(plotdata[,1])),
      y = min(plotdata[,2]),
      label = paste0("Cells:",format(nrow(plotdata),big.mark = ",")),
      size = num_size,
      fontface = "bold",
      family = family_theme
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = legend_size,alpha = 0.8))
    ) +
    coord_fixed() +
    theme(legend.position = legend_position)

  if (autolabel) {
    cell_centers <- plotdata %>%
      group_by(groupinfo) %>%
      summarise(
        umap1_center = mean(UMAP_1),
        umap2_center = mean(UMAP_2),
        .groups = "drop"
      )

    p <- p +
      shadowtext::geom_shadowtext(
        data = cell_centers[!(cell_centers$groupinfo %in% remove_name),],
        mapping = aes(x = umap1_center,y = umap2_center,label = groupinfo),
        size = label_text_size,
        bg.color = text_bg_color,
        fontface = "bold",
        family = family_theme,
        check_overlap = TRUE,
        show.legend = FALSE
      )


  }

  return(p)
}



#' 在ggplot图形上添加带背景阴影的文字标签
#'
#' 该函数用于在已有ggplot对象上叠加高可读性文字标签，特别适用于散点图、UMAP图等需要标记特定位置的应用场景。
#'
#' @param xpos 标签的X轴坐标（必需），单位与图形坐标系统一致
#' @param ypos 标签的Y轴坐标（必需），单位与图形坐标系统一致
#' @param labelname 标签文本内容（默认："label"），支持表达式
#' @param ggplotobj 目标ggplot对象（必需），需提前创建的基础图形
#' @param text_bg_color 文字背景色（默认："black"），增强标签与背景的对比度
#' @param label_text_size 标签文字大小（默认：3），单位pt
#' @param family_theme 字体类型（默认："Times New Roman"），推荐使用等宽字体提升可读性
#'
#' @return 返回修改后的ggplot对象，保留原始图形所有属性并叠加新标签
#' @export
add_label <- function(
    xpos = 0,
    ypos = 0,
    labelname = "label",
    ggplotobj = p,
    text_bg_color = "black",
    label_text_size = 3,
    family_theme = "Times New Roman"
) {
  plotdata = data.frame(
    umap1_center = xpos,
    umap2_center = ypos,
    groupinfo = labelname
  )
  p <- ggplotobj +
    shadowtext::geom_shadowtext(
      data = plotdata,
      mapping = aes(x = umap1_center,y = umap2_center,label = groupinfo),
      size = label_text_size,
      bg.color = text_bg_color,
      fontface = "bold",
      family = family_theme,
      check_overlap = TRUE,
      show.legend = FALSE
    )

  return(p)
}




#' 分组比例堆叠条形图/流行图生成器
#'
#' 该函数基于Seurat对象创建专业级的分组比例可视化图，支持两种模式：
#' - 经典堆叠百分比条形图（默认）
#' - 添加流行图(alluvial)效果展示组间流动关系
#' 特别适用于展示单细胞数据中样本与细胞类型的组成关系，或任何分组变量的比例分布。
#'
#' @param seurat_obj Seurat对象（必需），包含元数据
#' @param group1 主分组变量（默认：sample_ID），作为X轴类别
#' @param group2 次级分组变量（默认：cell_type），作为填充类别
#' @param groupcolor 颜色映射向量（必需），格式：`c("类型1"="#FF0000", ...)`
#' @param legend_title 图例标题（默认："Cell Type"）
#' @param ifflow 是否启用流行图效果（默认：FALSE），TRUE时添加流动曲线
#' @param legend_text_size 图例文本大小（默认：8）
#' @param legend_size 图例色块大小（默认：5），单位mm
#' @param legend_title_size 图例标题大小（默认：8）
#' @param axis_title_size 坐标轴标题大小（默认：10）
#' @param axis_text_size 坐标轴文本大小（默认：10），X轴文本默认旋转45°
#' @param plot_title_size 主标题大小（默认：12）
#' @param family_theme 字体族（默认："Times New Roman"），推荐"Arial"等无衬线字体
#' @param axix_line_width 坐标轴线宽（默认：0.8）
#' @param plot_title 主标题文本（默认：空）
#' @param x_title X轴标题（默认：空）
#' @param y_title Y轴标题（默认：空）
#'
#' @return 返回ggplot2对象，包含以下特征：
#' - 堆叠条形图高度总和为100%
#' - 自动计算组内比例
#' - 学术级绘图主题（白底+黑色轴线）
#' - 专业颜色编码系统
#'
#' @export
percent_bar <- function(
    seurat_obj = seuobj,
    group1 = sample_ID,
    group2 = cell_type,
    groupcolor = main_cell_type_color,
    legend_title = "Cell Type",
    ifflow = FALSE,
    legend_text_size = 8,
    legend_size = 5,
    legend_title_size = 8,
    axis_title_size = 10,
    axis_text_size = 10,
    plot_title_size = 12,
    family_theme = "Times New Roman",
    axix_line_width = 0.8,
    plot_title = "",
    x_title = "",
    y_title = ""
){
  group1names <- enquo(group1)
  group2names <- enquo(group2)

  df_summary <- compute_group_num(seurat_obj,!!group1names,!!group2names)

  colnames(df_summary) <- c("First_group", "Secondary_group","Count")

  # 计算比例
  df_prop <- df_summary %>%
    group_by(First_group) %>%  # 按样本分组
    mutate(proportion = Count / sum(Count))

  cell_type_color <- groupcolor

  plotdata <- df_prop

  p <- ggplot(
    data = plotdata,
    mapping = aes(
      x = First_group,
      y = proportion,
      stratum = Secondary_group,
      alluvium = Secondary_group,
      fill = Secondary_group)
  ) +
    ggalluvial::geom_stratum(width = 0.5,alpha = 1) +
    scale_fill_manual(
      values = cell_type_color
    ) +
    scale_y_continuous(expand = expansion(mult = c(0,0)))+
    labs(x = x_title, y = y_title, fill = legend_title,title = plot_title) +
    theme(
      text = element_text(color = "black",family = family_theme),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),

      axis.line.x = element_line(colour = "black",linewidth = axix_line_width),
      axis.line.y = element_line(color = "black",linewidth = axix_line_width),

      axis.text.x = element_text(color = "black",size = axis_text_size,face = "bold",family = family_theme,hjust = 1,angle = 45),
      axis.text.y = element_text(color = "black",size = axis_text_size,face = "bold",family = family_theme),
      plot.title = element_text(color = "black",size = plot_title_size,hjust = 0.5,face = "bold",family = family_theme),
      axis.title = element_text(color = "black",size = axis_title_size,face = "bold",family = family_theme),
      legend.title = element_text(size = legend_title_size,face = "bold",family = family_theme),
      legend.text = element_text(color = "black",size = legend_text_size,face = "bold",family = family_theme)
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = legend_size,alpha = 0.8))
    )

  if (ifflow) {
    message("添加流行图")
    p <- p +
      ggalluvial::geom_alluvium(alpha = 0.6)
  }

  return(p)
}




#' 分组比例趋势折线图生成器
#'
#' 该函数基于Seurat对象创建时间序列或顺序样本的比例趋势图，直观展示不同分组变量（如细胞类型）在主要分组维度（如时间点/样本批次）中的比例变化趋势。特别适用于动态过程分析（如发育时序、治疗响应等）。
#'
#' @param seurat_obj Seurat对象（必需），包含元数据
#' @param group1 主分组变量（默认：sample_ID），作为X轴变量，应为顺序变量（如时间序列）
#' @param group2 次级分组变量（默认：cell_type），作为折线分组依据
#' @param groupcolor 颜色映射向量（必需），格式：`c("类型1"="#FF0000", ...)`，长度需匹配次级分组数
#' @param legend_title 图例标题（默认："Cell Type"）
#'
#' @return 返回ggplot2对象，包含以下特征：
#' - 带数据点的折线图（点填充白色，增强辨识度）
#' - 学术级绘图主题（透明背景+加粗轴线）
#' - 自动比例计算与颜色编码
#'
#' @export
percent_line <- function(
    seurat_obj = seuobj,
    group1 = sample_ID,
    group2 = cell_type,
    groupcolor = main_cell_type_color,
    legend_title = "Cell Type"
) {
  group1names <- enquo(group1)
  group2names <- enquo(group2)

  df_summary <- compute_group_num(seurat_obj,!!group1names,!!group2names)

  colnames(df_summary) <- c("First_group", "Secondary_group","Count")

  # 计算比例
  df_prop <- df_summary %>%
    group_by(First_group) %>%  # 按样本分组
    mutate(proportion = Count / sum(Count))

  cell_type_color <- groupcolor

  plotdata <- df_prop

  p <- ggplot(
    data = plotdata,
    mapping = aes(
      x = First_group,
      y = proportion,
      group = Secondary_group,
      color = Secondary_group)
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 4,fill = "white",stroke = 1.2) +
    scale_color_manual(
      values = cell_type_color
    ) +
    scale_y_continuous(expand = expansion(mult = c(0,0)))+
    labs(x = "Sample", y = "Proportion", fill = legend_title) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.background = element_blank(),
      plot.title = element_text(),

      panel.background = element_blank(),

      axis.line.x = element_line(colour = "black",linewidth = 1),
      axis.line.y = element_line(color = "black",linewidth = 1),

      axis.text.x = element_text(size = 12,face = "bold",colour = "black",angle = 45,hjust = 1),
      axis.text.y = element_text(size = 12,face = "bold",colour = "black"),

      axis.title.x = element_text(size = 20,face = "bold",colour = "black"),
      axis.title.y = element_text(size = 20,face = "bold",colour = "black")
    )

  return(p)
}




#' UMAP基因表达热图可视化
#'
#' 在单细胞转录组UMAP降维空间中可视化特定基因的表达分布，通过双色渐变系统（零值/非零值分层着色）精确展示基因表达的空间定位特征。支持亚群数据筛选与细胞数量标注。
#'
#' @param seurat_obj Seurat对象（必需），需已完成标准预处理（含UMAP降维）
#' @param plotgene 目标基因名（默认："CHCHD10"），必须存在于`SCT`assay的`data`层
#' @param plot_title 图形主标题（默认与基因名相同）
#' @param plot_title_color 标题颜色（默认："black"）
#' @param point_size UMAP点大小（默认：0.2），推荐范围0.1-1.0
#' @param point_alpha 点透明度（默认：0.5），范围0（全透）~1（不透明）
#' @param legend_title 图例标题（默认："express"）
#' @param color_palette 非零表达点色阶（默认：c("#f4d58b","#800025")），需包含2个及以上颜色值
#' @param zero_point_color 零表达点颜色（默认："gray80"）
#' @param subcol 元数据列名（可选），用于筛选细胞亚群
#' @param subneed 亚群标识（可选），需与`subcol`指定的列中值匹配
#' @param if_label_cell_num 是否标注细胞数量（默认：TRUE），显示于左下角
#' @param x_extend X轴扩展比例（默认：c(0.15,0.1)），控制绘图边界
#' @param y_extend Y轴扩展比例（默认：c(0.15,0.1)），控制绘图边界
#'
#' @return 返回ggplot2对象，包含以下特性：
#' - 双图层点图：零表达点（灰色）与非零表达点（渐变色）分层渲染
#' - 固定纵横比（`coord_fixed()`）保持UMAP空间结构
#' - 学术级主题（无边框+黑色轴线+加粗文本）
#'
#' @export
umap_value <- function(
    seurat_obj = seuobj,
    plotgene = "CHCHD10",
    plot_title = "CHCHD10",
    plot_title_color = "black",
    point_size = 0.2,
    point_alpha = 0.5,
    legend_title = "express",
    color_palette = c("#f4d58b","#800025"),
    zero_point_color = "gray80",
    subcol,
    subneed,
    if_label_cell_num = TRUE,
    x_extend = c(0.15,0.1),
    y_extend = c(0.15,0.1)
) {
  if (!plotgene %in% rownames(seurat_obj@assays[["SCT"]]@data)) {
    stop(plotgene,"not found!!!")
  }
  gene_data <- SeuratObject::LayerData(
    object = seurat_obj,
    assay = "SCT",
    layer = "data",
    features = plotgene
  ) %>%
    as.data.frame() %>%
    t()

  subcolname <- enquo(subcol)

  umap_df <- Embeddings(
    object = seurat_obj,
    reduction = "umap"
  )
  seurat_meta <- seurat_obj@meta.data

  if (missing(subcol)) {
    plotdata <- merge(umap_df,gene_data,by = "row.names")
    plotdata <- column_to_rownames(plotdata,var = colnames(plotdata)[1])

  } else {
    groupinfo <- seurat_meta[,rlang::as_label(subcolname),drop = FALSE]

    colnames(groupinfo) <- c("subgroup")
    umapob <- merge(umap_df,groupinfo,by = "row.names")
    umapob <- column_to_rownames(umapob,var = colnames(umapob)[1])
    colnames(umapob) <- c("UMAP_1","UMAP_2","subgroup")

    umapob <- umapob[umapob$subgroup == subneed,]
    plotdata <- merge(umapob,gene_data,by = "row.names")
    plotdata <- column_to_rownames(plotdata,var = colnames(plotdata)[1])
    plotdata <- plotdata[,c(1,2,4)]
  }


  colnames(plotdata) <- c("X","Y","Value")
  plotdata <- as.data.frame(plotdata)
  p <- ggplot() +
    geom_point(
      data = plotdata[plotdata$Value == 0,],
      mapping = aes(x = X, y = Y),
      color = zero_point_color,
      size = point_size,
      alpha = point_alpha,
      shape = 16
    ) +
    ggnewscale::new_scale_color() +
    geom_point(
      data = plotdata[plotdata$Value != 0,],
      mapping = aes(x = X, y = Y,color = Value),
      size = point_size,
      alpha = point_alpha,
      shape = 16
    ) +
    scale_color_continuous(palette = color_palette) +
    theme_bw() +
    labs(
      x = "UMAP1",
      y = "UMAP2",
      color = legend_title,
      title = plot_title
    ) +
    theme(
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5,face = "bold",color = plot_title_color),
      axis.line = element_line(colour = "black",linewidth = 1),
      axis.ticks = element_line(colour = "black",linewidth = 1),
      axis.text = element_text(size = 10,face = "bold"),
      axis.title = element_text(size = 12,face = "bold"),
      legend.title = element_text(size = 12,face = "bold"),
      legend.text = element_text(size = 10,face = "bold")
    ) +
    scale_y_continuous(expand = expansion(mult = x_extend)) +
    scale_x_continuous(expand = expansion(mult = y_extend)) +
    coord_fixed()

  if (if_label_cell_num == TRUE) {
    p <- p + annotate(
      geom = "text",
      x = min(plotdata[,1]) + 0.1*(DescTools::Range(plotdata[,1])),
      y = min(plotdata[,2]),
      label = paste0("Cells:",format(nrow(plotdata),big.mark = ",")),
      size = 5,
      fontface = "bold"
    )
  }

  return(p)
}






#' 多基因UMAP集成可视化
#'
#' 批量生成多个基因在UMAP空间的表达分布图并整合为面板，适用于标记基因筛选和空间表达模式比对。通过精简坐标系统增强多图对比效率。
#'
#' @param marker_vec 基因向量（必需），需存在于Seurat对象中
#' @param seuratobj Seurat对象（默认：seuobj），需已完成UMAP降维
#' @param color_vec 颜色梯度（默认：c("#f4d58b","#800025")），长度≥2的连续色阶
#' @param colnum 面板列数（默认：4），控制布局密度
#' @param if_byrow 排列顺序（默认：TRUE），TRUE=按行填充，FALSE=按列填充
#'
#' @return 返回cowplot整合的ggplot2面板对象，包含：
#' - 每个基因独立的UMAP表达图（无坐标轴/图例以节省空间）
#' - 统一化的点参数（size=0.05, alpha=0.8）
#' - 进度实时反馈（标记基因名+完成提示）
#'
#' @export
multi_marker_umap <- function(
    marker_vec = main_markers_vec,
    seuratobj = seuobj,
    color_vec = c("#f4d58b","#800025"),
    colnum = 4,
    if_byrow = TRUE
) {
  ggplot_list <- list()

  for (marker in marker_vec) {
    p <- umap_value(
      seurat_obj = seuobj,
      plotgene = marker,
      plot_title = marker,
      plot_title_color = "black",
      point_size = 0.05,
      point_alpha = 0.8,
      legend_title = "express",
      color_palette = color_vec,
      zero_point_color = "gray80",
      if_label_cell_num = FALSE,
      x_extend = c(0,0),
      y_extend = c(0,0)
    ) +
      theme(
        text = element_text(family = "Times New Roman"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
      )

    ggplot_list[[marker]] <- p
    message(marker," plot completed!")
  }

  combined_plot <- cowplot::plot_grid(
    plotlist = ggplot_list,
    ncol = colnum,
    byrow = if_byrow,
    align = "hv"
  )

  return(combined_plot)
}






#' 细胞类型-基因表达气泡图
#'
#' 通过气泡图展示基因在细胞类型中的表达特征，整合平均表达量(Z-score)和表达细胞比例双维度信息。适用于细胞类型标记基因鉴定和差异表达模式分析。
#'
#' @param seurat_obj Seurat对象（默认：seuobj）
#' @param groupcol 分组列（必需），元数据中的细胞类型列
#' @param grouporder 细胞类型顺序向量（可选），指定Y轴顺序
#' @param groupcolor 细胞类型颜色映射（必需），命名向量（如c("T细胞"="#1f77b4")）
#' @param plotmarker 目标基因向量（必需）
#' @param x_text_size X轴文本大小（默认：10）
#' @param y_text_size Y轴文本大小（默认：10），自动应用groupcolor着色
#' @param fill_palette 填充色阶（默认：c("#fffec8","#fc913c","#800025")）
#'
#' @return 返回ggplot2气泡图对象，特性包括：
#' - 气泡大小：表达该基因的细胞百分比
#' - 气泡颜色：Z-score标准化后的平均表达量
#' - 学术级热图式布局（等距坐标+网格引导）
#'
#' @export
bubble_genes <- function(
    seurat_obj = seuobj,
    groupcol = cell_type,
    grouporder = "",
    groupcolor = main_cell_type_color,
    plotmarker = unique_marker_select,
    x_text_size = 10,
    y_text_size = 10,
    fill_palette = c("#fffec8","#fc913c","#800025")
){

  groupname <- enquo(groupcol)

  common_gene <- intersect(plotmarker,rownames(seurat_obj@assays[["SCT"]]@data))
  diff_gene <- setdiff(plotmarker,rownames(seurat_obj@assays[["SCT"]]@data))
  message("These genes not dound:",diff_gene)

  avg_by_subtype <- AverageExpression(
    object = seurat_obj,
    assays = "SCT",
    layer = "data",
    group.by = rlang::as_label(groupname),
    features = common_gene
  )

  avg_by_subtype_long <- avg_by_subtype$SCT %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "cell_type",
      values_to = "average"
    )

  genes <- unique(avg_by_subtype_long$gene)
  genes <- factor(genes,levels = genes)

  meta_data <- seurat_obj@meta.data

  results_list <- lapply(genes, function(gene) {
    expr <- FetchData(seurat_obj, vars = gene)[, 1]
    meta_data$is_expressing <- ifelse(expr > 0, "Positive", "Negative")
    meta_data %>%
      group_by(!!groupname) %>%
      summarise(
        gene = gene,
        proportion = mean(is_expressing == "Positive") * 100
      )
  })

  combined_results <- bind_rows(results_list)
  colnames(combined_results) <- c("cell_type","gene","proportion")

  merged_dot_df <- merge(
    avg_by_subtype_long,
    combined_results,
    by = c("gene","cell_type")
  ) %>%
    group_by(gene) %>%
    mutate(z_score = scale(average)) %>%
    ungroup()


  merged_dot_df$gene <- factor(merged_dot_df$gene,levels = genes)

  merged_dot_df$cell_type <- factor(merged_dot_df$cell_type,levels = grouporder)
  colnames(merged_dot_df)[5] <- "zscore"

  subtype_color <- groupcolor

  plotdata <- merged_dot_df

  p <- ggplot() +
    geom_point(
      data = plotdata,
      mapping = aes(
        x = gene,
        y = cell_type,
        size = proportion,
        fill = zscore
      ),
      shape = 21

    ) +
    labs(
      x = "",
      y = ""
    ) +
    theme_bw() +
    scale_fill_continuous(palette = fill_palette) +
    scale_size_continuous(breaks = c(20,40,60,80),range = c(0.5,4))+
    theme(
      text = element_text(family = "Times New Roman"),
      plot.background = element_rect(color = "white"),
      plot.title = element_text(),

      panel.background = element_rect(color = "white"),
      panel.border = element_rect(color = "black",linewidth = 1),
      panel.grid.major = element_line(colour = "grey80",linewidth = 0.5),
      panel.grid.minor = element_blank(),

      axis.line.x = element_line(colour = "black",linewidth = 0),
      axis.line.y = element_line(color = "black",linewidth = 0),

      axis.text.x = element_text(size = x_text_size,angle = 45,vjust = 1,hjust = 1,face = "bold",family = "Times New Roman"),
      axis.text.y = element_text(size = y_text_size,face = "bold",color = subtype_color,family = "Times New Roman"),

      axis.title.x = element_text(size = 10,face = "bold",colour = "black",family = "Times New Roman"),
      axis.title.y = element_text(size = 10,face = "bold",colour = "black",family = "Times New Roman"),
      legend.text = element_text(size = 8,face = "bold",colour = "black",family = "Times New Roman"),
      legend.title = element_text(size = 10,face = "bold",colour = "black",family = "Times New Roman")
    ) +
    coord_equal()


  return(p)
}





#' 多基因小提琴图面板
#'
#' 生成垂直堆叠的基因表达小提琴图面板，展示目标基因在细胞亚群中的表达分布特征。特别适用于亚型特异性表达模式比对。
#'
#' @param seurat_obj Seurat对象（默认：seuobj）
#' @param subtype_color 细胞亚型颜色映射（必需），命名向量（格式：c("亚型1"="#color")）
#' @param genes_vector 基因向量（必需），需存在于SCT assay
#' @param cell_subtype 细胞亚型列（必需），元数据中的分类列
#'
#' @return 返回垂直堆叠的ggplot2面板对象，包含：
#' - 每个基因独立小提琴图（Y轴=表达量，X轴=细胞亚型）
#' - 首图显示亚型标签（颜色同步subtype_color）
#' - 统一颜色编码与分布带宽（bw=0.1保证分布平滑）
#'
#' @export
violin_marker_plot <- function(
    seurat_obj = seuobj,
    subtype_color = setNames(c("#c0321a","#547bb4","#eab883"),c("Excitatory neuron","Inhibitory neuron","Non neuron")),
    genes_vector = main_markers_vec,
    cell_subtype = sub_cell_type
) {

  cell_subtype_enquo <- enquo(cell_subtype)
  cell_subtype_chr <- rlang::as_label(cell_subtype_enquo)

  plots <- list()
  for (i in 1:length(genes_vector)) {
    gene <- genes_vector[i]
    if (!gene %in% rownames(seurat_obj@assays[["SCT"]]@data)) {
      stop(gene,"not found!!!")
    }

    plotdf <- LayerData(
      object = seurat_obj,
      assay = "SCT",
      layer = "data",
    )[gene,] %>% as.vector()

    df <- data.frame(
      expression = plotdf,
      celltype = seuobj@meta.data[[cell_subtype_chr]]
    )

    group_color <- subtype_color

    p <- ggplot(df,aes(x = celltype,y = expression,fill = celltype)) +
      geom_violin(
        width = 0.8,
        trim = TRUE,
        scale = "width",
        alpha = 1,
        color = "black",
        linewidth = 0.2,
        linetype = "solid",
        bw = 0.1,
        adjust = 1
      ) +
      labs(
        y = gene
      ) +
      scale_fill_manual(values = group_color) +
      scale_y_continuous(breaks = c(0,floor(max(df$expression))),expand = expansion(mult = c(0.01,0.1))) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(color = "white",linewidth = 2),
        plot.title = element_text(color = "black",face = "bold",hjust = 0.5),
        axis.line.x = element_line(colour = "black",linewidth = 1),
        axis.line.y = element_line(color = "black",linewidth = 1),
        axis.text.y = element_text(size = 8,face = "bold",colour = "black"),

        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,face = "bold",colour = "black",angle = 0,vjust = 0.5),
        axis.ticks.length = unit(2,"mm"),
        axis.minor.ticks.y.left = element_blank(),
        axis.ticks = element_line(linewidth = 1)
      )

    if (i == 1) {
      p <- p +
        theme(
          axis.text.x = element_text(size = 8,face = "bold",colour = group_color,angle = 45,vjust = 1,hjust = 1)
        )
    } else {
      p <- p +
        theme(
          axis.text.x = element_blank()
        )
    }

    plots[[i]] <- p
  }

  combined_plot <- cowplot::plot_grid(
    plotlist = plots,
    ncol = 1,
    byrow = TRUE,
    align = "hv"
  )

  return(combined_plot)
}






