#' @Title bivariate_kernel_densityplot_plotting
#'
#' @param plotdata
#' @param x
#' @param y
#' @param group
#' @param grouporder
#' @param groupcolor
#' @param addcontour
#' @param addpoint
#' @param plottitle
#' @param plot_xname
#' @param plot_yname
#'
#' @returns
#' @export
#'
#' @examples
bivariate_kernel_densityplot_plotting <- function(
  plotdata = NULL,
  x = NULL,
  y = NULL,
  group = NULL,
  grouporder = NULL,
  groupcolor = NULL,
  addcontour = TRUE,
  addpoint = FALSE,

  plottitle = "",
  plot_xname = "X",
  plot_yname = "Y"
  ) {

  # logic preprocess
  colnames(plotdata)[colnames(plotdata) == x] <- "x"
  colnames(plotdata)[colnames(plotdata) == y] <- "y"
  colnames(plotdata)[colnames(plotdata) == group] <- "group"


  if (is.null(x)) {
    stop("Please specify the colname x whose type is string.", call. = FALSE)
  }

  if (is.null(y)) {
    stop("Please specify the colname y whose type is string.", call. = FALSE)
  }

  if (is.null(group)) {
    stop("Please specify the colname group whose type is string.", call. = FALSE)
  }

  if (is.null(grouporder)) {
    grouporder <- unique(plotdata$group)
  }

  groupnum <- length(grouporder)

  # assign the groupcolor if user do not specify the groupcolor.
  if (is.null(groupcolor)) {
    message("You do not specify the groupcolor,so use the default color.")
    groupcolor <- rainbow(10)[1:groupnum]
  }

  # Iterative plotting
  looplotobject <- ggplot(plotdata, aes(x = x, y = y))

  # add density layer
  for (num in 1:length(grouporder)) {
    looplotobject <- add_density_layer(
      alldata = plotdata,
      plotobject = looplotobject,
      subgroup = grouporder[num],
      colorvec = colorRampPalette(c("#FFFFFF", groupcolor[num]), space = "rgb")(5),
      densityalpha = 0.7
    )
  }

  # whether or not to add contour lines layer
  if (addcontour == TRUE) {
    for (num in 1:length(grouporder)) {
      looplotobject <- add_contour_layer(
        alldata = plotdata,
        plotobject = looplotobject,
        subgroup = grouporder[num],
        color = colorRampPalette(c("#FFFFFF", groupcolor[num]), space = "rgb")(5)[4],
        densityalpha = 0.15,
        contourbins = 20
      )
    }
  } else {
    message("Do not add the contour lines layer.")
  }


  # whether or not to add points layer
  if (addpoint == TRUE) {
    looplotobject <- looplotobject +
      ggnewscale::new_scale_fill() +
      geom_point(aes(fill = grouporder,color = grouporder), alpha = 0.5, size = 1.5,shape = 21) +
      scale_color_manual(values = groupcolor) +
      scale_fill_manual(values = groupcolor)
  }

  # add theme settings
  looplotobject <- looplotobject +
    labs(
      title = plottitle,
      x = plot_xname,
      y = plot_yname
      ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.background = element_blank(),
      panel.border  = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1.2,
        linetype = "solid"
      ),

      plot.title  = element_text(face = "bold", hjust = 0.5,size = 20),

      axis.title.x = element_text(margin = margin(t = 10),size = 15),
      axis.title.y = element_text(margin = margin(r = 10),size = 15),

      axis.line.x.bottom = element_line(color = "black",linewidth = 0.5),
      axis.line.x.top = element_line(color = "black",linewidth = 0.5),
      axis.line.y.left = element_line(color = "black",linewidth = 0.5),
      axis.line.y.right = element_line(color = "black",linewidth = 0.5),

      axis.text.x = element_text(color = "black",size = 15),
      axis.text.y = element_text(color = "black",size = 15),

      axis.ticks  = element_line(
        colour = "black",
        linewidth = 1,
        lineend = "round"
      ),
      axis.ticks.length  = unit(0.15, "cm"),
      legend.position  = "right",
      legend.key.height  = unit(0.5, "cm"),
    )


  return(looplotobject)

}


add_density_layer <- function(
    alldata = NULL,
    plotobject = baseplotobject,
    subgroup = "Group1",
    colorvec = c("#ffffff","#FFE5E5", "#FFB3B3", "#FF0000", "#B30000"),
    densityalpha = 0.5
) {
  added_density_layer <- plotobject +
    ggnewscale::new_scale_fill() +
    stat_density_2d(
        data = subset(alldata, group == subgroup),
        aes(fill = after_stat(density)),
        geom = "raster",
        contour = FALSE,
        alpha = densityalpha,
        adjust = c(1,1),
        interpolate = TRUE,
        n = 200
      ) +
    scale_fill_gradientn(
        colours = colorvec,
        name = paste0(subgroup," ","Density"),
        guide = guide_colorbar(order = 1)
      )

  return(added_density_layer)
}






add_contour_layer <- function(
    alldata = NULL,
    plotobject = baseplotobject,
    subgroup = "Group1",
    color = "#FFB3B3",
    densityalpha = 0.2,
    contourbins = 15
) {
  added_contour_layer <- plotobject +
    ggnewscale::new_scale_fill() +
    stat_density_2d(
      data = subset(alldata, group == subgroup),
      color = color,
      bins = contourbins,
      contour = TRUE,
      alpha = densityalpha,
      adjust = c(1,1),
      n = 200
    )

  return(added_contour_layer)
}
