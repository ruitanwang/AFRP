#' Title
#'
#' @param plotdata
#' @param x
#' @param group
#' @param grouporder
#' @param groupfill
#' @param grouplabel
#' @param plottitle
#' @param plot_xname
#' @param plot_yname
#' @param legendname
#'
#' @returns
#' @export
#'
#' @examples
univar_density_barplot_plotting <- function(
    plotdata = data,
    x = NULL,
    group = NULL,
    grouporder = NULL,
    groupfill = NULL,
    grouplabel = NULL,

    plottitle = "",
    plot_xname = "Value",
    plot_yname = "Density",

    legendname = "Groups"
) {

  # logic preprocess
  colnames(plotdata)[colnames(plotdata) == x] <- "x"
  colnames(plotdata)[colnames(plotdata) == group] <- "group"

  if (is.null(x)) {
    message("Please specify the colname x whose type is string.")
    break
  }

  if (is.null(group)) {
    message("Please specify the colname group whose type is string.")
    break
  }

  if (is.null(grouporder)) {
    message("You do not specify the order of group,the default order will be used.")

    grouporder <- unique(plotdata$group)
  }

  if (is.null(groupfill)) {
    message("You do not specify the groupfill,the default fill will be used.")
    groupnum <- length(grouporder)
    groupfill <- RColorBrewer::brewer.pal(groupnum, "Set2")
  }

  if (is.null(grouplabel)) {
    message("You do not specify the group labels,the default group names will be used.")
    grouplabel <- grouporder
  }



  # plotting
  plotobject <- ggplot(plotdata, aes(x = x, fill = factor(group,levels = grouporder))) +
    geom_histogram(
      mapping = aes(y = after_stat(density)),
      alpha = 0.6,
      linewidth = 0.8,
      position = "identity",
      bins = 50,
      color = "black",
      lineend = "round",
      linejoin = "mitre") +
    scale_fill_manual(
      name = legendname,
      values = groupfill,
      labels = grouplabel
    ) +
    labs(
      title = plottitle,
      x = plot_xname,
      y = plot_yname,
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), expand = expansion(add = 0)) +
    theme(
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.background = element_blank(),

      plot.title  = element_text(face = "bold", hjust = 0.5,size = 20),

      axis.title.x = element_text(margin = margin(t = 10),size = 15),
      axis.title.y = element_text(margin = margin(r = 10),size = 15),

      axis.line.x = element_line(color = "black",linewidth = 0.5),
      axis.line.y = element_line(color = "black",linewidth = 0.5),

      axis.text.x = element_text(color = "black",size = 15),
      axis.text.y = element_text(color = "black",size = 15),

      axis.minor.ticks.length.x = unit(2,"cm"),
      legend.position  = "right",
      legend.key.height  = unit(0.5, "cm"),
    )

  return(plotobject)

}
