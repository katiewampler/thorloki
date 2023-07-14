#' Add pretty stream layer to ggplot2
#'
#' Uses an NHDplus stream layer to add steams to a geospatial ggplot, uses
#' the visibility column in the stream layer to make bigger streams thicker on
#' the plot.
#' @importFrom displease seq_ease
#'
#' @param plot the ggplot2 object you want to add the streams to
#' @param stream the sf layer with streams (should be an NHDplus layer)
#' @param min the line size of the smallest stream
#' @param max the line size of the largest streams
#' @param color the color of the streams
#' @param min_vis the minimum visibility value to display (visibility column divided by 10000)
#'
#' @return a plot with the streams added
#' @export
#'
add_streams <- function(plot, stream, min=0.2, max=0.7, color="dodgerblue4",
                        min_vis = 10){
  #convert visibility to be more user friendly
  stream$visibility <- stream$visibility/10000

  #find levels
  levels <- unique(stream$visibility)
  levels <- levels[levels >=min_vis]

  #figure out widths
  linewidths <- displease::seq_ease(min, max, n=length(levels), type = 'quad-in')

  #create df of levels
  df <- data.frame(width = linewidths, level=levels)

  for(x in 1:length(levels)){
    data_add <- subset(stream, stream$visibility == df$level[x])
    plot <- plot + geom_sf(data=data_add, color=color, linewidth=df$width[x])
  }
  plot
}
