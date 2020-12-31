library(ggplot2)

# Define custom plot theme
StandardPlotTheme <- ggplot2::theme(
  text = ggplot2::element_text(
    family = 'Helvetica',
    face = 'plain',
    size = 8,
    color = '#141414'
  ),
  axis.title = ggplot2::element_text(
    family = 'Helvetica',
    face = 'plain',
    size = 10,
    hjust = 0.5,
    color = '#141414'
  ),
  axis.text.x = ggplot2::element_text(
    family = 'Helvetica',
    face = 'plain',
    size = 8,
    color = '#141414'
  ),
  axis.text.y = ggplot2::element_text(
    family = 'Helvetica',
    face = 'plain',
    size = 8,
    color = '#141414'
  ),
  axis.ticks = ggplot2::element_line(
    colour = '#d3d3d3',
    size = 0.2,
    lineend = 'round'
  ),
  axis.line = ggplot2::element_line(
    colour = '#d3d3d3',
    size = 0.2,
    lineend = 'round'
  ),
  plot.title = ggplot2::element_text(
    family = 'mono',
    face = 'plain',
    size = 12,
    hjust = 0.5,
    color = '#141414'
  )
)

#Define color palette
StandardFills = c('#89cff0', '#fff8dc', '#c9ffe5', '#fa6e79', 'f5f5dc', '#5d3f6a')
