suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
  library(tidyr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep('^--file=', args, value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(sub('^--file=', '', script_arg[1])) else '.'
root <- normalizePath(file.path(script_dir, '..'), mustWork = FALSE)
setwd(root)

dir.create('results/figures', recursive = TRUE, showWarnings = FALSE)

screen_colors <- list(
  essential = '#E64B35',
  enriched = '#3C5488',
  ns = 'grey80',
  synthetic_lethal = '#E64B35',
  resistance = '#3C5488',
  common_essential = '#F39B7F',
  crispri = '#4DBBD5',
  crispra = '#91D1C2'
)

theme_screen <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = 'grey20', linewidth = 0.4),
      axis.ticks = element_line(color = 'grey20', linewidth = 0.3),
      axis.text = element_text(color = 'grey20'),
      plot.title = element_text(face = 'bold', size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1),
      legend.title = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank()
    )
}

save_plot <- function(filename, plot, width, height) {
  ggsave(filename, plot, width = width, height = height, dpi = 300)
}
