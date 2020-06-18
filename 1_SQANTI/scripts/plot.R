library(cowplot)
library(ggtext)
library(ggplot2)
library(here)

# Import ------------------------------------------------------------------

a <- readRDS(here('..', 'plots', 'Fig1A.rds')) + theme(axis.title.y = element_text(size = 10),
                                                       strip.text.x = element_markdown(size = 8),
                                                       axis.text = element_text(size = 8))
b <- readRDS(here('..', 'plots', 'Fig1B.rds')) + theme(axis.title.y = element_text(size = 10),
                                                       strip.text.x = element_markdown(size = 8),
                                                       axis.text = element_text(size = 8))
c <- readRDS(here('..', 'plots', 'Fig1C.rds')) + theme(axis.title.y = element_text(size = 10),
                                                       strip.text.x = element_markdown(size = 8),
                                                       axis.text = element_text(size = 8))
d <- readRDS(here('..', '..', '2_TranscriptCoverage', 'plots', 'Fig1D.rds')) + theme(axis.title.y = element_text(size = 10),
                                                                                     axis.title.x = element_text(size = 8),
                                                                                     strip.text.x = element_markdown(size = 8),
                                                                                     axis.text = element_text(size = 8))
e <- readRDS(here('..', 'plots', 'Fig1E.rds')) + theme(axis.title = element_text(size = 10),
                                                       strip.text.x = element_markdown(size = 8),
                                                       axis.text = element_text(size = 10),
                                                       legend.text = element_text(size = 8)) +
  scale_y_continuous(limits = c(0, 22000), expand = c(0, 0))

##### Final figure

top_left <- plot_grid(a, b, c,
                 nrow = 1, align = 'h', axis = 'tb',
                 rel_widths = c(1, 1, 1.03),
                 labels = c('A', 'B', 'C'), label_size = 12)

t_legend <- get_legend(b + theme(legend.position = 'bottom', legend.text = element_text(size = 8)))

top_left_legend <- plot_grid(top_left, t_legend,
                        nrow = 2, rel_heights = c(1, 0.1))

top <- plot_grid(top_left_legend, d,
                 nrow = 1, rel_widths = c(2.75, 1),
                 labels = c('', 'D'), label_size = 12)

fig1 <- plot_grid(top, e,
                  nrow = 2, rel_heights = c(1, 1.5),
                  labels = c('', 'E'), label_size = 12)

save_plot(here('..', 'plots', 'Fig1.pdf'), fig1, base_width = 7.5, base_height = 5)
save_plot(here('..', 'plots', 'Fig1.png'), fig1, base_width = 7.5, base_height = 5)


