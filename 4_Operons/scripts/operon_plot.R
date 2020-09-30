library(cowplot)
library(ggplot)

A <- readRDS(here('..', 'plots', 'Fig5A.rds')) + theme(axis.title = element_text(size = 10),
                                                       axis.text = element_text(size = 8),
                                                       legend.position = 'none')

B <- readRDS(here('..', 'plots', 'Fig5B.rds')) + theme(axis.title = element_text(size = 10),
                                                       axis.text = element_text(size = 8),
                                                       legend.position = 'none')

C <- readRDS(here('..', 'plots', 'Fig5C.rds')) + theme(axis.title = element_text(size = 10),
                                                       axis.text = element_text(size = 8),
                                                       legend.position = 'none')

D <- readRDS(here('..', 'plots', 'Fig5D.rds')) + theme(axis.title = element_text(size = 10),
                                                       axis.text = element_text(size = 8))

E <- readRDS(here('..', 'plots', 'Fig5E.rds'))  + theme(axis.title = element_text(size = 10),
                                                        axis.text = element_text(size = 8))

F <- readRDS(here('..', 'plots', 'Fig5F.rds'))  + theme(axis.title = element_text(size = 10),
                                                        axis.text = element_text(size = 8))

G <- readRDS(here('..', 'plots', 'Fig5G.rds'))  + theme(axis.title = element_text(size = 10),
                                                        axis.text = element_text(size = 8))

top <- plot_grid(A, B, C,
                 nrow = 1, rel_widths = c(0.7, 1, 0.75),
                 labels = c('A', 'B', 'C'), axis = 'b', align = 'h',
                 scale = 0.95
)

top_legend <- plot_grid(get_legend(A + theme(legend.position = "top")), top,
                         nrow = 2, rel_heights = c(0.1, 1))

middle <- plot_grid(D, E, F,
                    nrow = 1, rel_widths = c(1, 1, 1),
                    labels = c('D', 'E', 'F'), vjust = 0.7, axis = 'b', align = 'h',
                    scale = 0.95
)

bottom <- plot_grid(G, get_legend(G + theme(legend.position = 'bottom')),
                    nrow = 2, rel_heights = c(1, 0.1), 
                    labels = c('G', ''),
                    scale = 0.95
)

final <- plot_grid(top_legend, NULL, middle, NULL, bottom,
                   nrow = 5, rel_heights = c(1, 0.05, 1, 0.05, 1.1))


save_plot(here('..', 'plots', 'Fig5.pdf'), final, base_height = 8.75, base_width = 7.5)
save_plot(here('..', 'plots', 'Fig5.png'), final, base_height = 8.75, base_width = 7.5)
