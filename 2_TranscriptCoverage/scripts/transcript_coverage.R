# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(patchwork)
library(cowplot)
library(ggtext)

# misc
library(conflicted)

# project managment
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")
conflict_prefer("here", "here")
conflict_prefer("summarize", "dplyr")

# Import ------------------------------------------------------------------

sample_n_groups = function(tbl, size, replace = FALSE, weight = NULL) {
  # regroup when done
  grps = tbl %>% groups %>% lapply(as.character) %>% unlist
  # check length of groups non-zero
  keep = tbl %>% summarise() %>% ungroup() %>% sample_n(size, replace, weight)
  # keep only selected groups, regroup because joins change count.
  # regrouping may be unnecessary but joins do something funky to grouping variable
  tbl %>% right_join(keep, by=grps) %>% group_by_(.dots = grps)
}

bm_coverage <- readRDS(here('..', 'data', 'bm_gene_coverage.rds')) %>%
  mutate(species = 'brugia_malayi')
di_coverage <- readRDS(here('..', 'data', 'di_gene_coverage.rds')) %>%
  mutate(species = 'dirofilaria_immitis')

coverage <- bind_rows(bm_coverage, di_coverage)

rm(bm_coverage, di_coverage)

# Tidy --------------------------------------------------------------------

high_depth_genes <- distinct(filter(coverage, depth > 40), gene_id)

low_depth_genes <- distinct(filter(group_by(coverage, gene_id), max(depth) < 5), gene_id)

fully_covered_genes <- group_by(coverage, gene_id) %>% 
  filter(index == max(index) | index == 1) %>%
  mutate(index = case_when(
    index == 1 ~ "start",
    TRUE ~ "stop"
  )) %>%
  pivot_wider(names_from = index, values_from = depth) %>%
  filter(start == stop) %>%
  distinct(gene_id)

filtered_genes <- distinct(bind_rows(high_depth_genes, low_depth_genes, fully_covered_genes))

mid_depth_genes <- coverage %>% 
  filter(!gene_id %in% filtered_genes$gene_id) 

gene_coverage <- mid_depth_genes %>%
  group_by(gene_id) %>%
  # flip index for genes on the negative strand
  mutate(index = case_when(
    strand == "-" ~ max(index) - index + 1,
    TRUE ~ index
  )) %>%
  # normalize index from 0 to 1
  mutate(index_norm = index / max(index)) %>%
  # normalize depth from 0 to 1
  mutate(depth_norm = depth / max(depth), depth_norm = replace_na(depth_norm, 0))

# remove all the duplicated data, get the start/stop of each depth segment
gene_reduced <- ungroup(gene_coverage) %>%
  group_by(gene_id, depth) %>%
  filter(index_norm == min(index_norm) | index_norm == max(index_norm)) %>%
  mutate(bin = cut(index_norm, breaks = seq(0, 1, 0.01), labels = FALSE))

# Plot --------------------------------------------------------------------

label <- tribble(~species, ~x, ~y, ~color,
                 "*B. malayi*", 0.75, 0.4, 'black',
                 "*D. immitis*", 0.5, 0.8, 'grey')

line_plot <- ggplot(gene_reduced, aes(x = index_norm, y = depth_norm)) +
  geom_smooth(aes(color = species), se = FALSE) +
  geom_richtext(data = label, aes(x = x, y = y, label = species, color = color),
                fill = NA, label.color = NA, size = 3) +
  scale_color_manual(values = c('black', 'black', 'grey', 'grey')) +
  lims(y = c(0, 1)) +
  labs(x = "Gene Position (5' to 3')", y = "Normalized Depth") +
  theme_minimal(base_family = "Helvetica", base_size = 16) +
  theme(
    axis.ticks = element_line(size = 0.25),
    axis.line = element_line(size = 0.3),
    legend.position = "none"
  ) +
  NULL

save_plot(here('..', 'plots', 'Fig1D.pdf'), line_plot, base_height = 4, base_width = 4)

saveRDS(line_plot, here('..', 'plots', 'Fig1D.rds'))
