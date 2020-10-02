# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(ggridges)
library(patchwork)
library(cowplot)
library(ggtext)
library(LaCroixColoR)

# stats 
library(broom)

# misc
library(conflicted)

# project management
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")
conflict_prefer("here", "here")
conflict_prefer("summarize", "dplyr")
conflict_prefer("select", "dplyr")

# function to choose groups rather than rows
sample_n_groups = function(tbl, size, replace = FALSE, weight = NULL) {
  # regroup when done
  grps = tbl %>% groups %>% lapply(as.character) %>% unlist
  # check length of groups non-zero
  keep = tbl %>% summarise() %>% ungroup() %>% sample_n(size, replace, weight)
  # keep only selected groups, regroup because joins change count.
  # regrouping may be unnecessary but joins do something funky to grouping variable
  tbl %>% right_join(keep, by=grps) %>% group_by(.dots = grps)
}

# Import data -------------------------------------------------------------

final_operonic_genes <- readRDS(here('..', 'data', "final_operonic_genes.rds")) %>%
  ungroup() %>%
  # for simplicity, only keep the first 2 cistrons of any operon
  filter(cistron <= 2) %>%
  group_by(operon_id) %>%
  mutate(distance = first(distance)) %>%
  filter(operon_type != 'Assembly')

rnaseq <- readRDS(here('..', 'data', 'staged_RNAseq.rds')) %>%
  clean_names() %>%
  filter(metric == "TPM_g")

gtf <- readRDS(here('..', '..', '1_SQANTI', 'data', 'bm_gtf.rds')) %>%
  filter(feature == "CDS") %>%
  filter(str_detect(attribute, 'protein_coding')) %>%
  mutate(gene_id = str_remove(attribute, ";.*")) %>%
  mutate(gene_id = str_remove_all(gene_id, "\"")) %>%
  mutate(gene_id = str_remove(gene_id, "gene_id ")) %>%
  mutate(transcript_id = str_remove(attribute, ".*transcript_id ")) %>%
  mutate(transcript_id = str_remove(transcript_id, ";.*")) %>%
  mutate(transcript_id = str_remove_all(transcript_id, "\"")) %>%
  mutate(transcript_id = str_remove(transcript_id, "transcript_id ")) %>%
  mutate(transcript_id = str_remove(transcript_id, ".[0-9]*$")) %>%
  # only keep principal isoform
  filter(!str_detect(transcript_id, "[b-z]$"))

# have to use CDS start/end because of overlapping UTRs 
cds_start_end <- gtf %>%
  group_by(scaffold, strand, gene_id, transcript_id) %>%
  summarise(cds_start = min(start), cds_end = max(end)) %>%
  arrange(scaffold, strand, cds_start)

set.seed(1234)
# create df of genes that occur next to each other within 5000 nt and don't include genes in the putative operon list
pseudo_operons <- filter(cds_start_end, !gene_id %in% final_operonic_genes$gene_id) %>%
  ungroup() %>%
  slice(2:n()) %>%
  mutate(pseudo_operon = rep(1:(length(gene_id) / 2), each = 2)) %>%
  mutate(pseudo_operon = paste0('BMPSO', pseudo_operon)) %>%
  group_by(scaffold) %>%
  # remove pseudo-operons at the beginning and end of each scaffold
  filter(pseudo_operon != max(pseudo_operon), pseudo_operon != min(pseudo_operon)) %>% 
  group_by(pseudo_operon) %>%
  mutate(distance = max(cds_start) - min(cds_end)) %>%
  filter(distance > 0, distance < 5000) %>%
  group_by(pseudo_operon) %>%
  mutate(cistron = row_number()) %>%
  # randomly select the same number of operons as are in filtered_operons
  sample_n_groups(size = length(unique(final_operonic_genes$operon_id))) 

# Tidy --------------------------------------------------------------------

operons_rnaseq <- left_join(final_operonic_genes, rnaseq, by = "gene_id") %>%
  select(gene_id, operon_id, sample_name, distance, cistron, expression) %>%
  group_by(gene_id, operon_id, distance, sample_name, cistron) %>%
  summarize(mean_expression = mean(expression)) %>%
  ungroup() %>%
  # convert to wide format where cistrons 1 and 2 are in separate columns
  pivot_wider(id_cols = c(operon_id, distance, sample_name), 
              names_from = cistron, 
              values_from = mean_expression, 
              names_prefix = "cistron_") %>%
  mutate(difference = cistron_1 - cistron_2, mean = (cistron_1 + cistron_2) / 2) %>%
  mutate(type = "operon") %>%
  rename(name = operon_id)

pseudo_operons_rnaseq <- left_join(pseudo_operons, rnaseq) %>%
  select(gene_id, pseudo_operon, sample_name, distance, cistron, expression) %>%
  group_by(gene_id, pseudo_operon, sample_name, distance, cistron) %>%
  summarize(mean_expression = mean(expression)) %>%
  ungroup() %>%
  # convert to wide format where cistrons 1 and 2 are in separate columns
  pivot_wider(id_cols = c(pseudo_operon, sample_name, distance), names_from = cistron, values_from = mean_expression, names_prefix = "cistron_") %>%
  mutate(difference = cistron_1 - cistron_2, mean = (cistron_1 + cistron_2) / 2) %>%
  mutate(type = "pseudooperon") %>%
  rename(name = pseudo_operon)

full <- bind_rows(operons_rnaseq, pseudo_operons_rnaseq) %>%
  filter(difference != 0)


# Model -------------------------------------------------------------------

# create linear model where cistron_1 expression at each stage is the predictor for cistron_2
lm_model <- group_by(full, name, type) %>%
  arrange(sample_name) %>%
  group_nest() %>%
  mutate(
    fit = map(data, ~ lm(cistron_2 ~ cistron_1, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance)) 

lm_glanced <- unnest(lm_model, glanced) %>%
  select(-data, -fit, -tidied) %>%
  left_join(., distinct(dplyr::select(full, type, name, distance))) 

lm_glanced <- lm_glanced %>%
  mutate(distance_bin = cut(lm_glanced$distance, breaks = c(0, 500, 1000, 1500, 2000, 7000)))

group_by(lm_glanced, type) %>%
  summarise(median = median(r.squared))

lm_tidied <- unnest(lm_model, tidied) %>%
  select(-data, -fit, -glanced) %>%
  left_join(., distinct(select(full, type, name, distance)))

# Plot --------------------------------------------------------------------

example <- group_by(full, name, type) %>%
  filter(name %in% c('BMOP2395', 'BMPSO3802')) %>%
  left_join(., distinct(rnaseq, sample_name, dev_stage_s)) %>%
  mutate(label = case_when(
    dev_stage_s == 'MF/L1' ~ 'L1',
    sample_name == 'AM' ~ 'AM',
    sample_name == 'AF' ~ 'AF',
    sample_name == 'E' ~ 'E',
    TRUE ~ as.character(dev_stage_s)
  )) 

example_line <- left_join(example, lm_tidied) %>%
  distinct(term, estimate, .keep_all = TRUE) %>%
  pivot_wider(id_cols = c('name', 'type'), names_from = term, values_from = estimate) %>%
  rename(intercept = `(Intercept)`) %>%
  left_join(., lm_glanced) %>%
  dplyr::select(name:r.squared) %>%
  ungroup() %>%
  mutate(x = 50, y = c(200, 150))

(example_plot <- ggplot(example, aes(x = cistron_1, y = cistron_2)) +
    geom_abline(data = example_line, aes(slope = cistron_1, intercept = intercept, color = type)) +
    geom_point(aes(color = type), size = 6) +
    geom_text(aes(label = label), color = 'white', fontface = 'bold', size = 2) +
    geom_richtext(data = example_line, aes(x = x, y = y, color = type,
                                           label = paste(expression(R^2), " = ", round(r.squared, 4))),
                  fontface = 'bold', size = 2.5, fill = NA) +
    scale_color_manual(values = c('black', 'grey')) +
    scale_fill_manual(values = c('black', 'grey')) +
    labs(x = "Expression of Cistron 1 (TPM)", y = "Expression of Cistron 2 (TPM)") +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      legend.position = "none",
      legend.title = element_blank()
    ) +
    NULL)

save_plot(here('..', 'plots', "Fig5E.pdf"), example_plot, base_height = 5, base_width = 5)

saveRDS(example_plot, here('..', 'plots', "Fig5E.rds"))

(rsquare_plot <- ggplot(filter(lm_glanced, !is.na(r.squared))) +
    geom_density_ridges(aes(x = r.squared, y = type, fill = type, height = ..density..)
    ) +
    scale_fill_manual(values = c('black', 'grey'), labels = c("Operon", "Pseudo-Operon")) +
    scale_y_discrete(limits = c("operon", "pseudooperon"), labels = c("Operon", "Pseudo-\noperon"), expand = expansion(c(0, 2))) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = expression(R^2), y = "Density") +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      axis.text.y = element_blank(),
      legend.position = "none",
      legend.title = element_blank()
    )
)

save_plot(here('..', 'plots', "Fig5F.pdf"), rsquare_plot, base_height = 3, base_width = 5)

saveRDS(rsquare_plot, here('..', 'plots', "Fig5F.rds"))

(rsquare_distance <- ggplot(filter(lm_glanced, !is.na(statistic), !is.na(distance)), aes(x = distance_bin, y = r.squared)) +
    gghalves::geom_half_boxplot(aes(fill = type), position = 'dodge', outlier.shape = NA) +
    gghalves::geom_half_point(aes(color = type), transformation = position_jitter(height = 0)) +
    scale_fill_manual(values = c('black', 'grey')) +
    scale_color_manual(values = c('black', 'grey')) +
    scale_x_discrete(labels = c('0-500', '501-1000', '1001-1500', '1501-2000', '>2000')) +
    labs(x = "Distance Between Genes", y = expression(R^2)) +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      legend.position = "none",
      legend.title = element_blank()
    ) +
    NULL)

save_plot(here('..', 'plots', "Fig5G.pdf"), rsquare_distance, base_height = 4, base_width = 8)

saveRDS(rsquare_distance, here('..', 'plots', "Fig5G.rds"))

