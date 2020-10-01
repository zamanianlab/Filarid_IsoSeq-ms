# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(ggbeeswarm)
library(ggridges)
library(patchwork)
library(cowplot)
library(ggtext)
library(LaCroixColoR)

# misc
library(conflicted)

# project managment
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")
conflict_prefer("here", "here")
conflict_prefer("summarize", "dplyr")
conflict_prefer("lag", "dplyr")

# Import ------------------------------------------------------------------

# structural classifications
bm_classifications <- read_delim(here("..", "data", "bm_all_transcripts_trim_classification.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "brugia_malayi") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  mutate(structural_category = str_to_title(str_replace_all(structural_category, "_", " "))) %>%
  filter(rts_stage == FALSE, structural_category != "Genic Intron", !str_detect(isoform, '_dup'))

di_classifications <- read_delim(here("..", "data", "di_all_transcripts_trim_classification.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "dirofilaria_immitis") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  mutate(structural_category = str_to_title(str_replace_all(structural_category, "_", " "))) %>%
  filter(rts_stage == FALSE, structural_category != "Genic Intron", !str_detect(isoform, '_dup'))

classifications <- bind_rows(bm_classifications, di_classifications) %>%
  filter(structural_category != "Genic Intron", !str_detect(isoform, '_dup'))

# junctions
bm_junctions <- read_delim(here("..", "data", "bm_all_transcripts_trim_junctions.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "brugia_malayi") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  group_by(isoform) %>%
  filter(!any(rts_junction == TRUE), !any(bite_junction == TRUE)) %>%
  ungroup()

di_junctions <- read_delim(here("..", "data", "di_all_transcripts_trim_junctions.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "dirofilaria_immitis") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  group_by(isoform) %>%
  filter(!any(rts_junction == TRUE), !any(bite_junction == TRUE)) %>%
  ungroup()

junctions <- bind_rows(bm_junctions, di_junctions)

# Colors ------------------------------------------------------------------

pamplemousse <- lacroix_palette("Pamplemousse", 9)
fsm_ism_palette <- pamplemousse[c(2,5)]

peach <- lacroix_palette("PeachPear", 3)

order <- ordered(c("Incomplete-Splice Match", "Full-Splice Match", "Novel Not In Catalog", "Intergenic", "Novel In Catalog", "Antisense", "Genic", "Fusion"))

# Structural categories ---------------------------------------------------

structural_summary <- group_by(classifications, species, structural_category) %>%
  summarise(category_total = n()) %>%
  group_by(species) %>%
  mutate(category_percent = category_total / sum(category_total)) %>%
  mutate(species = case_when(
    species == 'brugia_malayi' ~ "*B. malayi*",
    species == 'dirofilaria_immitis' ~ "*D. immitis*",
  ))

structural_plot <- ggplot(structural_summary, aes(x = structural_category)) +
  geom_col(aes(y = category_total, fill = structural_category)) +
  geom_text(aes(y = category_total + 800, label = round(category_percent, digits = 3), color = structural_category),
             size = 3) +
  scale_x_discrete(limits = order) +
  scale_y_continuous(limits = c(0, 21000), expand = c(0, 0)) +
  scale_color_manual(values = pamplemousse) +
  scale_fill_manual(values = pamplemousse) +
  facet_grid(cols = vars(species)) +
  labs(y = "Filtered CCSs") +
  theme_minimal(base_family = "Helvetica", base_size = 16) +
  theme(
    axis.line = element_line(size = 0.3),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_line(size = 0.25),
    strip.text = element_markdown(),
    legend.position = 'bottom',
    legend.justification = 'left',
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  ) +
  NULL

save_plot(here('..', 'plots', "Fig2E.pdf"), structural_plot, base_width = 8.5, base_height = 4)

saveRDS(structural_plot, here('..', 'plots', 'Fig2E.rds'))

# Tallies -----------------------------------------------------------------

total <- group_by(classifications, species, sex) %>%
  summarise(sum_af = sum(sex == "AF"), sum_am = sum(sex == "AM")) %>%
  mutate(total = sum_af + sum_am) %>%
  select(-sum_af, -sum_am) %>%
  mutate(species = case_when(
    species == 'brugia_malayi' ~ "*B. malayi*",
    species == 'dirofilaria_immitis' ~ "*D. immitis*",
  )) %>%
  ungroup() %>%
  mutate(y = c(45000, 18000, 50000, 18000))

ccs_by_sex <- ggplot(total) +
  geom_bar(aes(x = species, y = total, fill = sex), stat = "identity", alpha = 0.95) +
  geom_text(aes(x = species, y = y, label = total),
             color = "white", size = 3) +
  scale_fill_manual(limits = c("AF", "AM"), values = peach) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Sex", y = "Filtered CCSs", fill = "Sex") +
  facet_grid(cols = vars(species), scales = 'free') +
  theme_minimal(base_family = "Helvetica", base_size = 16) +
  theme(
    axis.line = element_line(size = 0.3),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_line(size = 0.25),
    panel.spacing.x = unit(0, "null"),
    strip.text = element_markdown(),
    legend.position = "none"
  ) +
  NULL

save_plot(here('..', 'plots', "Fig2A.pdf"), ccs_by_sex, base_width = 3, base_height = 3)

saveRDS(ccs_by_sex, here('..', 'plots', 'Fig2A.rds'))

# collapse into reads with identitical structure and summed by sex (1 means that isoform is found in that sex, 0 means it isn't or it's in both)
transcripts_mapped <- group_by(classifications, species, chrom, associated_gene, associated_transcript) %>%
  summarise(sum_af = sum(sex == "AF"), sum_am = sum(sex == "AM")) %>%
  mutate(isoform_category = case_when(
    sum_af > 0 & sum_am > 0 ~ "Both",
    sum_af > 0 & sum_am == 0 ~ "Female Only",
    sum_af == 0 & sum_am > 0 ~ "Male Only"
  )) %>%
  select(-sum_af, -sum_am)

# total number of collapsed reads (non-redundant isoforms) for each sex (or both)
total_transcripts_mapped <- ungroup(transcripts_mapped) %>%
  group_by(species, isoform_category) %>%
  summarise(total_transcripts_by_sex = n()) %>%
  mutate(species = case_when(
    species == 'brugia_malayi' ~ "*B. malayi*",
    species == 'dirofilaria_immitis' ~ "*D. immitis*",
  )) %>%
  ungroup() %>%
  mutate(y = c(14000, 9500, 3300, 17500, 10500, 2500))

# the number of distince reference transcripts 
transcripts_by_sex <- ggplot(total_transcripts_mapped) +
  geom_bar(aes(x = species, y = total_transcripts_by_sex, fill = isoform_category), stat = "identity") + 
  geom_text(aes(x = species, y = y, label = total_transcripts_by_sex), 
             color = "white", size = 3) +
  scale_fill_manual(limits = c("Female Only", "Male Only", "Both"), values = peach) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Sex", y = "Unique Transcripts") +
  facet_grid(cols = vars(species), scales = 'free') +
  theme_minimal(base_family = "Helvetica", base_size = 16) +
  theme(
    axis.line = element_line(size = 0.3),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_line(size = 0.25),
    panel.spacing.x = unit(0, "null"),
    strip.text = element_markdown(),
    legend.position = "none"
  ) +
  NULL

save_plot(here('..', 'plots', "Fig2C.pdf"), transcripts_by_sex, base_width = 3, base_height = 3)

saveRDS(transcripts_by_sex, here('..', 'plots', 'Fig2C.rds'))

#  reference genes did we mapped in each sex
gene_summary <- ungroup(classifications) %>%
  group_by(species, sex, associated_gene) %>%
  filter(!str_detect(associated_gene, "novel")) %>%
  summarise(ccs_per_gene = n()) %>%
  pivot_wider(id_cols = c(species, associated_gene), names_from = sex, values_from = ccs_per_gene) %>%
  mutate(isoform_category = case_when(
    AF > 0 & AM > 0 ~ "Both",
    AF > 0 & is.na(AM) ~ "Female Only",
    is.na(AF) & AM > 0 ~ "Male Only"
  )) %>%
  ungroup() %>%
  group_by(species, isoform_category) %>%
  summarise(genes_per_sex = n()) %>%
  mutate(species = case_when(
    species == 'brugia_malayi' ~ "*B. malayi*",
    species == 'dirofilaria_immitis' ~ "*D. immitis*",
  )) %>%
  ungroup() %>%
  mutate(y = c(4500, 2450, 700, 4400, 2300, 600))

genes_by_sex <- ggplot(gene_summary) +
  geom_bar(aes(x = species, y = genes_per_sex, fill = isoform_category), stat = "identity") + 
  geom_text(aes(x = species, y = y + 100, label = genes_per_sex),
             color = "white", size = 3) +
  # scale_x_discrete(limits = c("Female Only", "Male Only", "Both")) +
  scale_fill_manual(limits = c("Female Only", "Male Only", "Both"), values = peach) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Reference Genes") +
  facet_grid(cols = vars(species), scales = 'free') +
  theme_minimal(base_family = "Helvetica", base_size = 16) +
  theme(
    axis.line = element_line(size = 0.3),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_line(size = 0.25),
    strip.text = element_markdown(),
    panel.spacing.x = unit(0, "null"),
    legend.position = "none",
    legend.title = element_blank()
  ) +
  NULL

save_plot(here('..', 'plots', "Fig2B.pdf"), genes_by_sex, base_width = 3, base_height = 3)

saveRDS(genes_by_sex, here('..', 'plots', 'Fig2B.rds'))

# Gene model extension ----------------------------------------------------

# reference GTFs
gtf <- bind_rows(as_tibble(readRDS(here('..', 'data', 'bm_gtf.rds')) %>% mutate(species = 'brugia_malayi')),
                 as_tibble(readRDS(here('..', 'data', 'di_gtf.rds')) %>% mutate(species = 'dirofilaria_immitis'))
)

# alignment GTFs
alignment <- bind_rows(as_tibble(readRDS(here('..', 'data', 'bm_isoseq_gtf.rds')) %>% mutate(species = 'brugia_malayi')),
                       as_tibble(readRDS(here('..', 'data', 'di_isoseq_gtf.rds')) %>% mutate(species = 'dirofilaria_immitis'))
)

# get start and stop codons for every transcript
start_stop <- filter(gtf, feature %in% c("start_codon", "stop_codon"), frame == 0) %>%
  mutate(associated_gene = str_remove(attribute, ";.*"),
         associated_gene = str_remove_all(associated_gene, "\""),
         associated_gene = str_remove(associated_gene, "gene_id "),
         associated_transcript = str_remove(attribute, ".*transcript_id "),
         associated_transcript = str_remove(associated_transcript, ";.*"),
         associated_transcript = str_remove_all(associated_transcript, "\"")) %>%
  select(species, chrom = scaffold, associated_gene, associated_transcript, feature, start, end, -attribute, -source, -score, -strand, -frame) %>%
  pivot_wider(id_cols = species:associated_transcript, names_from = feature, values_from = start) %>%
  select(species:associated_transcript, start_codon, stop_codon)

# get beginning and end for every read
alignment_begin_end <- group_by(alignment, isoform) %>%
  summarise(start = min(start), end = max(end)) %>%
  select(isoform, alignment_begin = start, alignment_end = end)

# calculate distance to start and stop codons
classifications <- left_join(classifications, start_stop) %>%
  left_join(., alignment_begin_end) %>%
  mutate(diff_to_start_codon = case_when(
    strand == "+" ~ alignment_begin - start_codon, # negative includes 5' UTR
    strand == "-" ~ start_codon - alignment_end
  )) %>%
  mutate(diff_to_stop_codon = case_when(
    strand == "+" ~ stop_codon - alignment_end, # negative includes 3' UTR
    strand == "-" ~ alignment_begin - stop_codon # negative includes 3' UTR
  )) %>%
  mutate(species = case_when(
    species == 'brugia_malayi' ~ "*B. malayi*",
    species == 'dirofilaria_immitis' ~ "*D. immitis*",
  ))

totals <- group_by(classifications, species, structural_category) %>%
  summarise(total = n())

# What fraction of FSM + ISM had extended 5' UTR?
percent_true_tss <- group_by(filter(classifications, structural_category %in% c("Full-Splice Match", "Incomplete-Splice Match")), species, structural_category) %>%
  filter(diff_to_tss < 0) %>%
  summarise(number = n()) %>%
  left_join(., totals) %>%
  mutate(percent = number / total, x = -500, y = c(3500, 3000))

extended_five <- ggplot(classifications,
                        aes(x = diff_to_tss, fill = structural_category)) +
  geom_histogram(data = filter(classifications, structural_category == "Full-Splice Match"), alpha = 0.75, center = 0, bins = 50) +
  geom_histogram(data = filter(classifications, structural_category == "Incomplete-Splice Match"), alpha = 0.75, center = 0, bins = 50) +
  geom_vline(xintercept = 0, linetype = 18) + 
  geom_text(data = percent_true_tss, aes(x = x, y = y, color = structural_category,
                                         label = paste0(round(percent * 100, digits = 0), "%; ", number)),
            size = 3) +
  scale_color_manual(values = fsm_ism_palette) +
  scale_fill_manual(values = fsm_ism_palette) +
  scale_x_continuous(limits = c(-1000, 1000), breaks = c(-500, 0, 500), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 5000), expand = c(0, 0)) +
  labs(x = "Distance of CCS Start to TSS", y = "Filtered CCSs") +
  annotate("segment", x = -50, xend = -690, y = 4000, yend = 4000, color = "grey40", size = 0.5, arrow = arrow(type = 'closed', length = unit(0.30, "cm"))) +
  annotate("text", x = -550, y = 4600, label = "Extends 5' UTR", size = 3) +
  facet_grid(cols = vars(species)) +
  theme_minimal(base_family = "Helvetica", base_size = 12) +
  theme(strip.text = element_markdown(),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
  ) +
  NULL

save_plot(here('..', 'plots', "Fig3A.pdf"), extended_five, base_width = 5, base_height = 3)

saveRDS(extended_five, here('..', 'plots', 'Fig3A.rds'))

# What fraction of FSM + ISM had extended 3' UTR?
percent_true_tts <- group_by(filter(classifications, structural_category %in% c("Full-Splice Match", "Incomplete-Splice Match")), species,  structural_category) %>%
  filter(diff_to_tts < 0) %>%
  summarise(number = n()) %>%
  left_join(., totals) %>%
  mutate(percent = number / total, x = -500, y = c(3500, 3000))

extended_three <- ggplot(classifications,
                         aes(x = diff_to_tts, fill = structural_category)) +
  geom_histogram(data = filter(classifications, structural_category == "Full-Splice Match"), alpha = 0.75, center = 0) +
  geom_histogram(data = filter(classifications, structural_category == "Incomplete-Splice Match"), alpha = 0.75, center = 0) +
  geom_vline(xintercept = 0, linetype = 18) + 
  geom_text(data = percent_true_tts, aes(x = x, y = y, color = structural_category,
                                         label = paste0(round(percent * 100, digits = 0), "%; ", number)), 
            size = 3) +
  scale_fill_manual(values = fsm_ism_palette) +
  scale_color_manual(values = fsm_ism_palette) +
  scale_x_reverse(limits = c(1000, -1000), breaks = c(-500, 0, 500), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 5000), expand = c(0, 0)) +
  labs(x = "Distance of CCS End to TTS", y = "Filtered CCSs") +
  annotate("segment", x = -50, xend = -690, y = 4000, yend = 4000, color = "grey40", size = 0.5, arrow = arrow(type = 'closed', length = unit(0.30, "cm"))) +
  annotate("text", x = -550, y = 4600, label = "Extends 3' UTR", size = 3) +
  facet_grid(cols = vars(species)) +
  theme_minimal(base_family = "Helvetica", base_size = 12) +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        strip.text = element_markdown(),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
  ) +
  NULL

save_plot(here('..', 'plots', "Fig3B.pdf"), extended_three, base_width = 5, base_height = 3)

saveRDS(extended_three, here('..', 'plots', 'Fig3B.rds'))

legend <- get_legend(extended_five + theme(legend.position = "bottom", legend.title = element_blank()))

transcriptional_lengthiness <- plot_grid(extended_five, extended_three, nrow = 2, labels = "AUTO", label_size = 12)
transcriptional_lengthiness <- plot_grid(transcriptional_lengthiness, legend, nrow = 2, rel_heights = c(1, 0.1))

save_plot(here('..', 'plots', "Fig3.pdf"), transcriptional_lengthiness, base_width = 5, base_height = 5)
save_plot(here('..', 'plots', "Fig3.png"), transcriptional_lengthiness, base_width = 5, base_height = 5)


# Junctions ---------------------------------------------------------------


# structural classifications
old_bm_classifications <- read_delim(here("..", "data", "old_bm_all_transcripts_trim_classification.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "brugia_malayi") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  mutate(structural_category = str_to_title(str_replace_all(structural_category, "_", " "))) %>%
  filter(rts_stage == FALSE, bite == FALSE | is.na(bite))

old_di_classifications <- read_delim(here("..", "data", "old_di_all_transcripts_trim_classification.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "dirofilaria_immitis") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  mutate(structural_category = str_to_title(str_replace_all(structural_category, "_", " "))) %>%
  filter(rts_stage == FALSE, bite == FALSE | is.na(bite))

old_classifications <- bind_rows(old_bm_classifications, old_di_classifications) 

# junctions
old_bm_junctions <- read_delim(here("..", "data", "old_bm_all_transcripts_trim_junctions.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "brugia_malayi") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  group_by(isoform) %>%
  filter(!any(rts_junction == TRUE), !any(bite_junction == TRUE)) %>%
  ungroup()

old_di_junctions <- read_delim(here("..", "data", "old_di_all_transcripts_trim_junctions.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "dirofilaria_immitis") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  group_by(isoform) %>%
  filter(!any(rts_junction == TRUE), !any(bite_junction == TRUE)) %>%
  ungroup()

old_junctions <- bind_rows(old_bm_junctions, old_di_junctions)

# junction location on transcript (have to use SQANTI1 output)
junction_features <- old_junctions %>%
  left_join(., old_classifications) %>%
  select(species:junction_category, splice_site, canonical, read_length = length, exons, structural_category, associated_gene, associated_transcript, ref_length) %>%
  mutate(rel_transcript_coord = case_when(
    is.na(read_length) ~ as.double(NA),
    !is.na(read_length) ~ transcript_coord / read_length
  ))  %>%
  filter(!is.na(structural_category)) %>%
  mutate(site_category = case_when(
    transcript_coord > -1 & transcript_coord < 41 ~ "A",
    transcript_coord > 40 & transcript_coord < 81 ~ "B",
    transcript_coord > 80 & transcript_coord < 121 ~ "C",
    transcript_coord > 120 & transcript_coord < 161 ~ "D",
    transcript_coord > 160 & transcript_coord < 201 ~ "E",
    transcript_coord > 200 ~ "F"
  )) %>%
  mutate(structural_category = factor(structural_category, levels = order)) %>%
  mutate(species = case_when(
    species == 'brugia_malayi' ~ "*B. malayi*",
    species == 'dirofilaria_immitis' ~ "*D. immitis*",
  ))

rel_junction_stats <- filter(junction_features, 
                             species == '*B. malayi*',
                             structural_category %in% c("Full-Splice Match", "Incomplete-Splice Match", "Novel Not In Catalog", "Novel In Catalog"))
  
model <- lm(rel_transcript_coord ~ structural_category, data = rel_junction_stats)
aov <- aov(model)
tukey <- broom::tidy(TukeyHSD(aov))

rel_junctions_coords <- ggplot(filter(junction_features, 
                                      species == '*B. malayi*',
                                      structural_category %in% c("Full-Splice Match", "Incomplete-Splice Match", "Novel Not In Catalog", "Novel In Catalog"))) +
  stat_density_ridges(aes(x = rel_transcript_coord, y = structural_category, fill = structural_category),
                      quantile_lines = TRUE) +
  scale_fill_manual(values = pamplemousse[c(8, 9, 5, 2)],
                    limits = rev(c("Full-Splice Match", "Incomplete-Splice Match", "Novel Not In Catalog", "Novel In Catalog"))) +
  scale_y_discrete(limits = c("Full-Splice Match", "Incomplete-Splice Match", "Novel Not In Catalog", "Novel In Catalog"),
                   expand = expansion(add = c(0, 1.05))) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Splice Junction Postion on CCS", y = "Fraction Total", fill = "Structural Category") +
  theme_minimal(base_family = "Helvetica", base_size = 12) +
  theme(axis.line = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.25),
        axis.text.y = element_blank(),
        # legend.position = "none",
        plot.margin = margin(l = 15, r = 15, t = 15, b = 15)
  ) +
  NULL

save_plot(here('..', 'plots', "Fig4B.pdf"), rel_junctions_coords, base_width = 5, base_height = 5)

saveRDS(rel_junctions_coords, here('..', 'plots', 'Fig4B.rds'))

junction_sites <- group_by(junction_features, species, structural_category, site_category) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(species, structural_category) %>%
  mutate(total_junctions = sum(count))

junction_structural_plot <- ggplot(filter(junction_sites, 
                                          species == '*B. malayi*',
                                          structural_category %in% c("Full-Splice Match", "Incomplete-Splice Match", "Novel Not In Catalog", "Novel In Catalog"))) +
  geom_col(aes(x = structural_category, y = count / total_junctions, fill = site_category)) +
  scale_x_discrete(limits = c("Full-Splice Match", "Incomplete-Splice Match", "Novel Not In Catalog", "Novel In Catalog"),
                   labels = c("Full-Splice\nMatch", "Incomplete-\nSplice Match", "Novel Not\nIn Catalog", "Novel In\nCatalog")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = lacroix_palette("PassionFruit", n = 6), 
                    labels = c("0-40", "41-80", "81-120", "121-160", "161-200", ">200")) +
  labs(y = "Fraction Total", fill = "Splice Junction\nPosition on CCS (nt)") +
  theme_minimal(base_family = "Helvetica", base_size = 12) +
  theme(axis.line = element_line(size = 0.3),
        axis.ticks.y = element_line(size = 0.25),
        axis.title.y = element_blank()
  ) +
  coord_flip() +
  NULL

save_plot(here('..', 'plots', "Fig4C.pdf"), junction_structural_plot, base_width = 5, base_height = 5)

saveRDS(junction_structural_plot, here('..', 'plots', 'Fig4C.rds'))

# junction types
junction_summary <- left_join(junctions, classifications) %>%
  select(species, chrom, strand, genomic_start_coord, genomic_end_coord, junction_category, start_site_category, end_site_category, splice_site, canonical) %>%
  distinct() %>%
  group_by(species, junction_category, splice_site, canonical) %>%
  summarise(count = n()) %>%
  filter(!str_detect(splice_site, 'N')) %>%
  arrange(desc(count)) %>%
  mutate(junction_category = case_when(
    junction_category == 'known' ~ 'Known',
    junction_category == 'novel' ~ 'Novel'
  )) %>%
  mutate(canonical = case_when(
    canonical == 'canonical' ~ 'Canonical',
    canonical == 'non_canonical' ~ 'Non-Canonical'
  )) %>%
  mutate(splice_site = str_replace(splice_site, 'T', 'U')) %>%
  bind_rows(., tibble(species = 'brugia_malayi', junction_category = 'Known', splice_site = 'GUAA', canonical = 'Non-Canonical', count = 0))

plot_x <- filter(junction_summary, species == 'brugia_malayi') %>%
  group_by(splice_site) %>%
  summarise(total = sum(count)) %>%
  slice_max(total, n = 10)

plot_order <- factor(plot_x$splice_site, ordered = TRUE)

plot_data <- filter(junction_summary, species == 'brugia_malayi' & splice_site %in% plot_x$splice_site) %>%
  mutate(splice_site = factor(splice_site, levels = plot_order))

new_junctions <- ggplot(plot_data) +
  geom_col(aes(x = splice_site, y = log2(count), fill = junction_category),
           position = position_dodge()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = lacroix_palette("PeachPear", type = "discrete")[c(4, 6)]) +
  facet_grid(cols = vars(canonical), scales = 'free_x', space = 'free', drop = TRUE) +
  labs(x = "Splice Sequence", y = "log<sub>2</sub>(Count)") +
  theme_minimal(base_family = "Helvetica", base_size = 12) +
  theme(
    axis.line = element_line(size = 0.3),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    axis.ticks.y = element_line(size = 0.25),
    axis.title.y = element_markdown(),
    legend.position = "right",
    legend.title = element_blank()
  ) +
  NULL

save_plot(here('..', 'plots', "Fig4A.pdf"), new_junctions, base_width = 6, base_height = 4)

fig4 <- plot_grid(new_junctions + theme(legend.position = "right",
                                        legend.text = element_text(size = 8)),
                  rel_junctions_coords + theme(legend.position = "right",
                                               legend.text = element_text(size = 8),
                                               legend.title = element_text(size = 10)),
                  junction_structural_plot + theme(legend.position = "right",
                                                   legend.text = element_text(size = 8),
                                                   legend.title = element_text(size = 10)), 
                  nrow = 3, align = "v", axis = "lr", rel_heights = c(1, 1, 0.75), 
                  labels = "AUTO", label_size = 12)

save_plot(here('..', 'plots', "Fig4.pdf"), fig4, base_width = 6.5, base_height = 8)
save_plot(here('..', 'plots', "Fig4.png"), fig4, base_width = 6.5, base_height = 8)

saveRDS(junctions_by_structure, here('..', 'plots', "Fig4.rds"))
