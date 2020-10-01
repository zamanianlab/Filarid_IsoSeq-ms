# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
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

# Import ------------------------------------------------------------------

# SQANTI classifications
classifications <- read_delim(here('..', '..', "1_SQANTI", 'data', "bm_all_transcripts_trim_classification.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "brugia_malayi") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  mutate(structural_category = str_to_title(str_replace_all(structural_category, "_", " "))) %>%
  filter(structural_category != "Genic Intron", !str_detect(isoform, '_dup'))

# reference GTF
gtf <- readRDS(here('..', '..', "1_SQANTI", 'data', "bm_gtf.rds")) %>%
  mutate(gene_id = str_remove(attribute, ";.*"),
         gene_id = str_remove_all(gene_id, "\""),
         gene_id = str_remove(gene_id, "gene_id "),
         transcript_id = str_remove(attribute, ".*transcript_id "),
         transcript_id = str_remove(transcript_id, ";.*"),
         transcript_id = str_remove_all(transcript_id, "\"")) %>%
  mutate(transcript_id = case_when(
    feature == "gene" ~ as.character(NA),
    TRUE ~ transcript_id
  ))

# BED file of junctions for each read
bed <- readRDS(here('..', 'data', "brugia_malayi_gene_intersects.rds")) %>%
  select(!contains("null")) %>%
  mutate(start = start + 1) # change from base-0 to base-1

# CSV of target metadata
targets <- read_delim(here('..', 'data', "S4_Table.csv"),
                      delim = ",",
                      col_names = TRUE)

# Tidy -------------------------------------------------------------------

pamplemousse <- lacroix_palette("Pamplemousse", 9)

# get SQANTI data for each target
targets_tidy <- left_join(targets, classifications, by = c("gene_id" = "associated_gene")) %>%
  filter(!is.na(isoform))

targets_collapsed <- group_by(targets_tidy, gene_id, transcript_id, name, type, subtype, family, structural_category) %>%
  summarize(supporting_reads = n())


##### Figure 7

# slo-1 -------------------------------------------------------------------

slo_1 <- filter(targets_tidy, name == "slo-1")

slo_1_exons <- filter(gtf, gene_id == unique(slo_1$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  filter(str_detect(object, "[a,f].1")) %>%
  mutate(structural_category = "Reference Transcripts") %>%
  mutate(sex = case_when(
    str_detect(object, "a.1") ~ "AF",
    str_detect(object, "f.1") ~ "AM",
  )) %>%
  bind_rows(., filter(., str_detect(object, "f.1")) %>% mutate(sex = "AF"))

slo_1_bed <- filter(bed, isoform %in% slo_1$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) 

slo_1_features <- bind_rows(slo_1_exons, slo_1_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

slo_1_order <- reorder(slo_1_bed$object, slo_1_bed$object, function(x) - length(x))
slo_1_order <- as.vector(levels(slo_1_order))
slo_1_order <- c(unique(slo_1_exons$object), slo_1_order)

slo_1_features <- bind_rows(slo_1_exons, slo_1_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) %>%
  mutate(object = factor(object, levels = slo_1_order))

slo_1_as <- tribble(~sex, ~xmin, ~xmax, ~ymin, ~ymax,
                    "Male", 10092486, 10092598, 0, length(distinct(filter(slo_1_features, sex == "Male"), object)$object),
                    "Female", 10092486, 10092598, 0, length(distinct(filter(slo_1_features, sex == "Female"), object)$object),
                    "Male", 10091847, 10091960, 0, length(distinct(filter(slo_1_features, sex == "Male"), object)$object),
                    "Female", 10091847, 10091960, 0, length(distinct(filter(slo_1_features, sex == "Female"), object)$object))

slo_1_label <- tribble(~sex, ~x, ~y, ~label,
                       "Female", 10085000, 23.1, "Sex-specific<br>alternative splicing")

(slo_1_viewer <- ggplot(slo_1_features) +
    geom_rect(data = slo_1_as, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.6, fill = "orange") +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 4, alpha = 0.85) +
    labs(x = "Chr X Position (Mb)", color = "Structural Category", title = "*Bma-slo-1*") +
    scale_color_manual(values = c(pamplemousse[c(2, 5, 7, 8)], "grey60")) +
    scale_y_discrete() +
    scale_x_continuous(labels = function(x) (x / 10E6)) +
    facet_grid(rows = vars(sex), space = "free", scales = "free", drop = TRUE) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "bottom"
    ) +
    NULL
)

save_plot(here('..', 'plots', 'Fig7A.pdf'), slo_1_viewer, base_width = 4)
saveRDS(slo_1_viewer, here('..', 'plots', 'Fig7A.rds'))

# osm-9 ---------------------------------------------------------------

osm_9 <- filter(targets_tidy, name == "osm-9")

osm_9_exons <- filter(gtf, str_detect(transcript_id, "Bm1711[a-z].[0-9]")) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

osm_9_bed <- filter(bed, isoform %in% osm_9$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  arrange(object, desc(start)) 

osm_9_features <- bind_rows(osm_9_exons, osm_9_bed)

osm_9_order <- reorder(osm_9_bed$object, osm_9_bed$object, function(x) - length(x))
osm_9_order <- as.vector(levels(osm_9_order))
osm_9_order <- c(unique(osm_9_exons$object), osm_9_order)

osm_9_label <- tribble(~x, ~y, ~label,
                       9343590, 10.85, "Corrected<br>splice-sites")

(osm_9_viewer <- ggplot(osm_9_features) +
    annotate("rect", xmin = 9337065, xmax = 9337109, ymin = 0, ymax = 12, alpha = 0.6, fill = "orange") +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 4, alpha = 0.85) +
    labs(x = "Chr X Position (Mb)", color = "Structural Category", title = "*Bma-osm-9*") +
    scale_color_manual(values = c(pamplemousse[c(8)], "grey60")) +
    scale_y_discrete(limits = osm_9_order) +
    scale_x_reverse(limits = c(9344000, 9336000), minor_breaks = seq(9344000, 9336000, -500), labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

save_plot(here('..', 'plots', 'Fig7B.pdf'), osm_9_viewer, base_width = 4)
saveRDS(osm_9_viewer, here('..', 'plots', 'Fig7B.rds'))

# gar-3 ---------------------------------------------------------------

gar_3 <- filter(targets_tidy, name == "gar-3")

gar_3_exons <- filter(gtf, gene_id == "WBGene00233845") %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

gar_3_bed <- filter(bed, isoform %in% gar_3$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  arrange(object, desc(start)) 

gar_3_features <- bind_rows(gar_3_exons, gar_3_bed)

gar_3_order <- reorder(gar_3_bed$object, gar_3_bed$object, function(x) - length(x))
gar_3_order <- as.vector(levels(gar_3_order))
gar_3_order <- c(unique(gar_3_exons$object), gar_3_order)

gar_3_label <- tribble(~x, ~y, ~label,
                      1874630, 3.9, "Validated<br>isoforms")

(gar_3_viewer <- ggplot(gar_3_features) +
    annotate("rect", xmin = 1872500, xmax = 1864500, ymin = 3.75, ymax = 4.25, alpha = 0.4, fill = "steelblue") +
    annotate("rect", xmin = 1872500, xmax = 1864500, ymin = 1.75, ymax = 2.25, alpha = 0.4, fill = "steelblue") +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 4, alpha = 0.85) +
    labs(x = "Chr X Position (Mb)", color = "Structural Category", title = "*Bma-gar-3*") +
    scale_color_manual(values = c(pamplemousse[c(2, 5)], "grey60")) +
    scale_y_discrete(limits = gar_3_order) +
    scale_x_reverse(limits = c(1875000, 1864000), minor_breaks = seq(1875000, 1865000, -500), labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

save_plot(here('..', 'plots', 'Fig7C.pdf'), gar_3_viewer, base_width = 4)
saveRDS(gar_3_viewer, here('..', 'plots', 'Fig7C.rds'))

# Bm17517 ------------------------------------------------------------------

Bm17517 <- filter(targets_tidy, str_detect(transcript_id, "Bm17517"))

Bm17517_exons <- filter(gtf, gene_id == unique(Bm17517$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

Bm17517_bed <- filter(bed, isoform %in% Bm17517$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

Bm17517_features <- bind_rows(Bm17517_exons, Bm17517_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

Bm17517_order <- reorder(Bm17517_bed$object, Bm17517_bed$object, function(x) - length(x))
Bm17517_order <- as.vector(levels(Bm17517_order))
Bm17517_order <- c(unique(Bm17517_exons$object), Bm17517_order)

Bm17517_label <- tribble(~x, ~y, ~label,
                               7909100, 1.5, "Corrected/alternative<br>TSS and ATG")

(Bm17517_viewer <- ggplot(Bm17517_features) +
    # ggtext::geom_richtext(data = Bm17517_label, aes(x = x, y = y, label = label), color = "white", fill = "grey30", size = 3) +
    annotate("rect", xmin = 7908500, xmax = 7906800, ymin = 1.85, ymax = 2.15, alpha = 0.4, fill = "steelblue") +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 3, alpha = 0.85) +
    labs(x = "Chr 2 Position (Mb)", color = "Structural Category", title = "*Bm17517* (srt)") +
    scale_color_manual(values = c(pamplemousse[c(5)], "grey60")) +
    scale_y_discrete(limits = Bm17517_order) +
    scale_x_reverse(limits = c(7909500, 7906800), labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

save_plot(here('..', 'plots', 'Fig7D.pdf'), Bm17517_viewer, base_width = 4)
saveRDS(Bm17517_viewer, here('..', 'plots', 'Fig7D.rds'))


# Fig 7 -------------------------------------------------------------------

legend <- get_legend(slo_1_viewer + theme(legend.title = element_blank(),
                                          legend.text = element_text(size = 8)))

right <- plot_grid(osm_9_viewer + theme(legend.position = "none"), 
                   gar_3_viewer + theme(legend.position = "none"),
                   Bm17517_viewer  + theme(legend.position = "none"),
                   nrow = 3, labels = c("B", "C", "D"),
                   rel_heights = c(1.5, 1, 0.85), hjust = 0.25)

panel <- plot_grid(slo_1_viewer + theme(legend.position = "none"),
                   right,
                   # align = "h", axis = "b",
                   nrow = 1, labels = c("A", ""), hjust = 0.25)

panel_legend <- plot_grid(panel, legend, nrow = 2, rel_heights = c(1, 0.1), scale = 0.96)

save_plot(here('..', 'plots', "Fig7.pdf"), panel_legend, base_width = 7.5, base_height = 7)

##### Supplements

# btubs --------------------------------------------------------------------

# btub-1
btub_1 <- filter(targets_tidy, name == "btub-1")

btub_1_exons <- filter(gtf, gene_id == unique(btub_1$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

btub_1_bed <- filter(bed, isoform %in% btub_1$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) 

btub_1_features <- bind_rows(btub_1_exons, btub_1_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

btub_1_order <- reorder(btub_1_bed$object, btub_1_bed$object, function(x) - length(x))
btub_1_order <- as.vector(levels(btub_1_order))
btub_1_order <- c(unique(btub_1_exons$object), btub_1_order)

(btub_1_viewer <- ggplot(btub_1_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 2.5, alpha = 0.85) +
    labs(x = "Chr 2 Position (Mb)", color = "Structural Category", title = "*btub-1*") +
    scale_color_manual(values = c(pamplemousse[c(2, 5, 8)], "grey60")) +
    scale_x_continuous(labels = function(x) (x / 10E6)) +
    scale_y_discrete(limits = btub_1_order) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "grey90", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "bottom"
    ) +
    NULL
)

# btub-2
btub_2 <- filter(targets_tidy, name == "btub-2")

btub_2_exons <- filter(gtf, gene_id == unique(btub_2$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

btub_2_bed <- filter(bed, isoform %in% btub_2$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) 

btub_2_features <- bind_rows(btub_2_exons, btub_2_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

btub_2_order <- reorder(btub_2_bed$object, btub_2_bed$object, function(x) - length(x))
btub_2_order <- as.vector(levels(btub_2_order))
btub_2_order <- c(unique(btub_2_exons$object), btub_2_order)

(btub_2_viewer <- ggplot(btub_2_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 2.5, alpha = 0.85) +
    labs(x = "Chr 2 Position (Mb)", color = "Structural Category", title = "*btub-2*") +
    scale_color_manual(values = c(pamplemousse[c(2, 4, 5, 8)], "grey60")) +
    scale_x_continuous(labels = function(x) (x / 10E6)) +
    scale_y_discrete(limits = btub_2_order) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "grey90", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "bottom", 
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    ) +
    NULL
)

btub_legend <- get_legend(btub_1_viewer)

btubs <- plot_grid(btub_1_viewer + theme(legend.position = "none"), btub_2_viewer + theme(legend.position = "none"),
                   align = "h", axis = "tb", labels = "AUTO",
                   rel_widths = c(1, 0.6), hjust = -1)

btubs_final <- plot_grid(btubs, btub_legend,
                         nrow = 2, rel_heights = c(1, 0.1))

save_plot(here('..', 'plots', "S2_Fig.pdf"), btubs_final, base_width = 8.5, base_height = 6)

# acr-12 ------------------------------------------------------------------

acr_12 <- filter(targets_tidy, name == "acr-12")

acr_12_exons <- filter(gtf, gene_id == unique(acr_12$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

acr_12_bed <- filter(bed, isoform %in% acr_12$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

acr_12_features <- bind_rows(acr_12_exons, acr_12_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

acr_12_order <- reorder(acr_12_bed$object, acr_12_bed$object, function(x) - length(x))
acr_12_order <- as.vector(levels(acr_12_order))
acr_12_order <- c(unique(acr_12_exons$object), acr_12_order)

(acr_12_viewer <- ggplot(acr_12_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr X Position (Mb)", color = "Structural Category", title = "*Bma-acr-12*") +
    scale_color_manual(values = c(pamplemousse[c(8)], "grey60")) +
    scale_x_continuous(labels = function(x) (x / 10E6)) +
    scale_y_discrete(limits = acr_12_order) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# acr-8 -------------------------------------------------------------------

acr_8 <- filter(targets_tidy, str_detect(name, "acr-8"))

acr_8_exons <- filter(gtf, gene_id == unique(acr_8$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

acr_8_bed <- filter(bed, isoform %in% acr_8$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

acr_8_features <- bind_rows(acr_8_exons, acr_8_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

acr_8_order <- reorder(acr_8_bed$object, acr_8_bed$object, function(x) - length(x))
acr_8_order <- as.vector(levels(acr_8_order))
acr_8_order <- c(unique(acr_8_exons$object), acr_8_order)

(acr_8_viewer <- ggplot(acr_8_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr X Position (Mb)", color = "Structural Category", title = "*Bma-acr-8*") +
    scale_color_manual(values = c(pamplemousse[c(2, 5, 8)], "grey60")) +
    scale_y_discrete(limits = acr_8_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# avr-14 ------------------------------------------------------------------

avr_14 <- filter(targets_tidy, str_detect(name, "avr-14"))

avr_14_exons <- filter(gtf, gene_id == unique(avr_14$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

avr_14_bed <- filter(bed, isoform %in% avr_14$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

avr_14_features <- bind_rows(avr_14_exons, avr_14_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

avr_14_order <- reorder(avr_14_bed$object, avr_14_bed$object, function(x) - length(x))
avr_14_order <- as.vector(levels(avr_14_order))
avr_14_order <- c(unique(avr_14_exons$object), avr_14_order)

(avr_14_viewer <- ggplot(avr_14_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 3 Position (Mb)", color = "Structural Category", title = "*Bma-avr-14*") +
    scale_color_manual(values = c(pamplemousse[c(5, 7)], "grey60")) +
    scale_x_continuous(labels = function(x) (x / 10E6)) +
    scale_y_discrete(limits = avr_14_order) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# lgc-27 ------------------------------------------------------------------

lgc_27 <- filter(targets_tidy, str_detect(name, "lgc-27"))

lgc_27_exons <- filter(gtf, gene_id == unique(lgc_27$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

lgc_27_bed <- filter(bed, isoform %in% lgc_27$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

lgc_27_features <- bind_rows(lgc_27_exons, lgc_27_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

lgc_27_order <- reorder(lgc_27_bed$object, lgc_27_bed$object, function(x) - length(x))
lgc_27_order <- as.vector(levels(lgc_27_order))
lgc_27_order <- c(unique(lgc_27_exons$object), lgc_27_order)

(lgc_27_viewer <- ggplot(lgc_27_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 2 Position (Mb)", color = "Structural Category", title = "*Bma-lgc-27*") +
    scale_color_manual(values = c(pamplemousse[c(2, 5, 8)], "grey60")) +
    scale_y_discrete(limits = lgc_27_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# unc-38 ------------------------------------------------------------------

unc_38 <- filter(targets_tidy, str_detect(name, "unc-38"))

unc_38_exons <- filter(gtf, gene_id == unique(unc_38$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

unc_38_bed <- filter(bed, isoform %in% unc_38$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

unc_38_features <- bind_rows(unc_38_exons, unc_38_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

unc_38_order <- reorder(unc_38_bed$object, unc_38_bed$object, function(x) - length(x))
unc_38_order <- as.vector(levels(unc_38_order))
unc_38_order <- c(unique(unc_38_exons$object), unc_38_order)

(unc_38_viewer <- ggplot(unc_38_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 3 Position (Mb)", color = "Structural Category", title = "*Bma-unc-38*") +
    scale_color_manual(values = c(pamplemousse[c(5)], "grey60")) +
    scale_y_discrete(limits = unc_38_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# lgc-44 ------------------------------------------------------------------

lgc_44 <- filter(targets_tidy, str_detect(name, "lgc-44"))

lgc_44_exons <- filter(gtf, gene_id == unique(lgc_44$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

lgc_44_bed <- filter(bed, isoform %in% lgc_44$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

lgc_44_features <- bind_rows(lgc_44_exons, lgc_44_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

lgc_44_order <- reorder(lgc_44_bed$object, lgc_44_bed$object, function(x) - length(x))
lgc_44_order <- as.vector(levels(lgc_44_order))
lgc_44_order <- c(unique(lgc_44_exons$object), lgc_44_order)

(lgc_44_viewer <- ggplot(lgc_44_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 2 Position (Mb)", color = "Structural Category", title = "*Bma-lgc-44*") +
    scale_color_manual(values = c(pamplemousse[c(2)], "grey60")) +
    scale_y_discrete(limits = lgc_44_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# lgc-45 ------------------------------------------------------------------

lgc_45 <- filter(targets_tidy, str_detect(name, "lgc-45"))

lgc_45_exons <- filter(gtf, gene_id == unique(lgc_45$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

lgc_45_bed <- filter(bed, isoform %in% lgc_45$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

lgc_45_features <- bind_rows(lgc_45_exons, lgc_45_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

lgc_45_order <- reorder(lgc_45_bed$object, lgc_45_bed$object, function(x) - length(x))
lgc_45_order <- as.vector(levels(lgc_45_order))
lgc_45_order <- c(unique(lgc_45_exons$object), lgc_45_order)

(lgc_45_viewer <- ggplot(lgc_45_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 2 Position (Mb)", color = "Structural Category", title = "*Bma-lgc-45*") +
    scale_color_manual(values = c(pamplemousse[c(7)], "grey60")) +
    scale_y_discrete(limits = lgc_45_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# glc-4 ------------------------------------------------------------------

glc_4 <- filter(targets_tidy, str_detect(name, "glc-4"))

glc_4_exons <- filter(gtf, gene_id == unique(glc_4$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

glc_4_bed <- filter(bed, isoform %in% glc_4$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

glc_4_features <- bind_rows(glc_4_exons, glc_4_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

glc_4_order <- reorder(glc_4_bed$object, glc_4_bed$object, function(x) - length(x))
glc_4_order <- as.vector(levels(glc_4_order))
glc_4_order <- c(unique(glc_4_exons$object), glc_4_order)

(glc_4_viewer <- ggplot(glc_4_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 2 Position (Mb)", color = "Structural Category", title = "*Bma-glc-4*") +
    scale_color_manual(values = c(pamplemousse[c(5)], "grey60")) +
    scale_y_discrete(limits = glc_4_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# LGICs -------------------------------------------------------------------

lgic_legend <- get_legend(lgc_27_viewer + theme(legend.position = "right"))

top <-  plot_grid(lgc_27_viewer + theme(legend.position = "none"),
                  acr_12_viewer + theme(legend.position = "none"), 
                  avr_14_viewer + theme(legend.position = "none"),
                  labels = "AUTO",
                  nrow = 1)

middle <- plot_grid(lgc_44_viewer + theme(legend.position = "none"),
                    acr_8_viewer + theme(legend.position = "none"),
                    unc_38_viewer + theme(legend.position = "none"),
                    labels = c('D', 'E', 'F'),
                    nrow = 1)

bottom <- plot_grid(lgc_45_viewer + theme(legend.position = "none"),
                    glc_4_viewer + theme(legend.position = "none"),
                    lgic_legend,
                    labels = c('G', 'H', ''),
                    nrow = 1)

lgic_panel <- plot_grid(top,
                        middle,
                        bottom,
                        nrow = 3,
                        rel_heights = c(1, 0.6, 0.4))

save_plot(here('..', 'plots', "S3_Fig.pdf"), lgic_panel, base_width = 7, base_height = 8)

# tax-2 ------------------------------------------------------------------

tax_2 <- filter(targets_tidy, str_detect(name, "tax-2"))

tax_2_exons <- filter(gtf, gene_id == unique(tax_2$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

tax_2_bed <- filter(bed, isoform %in% tax_2$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

tax_2_features <- bind_rows(tax_2_exons, tax_2_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

tax_2_order <- reorder(tax_2_bed$object, tax_2_bed$object, function(x) - length(x))
tax_2_order <- as.vector(levels(tax_2_order))
tax_2_order <- c(unique(tax_2_exons$object), tax_2_order)

(tax_2_viewer <- ggplot(tax_2_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 3 Position (Mb)", color = "Structural Category", title = "*Bma-tax-2*") +
    scale_color_manual(values = c(pamplemousse[c(8)], "grey60")) +
    scale_y_discrete(limits = tax_2_order) +
    scale_x_continuous(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# Bm5246 ------------------------------------------------------------------

Bm5246 <- filter(targets_tidy, str_detect(transcript_id, "Bm5246"))

Bm5246_exons <- filter(gtf, gene_id == unique(Bm5246$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

Bm5246_bed <- filter(bed, isoform %in% Bm5246$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

Bm5246_features <- bind_rows(Bm5246_exons, Bm5246_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

Bm5246_order <- reorder(Bm5246_bed$object, Bm5246_bed$object, function(x) - length(x))
Bm5246_order <- as.vector(levels(Bm5246_order))
Bm5246_order <- c(unique(Bm5246_exons$object), Bm5246_order)

(Bm5246_viewer <- ggplot(Bm5246_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 1 (+ strand)", color = "Structural Category", title = "Bma-trp-2") +
    scale_color_manual(values = c(pamplemousse[c(2, 5)], "grey60")) +
    scale_y_discrete(limits = Bm5246_order) +
    # scale_x_reverse() +
    theme_bw(12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

misc_channel <- plot_grid(tax_2_viewer, Bm5246_viewer, nrow = 2)

save_plot(here('..', 'plots', 'S4_Fig.pdf'), misc_channel, base_height = 7, base_width = 5)

# ser-1 ------------------------------------------------------------------

ser_1 <- filter(targets_tidy, str_detect(name, "ser-1"))

ser_1_exons <- filter(gtf, gene_id == unique(ser_1$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

ser_1_bed <- filter(bed, isoform %in% ser_1$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

ser_1_features <- bind_rows(ser_1_exons, ser_1_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

ser_1_order <- reorder(ser_1_bed$object, ser_1_bed$object, function(x) - length(x))
ser_1_order <- as.vector(levels(ser_1_order))
ser_1_order <- c(unique(ser_1_exons$object), ser_1_order)

(ser_1_viewer <- ggplot(ser_1_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr X Position (Mb)", color = "Structural Category", title = "*Bma-ser-1*") +
    scale_color_manual(values = c(pamplemousse[c(8)], "grey60")) +
    scale_y_discrete(limits = ser_1_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# Bm635 ------------------------------------------------------------------

Bm635 <- filter(targets_tidy, str_detect(transcript_id, "Bm635"))

Bm635_exons <- filter(gtf, gene_id == unique(Bm635$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

Bm635_bed <- filter(bed, isoform %in% Bm635$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

Bm635_features <- bind_rows(Bm635_exons, Bm635_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

Bm635_order <- reorder(Bm635_bed$object, Bm635_bed$object, function(x) - length(x))
Bm635_order <- as.vector(levels(Bm635_order))
Bm635_order <- c(unique(Bm635_exons$object), Bm635_order)

(Bm635_viewer <- ggplot(Bm635_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 4 Position (Mb)", color = "Structural Category", title = "*Bm635* (srab)") +
    scale_color_manual(values = c(pamplemousse[c(2, 5)], "grey60")) +
    scale_y_discrete(limits = Bm635_order) +
    scale_x_reverse(n.breaks = 4, labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# Bm13207 ------------------------------------------------------------------

Bm13207 <- filter(targets_tidy, str_detect(transcript_id, "Bm13207"))

Bm13207_exons <- filter(gtf, gene_id == unique(Bm13207$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

Bm13207_bed <- filter(bed, isoform %in% Bm13207$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

Bm13207_features <- bind_rows(Bm13207_exons, Bm13207_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

Bm13207_order <- reorder(Bm13207_bed$object, Bm13207_bed$object, function(x) - length(x))
Bm13207_order <- as.vector(levels(Bm13207_order))
Bm13207_order <- c(unique(Bm13207_exons$object), Bm13207_order)

(Bm13207_viewer <- ggplot(Bm13207_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 4 Position (Mb)", color = "Structural Category", title = "*Bm13207* (srab)") +
    scale_color_manual(values = c(pamplemousse[c(5)], "grey60")) +
    scale_y_discrete(limits = Bm13207_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# Bm2601 ------------------------------------------------------------------

Bm2601 <- filter(targets_tidy, str_detect(transcript_id, "Bm2601"))

Bm2601_exons <- filter(gtf, gene_id == unique(Bm2601$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

Bm2601_bed <- filter(bed, isoform %in% Bm2601$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

Bm2601_features <- bind_rows(Bm2601_exons, Bm2601_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

Bm2601_order <- reorder(Bm2601_bed$object, Bm2601_bed$object, function(x) - length(x))
Bm2601_order <- as.vector(levels(Bm2601_order))
Bm2601_order <- c(unique(Bm2601_exons$object), Bm2601_order)

(Bm2601_viewer <- ggplot(Bm2601_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 4 Position (Mb)", color = "Structural Category", title = "*Bm2601* (srab)") +
    scale_color_manual(values = c(pamplemousse[c(5, 8)], "grey60")) +
    scale_y_discrete(limits = Bm2601_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# Bm6043 ------------------------------------------------------------------

Bm6043 <- filter(targets_tidy, str_detect(transcript_id, "Bm6043"))

Bm6043_exons <- filter(gtf, gene_id == unique(Bm6043$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

Bm6043_bed <- filter(bed, isoform %in% Bm6043$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

Bm6043_features <- bind_rows(Bm6043_exons, Bm6043_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

Bm6043_order <- reorder(Bm6043_bed$object, Bm6043_bed$object, function(x) - length(x))
Bm6043_order <- as.vector(levels(Bm6043_order))
Bm6043_order <- c(unique(Bm6043_exons$object), Bm6043_order)

(Bm6043_viewer <- ggplot(Bm6043_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 4 Position (Mb)", color = "Structural Category", title = "*Bm6043* (srw)") +
    scale_color_manual(values = c(pamplemousse[c(5, 8)], "grey60")) +
    scale_y_discrete(limits = Bm6043_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# Bm17479 ------------------------------------------------------------------

Bm17479 <- filter(targets_tidy, str_detect(transcript_id, "Bm17479"))

Bm17479_exons <- filter(gtf, gene_id == unique(Bm17479$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

Bm17479_bed <- filter(bed, isoform %in% Bm17479$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

Bm17479_features <- bind_rows(Bm17479_exons, Bm17479_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

Bm17479_order <- reorder(Bm17479_bed$object, Bm17479_bed$object, function(x) - length(x))
Bm17479_order <- as.vector(levels(Bm17479_order))
Bm17479_order <- c(unique(Bm17479_exons$object), Bm17479_order)

(Bm17479_viewer <- ggplot(Bm17479_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 2 Position (Mb)", color = "Structural Category", title = "*Bm17479* (srw)") +
    scale_color_manual(values = c(pamplemousse[c(5)], "grey60")) +
    scale_y_discrete(limits = Bm17479_order) +
    scale_x_continuous(labels = function(x) (x / 10E6), n.breaks = 4) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# Bm271 ------------------------------------------------------------------

Bm271 <- filter(targets_tidy, str_detect(transcript_id, "Bm271"))

Bm271_exons <- filter(gtf, gene_id == unique(Bm271$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

Bm271_bed <- filter(bed, isoform %in% Bm271$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

Bm271_features <- bind_rows(Bm271_exons, Bm271_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

Bm271_order <- reorder(Bm271_bed$object, Bm271_bed$object, function(x) - length(x))
Bm271_order <- as.vector(levels(Bm271_order))
Bm271_order <- c(unique(Bm271_exons$object), Bm271_order)

(Bm271_viewer <- ggplot(Bm271_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr 2 Position (Mb)", color = "Structural Category", title = "*Bm271* (srxa)") +
    scale_color_manual(values = c(pamplemousse[c(2, 8)], "grey60")) +
    scale_y_discrete(limits = Bm271_order) +
    scale_x_continuous(labels = function(x) (x / 10E6), n.breaks = 4) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# Bm18045 ------------------------------------------------------------------

Bm18045 <- filter(targets_tidy, str_detect(transcript_id, "Bm18045"))

Bm18045_exons <- filter(gtf, gene_id == unique(Bm18045$gene_id)) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

Bm18045_bed <- filter(bed, isoform %in% Bm18045$isoform) %>%
  left_join(., classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

Bm18045_features <- bind_rows(Bm18045_exons, Bm18045_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

Bm18045_order <- reorder(Bm18045_bed$object, Bm18045_bed$object, function(x) - length(x))
Bm18045_order <- as.vector(levels(Bm18045_order))
Bm18045_order <- c(unique(Bm18045_exons$object), Bm18045_order)

(Bm18045_viewer <- ggplot(Bm18045_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr X Position (Mb)", color = "Structural Category", title = "*Bm18045* (srab)") +
    scale_color_manual(values = c(pamplemousse[c(5)], "grey60")) +
    scale_y_discrete(limits = Bm18045_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

# GPCRs -------------------------------------------------------------------

gpcr_legend <- get_legend(lgc_27_viewer + theme(legend.position = "right"))

top <-  plot_grid(Bm635_viewer + theme(legend.position = "none"),
                  Bm2601_viewer + theme(legend.position = "none"), 
                  Bm6043_viewer + theme(legend.position = "none"),
                  labels = "AUTO", align = 'h', axis = 'tb',
                  nrow = 1)

middle <- plot_grid(Bm13207_viewer + theme(legend.position = "none"),
                    Bm271_viewer + theme(legend.position = "none"),
                    ser_1_viewer + theme(legend.position = "none"),
                    labels = c('D', 'E', 'F'), align = 'h', axis = 'tb',
                    nrow = 1)

bottom <-  plot_grid(Bm17479_viewer + theme(legend.position = "none"),
                     Bm18045_viewer + theme(legend.position = "none"),
                     gpcr_legend,
                     labels = c('G', 'H', ''),
                     nrow = 1)

gpcr_panel <- plot_grid(top,
                        middle,
                        bottom,
                        nrow = 3,
                        rel_heights = c(1.2, 1, 1))

save_plot(here('..', 'plots', "S5_Fig.pdf"), gpcr_panel, base_width = 7, base_height = 7)



# D. immitis --------------------------------------------------------------



# SQANTI classifications
di_classifications <- read_delim(here('..', '..', "1_SQANTI", 'data', "di_all_transcripts_trim_classification.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "brugia_malayi") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  mutate(structural_category = str_to_title(str_replace_all(structural_category, "_", " "))) %>%
  filter(structural_category != "Genic Intron", !str_detect(isoform, '_dup'))

# reference GTF
di_gtf <- readRDS(here('..', '..', "1_SQANTI", 'data', "di_gtf.rds")) %>%
  mutate(gene_id = str_remove(attribute, ";.*"),
         gene_id = str_remove_all(gene_id, "\""),
         gene_id = str_remove(gene_id, "gene_id "),
         transcript_id = str_remove(attribute, ".*transcript_id "),
         transcript_id = str_remove(transcript_id, ";.*"),
         transcript_id = str_remove_all(transcript_id, "\"")) %>%
  mutate(transcript_id = case_when(
    feature == "gene" ~ as.character(NA),
    TRUE ~ transcript_id
  ))

# BED file of junctions for each read
di_bed <- readRDS(here('..', 'data', "dirofilari_immitis_gene_intersects.rds")) %>%
  select(!contains("null")) %>%
  mutate(start = start + 1) # change from base-0 to base-1

di_tax4 <- filter(di_classifications, str_detect(associated_gene, 'nDi.2.2.2.g00427'))

di_tax4_exons <- filter(di_gtf, gene_id %in% c('nDi.2.2.2.g00427', 'nDi.2.2.2.g00428')) %>%
  filter(feature == "exon") %>%
  select(scaffold, object =  transcript_id, overlap_start = start, overlap_end = end) %>%
  mutate(structural_category = "Reference Transcripts")

di_tax4_bed <- filter(di_bed, isoform %in% di_tax4$isoform) %>%
  left_join(., di_classifications) %>%
  select(scaffold:isoform, structural_category) %>%
  rename(object = isoform) %>%
  separate(object, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  arrange(object, desc(start)) %>%
  filter(!is.na(structural_category))

di_tax4_features <- bind_rows(di_tax4_exons, di_tax4_bed) %>%
  mutate(sex = case_when(
    sex == "AM" ~ "Male",
    sex == "AF" ~ "Female"
  )) 

di_tax4_order <- reorder(di_tax4_bed$object, di_tax4_bed$object, function(x) - length(x))
di_tax4_order <- as.vector(levels(di_tax4_order))
di_tax4_order <- c(unique(di_tax4_exons$object), di_tax4_order)

(di_tax4_viewer <- ggplot(di_tax4_features) +
    geom_segment(aes(x = overlap_start, xend = overlap_end, y = object, yend = object, color = structural_category), size = 6, alpha = 0.85) +
    labs(x = "Chr X Position (Mb)", color = "Structural Category", title = "*Dim-tax-4*") +
    scale_color_manual(values = c(pamplemousse[c(3)], "grey60")) +
    scale_y_discrete(limits = di_tax4_order) +
    scale_x_reverse(labels = function(x) (x / 10E6)) +
    theme_bw(12) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    NULL
)

save_plot(here('..', 'plots', 'S1_Fig.pdf'), di_tax4_viewer, base_width = 4, base_height = 2)
