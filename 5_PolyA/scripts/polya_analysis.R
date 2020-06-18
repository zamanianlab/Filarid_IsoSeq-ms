# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(ggforce)
library(ggrepel)
library(ggbeeswarm)
library(ggridges)
library(patchwork)
library(cowplot)
library(ggtext)
library(LaCroixColoR)

# misc
library(conflicted)
library(seqinr)

# project managment
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")
conflict_prefer("here", "here")
conflict_prefer("summarize", "dplyr")

# poly(A) Length ----------------------------------------------------------

# df of poly(A) lengths
bm_length_data <- read_csv(here("..", 'data', "bm_all_transcripts_trim_polya.csv"), col_names = c("isoform", "polya_length"))
di_length_data <- read_csv(here("..", 'data', "di_all_transcripts_trim_polya.csv"), col_names = c("isoform", "polya_length"))
length_data <- bind_rows(bm_length_data %>% mutate(species = "*B. malayi*"), 
                         di_length_data %>% mutate(species = "*D. immitis*"))

(polya_plot <- ggplot(filter(length_data, species == "*B. malayi")) +
    geom_histogram(aes(x = polya_length, y = ..density..), center = median(length_data$polya_length), mutate(length_data, z = FALSE), bins = 250) +
    geom_histogram(aes(x = polya_length, y = ..density..), center = median(length_data$polya_length), mutate(length_data, z = TRUE), bins = 500) +
    geom_vline(aes(xintercept = median(length_data$polya_length)), mutate(length_data, z = TRUE), color = "red", linetype = "dashed") + # this add several seconds to the plotting
    scale_y_continuous(expand = c(0, 0)) +
    facet_zoom(xlim = c(15, 50), ylim = c(0, 0.15), horizontal = FALSE, zoom.data = z) +
    labs(x = "poly(A) Tail Length", y = "Density") +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      zoom.y = element_blank(),
      strip.background = element_rect(fill = "lightgrey", linetype = "dashed"),
      validate = FALSE) +
    NULL
)

save_plot(here('..', 'plots', 'Fig5A.pdf'), polya_plot, base_width = 3)
saveRDS(polya_plot, here('..', 'plots', 'Fig5A.rds'))

# kmers -------------------------------------------------------------------

# list of Bm scaffolds and lengths
bm_scaffolds <- read_tsv(here("..", 'data', "brugia_malayi.genome.bed"), col_names = c("scaffold", "length"))
di_scaffolds <- read_tsv(here("..", 'data', "dirofilaria_immitis.genome.bed"), col_names = c("scaffold", "length"))
scaffolds <- bind_rows(bm_scaffolds %>% mutate(species = "*B. malayi*"), 
                       di_scaffolds %>% mutate(species = "*D. immitis*"))

# BED file of all intergenic and intronic regions
bm_intergenic <- read_tsv(here("..", 'data', "brugia_malayi.intergenic.bed"), col_names = c("scaffold", "start", "stop")) 
di_intergenic <- read_tsv(here("..", 'data', "dirofilaria_immitis.intergenic.bed"), col_names = c("scaffold", "start", "stop")) 
intergenic <- bind_rows(bm_intergenic %>% mutate(species = "*B. malayi*"), 
                        di_intergenic %>% mutate(species = "*D. immitis*")) %>%
  mutate(length = stop - start)

# number of possible 6-mers in intergenic
n_6mer <- group_by(intergenic, species) %>%
  summarize(n_6mer = sum(length) - 6)

# tally of k-mers in the intergenic and intronic regions
bm_intergenic_kmer <- read_csv(here("..", 'data', "brugia_malayi.intergenic_6mer.csv")) 
di_intergenic_kmer <- read_csv(here("..", 'data', "dirofilaria_immitis.intergenic_6mer.csv")) 
intergenic_kmer <- bind_rows(bm_intergenic_kmer %>% mutate(species = "*B. malayi*"), 
                             di_intergenic_kmer %>% mutate(species = "*D. immitis*")) %>%
  clean_names() %>%
  left_join(., n_6mer) %>%
  mutate(freq = count / n_6mer) %>% # freq is the frequency of the kmer in the entire intergenic region
  mutate(kmer = str_replace_all(kmer, "T", "U")) %>%
  select(-n_6mer)

# df of most abundant k-mers
utr_kmer <- read_csv(here('..', 'data', 'bm_all_transcripts_trim_6mer.csv')) %>%
  mutate(species = '*B. malayi*') %>%
  bind_rows(., read_csv(here('..', 'data', 'bm_all_transcripts_trim_6mer.csv')) %>% mutate(species = '*D. immitis*')) %>%
  rename(count = Count) %>%
  arrange(desc(count)) %>%
  mutate(kmer = str_replace_all(kmer, "T", "U")) %>%
  left_join(., select(intergenic_kmer, species, kmer, freq)) %>%
  rename(exp_freq = freq) %>%
  mutate(exp_freq = exp_freq / sum(exp_freq)) %>%
  mutate(obs_freq = count / (length(length_data$isoform) * 50 - 6)) %>%
  mutate(exp_count = exp_freq * (length(length_data$isoform) * 50 - 6)) %>%
  mutate(ratio  = obs_freq / exp_freq) %>%
  mutate(z = (ratio - mean(ratio)) / sd(ratio)) %>%
  mutate(weighted_z = z * obs_freq) %>%
  arrange(-weighted_z)

# write out the reordered k-mer list by species

######### rerun polyAudit with -PAS and file generated above  

# df of the putative PAS for each transcript
pas_index <- read_csv(here("..", 'data', "bm_all_transcripts_trim_PAS.csv"), 
                      col_names = c("null", "isoform", "3'_sequence", "kmer", "index"),
                      skip = 1) %>%
  mutate(species = "*B. malayi*") %>%
  select(-null) %>%
  mutate(kmer = str_replace_all(kmer, "T", "U")) %>%
  select(species, isoform, kmer) %>%
  group_by(species, kmer) %>%
  summarise(genes_with_kmer = n()) %>%
  mutate(kmer = fct_reorder(kmer, -genes_with_kmer)) %>%
  arrange(-genes_with_kmer)

# take the top 10 most abundant and add together all others
pas_plot_data <- group_by(pas_index, species) %>%
  top_n(10, genes_with_kmer) 
  
# color palette
palette <- viridisLite::viridis(length(pas_plot_data$kmer))

# plot the number of times a given PAS appears in the CCSs
(isoform_plot <- ggplot(filter(pas_plot_data), 
                        aes(x = kmer, y = genes_with_kmer / sum(filter(pas_index, species == "*B. malayi*")$genes_with_kmer))) +
    geom_col(aes(fill = kmer)) +
    labs(x = "", y = "Fraction CCSs with PAS") + 
    scale_fill_manual(values = palette) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      legend.position = "none") +
    NULL
)

save_plot(here('..', 'plots', 'Fig5B.pdf'), isoform_plot, base_width = 3)
saveRDS(isoform_plot, here('..', 'plots', 'Fig5B.rds'))

# plot location of the PAS with respect to the 3' end of the PAS
pas_index_summary <- read_csv(here("..", 'data', "bm_all_transcripts_trim_PAS.csv"), 
                              col_names = c("null", "isoform", "3'_sequence", "kmer", "index"),
                              skip = 1) %>%
  mutate(species = "*B. malayi*") %>%
  select(-null) %>%
  mutate(kmer = str_replace_all(kmer, "T", "U")) %>%
  select(species, kmer, index) %>%
  filter(kmer %in% filter(pas_plot_data, species == '*B. malayi*')$kmer) %>%
  mutate(kmer = factor(kmer, levels = levels(pas_plot_data$kmer)))

(index_plot <- ggplot(filter(pas_index_summary, !is.na(kmer))) +
    geom_density_ridges(aes(x = index, y = kmer, fill = kmer), scale = 3, rel_min_height = 0.01) +
    scale_x_reverse(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(droplevels(pas_index_summary$kmer))),
                     expand = c(0, 0)) +
    scale_fill_manual(values = palette) +
    labs(y = "", x = "PAS distance from poly(A)") +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      legend.position = "none") +
    NULL
  )

save_plot(here('..', 'plots', 'Fig5C.pdf'), index_plot, base_width = 3)
saveRDS(index_plot, here('..', 'plots', 'Fig5C.rds'))

pas_nest <- read_csv(here("..", 'data', "bm_all_transcripts_trim_PAS.csv"), 
                     col_names = c("null", "isoform", "3'_sequence", "kmer", "index"),
                     skip = 1) %>%
  mutate(species = "*B. malayi*") %>%
  select(-null) %>%
  mutate(kmer = str_replace_all(kmer, "T", "U")) %>%
  filter(kmer %in% pas_plot_data$kmer) %>%
  group_nest(kmer)

######### rerun polyAudit with -PSSM for each PAS

pssm_files <- tibble(file = list.files(path = here('..', 'data', 'pssm', 'brugia_malayi'), pattern = ".*_pssm.txt$", recursive = TRUE)) %>%
  separate(file, c("kmer", "txt"), sep = "_", remove = FALSE) %>%
  select(file, kmer)

import_pssm <- function(file) {
  
  data <- read_table2(here('..', 'data', 'pssm', 'brugia_malayi', file[1]), skip = 1, col_names = c("X", "A", "C", "G", "T")) %>%
    select(-1) %>%
    mutate(kmer = file[2])
  
}

pssm <- apply(pssm_files, 1, import_pssm) %>%
  bind_rows() %>%
  group_by(kmer) %>%
  mutate(index = seq(1, 65, 1)) %>%
  select(index, kmer, everything()) %>%
  rename(U = T) %>%
  pivot_longer(A:U, names_to = "nucleotide", values_to = "count") %>%
  ungroup() %>%
  mutate(kmer = factor(kmer, levels = levels(droplevels(pas_index_summary$kmer))))

pssm_summary <- group_by(pssm, kmer, index) %>%
  summarise(sum = sum(count))

pssm <- left_join(pssm, pssm_summary) %>%
  mutate(frequency = count / sum)

(freq <- ggplot(pssm) +
    geom_line(aes(x = index - 50, y = frequency, color = nucleotide), size = 1) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    scale_color_manual(values = lacroix_palette("PassionFruit", n = 4, type = "discrete")) +
    scale_x_continuous(limits = c(-50, 15), breaks = seq(-45, 15, 15), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    facet_grid(rows = vars(kmer)) +
    labs(x = "Index", y = "Frequency") +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      strip.background = element_rect(fill = "lightgrey", linetype = "dashed"),
      strip.text = element_text(size = 7),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 8)) +
    NULL
  )

save_plot(here('..', 'plots', 'Fig5D.pdf'), freq, base_width = 3)
saveRDS(freq, here('..', 'plots', 'Fig5D.rds'))

# Final plotting ----------------------------------------------------------

isoform_plot <- isoform_plot + theme(axis.text.x = element_text(size = 8))

middle <- isoform_plot / index_plot
left <- polya_plot + middle
final <- left + freq + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold'))

save_plot(here('..', 'plots', "Fig5.pdf"), final, base_height = 7.5, base_width = 7.5)
save_plot(here('..', 'plots', "Fig5.png"), final, base_height = 7.5, base_width = 7.5)





# D. immitis --------------------------------------------------------------

(polya_plot <- ggplot(filter(length_data, species == "*B. malayi")) +
   geom_histogram(aes(x = polya_length, y = ..density..), center = median(length_data$polya_length), mutate(length_data, z = FALSE), bins = 250) +
   geom_histogram(aes(x = polya_length, y = ..density..), center = median(length_data$polya_length), mutate(length_data, z = TRUE), bins = 500) +
   geom_vline(aes(xintercept = median(length_data$polya_length)), mutate(length_data, z = TRUE), color = "red", linetype = "dashed") + # this add several seconds to the plotting
   facet_zoom(xlim = c(15, 50), ylim = c(0, 0.15), horizontal = FALSE, zoom.data = z) +
   labs(x = "poly(A) Tail Length", y = "Density") +
   theme_minimal() +
   theme(zoom.y = element_blank(),
         strip.background = element_rect(fill = "lightgrey", linetype = "dashed"),
         validate = FALSE) +
   NULL
)

save_plot(here('..', 'plots', 'Fig5A.pdf'), polya_plot, base_width = 3)
saveRDS(polya_plot, here('..', 'plots', 'Fig5A.rds'))

(isoform_plot <- ggplot(filter(pas_plot_data, species == "*B. malayi*"), 
                        aes(x = kmer, y = genes_with_kmer / sum(filter(pas_index, species == "*B. malayi*")$genes_with_kmer))) +
    geom_col(aes(fill = kmer)) +
    labs(x = "", y = "Fraction CCSs with PAS") + 
    scale_fill_manual(values = palette) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          legend.position = "none") +
    NULL
)