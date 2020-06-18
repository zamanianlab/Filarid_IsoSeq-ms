# data wrangling/plotting
library(tidyverse)
library(janitor)

# misc
library(conflicted)

# project managment
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")
conflict_prefer("here", "here")
conflict_prefer("summarize", "dplyr")

# Import ------------------------------------------------------------------

# structural classifications
bm_classifications <- read_delim(here("..", "..", "1_SQANTI", "data", "bm_all_transcripts_trim_classification.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "brugia_malayi") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  mutate(structural_category = str_to_title(str_replace_all(structural_category, "_", " "))) %>%
  filter(rts_stage == FALSE, structural_category != "Genic Intron", !str_detect(isoform, '_dup'))

di_classifications <- read_delim(here("..", "..", "1_SQANTI", "data", "di_all_transcripts_trim_classification.txt"), delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  mutate(species = "dirofilaria_immitis") %>%
  separate(isoform, into = c("sex", "quality"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(species, sex, quality, isoform, everything()) %>%
  mutate(structural_category = str_to_title(str_replace_all(structural_category, "_", " "))) %>%
  filter(rts_stage == FALSE, structural_category != "Genic Intron", !str_detect(isoform, '_dup'))

classifications <- bind_rows(bm_classifications, di_classifications) %>%
  filter(structural_category != "Genic Intron", !str_detect(isoform, '_dup'))

# reference GTFs
gtf <- bind_rows(as_tibble(readRDS(here("..", "..", "1_SQANTI", 'data', 'bm_gtf.rds')) %>% mutate(species = 'brugia_malayi')),
                 as_tibble(readRDS(here("..", "..", "1_SQANTI", 'data', 'di_gtf.rds')) %>% mutate(species = 'dirofilaria_immitis'))
)

# ncRNA -------------------------------------------------------------------

ncRNA <- filter(gtf, str_detect(attribute, "ncRNA")) %>%
  filter(feature == 'transcript') %>%
  mutate(associated_gene = str_remove(attribute, ";.*"),
         associated_gene = str_remove_all(associated_gene, "\""),
         associated_gene = str_remove(associated_gene, "gene_id "),
         associated_transcript = str_remove(attribute, ".*transcript_id "),
         associated_transcript = str_remove(associated_transcript, ";.*"),
         associated_transcript = str_remove_all(associated_transcript, "\""),
         associated_transcript = str_remove(associated_transcript, "gene_id "))

ncRNA_sqanti <- filter(classifications, associated_gene %in% ncRNA$associated_gene)

# Novel protein coding ----------------------------------------------------

novel_protein_coding <- filter(classifications, str_detect(associated_gene, 'novel'), exons > 1, perc_a_downstream_tts > 0)

write_csv(select(novel_protein_coding, isoform), here('..', 'data', "novel_genes_coding.txt"), col_names = FALSE)

novel_protein_coding <- filter(classifications, str_detect(associated_gene, 'novel'), exons > 1) %>%
  group_by(species) %>%
  tally()

# Find number of loci mapped ----------------------------------------------

bm_novel_ccs <- read_tsv(here("..", "data", "bm_novel_genes_coding_corrected_blastn_filter_3.out"),
                 col_names = c("isoform", 'chrom', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')) %>%
  group_by(isoform, chrom) %>%
  summarise(min_sstart = min(sstart), max_send = max(send)) %>%
  arrange(chrom, min_sstart) %>%
  ungroup()

bm_novel_loci <- bm_novel_ccs %>%
  select(-isoform) %>%
  distinct() %>%
  filter(!(chrom == lead(chrom) & max_send >= lead(min_sstart)))

di_novel_ccs <- read_tsv(here("..", "data", "di_novel_genes_coding_corrected_blastn_filter_3.out"),
                          col_names = c("isoform", 'chrom', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')) %>%
  group_by(isoform, chrom) %>%
  summarise(min_sstart = min(sstart), max_send = max(send)) %>%
  arrange(chrom, min_sstart) %>%
  ungroup() 

di_novel_loci <- di_novel_ccs %>%
  select(-isoform) %>%
  distinct() %>%
  filter(!(chrom == lead(chrom) & max_send >= lead(min_sstart)))

novel_protein_coding_loci <- bind_rows(bm_novel_loci, di_novel_loci) %>%
  rename(scaffold = chrom, start = min_sstart, end = max_send)

write_csv(novel_protein_coding_loci, here('..', 'supp_data', 'S2_Table.csv'))
  
