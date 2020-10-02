# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(ggridges)
library(patchwork)
library(cowplot)
library(ggtext)
library(LaCroixColoR)

# misc
library(conflicted)

# project management
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("here", "here")
conflict_prefer("summarize", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("lag", "dplyr")


# Import data -------------------------------------------------------------

# read in GTF and parse to keep relevant metadata
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

# get the 5' and 3' termini for all CDSs
cds_start_end <- gtf %>%
  group_by(scaffold, strand, gene_id, transcript_id) %>%
  summarise(cds_start = min(start), cds_end = max(end)) 

# read in bed file that describes intersections between isoforms and genes
bed <- readRDS(here('..', 'data', "brugia_malayi_gene_intersects.rds")) %>%
  select(scaffold, gene_start, gene_end, gene_id, strand, isoform) %>%
  distinct()

# Operon identification ---------------------------------------------------

# an operon is any isoform that maps to multiple genes
operons <- group_by(bed, isoform) %>%
  count(name = "count") %>%
  filter(count > 1)

# calculate the distance between the end of one cistron and the start of the next
bed_operons <- filter(bed, isoform %in% operons$isoform) %>%
  filter(scaffold != "Bm_029") %>% # remove a scaffold that was causing a lot of issues due to repetitive regions
  select(isoform, gene_id) %>%
  left_join(., cds_start_end) %>%
  arrange(isoform, scaffold, strand, cds_start) %>%
  # calculate the distance between overlapped CDSs
  group_by(isoform) %>%
  mutate(distance = case_when(
    scaffold == lead(scaffold) & strand == lead(strand) ~ lead(cds_start) - cds_end, # when they're on the same scaffold/ & strand, calculate the distance
    scaffold == lead(scaffold) & strand != lead(strand) ~ NA_real_, # when they're not on the same strand, mark for filtration
    scaffold != lead(scaffold) ~ NA_real_  # when they're not on the same scaffold, mark for filtration
  )) %>%
  group_by(isoform, scaffold, strand) %>%
  # keep any operon that has a distance calculation
  filter(any(!is.na(distance))) %>%
  # number cistrons
  mutate(cistron = row_number()) %>%
  # reverse cistron order for - strand
  mutate(cistron = case_when(
    strand == '-' ~ max(cistron) - cistron + 1,
    TRUE ~ as.double(cistron)
  )) %>%
  arrange(isoform, scaffold, strand, cistron) %>%
  # rearrange the distance calculations for - strand operons so that cistron 1 contains the distance
  mutate(distance = case_when(
    strand == '-' & is.na(distance) ~ lead(distance),
    strand == '-' & cistron == max(cistron) ~ NA_real_,
    strand == '-' & cistron != max(cistron) ~ lead(distance),
    TRUE ~ distance
  )) %>%
  select(isoform, gene_id, distance, everything()) %>%
  group_by(isoform) %>%
  arrange(distance)

# give each operon a name
operon_names <- filter(bed_operons, cistron == 1) %>%
  group_by(scaffold, cds_start, strand) %>%
  mutate(operon_id = str_glue('BMISO', cur_group_id())) %>% # name by group, BMISO = B. malayi Iso-Seq Operon
  ungroup() %>%
  select(isoform, scaffold, strand, operon_id)

bed_operons <- left_join(bed_operons, operon_names) %>%
  ungroup() %>%
  group_by(operon_id) %>%
  arrange(distance, operon_id) %>%
  # calculate how many reads map to each operon
  mutate(support = sum(cistron == 1))

# remove isoform info
collapsed_operons <- select(bed_operons, -isoform) %>%
  distinct()

# filter by distance
filtered_operons <- group_by(collapsed_operons, operon_id) %>%
  # the vast majority of assembly-annotated operons have intercistronic distances of < 5000
  filter(max(distance, na.rm = TRUE) < 5000) %>%
  group_by(operon_id) %>%
  arrange(operon_id, cistron) %>%
  select(operon_id, distance, cistron, gene_id, transcript_id, scaffold, strand, cds_start, cds_end, support)

# Concatenate cistron sequences to blast against Ce -----------------------

# Identify fragmented models ----------------------------------------------

# write out all the transcript_ids that are in putative operons
write_csv(select(ungroup(filtered_operons), transcript_id), here('..', 'data', "polycistron_transcripts.txt"), col_names = FALSE)

fasta <- seqinr::read.fasta(here('..', 'data', 'polycistron_transcripts.fasta'), as.string = TRUE)

# concatenate transcript sequences to make near-complete operonic sequences
fasta_tab <- tibble(id = sapply(fasta, attr, 'Annot'), seq = unlist(fasta)) %>%
  mutate(transcript_id = str_remove_all(id, '>')) %>%
  mutate(transcript_id = str_remove(transcript_id, ' wormpep.*')) %>%
  select(transcript_id, seq) %>%
  left_join(., select(filtered_operons, transcript_id, operon_id)) %>%
  distinct() %>%
  group_by(operon_id) %>%
  summarise(seq_cat = str_c(seq, collapse = ""))

# write out the concatenated sequences
seqinr::write.fasta(as.list(fasta_tab$seq_cat), as.list(fasta_tab$operon_id), here('..', 'data',"operons.fasta"))

# read in the blastp output of operons.fasta vs the Ce proteome
operon_blast <- read_tsv(here('..', 'data', 'operons_blastp_1.out'), col_names = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')) %>%
  mutate(sseqid = str_remove(sseqid, '[a-z]$')) %>%
  group_by(qseqid, sseqid) %>%
  # keep the most similar region of each hit
  slice_min(n = 1, order_by = 'evalue') %>%
  # keep only significant hits
  filter(evalue < 1e-1)

# read in the blastp output of polycistron_transcripts.fasta vs the Ce proteome
ce_blast <- read_tsv(here('..', 'data', "polycistron_transcripts_blastp_1.out"), col_names = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')) %>%
  mutate(sseqid = str_remove(sseqid, '[a-z]$')) %>%
  group_by(qseqid) %>%
  arrange(qseqid, evalue) %>%
  # keep the top hit
  slice_head(n = 1) %>%
  select(transcript_id = qseqid, blast_hit = sseqid)

filtered_operons_blast <- left_join(filtered_operons, ce_blast)

# Manual curation ---------------------------------------------------------

# write to a Google Sheet for manual curation (please do not rerun)
# sheet <- googlesheets4::gs4_create('operon_curation')
# 
# filtered_operons_blast %>% googlesheets4::sheet_append(sheet)
# operon_blast %>% googlesheets4::sheet_append(sheet, sheet = 'operon_blast')

# read in the Google Sheet with manual curation and the official operon_id's
curations <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1wwIxXf0ELtPqPCZo5-lra7rEIT2fZez6rBk9I-qQg3U/edit#gid=316400157") %>%
  janitor::clean_names()

curated_operons <- curations %>% 
  filter(str_detect(operon_id, "BM*") == TRUE) 

operon_curation_summary <- curated_operons %>% 
  group_by(operon_id) %>% 
  filter(curation_type != 'Killed') %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  group_by(curation_type) %>% 
  distinct(operon_id) %>% 
  summarise(Count = n())

merge_curation_summary <- curations %>% 
  filter(str_detect(operon_id, "MERGE") == TRUE) %>% 
  group_by(operon_id) %>% 
  filter(curation_type != "Killed", row_number() == 1) %>% 
  ungroup() %>% 
  group_by(curation_type) %>% 
  distinct(operon_id) %>% 
  summarise(Count = n())

# transfer over official id's
bed_operons <- ungroup(bed_operons) %>% 
  select(-operon_id) %>% 
  left_join(., select(curated_operons, operon_id, gene_id, transcript_id)) %>% 
  group_by(operon_id, gene_id, distance, scaffold, strand, transcript_id, cds_start, cds_end, support) %>% 
  summarise(cistron = min(cistron))

# Get operons from Bmal-4.0 -----------------------------------------------

# operons grep'ed from the GFF and converted to BED
assembly_operons <- read_delim(here('..', 'data', "assembly_operons.bed"),
                               delim = "\t",
                               comment = "#",
                               col_names = c("scaffold", "gene_start", "gene_end", "NULL", "frame", "strand", "source", "type", "NULL", "attribute")) %>%
  select(scaffold, gene_start, gene_end, strand, attribute) %>%
  separate(attribute, sep = ";", into = c("operon_id", "gene_ids")) %>%
  mutate(operon_id = str_replace_all(operon_id, "Name=", "")) %>%
  mutate(gene_ids = str_replace_all(gene_ids, "genes=", "")) %>%
  mutate(gene_ids = str_split(gene_ids, ",")) %>%
  unnest(gene_ids) %>%
  rename(gene_id = gene_ids) %>%
  select(-gene_start, -gene_end) %>%
  # merge with the GTF CDSs because start/end in the BED is for the operons, not the cistrons
  left_join(., cds_start_end) %>%
  arrange(operon_id, scaffold, strand, cds_start) %>%
  group_by(operon_id) %>%
  mutate(distance = case_when(
    scaffold == lead(scaffold) & strand == lead(strand) ~ lead(cds_start) - cds_end, # when they're on the same scaffold/ & strand, calculate the distance
    scaffold == lead(scaffold) & strand != lead(strand) ~ NA_real_, # when they're not on the same strand, mark for filtration
    scaffold != lead(scaffold) ~ NA_real_  # when they're not on the same scaffold, mark for filtration
  )) %>%
  filter(any(!is.na(distance))) %>%
  mutate(cistron = row_number()) %>%
  # reverse cistron order for - strand
  mutate(cistron = case_when(
    strand == '-' ~ max(cistron) - cistron + 1,
    TRUE ~ as.double(cistron)
  )) %>%
  arrange(scaffold, strand, cistron) %>%
  # rearrange the distance calculations for - strand operons to max the order for + strand
  mutate(distance = case_when(
    strand == '-' & is.na(distance) ~ lead(distance),
    strand == '-' & cistron == max(cistron) ~ NA_real_,
    strand == '-' & cistron != max(cistron) ~ lead(distance),
    TRUE ~ distance
  )) %>%
  select(gene_id, distance, everything()) %>%
  arrange(distance)

# create final list of B. malayi operons
final_operons <- curated_operons %>% 
  group_by(operon_id) %>% 
  mutate(operon_type = case_when(
    curation_type %in% 'Novel operon' ~ 'Novel',
    curation_type %in% 'No change' ~ 'Shared',
    curation_type %in% 'Merge & Keep' ~ 'Shared',
    curation_type %in% 'Gene added' ~ 'Shared',
    curation_type %in% 'Model correction' ~ 'Shared'
  )) %>% 
  filter(!is.na(operon_type)) %>% 
  select(operon_id, operon_type) %>%
  distinct()

final_operons <- bind_rows(final_operons, 
            distinct(select(filter(ungroup(assembly_operons), !operon_id %in% final_operons$operon_id),
                   operon_id) %>% mutate(operon_type = 'Assembly'))) %>% 
  mutate(species = 'B. malayi')

final_operonic_genes <- left_join(final_operons, bind_rows(assembly_operons, select(filter(bed_operons, !operon_id %in% assembly_operons$operon_id), -support))) %>%
  mutate(operon_id = as.character(operon_id)) %>% 
  distinct()

# saveRDS(final_operonic_genes, here('..', 'data', "final_operonic_genes.rds"))

# get all the genes that aren't in operons
non_operonic_genes <- filter(cds_start_end, !gene_id %in% final_operonic_genes$gene_id) %>%
  arrange(scaffold, strand, cds_start) %>%
  group_by(scaffold, strand) %>%
  mutate(distance = case_when(
    scaffold == lead(scaffold) & strand == lead(strand) ~ lead(cds_start) - cds_end, # when they're on the same scaffold/ & strand, calculate the distance
    scaffold == lead(scaffold) & strand != lead(strand) ~ NA_real_, # when they're not on the same strand, mark for filtration
    scaffold != lead(scaffold) ~ NA_real_  # when they're not on the same scaffold, mark for filtration
  )) %>%
  mutate(operon_id = factor(NA), operon_type = "Monocistron", cistron = NA, species = "B. malayi")

# Get operons from Ce WS276 -----------------------------------------------

ce_operons <- read_delim(here('..', 'data', "ce_operons.bed"),
                         delim = "\t",
                         comment = "#",
                         col_names = c("scaffold", "gene_start", "gene_end", "NULL", "frame", "strand", "subtype", "type", "NULL", "attribute")) %>%
  filter(subtype != "deprecated_operon") %>%
  separate(attribute, sep = ";", into = c("operon_id", "gene_ids")) %>%
  mutate(operon_id = str_replace_all(operon_id, "Name=", "")) %>%
  mutate(gene_ids = str_replace_all(gene_ids, "genes=", "")) %>%
  mutate(gene_ids = str_split(gene_ids, ",")) %>%
  unnest(gene_ids) %>%
  rename(gene_id = gene_ids) %>%
  select(scaffold, strand, operon_id, gene_id)

ce_gtf <- readRDS(here('..', 'data', 'ce_gtf.rds')) %>%
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
ce_cds_start_end <- ce_gtf %>%
  group_by(scaffold, strand, gene_id, transcript_id) %>%
  summarise(cds_start = min(start), cds_end = max(end)) 

ce_operons_cds <- left_join(ce_operons, ce_cds_start_end) %>%
  distinct() %>%
  arrange(operon_id, scaffold, strand, cds_start) %>%
  group_by(operon_id) %>%
  mutate(distance = case_when(
    scaffold == lead(scaffold) & strand == lead(strand) ~ lead(cds_start) - cds_end, # when they're on the same scaffold/ & strand, calculate the distance
    scaffold == lead(scaffold) & strand != lead(strand) ~ NA_real_, # when they're not on the same strand, mark for filtration
    scaffold != lead(scaffold) ~ NA_real_  # when they're not on the same scaffold, mark for filtration
  )) %>%
  group_by(operon_id, scaffold, strand) %>%
  filter(any(!is.na(distance))) %>%
  mutate(cistron = row_number()) %>%
  # reverse cistron order for - strand
  mutate(cistron = case_when(
    strand == '-' ~ max(cistron) - cistron + 1,
    TRUE ~ as.double(cistron)
  )) %>%
  arrange(scaffold, strand, cistron) %>%
  # rearrange the distance calculations for - strand operons to max the order for + strand
  mutate(distance = case_when(
    strand == '-' & is.na(distance) ~ lead(distance),
    strand == '-' & cistron == max(cistron) ~ as.numeric(NA),
    strand == '-' & cistron != max(cistron) ~ lead(distance),
    TRUE ~ distance
  )) %>%
  select(gene_id, distance, everything()) %>%
  mutate(operon_type = "Assembly", species = "C. elegans") 

ce_non_operons <- filter(ce_cds_start_end, !gene_id %in% ce_operons_cds$gene_id) %>%
  arrange(scaffold, strand, cds_start) %>%
  group_by(scaffold, strand) %>%
  mutate(distance = case_when(
    scaffold == lead(scaffold) & strand == lead(strand) ~ lead(cds_start) - cds_end, # when they're on the same scaffold/ & strand, calculate the distance
    scaffold == lead(scaffold) & strand != lead(strand) ~ NA_real_, # when they're not on the same strand, mark for filtration
    scaffold != lead(scaffold) ~ NA_real_  # when they're not on the same scaffold, mark for filtration
  )) %>%
  mutate(operon_id = factor(NA), operon_type = "Monocistron", cistron = NA, species = "C. elegans")

# Orthologs ---------------------------------------------------------------

# use BioMart from WormBase ParaSite to get Ce orthologs of Bm genes in operons
library(biomaRt)

mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

orthologs <- getBM(mart = mart, 
                   filters = c("wbps_gene_id"),
                   value = list(final_operonic_genes$gene_id),
                   attributes = c("production_name_1010", "wbps_gene_id", "caelegprjna13758_gene")) %>%
  clean_names() %>%
  dplyr::select(gene_id = wbps_gene_id, ce_ortholog = caelegprjna13758_gene) %>%
  filter(ce_ortholog != "")

ortholog_operons <- left_join(final_operonic_genes, orthologs) %>%
  mutate(ortholog_type = case_when(
    is.na(ce_ortholog) ~ "No ortholog",
    !is.na(ce_ortholog) & ce_ortholog %in% ce_operons$gene_id ~ "Ortholog is in an operon",
    !is.na(ce_ortholog) & !ce_ortholog %in% ce_operons$gene_id ~ "Ortholog is not in an operon"
  ))

ortholog_summary <- group_by(ortholog_operons, operon_type, ortholog_type) %>%
  summarise(total_genes = n())

# detach biomaRt because it has its own select() command that interferes with dplyr and doesn't play well with conflicted
unloadNamespace('biomaRt') 

# Plotting ----------------------------------------------------------------

plot_operons <- bind_rows(final_operonic_genes, non_operonic_genes, ce_operons_cds, ce_non_operons) %>%
  mutate(operon_type = factor(operon_type, levels = c('Novel', 'Shared', 'Assembly', 'Monocistron')))

# stacked bar comparing operons already annotated in assembly and new/shared operons from Iso-Seq
(comparison <- ggplot(plot_operons, aes(x = species, fill = operon_type)) +
    geom_bar() +
    scale_fill_manual(values = c(lacroix_palette("PassionFruit", 3), 'grey'), labels = c('Novel', 'Iso-Seq confirmed', 'Assembly annotated', 'Non-operonic')) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = "Total Genes") +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      axis.title.x = element_blank(),
      axis.text.x = element_text(face = 'italic'),
      legend.title = element_blank()
    ) +
    NULL
)

save_plot(here('..', 'plots', "Fig5A.pdf"), comparison, base_width = 4)

saveRDS(comparison, here('..', 'plots', "Fig5A.rds"))

# stacked bar of which cistrons have orthologs in Ce
(ortholog_plot <- ggplot(ortholog_summary) +
    geom_col(aes(x = ortholog_type, y = total_genes, fill = factor(operon_type, levels = c('Novel', 'Shared', 'Assembly')))) +
    labs(y = "Total *B. malayi* Genes") +
    scale_fill_manual(values = c(lacroix_palette("PassionFruit", 3)),
                      labels = c('New', 'Iso-Seq confirmed', 'Assembly annotated')) +
    scale_x_discrete(labels = c("No\nOrthologs", "Orthologous;\nOperonic", "Orthologous;\nNon-operonic")) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      axis.title.y = element_markdown(),
      axis.title.x = element_blank(),
      legend.text = element_markdown(),
      legend.position = "top",
      legend.title = element_blank()
    ) +
    NULL
)

save_plot(here('..', 'plots', 'Fig5B.pdf'), ortholog_plot, base_width = 4)

saveRDS(ortholog_plot, here('..', 'plots', 'Fig5B.rds'))

# median intergenic distances by type
type_median <- filter(plot_operons, species == "B. malayi") %>% 
  group_by(operon_type) %>%
  summarize(median = median(distance, na.rm = TRUE))

# comparison of intergenic distances between operonic cistrons and monocistronic genes
(type_histogram <- ggplot(data = (filter(plot_operons, species == "B. malayi", cistron %in% c(NA, 1)) %>%
                                    mutate(operon_type = fct_rev(operon_type)))) +
    geom_histogram(aes(x = distance, fill = operon_type), bins = 25) +
    geom_vline(data = type_median, aes(xintercept = median, color = operon_type), linetype = "dashed", size = 1) +
    scale_x_continuous(limits = c(0, 11000), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = rev(c(lacroix_palette("PassionFruit", 3), 'grey')),
                      labels = rev(c('New', 'Iso-Seq confirmed', 'Assembly annotated', 'Non-operonic'))) +
    scale_color_manual(values = c(lacroix_palette("PassionFruit", 3), 'grey'),
                       labels = NULL, guide = "none") +
    labs(x = "Distance Between Genes", y = "Count", fill = "type") +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      legend.text = element_markdown(),
      legend.position = "top",
      legend.title = element_blank()
    ) +
    NULL)

save_plot(here('..', 'plots', 'Fig5C.pdf'), type_histogram, base_width = 6)

saveRDS(type_histogram, here('..', 'plots', 'Fig5C.rds'))

# median intercistronic distances by species
species_median <- filter(plot_operons, operon_type != "Monocistron") %>%
  group_by(species) %>%
  summarize(median = median(distance, na.rm = TRUE))

plot_operons <- plot_operons %>%
  mutate(species = fct_rev(species))

(species_histogram <- ggplot(filter(plot_operons, operon_type != "Monocistron", cistron %in% c(NA, 1))) +
    geom_density(aes(x = distance, y = ..density.., fill = species)) +
    geom_vline(data = species_median, aes(xintercept = median, color = species), linetype = "dashed", size = 1) +
    scale_x_continuous(limits = c(0, 5000), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c('orange', 'purple'),
                      labels = c("*C. elegans*", "*B. malayi*")
    ) +
    scale_color_manual(values = rev(c('orange', 'purple')),
                       labels = NULL, guide = "none") +
    labs(x = "Distance Between Genes", y = "Density", fill = "type") +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.3),
      legend.text = element_markdown(),
      legend.position = "top",
      legend.title = element_blank()
    ) +
    NULL
)

save_plot(here('..', 'plots', 'Fig5D.pdf'), species_histogram, base_width = 4)

saveRDS(species_histogram, here('..', 'plots', 'Fig5D.rds'))
