#!bin/bash

### Define project directories
proj="Filarid_IsoSeq-ms"

local_dir="${GIT_DATA}/${proj}"

genome_dir="${local_dir}"/genomes
operon_dir="${local_dir}"/operon
alignment_dir="${local_dir}"/alignments
receptor_dir="${local_dir}"/receptor

# Operons ----------------------------------------------------------------------

# convert from alignment SAM to BED, remove secondary alignments
samtools view -b -F 256 "${alignment_dir}"/brugia_malayi/all_sorted.bam | bedtools bamtobed > "${operon_dir}"/brugia_malayi/all_sorted.bed

# get regions where alignments intersect with genes
bedtools intersect -a "${receptor_dir}"/brugia_malayi/brugia_malayi.gene.bed -b "${operon_dir}"/brugia_malayi/all_sorted.bed -wa -wb > "${operon_dir}"/brugia_malayi/brugia_malayi_gene_intersects.bed

# get sequences from Bm transcripts that may be in operons
seqtk subseq "${genome_dir}"/brugia_malayi/brugia_malayi.PRJNA10729.WBPS14.protein.fa "${operon_dir}"/brugia_malayi/polycistron_transcripts.txt > "${operon_dir}"/brugia_malayi/polycistron_transcripts.fasta

# blast against the Ce genome to see if the non-overlapping segments (i.e., cistrons) hit to the same protein
blastp -query "${operon_dir}"/brugia_malayi/polycistron_transcripts.fasta -db "${genome_dir}"/other/caenorhabditis_elegans.PRJNA13758.WBPS14.protein.fa -num_threads 4 -outfmt 6 > "${operon_dir}"/brugia_malayi/polycistron_transcripts_blastp_1.out

# blast concatenated transcripts against the Ce genome to see if an operon hits to one protein
blastp -query "${operon_dir}"/brugia_malayi/operons.fasta -db "${genome_dir}"/other/caenorhabditis_elegans.PRJNA13758.WBPS14.protein.fa -num_threads 4 -outfmt 6 > "${operon_dir}"/brugia_malayi/operons_blastp_1.out
