#!bin/bash

### Define project directories
proj="Filarid_IsoSeq-ms"

local_dir="${GIT_DATA}/${proj}"

genome_dir="${local_dir}"/genomes
alignment_dir="${local_dir}"/alignments
target_dir="${local_dir}"/targets


# Targets ----------------------------------------------------------------------

# convert GTF to BED
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id "";"; }' "${genome_dir}"/brugia_malayi/brugia_malayi.PRJNA10729.WBPS14.canonical_geneset.gtf | gtf2bed - > "${target_dir}"/brugia_malayi/brugia_malayi.gtf.bed
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id "";"; }' "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.PRJEB1797.WBPS14.canonical_geneset.gtf | gtf2bed - > "${target_dir}"/dirofilaria_immitis/dirofilaria_immitis.gtf.bed

# only keep genes
cat "${target_dir}"/brugia_malayi.gtf.bed | sed '/\tgene\t/!d' | sed '/protein_coding/!d' | awk -v OFS='\t' '{print $1, $2, $3, $4, $6}' > "${target_dir}"/brugia_malayi/brugia_malayi.gene.bed
cat "${target_dir}"/dirofilaria_immitis/dirofilaria_immitis.gtf.bed | sed '/\tgene\t/!d' | sed '/protein_coding/!d' | awk -v OFS='\t' '{print $1, $2, $3, $4, $6}' > "${target_dir}"/dirofilaria_immitis/dirofilaria_immitis.gene.bed

# convert SAM to BAM, keep secondary alignments
samtools view -S -b "${alignment_dir}"/brugia_malayi/all_sorted.sam  > "${alignment_dir}"/brugia_malayi/all_sorted.bam
samtools view -S -b "${alignment_dir}"/dirofilaria_immitis/all_sorted.sam  > "${alignment_dir}"/dirofilaria_immitis/all_sorted.bam

# convert BAM to BED (split keeps splicing/exon info)
bedtools bamtobed -split -i "${alignment_dir}"/brugia_malayi/all_sorted.bam > "${target_dir}"/all_sorted.bed
bedtools bamtobed -split -i "${alignment_dir}"/dirofilaria_immitis/all_sorted.bam > "${target_dir}"/dirofilaria_immitis/all_sorted.bed

# get gene intersects
bedtools intersect -a "${target_dir}"/brugia_malayi/brugia_malayi.gene.bed -b "${target_dir}"/brugia_malayi/all_sorted.bed -wa -wb > "${target_dir}"/brugia_malayi/brugia_malayi_gene_intersects.bed
bedtools intersect -a "${target_dir}"/dirofilaria_immitis/dirofilaria_immitis.gene.bed -b "${target_dir}"/dirofilaria_immitis/all_sorted.bed -wa -wb > "${target_dir}"/dirofilaria_immitis/dirofilaria_immitis_gene_intersects.bed
