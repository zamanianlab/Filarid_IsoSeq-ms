#!bin/bash

### Define project directories
proj="Filarid_IsoSeq-ms"

gh_dir="${GIT_PATH}/${proj}"
local_dir="${GIT_DATA}/${proj}"

genome_dir="${local_dir}"/genomes
sqanti_dir="${local_dir}"/sqanti_out
novel_gene_dir="${local_dir}"/novel_genes


# Novel genes ------------------------------------------------------------------

# CCSs mapping to potentially novel genes were identified from SQANTI output (intergenic or antisense categories)
# extract sequences from CCSs mapping to putative novel protein coding genes
# use the reference-corrected sequences
seqtk subseq "${sqanti_dir}"/brugia_malayi/bm_all_transcripts_trim_corrected.fasta "${novel_gene_dir}"/brugia_malayi/novel_genes_coding.txt > "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected.fasta
seqtk subseq "${sqanti_dir}"/dirofilaria_immitis/di_all_transcripts_trim_corrected.fasta "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding.txt > "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected.fasta

# blast translated CCSs against Bm proteins, remove anything with a significant hit
blastx -query "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected.fasta -db "${genome_dir}"/brugia_malayi/brugia_malayi.PRJNA10729.WBPS14.protein.fa -num_threads 4 -outfmt 6 > "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastx.out
cat "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastx.out | \
  # sort by isoform, bitscore (descending), evalue, and percent similarity
sort -k1,1 -k12,12gr -k11,11g -k3,3gr | \
  # keep top hit per isoform
sort -u -k1,1 --merge | \
  # remove any isoform that had a significant hit
gawk '{ if ($11 > 1E-3) { print } }' | \
  # sort by blast hit
sort -k2 >  \
  "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastx_filter.out

blastx -query "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected.fasta -db "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.PRJEB1797.WBPS14.protein.fa -num_threads 4 -outfmt 6 > "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastx.out
cat "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastx.out | \
  sort -k1,1 -k12,12gr -k11,11g -k3,3gr | \
  sort -u -k1,1 --merge | \
  gawk '{ if ($11 > 1E-3) { print } }' | \
  sort -k2 >  \
  "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastx_filter.out

# extract surviving isoforms
cat "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastx_filter.out | gawk '{print $1}' > \
  "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_2.txt

cat "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastx_filter.out | gawk '{print $1}' > \
  "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_2.txt

# extract sequences of surviving isoforms
seqtk subseq "${sqanti_dir}"/brugia_malayi/bm_all_transcripts_trim_corrected.fasta "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_2.txt > "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_2.fasta
seqtk subseq "${sqanti_dir}"/dirofilaria_immitis/di_all_transcripts_trim_corrected.fasta "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_2.txt > "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_2.fasta

# blastx against outgroup nematodes
blastx -query "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_2.fasta -db "${genome_dir}"/other/all_outgroup.protein.fa -num_threads 4 -outfmt 6 > "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastx_2.out
cat "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastx_2.out | \
  # sort by isoform, bitscore (descending), evalue, and percent similarity
sort -k1,1 -k12,12gr -k11,11g -k3,3gr | \
  # keep top hit per isoform
sort -u -k1,1 --merge | \
  # keep any isoform that had a significant hit
gawk '{ if ($11 < 1E-3) { print } }' | \
  # sort by blast hit
sort -k2,2 >  \
  "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastx_filter_2.out

blastx -query "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_2.fasta -db "${genome_dir}"/other/all_outgroup.protein.fa -num_threads 4 -outfmt 6 > "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastx_2.out
cat "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastx_2.out | \
  sort -k1,1 -k12,12gr -k11,11g -k3,3gr | \
  sort -u -k1,1 --merge | \
  gawk '{ if ($11 < 1E-3) { print } }' | \
  sort -k2,2 >  \
  "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastx_filter_2.out

# extract surviving isoforms
cat "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastx_filter_2.out | gawk '{print $1}' > \
  "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_3.txt

cat "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastx_filter_2.out | gawk '{print $1}' > \
  "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_3.txt

# extract sequences of surviving isoforms
seqtk subseq "${sqanti_dir}"/brugia_malayi/bm_all_transcripts_trim_corrected.fasta "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_3.txt > "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_3.fasta
seqtk subseq "${sqanti_dir}"/dirofilaria_immitis/di_all_transcripts_trim_corrected.fasta "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_3.txt > "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_3.fasta

# blast against the Bm genome to see how many new loci were identified
blastn -query "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_3.fasta -db "${genome_dir}"/brugia_malayi/brugia_malayi.PRJNA10729.WBPS14.genomic.fa -num_threads 4 -outfmt 6 > "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastn_3.out
cat "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastn_3.out |\
  sort -k1,1 -k2,2 -k9,9 > \
  "${novel_gene_dir}"/brugia_malayi/novel_genes_coding_corrected_blastn_filter_3.out

blastn -query "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_3.fasta -db "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.PRJEB1797.WBPS14.genomic.fa -num_threads 4 -outfmt 6 > "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastn_3.out
cat "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastn_3.out |\
  # sort by isoform, blast hit, and scaffold start site
sort -k1,1 -k2,2 -k9,9 > \
  "${novel_gene_dir}"/dirofilaria_immitis/novel_genes_coding_corrected_blastn_filter_3.out
