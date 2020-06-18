#!bin/bash

### Define project directories
proj="Filarid_IsoSeq-ms"

local_dir="${GIT_DATA}/${proj}"

genome_dir="${local_dir}"/genomes
isoform_dir="${local_dir}"/isoforms
polya_dir="${local_dir}"/polya

# poly(A) prep -----------------------------------------------------------------

# convert GTF to BED
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id "";"; }' "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.PRJNA10729.WBPS14.canonical_geneset.gtf | gtf2bed - > "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.gtf.bed
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id "";"; }' "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.PRJEB1797.WBPS14.canonical_geneset.gtf | gtf2bed - > "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.gtf.bed

# get list of scaffolds/lengths
seqkit fx2tab -nl "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.PRJNA10729.WBPS14.genomic.fa | awk -v OFS='\t' '{print $1, 0, $2}' > "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.genome.bed
seqkit fx2tab -nl "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.PRJEB1797.WBPS14.genomic.fa | awk -v OFS='\t' '{print $1, 0, $3}' > "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.genome.bed

# only keep exons or CDS from BED file
cat "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.gtf.bed | sed -E '/\texon\t/!d' | sed '/protein_coding/!d' | awk -v OFS='\t' '{print $1, $2, $3, $4, $6}'| uniq -u > "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.exon.bed
cat "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.gtf.bed | sed -E '/\tCDS|exon\t/!d' | sed '/protein_coding/!d' | awk -v OFS='\t' '{print $1, $2, $3, $4, $6}'| uniq -u > "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.exon.bed

# subtract exons/CDS from genome to get the intergenic regions
bedtools subtract -a "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.genome.bed -b "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.exon.bed > "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.intergenic.bed
bedtools subtract -a "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.genome.bed -b "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.exon.bed > "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.intergenic.bed

# get sequences from the intergenic regions
bedtools getfasta -fi "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.PRJNA10729.WBPS14.genomic.fa -bed "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.intergenic.bed > "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.intergenic.fa
bedtools getfasta -fi "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.PRJEB1797.WBPS14.genomic.fa -bed "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.intergenic.bed > "${genome_dir}"/dirofilaria_immitis/dirofilaria_immitis.intergenic.fa

# run polyAudit ----------------------------------------------------------------

# activate conda environment
source activate polyaudit

### B. malayi

# get poly(A) lengths
python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/polyaudit_output/ -polyA

# re-generate the list of most likely PAS
python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -PAS "${polya_dir}"/dirofilaria_immitis/polyaudit_output/all_transcripts_trim_6mer_filter.csv

# generate PSSM files for each PAS
python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AAAAAA.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AAAUAA.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AAUAAA.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AUAAAA.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AUAAAC.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AUAAAG.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AUAAAU.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/UAAAAU.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/UAAUAA.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/UUUUGU.csv

### D. immitis

# get poly(A) lengths
python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/polyaudit_output/ -polyA

# re-generate the list of most likely PAS
python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -PAS "${polya_dir}"/dirofilaria_immitis/polyaudit_output/all_transcripts_trim_6mer_filter.csv

# generate PSSM files for each PAS
python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AAAAAA.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AAAUAA.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/UUUUUG.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AAAAUU.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AUAAAA.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AAUAAA.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AUAAAU.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/AUAAAG.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/UAAAAU.csv

python polyAudit.py "${isoform_dir}"/dirofilaria_immitis/all_transcripts_trim.fasta "${polya_dir}"/dirofilaria_immitis/ \
  -pssm "${polya_dir}"/dirofilaria_immitis/polyaudit_output/kmer_pssm/UAAAUU.csv
