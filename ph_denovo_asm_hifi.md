# Pacific Herring Genome Assembly

Using PacBio HiFi reads, we will attempt to generate a de novo assembly of the Pacific Herring genome

There is approx. 2 million years of divergence between Pacific and Atlantic Herring

The Altantic Herring has 26 autosomes with a total size of 726 Mb

## 0. directory/environment set up
```
#These are the pacific herring HiFi reads
/share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021/m64069_210418_020829.hifi_reads.bam
/share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021/m64202e_210514_194300.hifi_reads.bam

#These are the pacific herring Hi-C reads
/share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L4_1.fq.gz
/share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L4_2.fq.gz
/share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L5_1.fq.gz
/share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L5_2.fq.gz

# setting up the conda environment
cd /share/dennislab/projects/pacific_herring/
source /share/dennislab/programs/dennis-miniconda/etc/profile.d/conda.sh
```
## 1. HiFi and Hi-C reads QC
### 1.1 k-mer distance
```
mkdir /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/atlantic_herring
cd /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/atlantic_herring
## retrieve atlantic herring genome
curl -JLO [https://www.ncbi.nlm.nih.gov/assembly/GCF_000966335.1/#](https://www.ncbi.nlm.nih.gov/assembly/GCF_900700415.2#) -o GCF_000966335.fa.gz

mkdir /share/dennislab/projects/pacific_herring/denovo_asm
mkdir /share/dennislab/projects/pacific_herring/denovo_asm/samples_mash
cd /share/dennislab/projects/pacific_herring/samples_mash

conda activate chifi
bam2fasta -o m64069_210418_020829.hifi_reads /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021/m64069_210418_020829.hifi_reads.bam
mash sketch -k 21 -s 10000 -r -m 1 -o m64069_210418_020829.hifi_reads m64069_210418_020829.hifi_reads.fasta.gz 

mash sketch -k 21 -s 10000 -r -m 1 -o m64202e_210514_194300 /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021/m64202e_210514_194300.hifi_reads.bam
mash sketch -k 21 -s 10000 -r -m 1 -o m64202e_210514_194300

mash paste combined.msh *.msh
mash dist -t combined.msh combined.msh > combined.tbl

module load R/4.0.1
Rscript plot.R $out
```

### 1.2 Bam to Fastq
Convert HiFi bam files to fastq and evaluate for overall quality
```
mkdir /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifi
mkdir /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifi/fastq
cd /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifi/fastq


conda activate chifi
bam2fastq -o herring_DNA \
/share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021/m64069_210418_020829.hifi_reads.bam \
/share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021/m64202e_210514_194300.hifi_reads.bam

conda activate nanoplot
NanoPlot -t 10 --fastq herring_DNA.fastq.gz --N50 -f png -o nanoplot_herring
```
### 1.3 Fastqc
Evaluate Hi-C reads for overall quality
```
mkdir /share/dennislab/projects/pacific_herring/denovo_asm/herring_hic
mkdir /share/dennislab/projects/pacific_herring/denovo_asm/herring_hic/fastqc
cd /share/dennislab/projects/pacific_herring/denovo_asm/herring_hic/fastqc

conda activate fastqc #v0.12.1

fastqc /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L4_1.fq.gz
fastqc /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L4_2.fq.gz
fastqc /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L5_1.fq.gz
fastqc /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L5_2.fq.gz
```
## 2. Genome profiling
```
conda activate jellyfish

cd /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifi

jellyfish count -C -m 21 -s 1000000000 -t 10 <(gunzip -c /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifi/fastq/herring_DNA.fastq.gz) -o herring_DNA.jf

jellyfish histo -t 10 herring_DNA.jf > herring_DNA.histo
```
Histograms uploaded here: http://qb.cshl.edu/genomescope/genomescope2.0/

Output here: http://genomescope.org/genomescope2.0/analysis.php?code=mBBh0hGmHtxqWwFp5hkP

## 3. De novo assembly

### 3.1 Hifiasm
```
mkdir /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm
cd /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm

cat /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L4_1.fq.gz /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L5_1.fq.gz > /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L4_L5_1.fq.gz

cat /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L4_2.fq.gz /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L5_2.fq.gz > /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L4_L5_2.fq.gz

conda activate hifiasm-latest # hifiasm 0.18.5-r499

hifiasm \
-o herring_DNA.asm \
--h1 /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L4_L5_1.fq.gz \
--h2 /share/dennislab-backedup/illumina/pacific_herring/Undetermined_Undetermined_H7Y75CCX2_L4_L5_2.fq.gz \
-t 64 \
-l 1 \
../herring_hifi/fastq/herring_DNA.fastq.gz
# --purge-max 

awk '/^S/{print ">"$2;print $3}' herring_DNA.asm.hic.p_ctg.gfa > herring_DNA.asm.hic.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' herring_DNA.asm.hic.hap1.p_ctg.gfa > herring_DNA.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' herring_DNA.asm.hic.hap2.p_ctg.gfa > herring_DNA.asm.hic.hap2.p_ctg.fa
```
### 3.2 Draft assembly evaluation
contiguity
```
conda activate quast

cd /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm

quast --threads 32 herring_DNA.asm.hic.p_ctg.fa --est-ref-size 2051817670 -o output_quast_primary
quast --threads 32 herring_DNA.asm.hic.hap1.p_ctg.fa --est-ref-size 2051817670 -o output_quast_hap1
quast --threads 32 herring_DNA.asm.hic.hap2.p_ctg.fa --est-ref-size 2051817670 -o output_quast_hap2
```
completeness
```
conda activate busco #v5.4.4

cd /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm

busco -c 32 -i herring_DNA.asm.hic.p_ctg.fa -l eukaryota_odb10 -m geno --out output_busco_primary
busco -c 32 -i herring_DNA.asm.hic.hap1.p_ctg.fa -l eukaryota_odb10 -m geno --out output_busco_hap1
busco -c 32 -i herring_DNA.asm.hic.hap2.p_ctg.fa -l eukaryota_odb10 -m geno --out output_busco_hap2

busco -c 32 -i herring_DNA.asm.hic.p_ctg.fa -l actinopterygii_odb10 -m geno --out output_busco_actinopterygii_primary
busco -c 32 -i herring_DNA.asm.hic.hap1.p_ctg.fa -l actinopterygii_odb10 -m geno --out output_actinopterygii_busco_hap1
busco -c 32 -i herring_DNA.asm.hic.hap2.p_ctg.fa -l actinopterygii_odb10 -m geno --out output_actinopterygii_busco_hap2
```
correctness
```
export PATH=/share/dennislab/programs/meryl-1.0/Linux-amd64/bin:$PATH
export MERQURY=/share/dennislab/programs/merqury
sh $MERQURY/best_k.sh 2000000000 # Best k-mer size ~20

cd /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/output_meryl_primary

# we also have assembly from 10x assembly: phgenome_ncbi.fasta

# Illumina fastqs: /share/dennislab/users/jagill/shortread2020

meryl k=$k count *.fastq.gz output .meryl

meryl k=20 count output read1.meryl /share/dennislab/users/jagill/shortread2020/PH_BR_CKDL190143414-1a-D707-AK1545_H725HCCX2_L4_R1_paired.fq.gz
meryl k=20 count output read2.meryl /share/dennislab/users/jagill/shortread2020/PH_BR_CKDL190143414-1a-D707-AK1545_H725HCCX2_L4_R2_paired.fq.gz
meryl union-sum output illumina.meryl read*.meryl

ln -s $MERQURY/merqury.sh # Link merqury
./merqury.sh /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/output_meryl_primary/illumina.meryl /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/herring_DNA.asm.hic.p_ctg.fa meryloutput

## Use nucmer and dnadiff to do an assembly-to-assembly comparison

module load mummer/4.0.0rc1
nucmer -t 60 -l 100 -c 500 -p output /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/herring_DNA.asm.hic.p_ctg.fa /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/nucmer/phgenome_ncbi.fasta
dnadiff -d output.delta -p output.dnadiff

## Check for contamination
cd /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/kraken2

module load kraken2/2.1.2

kraken2 \
--db db \
--threads 10 \
--output assembly-kraken2-out.txt \
--report assembly-kraken2-report.txt \
/share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/herring_DNA.asm.hic.p_ctg.fa 

# create atlantic herring database

# mkdir atlanticDB

kraken2-build --download-taxonomy --db atlanticDB
kraken2-build --download-library archaea --db atlanticDB
kraken2-build --download-library bacteria --db atlanticDB
kraken2-build --download-library plasmid --db atlanticDB
kraken2-build --download-library viral --db atlanticDB
kraken2-build --download-library human --db atlanticDB
kraken2-build --download-library fungi --db atlanticDB
kraken2-build --download-library plant --db atlanticDB
kraken2-build --download-library protozoa --db atlanticDB
kraken2-build --download-library UniVec_Core --db atlanticDB
kraken2-build --add-to-library --db atlanticherring/GCF_900700415.2_Ch_v2.0.2_genomic.fna.gz atlanticDB
kraken2-build --build --db ./atlanticDB

module load kraken2/2.1.2

kraken2 \
--db atlanticDB \
--threads 10 \
--output assembly-ah-kraken2-out.txt \
--report assembly-ah-kraken2-report.txt \
/share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/herring_DNA.asm.hic.p_ctg.fa

7950
```

## 4. Duplication Purging

```
cd /share/dennislab/projects/pacific_herring/denovo_asm/purge_dups

## Primary assembly 

cd /share/dennislab/projects/pacific_herring/denovo_asm/purge_dups/primary

# Mapping HiFi reads to primary assembly using Minimap2 # version 2.17
minimap2 -t 32 -xasm20 /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/herring_DNA.asm.hic.p_ctg.fa /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifi/fastq/herring_DNA.fastq.gz | gzip -c - > herring_DNA.toPrimary.paf.gz

# Calculate coverage cutoff, base-level read depth, and create read depth histogram for PacBio data (calcuts+pbcstats)
# PBCSTAT base coverage: contains the base-level coverage information
# Calcuts-cutoff: incldues the thresholds calculated by purge dups
# Histogram plot
pbcstat herring_DNA.toPrimary.paf.gz # produces PB.base.cov and PB.stat files calcuts.log

# Segment draft assembly into contigs by cutting at blocks of Ns, and use minimap2 to generate an all by all self-alignment
split_fa /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/herring_DNA.asm.hic.p_ctg.fa > herring_DNA.asm.hic.p_ctg.fa.split
minimap2 -t 32 -xasm5 -DP herring_DNA.asm.hic.p_ctg.fa.split herring_DNA.asm.hic.p_ctg.fa.split | gzip -c - > herring_DNA.asm.hic.p_ctg.fa.split.self.paf.gz

# Use the self alignments and the cutoffs for identifying the haplotypic duplications. 
# This command will only remove haplotypic duplicaitons at the ends of the contigs.
# If you also want to remove the duplications in the middle, please remove -e option at your own risk,
# It may delete false duplicaitons
purge_dups -2 -T cutoffs -c PB.base.cov herring_DNA.asm.hic.p_ctg.fa.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/herring_DNA.asm.hic.p_ctg.fa

# Evaluation
conda activate quast
quast --threads 32 purged.fa --est-ref-size 506212431 -o output_quast_purged
conda activate busco
busco -c 32 -i purged.fa -l eukaryota_odb10 -m geno --out output_busco_eukaryota_purged
busco -c 32 -i purged.fa -l actinopterygii_odb10 -m geno --out output_busco_actinopterygii_purged 

## Haplotype 1

cd /share/dennislab/projects/pacific_herring/denovo_asm/purge_dups/hap1

# Mapping HiFi reads to primary assembly using Minimap2 # version 2.17
minimap2 -t 32 -xasm20 /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/herring_DNA.asm.hic.hap1.p_ctg.fa /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifi/fastq/herring_DNA.fastq.gz | gzip -c - > herring_DNA.toHap1.paf.gz

# Calculate coverage cutoff, base-level read depth, and create read depth histogram for PacBio data (calcuts+pbcstats)
# PBCSTAT base coverage: contains the base-level coverage information
# Calcuts-cutoff: incldues the thresholds calculated by purge dups
# Histogram plot
pbcstat herring_DNA.toHap1.paf.gz # produces PB.base.cov and PB.stat files calcuts.log

# Segment draft assembly into contigs by cutting at blocks of Ns, and use minimap2 to generate an all by all self-alignment
split_fa /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/herring_DNA.asm.hic.hap1.p_ctg.fa > herring_DNA.asm.hic.hap1.p_ctg.fa.split
minimap2 -t 32 -xasm5 -DP herring_DNA.asm.hic.hap1.p_ctg.fa.split herring_DNA.asm.hic.hap1.p_ctg.fa.split | gzip -c - > herring_DNA.asm.hic.hap1.p_ctg.fa.split.self.paf.gz

# Use the self alignments and the cutoffs for identifying the haplotypic duplications. 
# This command will only remove haplotypic duplicaitons at the ends of the contigs.
# If you also want to remove the duplications in the middle, please remove -e option at your own risk,
# It may delete false duplicaitons
purge_dups -2 -T cutoffs -c PB.base.cov herring_DNA.asm.hic.p_ctg.fa.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifiasm/herring_DNA.asm.hic.p_ctg.fa

# Evaluation
conda activate quast
quast --threads 32 purged.fa --est-ref-size 506212431 -o output_quast_purged
conda activate busco
busco -c 32 -i purged.fa -l eukaryota_odb10 -m geno --out output_busco_eukaryota_purged
busco -c 32 -i purged.fa -l actinopterygii_odb10 -m geno --out output_busco_actinopterygii_purged 
```
