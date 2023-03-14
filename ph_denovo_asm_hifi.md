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
## 1. HiFi reads QC
### 1.1 k-mer distance
```
mkdir /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/atlantic_herring
cd /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/atlantic_herring
## retrieve atlantic herring genome
curl -JLO [https://www.ncbi.nlm.nih.gov/assembly/GCF_000966335.1/#](https://www.ncbi.nlm.nih.gov/assembly/GCF_900700415.2#) -o GCF_000966335.fa.gz

mkdir /share/dennislab/projects/pacfic_herring/denovo_asm
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

conda activate nanopore
NanoPlot -t 10 --fastq herring_DNA.fastq.gz --N50 -f png -o nanoplot_herring
```

## 2. Genome profiling
```
conda activate jellyfish
mkdir /share/dennislab/projects/pacfic_herring/denovo_asm/herring_hifi
cd /share/dennislab/projects/pacfic_herring/denovo_asm/herring_hifi

jellyfish count -C -m 21 -s 1000000000 -t 10 <(gunzip -c /share/dennislab/projects/pacific_herring/denovo_asm/herring_hifi/fastq/herring_DNA.fastq.gz) -o herring_DNA.jf

jellyfish histo -t 10 herring_DNA.jf > herring_DNA.histo
```
### 3. De novo assembly

```
mkdir /share/dennislab/projects/pacfic_herring/denovo_asm/herring_hifiasm
cd /share/dennislab/projects/pacfic_herring/denovo_asm/herring_hifiasm

conda activate voles # hifiasm 0.18.5-r499

hifiasm \
-o herring_DNA.asm \
--h1 /Undetermined_Undetermined_H7Y75CCX2_L4_1.fq.gz \
--h2 /Undetermined_Undetermined_H7Y75CCX2_L4_2.fq.gz \
-t 64 \
-l 1 \
../herring_hifi/fastq/herring_DNA.fastq.gz
# --purge-max 

awk '/^S/{print ">"$2;print $3}' herring_DNA.asm.hic.p_ctg.gfa > herring_DNA.asm.hic.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' herring_DNA.asm.hic.hap1.p_ctg.gfa > herring_DNA.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' herring_DNA.asm.hic.hap2.p_ctg.gfa > herring_DNA.asm.hic.hap2.p_ctg.fa
```
