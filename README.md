# Pacific-Herring
Pacific Herring Genome Assembly
Using PacBio HiFi reads, we will attempt to generate a de novo assembly of the pacific herring genome
# 0. directory set up
```
cd /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021
cd /share/dennislab/projects/
mkdir /share/dennislab/projects/pacific_herring
mkdir /share/dennislab/projects/pacific_herring/fastq
```
# 1. HiFi reads QC
## 1.1 k-mer distance
```
cd /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/
mkdir /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/atlantic_herring
cd /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/atlantic_herring
## retrieve atlantic herring genome
curl -JLO https://www.ncbi.nlm.nih.gov/assembly/GCF_000966335.1/# -o GCF_000966335.fa.gz

mkdir /share/dennislab/projects/pacific_herring/samples_mash
```

## 1.2 Bam to Fastq
Convert HiFi bam files to fastq and evaluate for overall quality
```
mkdir /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021/fastq
cd /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021/fastq

conda activate chifi
bam2fasta -o m64069_210418_020829.hifi_reads. /share/dennislab/projects/vole/data/cshl/m54333U_210606_014450.hifi_reads.bam
mash sketch -k 21 -s 10000 -r -m 1 -o m54333U_210606_014450.hifi_reads m54333U_210606_014450.hifi_reads.fasta.gz
