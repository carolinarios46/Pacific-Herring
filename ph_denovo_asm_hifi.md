# Pacific Herring Genome Assembly

Using PacBio HiFi reads, we will attempt to generate a de novo assembly of the Pacific Herring genome

There is approx. 2 million years of divergence between Pacific and Atlantic Herring
The Altantic Herring has 26 autosomes with a total size of 726 Mb

## 0. directory set up
```
mkdir /share/dennislab/projects/pacific_herring/fastq
cd /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021
```
## 1. HiFi reads QC
### 1.1 k-mer distance
```
mkdir /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/atlantic_herring
cd /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/atlantic_herring
## retrieve atlantic herring genome
curl -JLO [https://www.ncbi.nlm.nih.gov/assembly/GCF_000966335.1/#](https://www.ncbi.nlm.nih.gov/assembly/GCF_900700415.2#) -o GCF_000966335.fa.gz

mkdir /share/dennislab/projects/pacific_herring/samples_mash
cd /share/dennislab/projects/pacific_herring/samples_mash

```

### 1.2 Bam to Fastq
Convert HiFi bam files to fastq and evaluate for overall quality
```
mkdir /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021/fastq
cd /share/dennislab-backedup/pacbio/pacific_herring/ph_genome/april19_2021/fastq

conda activate chifi
bam2fasta -o m64069_210418_020829.hifi_reads. /share/dennislab/projects/vole/data/cshl/m54333U_210606_014450.hifi_reads.bam
mash sketch -k 21 -s 10000 -r -m 1 -o m54333U_210606_014450.hifi_reads m54333U_210606_014450.hifi_reads.fasta.gz

