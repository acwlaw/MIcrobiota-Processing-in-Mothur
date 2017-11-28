# Microbiota Processing in Mothur

Project members: Alex Law, Anthony Yan, Fan Tina, Max Faz, Ryan Lou, Tiffany Leung
Link to the report can be found here: ____

## Method workflow

### Preparing Data using FASTQC
Datasets were run under FASTQC to determine the overall quality of the sequences. It was observed that majority of bases are above a quality score of 6. The tool *Trimmomatic* was used to eliminate the poor qualities at the end of the reads.

```
trimmomatic PE -phred33 -trimlog SI072_S3_100_log -threads 4 \
/home/micb405/data/project_3/SI072_S3_100_1.fastq \
/home/micb405/data/project_3/SI072_S3_100_2.fastq \
/home/micb405/Group10/Project3/Trimmed/SI072_S3_100_pe.1.fq \
/home/micb405/Group10/Project3/Trimmed/SI072_S3_100_se.1.fq \
/home/micb405/Group10/Project3/Trimmed/SI072_S3_100_pe.2.fq \
/home/micb405/Group10/Project3/Trimmed/SI072_S3_100_se.2.fq \
ILLUMINACLIP:/home/linuxbrew/.linuxbrew/share/trimmomatic/adapters/TruSeq3-PE.fa:2:3:10 \
LEADING:6 TRAILING:6 SLIDINGWINDOW:4:15 MINLEN:100 \
1>/home/micb405/Group10/Project3/Trimmed/log/trim.stdout \
2>/home/micb405/Group10/Project3/Trimmed/log/trim.stderr &
```

### Combining paired end reads
It is necessary to combine the two paired-end runs from an Illumina MiSeq into a single, longer read for every sequence to yield a contig. 

To do this, we must first prepare our sequences that create a tab-delimitted file of our trimmed FastQC sequences. Files were renamed to the appropriate depth.
```
make.file(inputdir=/home/micb405/data/project_3, prefix=Saanich)
```
Contigs were then created using the newly created files. `make.contigs()` utilizes a simple algorithim that aligns the pairs of sequences and identifies pairs that do not match. Depending on the quality score of the paired bases, the algorithm may assign the consensus base to an N.
```
make.contigs(file=Saanich.files)
```



