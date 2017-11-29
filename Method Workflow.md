# Microbiota Processing in Mothur

Project members: Alex Law, Anthony Yan, Fan Tina, Max Faz, Ryan Lou, Tiffany Leung

Link to the report can be found here: ____

## Method workflow

The output of all of our results have been set to: `set.dir(output=/home/micb405/Group10/Project3_2)`

### Preparing Data using FASTQC and Trimmomatic
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
The output is then summarized using `summary.seqs()`. 
```
summary.seqs(fasta=Saanich.trim.contigs.fasta)
```
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       121     121     0       3       1
2.5%-tile:      1       296     296     0       4       8319
25%-tile:       1       297     297     0       4       83182
Median:         1       297     297     0       4       166364
75%-tile:       1       299     299     0       5       249546
97.5%-tile:     1       301     301     1       6       324409
Maximum:        1       602     602     92      301     332727
Mean:   1       297.673 297.673 0.124613        4.48462
# of Seqs:      332727
```
In our newly formed contigs, it is determined that there are 332727 number of sequences that are generally ~297 bases long. 

### Removing poor quality sequences
Only the reads that are completely overlapped are used for downstream analyses which is the reason for setting a minimum and maximum length. Reads with long homopolymers of a certain length are also removed as these are more difficult to detect due to the limitations of Mi-Seq sequencing. The following the steps are used to clean up our sequences for a more careful analysis.

Using `screen.seq()`, sequences were filtered using these parameters:
* Minimum length: 200
* Maximum length: 350
* Maximum homopolymer length: 8
* Maximum number of ambiguous bases: 0 (removal of all ambiguous bases)

The minimum and maximum lengths were chosen based off the desired length of ~297 bases plus an additional 50 bases to accomodate for any adaptors. 
```
screen.seqs(fasta=Saanich.trim.contigs.fasta, group=Saanich.contigs.groups, maxambig=0, maxhomop=8, minlength=200, maxlength=350)
```
Files were outputted with the following file names to help distinguish editted sequences (both trimmed and untrimmed sequences):

* `/home/micb405/Group10/Project3_2/Saanich.trim.contigs.good.fasta`
* `/home/micb405/Group10/Project3_2/Saanich.trim.contigs.bad.accnos`
* `/home/micb405/Group10/Project3_2/Saanich.contigs.good.groups`

### Isolating unique sequences
Having removed all sequences that are poor quality, there is still the likelihood of having duplicate reads present in our data. These duplicate reads are removed to ensure a more computational efficient analysis in subsequent steps.

```
unique.seqs(fasta=Saanich.trim.contigs.good.fasta)
```

Outputted file names to distinguish between unique and duplicate containing reads:
* `/home/micb405/Group10/Project3_2/Saanich.trim.contigs.good.names`
* `/home/micb405/Group10/Project3_2/Saanich.trim.contigs.good.unique.fasta`

### Converting name and group files to count files 
Both name files and group files are converted to a count file to reduce the size of the commands in downstream analyses.

```
count.seqs(name=Saanich.trim.contigs.good.names, group=Saanich.contigs.good.groups)
```
The file was saved as:
```
/home/micb405/Group10/Project3_2/Saanich.trim.contigs.good.count_table
```
To see the result of the removing poor quality sequences as well as isolating duplicates, we run `summary.seqs()` to view the summary table of the newly created count table:

```
summary.seqs(fasta=Saanich.trim.contigs.good.unique.fasta, count=Saanich.trim.contigs.good.count_table)
```
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       204     204     0       3       1
2.5%-tile:      1       296     296     0       4       7693
25%-tile:       1       297     297     0       4       76929
Median:         1       297     297     0       4       153858
75%-tile:       1       299     299     0       5       230786
97.5%-tile:     1       301     301     0       6       300022
Maximum:        1       348     348     0       8       307714
Mean:   1       297.993 297.993 0       4.47064
# of unique seqs:       302601
total # of seqs:        307714
```
The data now contains reads that are between 200-350 bases long, have no ambiguous bases and a maximum homopolymer length of 8 bases. 





