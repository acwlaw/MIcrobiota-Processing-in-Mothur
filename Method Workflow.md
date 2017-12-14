# Microbiota Processing in Mothur Method Workflow

Project members: Max Fan, Tina Fan, Alex Law, Tiffany Leung, Ryan Lou, Anthony Yan

Link to the report can be found [here](https://drive.google.com/open?id=1zEXALDYhzRuwM7ROkxOkqXfr0d1kZPTi4xGyXSrXdWo)

## Preparation of data

The output of all of the results have been set to: `set.dir(output=/home/micb405/Group10/Project3_2)`

### Analyzing quality of sequences + trimming using FASTQC and Trimmomatic
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
Mean:           1       297.673 297.673 0.124613        4.48462
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

The minimum and maximum lengths were chosen based off the desired length of ~297 bases plus an additional 50 bases to accomodate for any adaptors. Maximum homopolymer lengths were capped at 8 as anything higher would be indicative of a sequencing error.
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
Mean:           1       297.993 297.993 0       4.47064
# of unique seqs:       302601
total # of seqs:        307714
```
The data now contains reads that are between 200-350 bases long, have no ambiguous bases and a maximum homopolymer length of 8 bases. The number of sequences have been reduced from 332727 to 307714 with a total of 302601 unique sequences.

### Aligning sequences against a database
Sequences were aligned against a Silva database of 16S sequences. Though GreenBase is also a valid database for alignment, it is often of lower quality and contains incomplete 16S sequences making alignment more difficult and error-prone. An additional parameter is given where `flip=T` indicates an  attempt to run the reverse complement of the sequence if the forward sequence falls below a certain threshold. The sequence that had a better alignment is returned.

```
align.seqs(fasta=Saanich.trim.contigs.good.unique.fasta, reference=/home/micb405/data/project_3/databases/silva.nr_v128.align, flip=T)
```
Output of the alignment:
```
Reading in the /home/micb405/data/project_3/databases/silva.nr_v128.align template sequences... DONE.
It took 290 to read  190661 sequences.
Aligning sequences from Saanich.trim.contigs.good.unique.fasta ...
[WARNING]: Some of your sequences generated alignments that eliminated too many bases, a list is provided in
/home/micb405/Group10/Project3_2/Saanich.trim.contigs.good.unique.flip.accnos.
If the reverse compliment proved to be better it was reported.

It took 7000 secs to align 302601 sequences.

Output File Names:
/home/micb405/Group10/Project3_2/Saanich.trim.contigs.good.unique.align
/home/micb405/Group10/Project3_2/Saanich.trim.contigs.good.unique.align.report
/home/micb405/Group10/Project3_2/Saanich.trim.contigs.good.unique.flip.accnos
```
The summary of the sequences following alignment was reported.
```
summary.seqs(fasta=Saanich.trim.contigs.good.unique.align, count=Saanich.trim.contigs.good.count_table)
```
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        6237    16119   8       0       2       1
2.5%-tile:      10366   25318   296     0       4       7693
25%-tile:       10370   25318   297     0       4       76929
Median:         10370   25318   297     0       4       153858
75%-tile:       10370   25432   299     0       5       230786
97.5%-tile:     11886   25436   301     0       6       300022
Maximum:        43102   43116   324     0       8       307714
Mean:           10627.6 25345.5 297.99  0       4.47062
# of unique seqs:       302601
total # of seqs:        307714
```
The data now shows a new starting point for many of the sequences at ~10370 bases. (??)

### Rescreening 
The purpose of rescreening our data is to verify that the sequences overlap the same region using the same start and end positions found in the most recent sequence summary.

```
screen.seqs(fasta=Saanich.trim.contigs.good.unique.align, count=Saanich.trim.contigs.good.count_table, summary=Saanich.trim.contigs.good.unique.summary, start=10370, end=25318)
```
Output:
```
Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.summary
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.align
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.bad.accnos
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.good.count_table
```
Another summary is produced:
```
summary.seqs(fasta=Saanich.trim.contigs.good.unique.good.align, count=Saanich.trim.contigs.good.good.count_table)
```
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        10255   25318   268     0       3       1
2.5%-tile:      10366   25318   297     0       4       6369
25%-tile:       10370   25318   297     0       4       63682
Median:         10370   25318   298     0       4       127364
75%-tile:       10370   25432   299     0       5       191046
97.5%-tile:     10370   25436   301     0       6       248359
Maximum:        10370   26158   324     0       8       254727
Mean:           10369.8 25350.4 298.127 0       4.37627
# of unique seqs:       250706
total # of seqs:        254727

```
### Filtering out overhangs from the sequences
This step filters out sequences that extend the specified regions of alignment. To ensure that overhangs on both sides of the sequences are filtered out, `vertical=T` is set to true. Characters such as '.' and '-' are inserted into the sequence to represent insertions and deletions relative to the alignment database, but do not contribute to downstream analyses and utilize unnecessary computational processing are therefore are removed. '-' is removed by default and `trump=.` removes '.' from the sequences.
```
filter.seqs(fasta=Saanich.trim.contigs.good.unique.good.align, vertical=T, trump=.)
```
```
Length of filtered alignment: 564
Number of columns removed: 49436
Length of the original alignment: 50000
Number of sequences used to construct filter: 250706


Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.filter
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.fasta
```
Following alignment and filtering out sequence overhangs, there is the need to de-replicate again to remove any additional identical sequences that may have formed.

```
unique.seqs(fasta=Saanich.trim.contigs.good.unique.good.filter.fasta, count=Saanich.trim.contigs.good.good.count_table)
```
```
Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.count_table
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.fasta
```

### Organizing sequences into pre-clusters
By organizing similar sequences into pre-clusters, we can further remove unnecessary sequences. `diff=3` indicates sequences that vary up to 3 sequences will be grouped into the same cluster to account for PCR errors. In addition to the grouping of sequences, the groups will be ordered in decreasing abundance.
```
pre.cluster(fasta=Saanich.trim.contigs.good.unique.good.filter.unique.fasta, count=Saanich.trim.contigs.good.unique.good.filter.count_table, diffs=3, processors=10)
```
```
Total number of sequences before pre.cluster was 242619.
pre.cluster removed 147377 sequences.


It took 2891 secs to cluster 242619 sequences.
It took 2905 secs to run pre.cluster.


Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.fasta
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.count_table
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.100m.map
```
Again, a `summary.seqs()` is called to give us a picture of the sequences thus far.
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       562     268     0       3       1
2.5%-tile:      1       564     296     0       4       6369
25%-tile:       1       564     297     0       4       63682
Median:         1       564     297     0       4       127364
75%-tile:       1       564     297     0       5       191046
97.5%-tile:     1       564     298     0       6       248359
Maximum:        5       564     310     0       8       254727
Mean:           1.00016 564     297.02  0       4.41622
# of unique seqs:       95242
total # of seqs:        254727


Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.summary


It took 1 secs to summarize 254727 sequences.
```
### Removing chimeras
The presence of chimeras in our data may inflate the actual number of unique sequences. This is because chimeras are essentially a combination of two separate unique sequences which make the sequence appear as a third unique sequence. For this reason, chimeras must be removed using `chimera.uchime()`. Sequences that are marked as chimeras in one sample are not removed in other samples simply because there is only one sample used and as such, we must pass in `dereplicate=t` as an option. Ten threads were used to attempt to speed up the removal process.
```
chimera.uchime(fasta=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.fasta,count=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.count_table,dereplicate=t,processors=10)
```
```
00:00  67Mb  100.0% Reading Saanich.trim.contigs.good.unique.good.filter.unique.precluster.temp
00:00  67Mb 95.2k sequences
44:53:16  82Mb  100.0% 19887/95241 chimeras found (20.9%)
uchime v4.2.40
by Robert C. Edgar
http://drive5.com/uchime
This code is donated to the public domain.

It took 161598 secs to check 95242 sequences from group 100m.

Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.chimeras
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos
```
`chimera.uchime()` only identifies chimeras. `remove.seqs()` takes care of removing the chimeras.
```
remove.seqs(fasta=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.fasta,count=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table,accnos=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos)
```
```
[NOTE]: The count file should contain only unique names, so mothur assumes your fasta, list and taxonomy files also contain only uniques.


Removed 19887 sequences from your fasta file.
Removed 0 sequences from your count file.


Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table
```
The summary table of the dataset following chimera identification and removal:
```
summary.seqs(fasta=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table)
```
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       562     268     0       3       1
2.5%-tile:      1       564     296     0       4       5791
25%-tile:       1       564     297     0       4       57902
Median:         1       564     297     0       4       115804
75%-tile:       1       564     297     0       5       173705
97.5%-tile:     1       564     298     0       6       225816
Maximum:        5       564     310     0       8       231606
Mean:   1.00002 564     297.022 0       4.41318
# of unique seqs:       75355
total # of seqs:        231606


Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.summary


It took 3 secs to summarize 231606 sequences.
```
### Removal of singletons
Sequences that appear only once in the dataset (singletons) will have minimal impact on downstream analyses and may make for a less efficient processing and should therefore be removed from the current dataset. These sequences are placed in the "rare" output file. Having `cutoff=1` keeps sequences that appear more than one time in the dataset. 
```
mothur > split.abund(fasta=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, cutoff=1)
```
```
Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.rare.count_table
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.abund.count_table
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.rare.fasta
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.abund.fasta
```
Two summaries are produced, the first is the summary of the removed singletons (rare) and the second of the dataset without singletons (abund).

**Removed Singletons (rare)**
```
summary.seqs(fasta=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.rare.fasta, count=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.rare.count_table)
```
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       562     268     0       3       1
2.5%-tile:      1       564     296     0       4       1387
25%-tile:       1       564     297     0       4       13866
Median:         1       564     297     0       4       27731
75%-tile:       1       564     297     0       5       41596
97.5%-tile:     1       564     298     0       6       54075
Maximum:        5       564     310     0       8       55461
Mean:           1.00009 563.999 296.978 0       4.46164
# of unique seqs:       55461
total # of seqs:        55461

Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.rare.summary
```
**Dataset (abund)**
```
summary.seqs(fasta=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.abund.fasta, count=Saanich.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.abund.count_table)
```
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       562     294     0       3       1
2.5%-tile:      1       564     296     0       4       4404
25%-tile:       1       564     297     0       4       44037
Median:         1       564     297     0       4       88073
75%-tile:       1       564     297     0       5       132109
97.5%-tile:     1       564     298     0       6       171742
Maximum:        1       564     301     0       8       176145
Mean:           1       564     297.036 0       4.39793
# of unique seqs:       19894
total # of seqs:        176145


Output File Names:
/home/micb405/Group10/Project3_2/screenSeq_+1/Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.abund.summary
```
The dataset has now been fully cleared of poor quality/undesirable sequences. To make for a more organized data directory, files were renamed to our appropriate depths.
```
cp Saanich.trim.contigs.good.unique.good.filter.unique.precluster.pick.abund.fasta ../OTUFiles/Saanich.100m.final.fasta
```

## Clustering OTU 
A new output directory is specified for clustering data: `set.dir(output=/home/micb405/Group10/Project3_2/OTUFiles)`

### Making OTU
Following sequence cleanup, it is now possible to make our operational taxonomic units (OTUs) for our dataset. The distance between the sequences are calculated and are then clustered according to their distances. `dist.seqs()` calculates a distance matrix of our sequences.
```
dist.seqs(fasta=Saanich.100m.final.fasta, output=lt)
```
```
Output File Names:
/home/micb405/Group10/Project3_2/OTUFiles/Saanich.100m.final.phylip.dist
```
The sequences are then clustered together using the default cutoff of 0.3.
```
cluster(phylip=Saanich.100m.final.phylip.dist, count=Saanich.100m.final.count, cutoff=0.03)
```
```
iter    time    label   num_otus        cutoff  tp      tn      fp      fn      sensitivity     specificity     ppv     npv     fdr     accuracy        mcc    f1score
0       0       0.03    19894   0.03    0       196243932       0       1631739 0       1       0       0.991754        0       0.991754        0       0
1       2       0.03    1069    0.03    1392079 195527170       716762  239660  0.853126        0.996348        0.660116        0.998776        0.339884       0.995167 0.748138        0.744312
2       2       0.03    811     0.03    1462878 195449481       794451  168861  0.896515        0.995952        0.648057        0.999137        0.351943       0.995132 0.759999        0.752303
3       2       0.03    800     0.03    1463610 195449099       794833  168129  0.896963        0.99595 0.648062        0.999141        0.351938        0.9951330.760194        0.752464
4       2       0.03    797     0.03    1463280 195449876       794056  168459  0.896761        0.995954        0.648233        0.999139        0.351767       0.995136 0.76021 0.752508

It took 199 seconds to cluster


Output File Names:
/home/micb405/Group10/Project3_2/OTUFiles/Saanich.100m.final.phylip.opti_mcc.list
/home/micb405/Group10/Project3_2/OTUFiles/Saanich.100m.final.phylip.opti_mcc.steps
/home/micb405/Group10/Project3_2/OTUFiles/Saanich.100m.final.phylip.opti_mcc.sensspec
```
The clusters for our data were then formatted to an OTU table. `label=0.3` specifies the correct species-level of OTUs.

```
make.shared(list=Saanich.100m.final.phylip.opti_mcc.list, count=Saanich.100m.final.count,label=0.03)
```

### Classifying OTU
Once our OTUs have been sorted based off our calculated distance matrix, the following step allows us to classify the OTU and determine what they are. 

`classify.seqs()` allows us to classify the sequences based on a database. Rather than using SILVA as before, the GreenGenes database is used to classify our OTU as the database have a low number of unsclassifieds compared to SILVA and can yield for better classifications of our OTU. We use the cutoff of 80 since we are 80% confident in the classification.

```
classify.seqs(fasta=Saanich.100m.final.fasta, count=Saanich.100m.final.count, template=/home/micb405/data/project_3/databases/gg_13_8_99.fasta, taxonomy=/home/micb405/data/project_3/databases/gg_13_8_99.gg.tax, cutoff=80)
```
```
Output File Names:
/home/micb405/Group10/Project3_2/OTUFiles/Saanich.100m.final.gg.wang.taxonomy
/home/micb405/Group10/Project3_2/OTUFiles/Saanich.100m.final.gg.wang.tax.summary

```
`classify.otu()` uses the referenced classification to get a consensus taxonomy for our OTU.
```
classify.otu(list=Saanich.100m.final.phylip.opti_mcc.list, taxonomy=Saanich.100m.final.gg.wang.taxonomy, count=Saanich.100m.final.count, label=0.03, cutoff=80, basis=otu)
```
```
0.03    797

Output File Names:
/home/micb405/Group10/Project3_2/OTUFiles/Saanich.100m.final.phylip.opti_mcc.0.03.cons.taxonomy
/home/micb405/Group10/Project3_2/OTUFiles/Saanich.100m.final.phylip.opti_mcc.0.03.cons.tax.summary
```
Finally, the OTU data is summarized using `summary.single()`
```
summary.single(shared=Saanich.100m.final.phylip.opti_mcc.shared, label=0.03, calc=nseqs-sobs-coverage)
```
```
Output File Names:
/home/micb405/Group10/Project3_2/OTUFiles/Saanich.100m.final.phylip.opti_mcc.groups.summary
```
`summary.single()` is also able to process the data differently using different calculators such as a Shannon diversity index to get a different perspective on the sample. 
```
summary.single(list=Saanich.100m.final.phylip.an.unique_list.list, calc=shannon)
```
```
Output File Names: 
/home/micb405/Group10/Assignment3_2/Saanich.100m.final.phylip.an.unique_list.summary
```

