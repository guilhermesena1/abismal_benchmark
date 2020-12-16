# abismal benchmark

This repository contains a set of programs used to compare mapping
results of abismal and multiple WGBS mappers. Here we list what the
programs do and why they are needed

### hash_counter

This is used to count the real and theoretical ratios in a reference
genome given a number of bits. It estimates the probabilities of weak
and strong bases from a genome, and it counts the number of k-mers
under two- and three-letter encodings. It then compares the two- and
three-letter encoding ratios, both in practice (using the k-mer counts
from the genome) and in the developed i.i.d. theory. For instance, to
summarize the hit ratios for 26 bits, run

```
$ hash_counter -b 26 -o stats.tsv /path/to/hg38.fa
```

### mr-to-sam

Converts walt output from mr to sam. This is equivalent to
`format_reads` in methpipe, where the resulting sam file is
single-ended.

```
$ mr-to-sam -o <sam-output> <mr-input>
```
### grepenc

This program is supposed to act like grep but for fasta files and allowing us
to look for two-  or three-letter encodings in the sequence. This is a
debugging tool to search for read locations when reads are mapped ambiguously
by one mapper but uniquely by the other. The `-a` flag allows both two- and
three-letter encodings to be queried. In three-letter encodings, As are 0s, Gs
are 1s and C/Ts are 2s.  For instance, find the number of occurrences of the
two-letter encoding `00110010` on the genome `genome.fa`:

```
$ grepenc -a 2 00110010 /path/to/genome.fa
```

### compare-sam

This program compares the result of a mapping algorithm with the true
status of reads. It takes two inputs: `<truth-sam>` and `<input-sam>`.
If `<truth-sam>` has ambiguously mapped reads, with the ambiguous flag
(0x100) set to 1, the program is capable of reporting false uniques in
the dataset, defined as reads on truth that are ambiguous and at least
as good as the one on input.

```
$ compare-sam <ground_truth_dataset.sam> <input_dataset.sam>
```

### fix-bismark-sam

This program takes the bismark input and corrects some things to make it
suitable for accuracy comparison: It uses the XM tag to subtract from the
edit distance the bases that were counted as mismatches due to bisulfite
conversion. It also formats the read name adding the .1 and .2 to read
names as it was standardized for other mappers. 

```
fix-bismark-sam -o bismark_fixed.sam bismark_original.sam
```

# R script to parse output files

To generate the figures and tables from the manuscript, run the following
command:

```
$ Rscript R/benchmark.r
```

Or alternatively, in an R environment
```
> source('R/benchmark.r')
```

This will generate a `tbl` object, with rows as the tests provided in the
`metadata/tests.txt` file, and columns are the features displayed in the
supplementary tables. The script will also generate a `times` table and a
`mems` table that parses the snakemake outputs showing the total wall time and
maximum resident set size used to run the job.
