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
hash_counter -b 26 -o stats.tsv /path/to/hg38.fa
```

### mr-to-sam

Converts walt output from mr to sam. This is equivalent to
`format_reads` in methpipe, where the resulting sam file is
single-ended.

```
mr-to-sam -o <sam-output> <mr-input>
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
grepenc -a 2 00110010 /path/to/genome.fa
```

### compare-sam

This program compares the result of a mapping algorithm with the true
status of reads. It takes two inputs: `<truth-sam>` and `<inpu-sam>`.
If `<truth-sam>` has ambiguously mapped reads, with the ambiguous flag
(0x100) set to 1, the program is capable of reporting false uniques in
the dataset, defined as reads on truth that are ambiguous and at least
as good as the one on input.

### fix-bismark-sam

This program takes the bismark input and corrects some things to make it
suitable for accuracy comparison: It uses the XM tag to subtract from the
edit distance the bases that were counted as mismatches due to bisulfite
conversion. It also formats the read name adding the .1 and .2 to read
names as it was standardized for other mappers. 

## Possible accuracy issues
1. [low sensitivity paradox] Mappers may report a higher number of
   uniquely mapped reads when they have lower sensitivity because
   among reads for which they find mapping locations, they are unable
   to detect the other locations where the read maps equally well. In
   order to identify and quantify this problem, we need a way to map a
   data set with as high sensitivity as possible, which will generally
   take a long time. We can do this because we have the resources.
   - *abismal* Doing this for abismal means not making any changes in
     the index (it will already be the full size), but does mean
     dramatically increasing the `max_candidates` variable and the
     number of seed shifts. We cannot guarantee full sensitivity to
     any number of mismatches, and we would not be able to as in real
     data this can climp to 10% and still be reasonable. Setting
     `max_candidates` to something like 100k, and using up to 8-10
     shifts should be sufficient. We already allow enough total
     indels.
   - *bsmap* The parameter for `index_interval` must be set to 1 for
     full sensitivity. The cost is a large index data structure. The
     parameter for `-k` which sets `max_kmer_ratio` must be set to 0.
2. [discrepancy in scoring] Some mappers may have more precise scoring
   schemes that better discriminate reads as uniquely mapping. This
   should only affect false reporting of unique mapping reads when
   they are outside of deadzones.
3. [Polypuring/pyrimidine tracts] These will be problematic for
   abismal. The question is how problematic. Most tracts like this are
   deadzones, but some are not. We need a quantification for this.
4. [paired-end] The way paired-ends are scored and reported may
   differ. If the pairs are scored using their distance, we have more
   criteria on which to identify uniquely mapping reads. For each
   mapper, we need to determine how it handles reporting each end of
   an ambiguous mapping paired-end read, or discordant mappings for
   ends.
5. [skipping doing alignment] If the mismatches are few enough, it
   makes sense to skip doing the local alignment. There is a question
   about how few this should be. It seems like 2-3 might be
   reasonable, but otherwise we risk missing a mapping or falsely
   reporting a unique mapping by not doing the alignment step.
