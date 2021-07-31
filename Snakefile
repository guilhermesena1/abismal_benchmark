# tsv with 1st column = fq srr, 2nd column = species to map to
input_dataset_tsv = "metadata/tests.txt"

# map results go to local disk
DIR_OUTPUT_MAP = \
  "/panfs/qcb-panasas/desenabr/ab/results/benchmark-files"

# everything else (that does not depend on speed and heavy IO) to cmb-06
DIR_OUTPUT_REST = "/project/andrewds_103/desenabr/abismal_devel/tests"

# scratch space for heavy temp files
DIR_DUMP = "/scratch2/desenabr/snakemake_dump"

# root directory for reference genome files
genome_dir = "/panfs/qcb-panasas/desenabr/ref_genomes/{species}/"
genome_dir_bismark = \
  "/panfs/qcb-panasas/desenabr/ref_genomes/{wildcards.species}/"

# reference genome fasta files
fa = genome_dir + "genome.fa"
fa_minimap2 = genome_dir + "genome_bisulfite.fa"
fai = genome_dir + "genome.fa.fai"

############### FASTQ FILES #################
# directory to fastq files
fa_dir = "/panfs/qcb-panasas/desenabr/ab/reads"

fa_single = fa_dir + "/single/{sample}_{species}.fq"
fa_single_minimap2 = fa_dir + "/single/{sample}_{species}_C_to_T.fq"

# pattern of r1 and r2 filenames
fa_paired_r1 = fa_dir + "/paired/{sample}_{species}_1.fq"
fa_paired_r2 = fa_dir + "/paired/{sample}_{species}_2.fq"

# pbat paths
fa_rpbat_r1 = fa_dir + "/rpbat/{sample}_{species}_1.fq"
fa_rpbat_r2 = fa_dir + "/rpbat/{sample}_{species}_2.fq"

############ INDEX PATHS ###############
# directory where Bisulfite_Genome lies for bismark
index = genome_dir + "index"
index_bismark = genome_dir_bismark + "index"

# walt index
index_walt = index + "/walt.dbindex"

# abismal index file
index_abismal = index + "/abismal_min.idx"

# minimap2 index file
index_minimap2 = index + "/minimap2.idx"

# hisat-3n
index_hisat = index + "/hisat"

########### END IND #############
species = ["hg38", "mm10", "danre11", "galgal6", "tair10"]
protocols = ["wgbs_single", "wgbs_paired", "wgbs_rpbat"]

# init the hash map
SAMPLES = {}
for i in species:
  SAMPLES[i] = {}
  for j in protocols:
    SAMPLES[i][j] = []

# populate the hash map
for x in open(input_dataset_tsv):
  y = x.strip().split()
  if y[0] != "srr": # the header
    species = y[1]
    protocol = y[2]
    SAMPLES[species][protocol].append(y[0])

# mappers used in benchmarking
MAPPERS_SINGLE = ["abismal", "bismark", "bsmap", "walt", "bwa", "hisat_3n", "minimap2"]

MAPPERS_PAIRED = ["abismal", "bismark", "bsmap", "walt", "bwa", "hisat_3n"]
MAPPERS_PAIRED = ["abismal", "bismark", "walt", "hisat_3n"]

MAPPERS_RPBAT = ["abismal", "bismark", "bsmap"]
MAPPERS_RPBAT = []

# output files
OUT_FILES = ["samstats", "bsrate", "levels"]

################### BEGIN RULE ALL ####################
#
rule all:
  input:
    # SE outputs
    expand("{dirr}/{mapper}/{sample}_hg38.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_SINGLE, sample = SAMPLES["hg38"]["wgbs_single"], outfile = OUT_FILES),
    expand("{dirr}/{mapper}/{sample}_mm10.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_SINGLE, sample = SAMPLES["mm10"]["wgbs_single"], outfile = OUT_FILES),
    expand("{dirr}/{mapper}/{sample}_danre11.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_SINGLE, sample = SAMPLES["danre11"]["wgbs_single"], outfile = OUT_FILES),
    expand("{dirr}/{mapper}/{sample}_galgal6.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_SINGLE, sample = SAMPLES["galgal6"]["wgbs_single"], outfile = OUT_FILES),
    expand("{dirr}/{mapper}/{sample}_tair10.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_SINGLE, sample = SAMPLES["tair10"]["wgbs_single"], outfile = OUT_FILES),
    # PE outputs
    expand("{dirr}/{mapper}/{sample}_hg38.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_PAIRED, sample = SAMPLES["hg38"]["wgbs_paired"], outfile = OUT_FILES),
    expand("{dirr}/{mapper}/{sample}_mm10.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_PAIRED, sample = SAMPLES["mm10"]["wgbs_paired"], outfile = OUT_FILES),
    expand("{dirr}/{mapper}/{sample}_danre11.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_PAIRED, sample = SAMPLES["danre11"]["wgbs_paired"], outfile = OUT_FILES),
    expand("{dirr}/{mapper}/{sample}_galgal6.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_PAIRED, sample = SAMPLES["galgal6"]["wgbs_paired"], outfile = OUT_FILES),
    expand("{dirr}/{mapper}/{sample}_tair10.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_PAIRED, sample = SAMPLES["tair10"]["wgbs_paired"], outfile = OUT_FILES),
    # RPBAT outputs
    expand("{dirr}/{mapper}/{sample}_hg38.{outfile}", dirr = DIR_OUTPUT_MAP, mapper = MAPPERS_RPBAT, sample = SAMPLES["hg38"]["wgbs_rpbat"], outfile = OUT_FILES),

############### ABISMAL #####################
rule abismal_map_single:
  input:
    r = fa_single,
    fa = fa
  output:
    sam = temp(DIR_OUTPUT_MAP + "/abismal/{sample}_{species}.sam"),
    mapstats = DIR_OUTPUT_MAP + "/abismal/{sample}_{species}.sam.mapstats",
  benchmark:
    DIR_OUTPUT_MAP + "/abismal/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    abismal -v -t 16 -g {input.fa} -o {output.sam} -s {output.mapstats} {input.r}
    """

rule abismal_map_paired:
  input:
    r1 = fa_paired_r1,
    r2 = fa_paired_r2,
    fa = fa
  output:
    sam = temp(DIR_OUTPUT_MAP + "/abismal/{sample}_{species}.sam"),
    mapstats = DIR_OUTPUT_MAP + "/abismal/{sample}_{species}.sam.mapstats",
  benchmark:
    DIR_OUTPUT_MAP + "/abismal/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    abismal -v -t 16 -g {input.fa} -o {output.sam} -s {output.mapstats} {input.r1} {input.r2}
    """

rule abismal_map_rpbat:
  input:
    r1 = fa_rpbat_r1,
    r2 = fa_rpbat_r2,
    fa = fa
  output:
    sam = temp(DIR_OUTPUT_MAP + "/abismal/{sample}_{species}.sam"),
    mapstats = DIR_OUTPUT_MAP + "/abismal/{sample}_{species}.sam.mapstats",
  benchmark:
    DIR_OUTPUT_MAP + "/abismal/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    abismal -v -t 16 -R -g {input.fa} -o {output.sam} -s {output.mapstats} {input.r1} {input.r2}
    """

############### BISMARK #################
rule bismark_map_single:
  input:
    r = fa_single
  output:
    bam = DIR_OUTPUT_MAP + \
      "/bismark/{sample}_{species}_bismark_bt2.bam",
    rep = temp(DIR_OUTPUT_MAP + \
      "/bismark/{sample}_{species}_bismark_bt2_SE_report.txt")
  benchmark:
    DIR_OUTPUT_MAP + "/bismark/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    bismark --temp_dir /scratch2/desenabr/bismark_dump --parallel 8 \
    -o %s/bismark --local --bowtie2 --icpc %s {input.r}
    """ % (DIR_OUTPUT_MAP, index_bismark)


rule bismark_map_paired:
  input:
    r1 = fa_paired_r1,
    r2 = fa_paired_r2
  output:
    bam = DIR_OUTPUT_MAP + \
      "/bismark/{sample}_{species}_1_bismark_bt2_pe.bam",
    rep = temp(DIR_OUTPUT_MAP + \
      "/bismark/{sample}_{species}_1_bismark_bt2_PE_report.txt")
  benchmark:
    DIR_OUTPUT_MAP + "/bismark/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    bismark --temp_dir /scratch2/desenabr/bismark_dump --parallel 8 \
    -o %s/bismark --local --bowtie2 --icpc -I 32 -X 3000 \
    -1 {input.r1} -2 {input.r2} %s
    """ % (DIR_OUTPUT_MAP, index_bismark)

rule bismark_map_rpbat:
  input:
    r1 = fa_rpbat_r1,
    r2 = fa_rpbat_r2
  output:
    bam = DIR_OUTPUT_MAP + \
      "/bismark/{sample}_{species}_1_bismark_bt2_pe.bam",
    rep = temp(DIR_OUTPUT_MAP + \
          "/bismark/{sample}_{species}_1_bismark_bt2_PE_report.txt")
  benchmark:
    DIR_OUTPUT_MAP + "/bismark/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    bismark --temp_dir /scratch2/desenabr/bismark_dump --bowtie2 \
            --non_directional --parallel 4 \
            -o %s/bismark --local --icpc -I 32 -X 3000 -1 {input.r1} -2 {input.r2} %s
    """ % (DIR_OUTPUT_MAP, index_bismark)

rule bismark_to_sam_se:
  input:
    bam = DIR_OUTPUT_MAP + "/bismark/{sample}_{species}_bismark_bt2.bam",
    rep = DIR_OUTPUT_MAP + "/bismark/{sample}_{species}_bismark_bt2_SE_report.txt"
  output:
    sam = temp(DIR_OUTPUT_MAP + "/bismark/{sample}_{species}.sam"),
    rep = DIR_OUTPUT_MAP + "/bismark/{sample}_{species}.mapstats"
  shell:
    """
    samtools view -h {input.bam} -o {output.sam}temp && \
    fix-bismark-sam -s -o {output.sam} {output.sam}temp && \
    cp {input.rep} {output.rep} && \
    rm {output.sam}temp
    """

rule bismark_to_sam_pe:
  input:
    bam = DIR_OUTPUT_MAP + "/bismark/{sample}_{species}_1_bismark_bt2_pe.bam",
    rep = DIR_OUTPUT_MAP + "/bismark/{sample}_{species}_1_bismark_bt2_PE_report.txt"
  output:
    sam = temp(DIR_OUTPUT_MAP + "/bismark/{sample}_{species}.sam"),
    rep = DIR_OUTPUT_MAP + "/bismark/{sample}_{species}.mapstats"
  shell:
    """
    samtools view -h {input.bam} -o {output.sam}temp && \
    fix-bismark-sam -o {output.sam} {output.sam}temp && \
    cp {input.rep} {output.rep} && \
    rm {output.sam}temp
    """

############### BSMAP #####################
rule bsmap_map_single:
  input:
    r = fa_single,
    fa = fa
  output:
    sam = temp(DIR_OUTPUT_MAP + "/bsmap/{sample}_{species}.sam"),
    out = DIR_OUTPUT_MAP + "/bsmap/{sample}_{species}.out",
  benchmark:
    DIR_OUTPUT_MAP + "/bsmap/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    bsmap -g 3 -v 0.1 -r 0 -p 16 -V 2 -a {input.r} \
    -d {input.fa} -o {output.sam} 2>{output.out}
    """

rule bsmap_map_paired:
  input:
    r1 = fa_paired_r1,
    r2 = fa_paired_r2,
    fa = fa
  output:
    sam = temp(DIR_OUTPUT_MAP + "/bsmap/{sample}_{species}.sam"),
    out = DIR_OUTPUT_MAP + "/bsmap/{sample}_{species}.out",
  benchmark:
    DIR_OUTPUT_MAP + "/bsmap/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    bsmap -g 3 -v 0.1 -m 32 -x 3000 -r 0 -p 16 -V 2 -a {input.r1} \
    -b {input.r2} -d {input.fa} -o {output.sam} 2>{output.out}
    """

rule bsmap_map_rpbat:
  input:
    r1 = fa_rpbat_r1,
    r2 = fa_rpbat_r2,
    fa = fa
  output:
    sam = temp(DIR_OUTPUT_MAP + "/bsmap/{sample}_{species}.sam"),
    out = DIR_OUTPUT_MAP + "/bsmap/{sample}_{species}.out",
  benchmark:
    DIR_OUTPUT_MAP + "/bsmap/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    bsmap -g 3 -v 0.1 -m 32 -x 3000 -r 0 -n 1 -p 16 -V 2 \
    -a {input.r1} -b {input.r2} -d {input.fa} -o {output.sam} \
    2>{output.out}
    """

############### MINIMAP2 ################
rule minimap2_map_single:
  input:
    fa = fa_single_minimap2,
    ref = index_minimap2
  output:
    temp(DIR_OUTPUT_MAP + "/minimap2/{sample}_{species}.samtemp")
  benchmark:
    DIR_OUTPUT_MAP + "/minimap2/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    minimap2 -x sr  -A 2 -B 2 -O 2,2 -E 2,2 -f 0.000001 -N 1 -t 16 -a {input.ref} {input.fa} >{output}
    """

rule fix_minimap2:
  input:
    DIR_OUTPUT_MAP + "/minimap2/{sample}_{species}.samtemp"
  output:
    temp(DIR_OUTPUT_MAP + "/minimap2/{sample}_{species}.sam")
  shell:
    """
    fix-minimap2-sam -o {output} {input}
    """

############### BWA-METH################
rule bwa_map_single:
  input:
    r = fa_single,
    ref = fa
  output:
    temp(DIR_OUTPUT_MAP + "/bwa/{sample}_{species}.samtemp")
  benchmark:
    DIR_OUTPUT_MAP + "/bwa/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    python /panfs/qcb-panasas/desenabr/software/bwa-meth/bwameth.py --threads 16 --reference \
    {input.ref} {input.r} >{output}
    """

rule bwa_map_paired:
  input:
    r1 = fa_paired_r1,
    r2 = fa_paired_r2,
    ref = fa
  output:
    temp(DIR_OUTPUT_MAP + "/bwa/{sample}_{species}.samtemp")
  benchmark:
    DIR_OUTPUT_MAP + "/bwa/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    python /panfs/qcb-panasas/desenabr/software/bwa-meth/bwameth.py --threads 16 --reference \
    {input.ref} {input.r1} {input.r2} >{output}
    """

rule fix_bwa:
  input:
    DIR_OUTPUT_MAP + "/bwa/{sample}_{species}.samtemp"
  output:
    temp(DIR_OUTPUT_MAP + "/bwa/{sample}_{species}.sam")
  shell:
    """
    fix-bwa-sam -o {output} {input}
    """

############### WALT #####################
rule walt_map_single:
  input:
    r = fa_single,
    index = index_walt
  output:
    mr = temp(DIR_OUTPUT_MAP + "/walt/{sample}_{species}.mr"),
    mapstats = DIR_OUTPUT_MAP + "/walt/{sample}_{species}.mr.mapstats"
  benchmark:
    DIR_OUTPUT_MAP + "/walt/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    walt -N 1000000 -m 20 -t 16 -i {input.index} -r {input.r} -o {output.mr}
    """

rule walt_map_paired:
  input:
    r1 = fa_paired_r1,
    r2 = fa_paired_r2,
    index = index_walt
  output:
    mr = temp(DIR_OUTPUT_MAP + "/walt/{sample}_{species}.mr"),
    mapstats = DIR_OUTPUT_MAP + "/walt/{sample}_{species}.mr.mapstats"
  benchmark:
    DIR_OUTPUT_MAP + "/walt/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    walt -N 1000000 -k 100 -m 20 -L 3000 -t 16 -i {input.index} -1 {input.r1} \
    -2 {input.r2} -o {output.mr}
    """

# converts all walt outputs to sam
rule walt_to_sam:
  input:
    mr = DIR_OUTPUT_MAP + "/walt/{sample}_{species}.mr",
    fai = fai
  output:
    DIR_OUTPUT_MAP + "/walt/{sample}_{species}.sam",
  shell:
    "mr-to-sam -h {input.fai} -o {output} {input.mr}"

############## HISAT-3N ########################

# hisat-3n --summary-file {output.stats} --index {input.index} -p 16 -I 32 -X 3000 -t -q -1 {input.r1} -2 {input.r2} -S {output.sam} --base-change C,T --no-spliced-alignment
rule hisat_map_single:
  input:
    r = fa_single,
    index = index_hisat
  output:
    sam = DIR_OUTPUT_MAP + "/hisat_3n/{sample}_{species}.sam",
    stats  = DIR_OUTPUT_MAP + "/hisat_3n/{sample}_{species}.mapstats",
  benchmark:
    DIR_OUTPUT_MAP + "/hisat_3n/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    hisat-3n --summary-file {output.stats} --index {input.index} \
    -p 16 -I 32 -X 3000 -t -q -U {input.r} -S {output.sam} --unique-only \
    --ignore-quals --mp 1,1 --sp 1,1 --np 1 --rdg 1,1 --rfg 1,1 --score-min L,0.0,-0.1 \
    --base-change C,T --no-spliced-alignment
    """

rule hisat_map_paired:
  input:
    r1 = fa_paired_r1,
    r2 = fa_paired_r2,
    index = index_hisat
  output:
    sam = DIR_OUTPUT_MAP + "/hisat_3n/{sample}_{species}.sam",
    stats  = DIR_OUTPUT_MAP + "/hisat_3n/{sample}_{species}.mapstats",
  benchmark:
    DIR_OUTPUT_MAP + "/hisat_3n/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    hisat-3n --summary-file {output.stats} --index {input.index} \
    -p 16 -I 32 -X 3000 -t -q -1 {input.r1} -2 {input.r2} -S {output.sam} --unique-only \
    --ignore-quals --mp 1,1 --sp 1,1 --np 1 --rdg 1,1 --rfg 1,1 --score-min L,0.0,-0.1 \
    --base-change C,T --no-spliced-alignment
    """

############## OUTPUT FORMATTING ################
# this perform the multiple formatting tasks to obtain
# a sorted, duplicate-removed single-end sam

rule format_sam:
  input:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.sam"
  output:
    temp(DIR_DUMP + "/{mapper}/{sample}_{species}.samf")
  shell:
    """
    format_reads -f {wildcards.mapper} -s 2 -o {output} {input}
    """

rule sort_formatted_sam:
  input:
    DIR_DUMP + "/{mapper}/{sample}_{species}.samf"
  output:
    temp(DIR_DUMP + "/{mapper}/{sample}_{species}.samfs")
  shell:
    """
    samtools sort -O sam -o {output} {input}
    """

rule duplicate_remove:
  input:
    DIR_DUMP + "/{mapper}/{sample}_{species}.samfs"
  output:
    o = temp(DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samd"),
    s = DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.drstats"
  shell:
    """
    duplicate-remover -D -S {output.s} {input} {output.o}
    """

############## METHPIPE ANALYSIS ################
# we will use cmb-06 to keep .meth files and symmetric .meth
# files because it is not heavy IO but take a lot of space

rule samtools_stats:
  input:
    sam = DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.sam",
    fa = fa
  output:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samstats",
  shell:
    """
    samtools stats -r {input.fa} -i 3000 {input.sam} >{output}
    """

rule samtools_stats_minimap2:
  input:
    sam = DIR_OUTPUT_MAP + "/minimap2/{sample}_{species}.sam",
    fa = fa_minimap2
  output:
    DIR_OUTPUT_MAP + "/minimap2/{sample}_{species}.samstats",
  shell:
    """
    samtools stats -r {input.fa} -i 3000 {input.sam} >{output}
    """
ruleorder: samtools_stats_minimap2 > samtools_stats

rule bsrate:
  input:
    mr = DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samd",
    fa = fa
  output:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.bsrate"
  shell:
    "bsrate -c {input.fa} -o {output} {input.mr}"

rule methcounts:
  input:
    mr = DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samd",
    fa = fa
  output:
    DIR_OUTPUT_REST + "/{mapper}/{sample}_{species}.meth"
  shell:
    "methcounts -c {input.fa} -o {output} {input.mr}"

rule levels:
  input:
    DIR_OUTPUT_REST + "/{mapper}/{sample}_{species}.meth"
  output:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.levels"
  shell:
    "levels -o {output} {input}"
