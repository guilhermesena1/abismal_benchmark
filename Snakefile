# tsv with 1st column = fq srr, 2nd column = species to map to
input_dataset_tsv = "metadata/tests.txt"

# map results go to local disk
DIR_OUTPUT_MAP = \
  "/panfs/qcb-panasas/desenabr/ab/results/benchmark-files"

# everything else (that does not depend on speed and heavy IO) to cmb-06
DIR_OUTPUT_REST = "/home/cmb-06/as/desenabr/abismal_devel/tests"

# root directory for reference genome files
genome_dir = "/panfs/qcb-panasas/desenabr/ref_genomes/{species}/"
genome_dir_bismark = \
  "/panfs/qcb-panasas/desenabr/ref_genomes/{wildcards.species}/"

# reference genome fasta files
fa = genome_dir + "genome.fa"
fai = genome_dir + "genome.fa.fai"

############### FASTQ FILES #################
# directory to fastq files
fa_dir = "/scratch2/desenabr"

# pattern of r1 and r2 filenames
fa_paired_r1 = fa_dir + "/paired_q0/{sample}_{species}_1.fq"
fa_paired_r2 = fa_dir + "/paired_q0/{sample}_{species}_2.fq"

# pbat paths
fa_rpbat_r1 = fa_dir + "/rpbat_q0/{sample}_{species}_1.fq"
fa_rpbat_r2 = fa_dir + "/rpbat_q0/{sample}_{species}_2.fq"

# SE sample paths
fa_single_r1 = fa_dir + "/single/{sample}_{species}.fq"

# SE truth
sam_truth = fa_dir + "/truth_sam/{sample}.sam"
sam_truth_walt = fa_dir + "/truth_sam/{sample}.samf"

############ INDEX PATHS ###############
# directory where Bisulfite_Genome lies for bismark
index = genome_dir + "index"
index_bismark= genome_dir_bismark + "index"

# walt index
index_walt = index + "/walt.dbindex"

# abismal index file
index_abismal = index + "/abismal_compressed.idx"
########### END INDEX PATHS #############

# samples to map
HUMAN_SAMPLES = []
MOUSE_SAMPLES = []
ZEBRAFISH_SAMPLES = []
CHICKEN_SAMPLES = []
ARABIDOPSIS_SAMPLES = []

for x in open(input_dataset_tsv):
  y = x.strip().split()
  if(y[7] == "yes"):
    species = y[1]
    if(species == "hg38"):
      HUMAN_SAMPLES.append(y[0])
    if(species == "mm10"):
      MOUSE_SAMPLES.append(y[0])
    if(species == "danre11"):
      ZEBRAFISH_SAMPLES.append(y[0])
    if(species == "galgal6"):
      CHICKEN_SAMPLES.append(y[0])
    if(species == "tair10"):
      ARABIDOPSIS_SAMPLES.append(y[0])

# mappers used in benchmarking
MAPPERS = ["abismal", "bismark", "bsmap"]

# output files
OUT_FILES = ["bsrate", "levels"]


################### END REPEAT CONFIG ####################
#
rule all:
  input:
   # what to keep from methpipe
   expand("{dirr}/{mapper}/{sample}_hg38.{methpipeout}",
           dirr = DIR_OUTPUT_MAP,
           mapper = MAPPERS,
           sample = HUMAN_SAMPLES,
           methpipeout = OUT_FILES),
   expand("{dirr}/{mapper}/{sample}_mm10.{methpipeout}",
           dirr = DIR_OUTPUT_MAP,
           mapper = MAPPERS,
           sample = MOUSE_SAMPLES,
           methpipeout = OUT_FILES),
   expand("{dirr}/{mapper}/{sample}_danre11.{methpipeout}",
           dirr = DIR_OUTPUT_MAP,
           mapper = MAPPERS,
           sample = ZEBRAFISH_SAMPLES,
           methpipeout = OUT_FILES),
   expand("{dirr}/{mapper}/{sample}_galgal6.{methpipeout}",
           dirr = DIR_OUTPUT_MAP,
           mapper = MAPPERS,
           sample = CHICKEN_SAMPLES,
           methpipeout = OUT_FILES),
   expand("{dirr}/{mapper}/{sample}_tair10.{methpipeout}",
           dirr = DIR_OUTPUT_MAP,
           mapper = MAPPERS,
           sample = ARABIDOPSIS_SAMPLES,
           methpipeout = OUT_FILES)

############### BSMAP #####################
rule bsmap_map_paired:
  input:
    r1 = fa_paired_r1,
    r2 = fa_paired_r2,
    fa = fa
  output:
    sam = DIR_OUTPUT_MAP + "/bsmap/{sample}_{species}.sam",
    out = DIR_OUTPUT_MAP + "/bsmap/{sample}_{species}.out",
  benchmark:
    DIR_OUTPUT_MAP + "/bsmap/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    bsmap -v 0.1 -m 32 -x 3000 -f 0 -r 0 -p 16 -V 2 -a {input.r1} \
    -b {input.r2} -d {input.fa} -o {output.sam} 2>{output.out}
    """

rule bsmap_map_pbat:
  input:
    r1 = fa_rpbat_r1,
    r2 = fa_rpbat_r2,
    fa = fa
  output:
    sam = DIR_OUTPUT_MAP + "/bsmap/{sample}_{species}.sam",
    out = DIR_OUTPUT_MAP + "/bsmap/{sample}_{species}.out",
  benchmark:
    DIR_OUTPUT_MAP + "/bsmap/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    bsmap -v 0.1 -m 32 -x 3000 -f 0 -r 0 -n 1 -p 16 -V 2 \
    -a {input.r1} -b {input.r2} -d {input.fa} -o {output.sam} \
    2>{output.out}
    """

############### WALT #####################
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
    walt -m 15 -L 3000 -t 16 -i {input.index} -1 {input.r1} \
    -2 {input.r2} -o {output.mr}
    """

rule walt_map_pbat:
  input:
    r1 = fa_rpbat_r1,
    r2 = fa_rpbat_r2,
    index = index_walt
  output:
    mr = temp(DIR_OUTPUT_MAP + "/walt/{sample}_{species}.mr"),
    mapstats = DIR_OUTPUT_MAP + "/walt/{sample}_{species}.mr.mapstats"
  benchmark:
    DIR_OUTPUT_MAP + "/walt/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    walt -P -m 15 -L 3000 -t 16 -i {input.index} -1 {input.r1} \
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

############### ABISMAL #####################
rule abismal_map_paired:
  input:
    r1 = fa_paired_r1,
    r2 = fa_paired_r2,
    index = index_abismal
  output:
    sam = DIR_OUTPUT_MAP + "/abismal/{sample}_{species}.sam",
    mapstats = DIR_OUTPUT_MAP + "/abismal/{sample}_{species}.sam.mapstats",
  benchmark:
    DIR_OUTPUT_MAP + "/abismal/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    abismal -v -t 16 -i {input.index} -o {output.sam} {input.r1} {input.r2}
    """

rule abismal_map_pbat:
  input:
    r1 = fa_rpbat_r1,
    r2 = fa_rpbat_r2,
    index = index_abismal
  output:
    sam = DIR_OUTPUT_MAP + "/abismal/{sample}_{species}.sam",
    mapstats = DIR_OUTPUT_MAP + "/abismal/{sample}_{species}.sam.mapstats",
  benchmark:
    DIR_OUTPUT_MAP + "/abismal/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    abismal -v -t 16 -R -i {input.index} -o {output.sam} {input.r1} {input.r2}
    """

############### BISMARK #################
rule bismark_map_paired:
  input:
    r1 = fa_paired_r1,
    r2 = fa_paired_r2,
  output:
    bam = DIR_OUTPUT_MAP + \
      "/bismark/{sample}_{species}_1_bismark_bt2_pe.bam",
    rep = temp(DIR_OUTPUT_MAP + \
      "/bismark/{sample}_{species}_1_bismark_bt2_PE_report.txt")
  benchmark:
    DIR_OUTPUT_MAP + "/bismark/snakemake_time_{sample}_{species}.txt"
  shell:
    """
    bismark --temp_dir /scratch2/desenabr/bismark_dump --parallel 8  \
    -o %s/bismark --bowtie2 --icpc -I 32 -X 3000 \
    -1 {input.r1} -2 {input.r2} %s
    """ % (DIR_OUTPUT_MAP, index_bismark)

rule bismark_map_pbat:
  input:
    r1 = fa_rpbat_r1,
    r2 = fa_rpbat_r2,
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
            --non_directional --parallel 8 \
            -o %s/bismark --icpc -I 32 -X 3000 -1 {input.r1} -2 {input.r2} %s
    """ % (DIR_OUTPUT_MAP, index_bismark)

rule bismark_to_sam_pe:
  input:
    bam = DIR_OUTPUT_MAP + \
      "/bismark/{sample}_{species}_1_bismark_bt2_pe.bam",
    rep = DIR_OUTPUT_MAP + "/bismark/{sample}_{species}_1_bismark_bt2_PE_report.txt"
  output:
    sam = DIR_OUTPUT_MAP + "/bismark/{sample}_{species}.sam",
    rep = DIR_OUTPUT_MAP + "/bismark/{sample}_{species}.mapstats",
  shell:
    """
    samtools view -h {input.bam} -o {output.sam}temp &&
    fix-bismark-sam -o {output.sam} {output.sam}temp &&
    cp {input.rep} {output.rep} &&
    rm {output.sam}temp
    """

############## GROUND TRUTH COMPARISON ################
rule sort_for_truth_comp:
  input:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.sam"
  output:
    temp(DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samt")
  shell:
    """
    source ~/.samsort && samsort {input} {output}
    """

rule compare_with_truth:
  input:
    truth = sam_truth,
    inp = DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samt"
  output:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.truth_stats"
  shell:
    """
    compare-sam {input.truth} {input.inp} >{output}
    """

rule compare_with_truth_walt:
  input:
    truth = sam_truth_walt,
    inp = DIR_OUTPUT_MAP + "/walt/{sample}_{species}.samt"
  output:
    DIR_OUTPUT_MAP + "/walt/{sample}_{species}.truth_stats"
  shell:
    """
    compare-sam {input.truth} {input.inp} >{output}
    """
ruleorder: compare_with_truth_walt > compare_with_truth

############## OUTPUT FORMATTING ################
# this perform the multiple formatting tasks to obtain
# a sorted, duplicate-removed single-end sam

rule format_sam:
  input:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.sam"
  output:
    temp(DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samf")
  shell:
    """
    format_reads -f {wildcards.mapper} -s 2 -o {output} {input}
    """

rule sort_formatted_sam:
  input:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samf"
  output:
    temp(DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samfs")
  shell:
    """
    samtools sort -O sam -o {output} {input}
    """

rule duplicate_remove:
  input:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samfs"
  output:
    o = temp(DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samd"),
    s = DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.drstats"
  shell:
    """
    duplicate-remover_sam -D -S {output.s} {input} {output.o}
    """

############## METHPIPE ANALYSIS ################
# we will use cmb-06 to keep .meth files and symmetric .meth
# files because it is not heavy IO but take a lot of space

rule bsrate:
  input:
    mr = DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samd",
    fa = fa
  output:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.bsrate"
  shell:
    "bsrate_sam -c {input.fa} -o {output} {input.mr}"

rule methcounts:
  input:
    mr = DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.samd",
    fa = fa
  output:
    DIR_OUTPUT_REST + "/{mapper}/{sample}_{species}.meth"
  shell:
    "methcounts_sam -c {input.fa} -o {output} {input.mr}"

rule sym_cpg:
  input:
     DIR_OUTPUT_REST + "/{mapper}/{sample}_{species}.meth"
  output:
     DIR_OUTPUT_REST + "/{mapper}/{sample}_{species}.meth.sym"
  shell:
    "symmetric-cpgs -o {output} {input}"

rule levels:
  input:
    DIR_OUTPUT_REST + "/{mapper}/{sample}_{species}.meth"
  output:
    DIR_OUTPUT_MAP + "/{mapper}/{sample}_{species}.levels"
  shell:
    "levels -o {output} {input}"
