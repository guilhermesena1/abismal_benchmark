source ~/.bashrc;
num_reads=1000000;
for err in 0 1 2 3 4 5
do
  for len in 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190;
  do
    the_dir="reads/err_${err}_len_${len}"
    in_r1="$PWD/${the_dir}/simulated_1.fastq";
    in_r2="$PWD/${the_dir}/simulated_2.fastq";
    genome_dir=/panfs/qcb-panasas/desenabr/ref_genomes/hg38
    genome_fa=${genome_dir}/genome.fa
    genome_abismal=~/abdev/hg38.idx
    genome_bismark=${genome_dir}/index
    genome_walt=${genome_dir}/index/walt.dbindex
    for mapper in abismal bismark bsmap walt
    do
      out_sam="$PWD/sam/${mapper}-err-${err}-len-${len}.sam"
      out_walt_temp="$PWD/sam/temp_${mapper}-err-${err}-len-${len}.sam"
      out_acc="$PWD/accuracy-out/${mapper}-err-${err}-len-${len}.txt"
      
      # commands
      abismal_command="~/abdev/src/abismal -v -t 16 -s ${out_sam}.mapstats -i ${genome_abismal} -o ${out_sam} ${in_r1} ${in_r2}"
      bsmap_command="bsmap -g 3 -v 0.1  -m 32 -x 3000 -r 0 -p 16 -V 2 -a ${in_r1} -b ${in_r2} -o ${out_sam} -d ${genome_fa}"
      walt_command="walt -b 200000 -k 100 -m 15 -L 3000 -t 16 -i ${genome_walt} -1 ${in_r1} -2 ${in_r2} -o ${out_sam}"
      bismark_command="cd ${the_dir}; bismark --temp_dir . --parallel 8 -o .  --local --bowtie2 --icpc -I 32 -X 3000 -1 ${in_r1} -2 ${in_r2} ${genome_bismark} && samtools view -h -o temp.sam simulated_1_bismark_bt2_pe.bam && fix-bismark-sam -o ${out_sam} temp.sam && rm temp.sam";
      if [ ! -f ${out_acc} ]
      then
        printf "mapper = ${mapper} err = ${err} len = ${len} ";
        command=""
        mem=""
        if [ ${mapper} == "abismal" ]
        then
          command=${abismal_command}
          mem="8GB"
        fi
        if [ ${mapper} == "bsmap" ]
        then
          command=${bsmap_command}
          mem="16GB";
        fi
        if [ ${mapper} == "walt" ]
        then
          command=${walt_command}
          mem="32GB";
        fi
        if [ ${mapper} == "bismark" ]
        then
          command=${bismark_command}
          mem="32GB"
        fi
        lazy_qsub "${command} && sherman-accuracy ${out_sam} >${out_acc}" ${mapper}-err-${err}-len-${len} 16 ${mem} 24:00:00
      fi
    done
  done
done
