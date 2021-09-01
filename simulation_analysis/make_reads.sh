source ~/.bashrc;
num_reads=1000000;
for err in 0 1 2 3 4 5
do
  for len in 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190;
  do
    if [ ! -f reads_rpbat/err_${err}_len_${len}/simulated_1.fastq ]
    then
      printf "sherman err = ${err} len = ${len} ";
      lazy_qsub "
      cd /panfs/qcb-panasas/desenabr/ab/simulation_analysis/reads_rpbat/err_${err}_len_${len};
      /panfs/qcb-panasas/desenabr/software/sherman/sherman --non_directional -l ${len} -n  ${num_reads} -pe -I 32 -X 500 -CH 99 -CG 30 -e ${err} --genome_folder /panfs/qcb-panasas/desenabr/ref_genomes/hg38/_tmp;
      " sherman_err_${err}_len_${len} 1 8GB
    fi
  done
done
