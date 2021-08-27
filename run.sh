#!/bin/bash
# run snakemake script on a SLURM machine
snakemake --printshellcmds --keep-going --ignore-incomplete \
          --cluster-config cluster.yaml \
          --latency-wait 60 --cluster \
          "sbatch --partition=qcb --time=200:00:00 --mem={cluster.mem} \
          --ntasks={cluster.threads} \
          --ntasks-per-node={cluster.threads} \
          --job-name=sm_{cluster.name} \
          --error=sm_{cluster.name}.error \
          --output=sm_{cluster.name}.output" \
          --jobs 100
