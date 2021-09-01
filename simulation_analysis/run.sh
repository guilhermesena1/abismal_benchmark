snakemake --printshellcmds --keep-going --ignore-incomplete \
          --cluster-config cluster.yaml \
          --latency-wait 300 --cluster \
          "sbatch --partition=qcb --time=10:00:00 --mem={cluster.mem} \
          --ntasks={cluster.threads} \
          --ntasks-per-node={cluster.threads} \
          --job-name=sm_sim_{cluster.name} \
          --error=sm_{cluster.name}.error \
          --output=sm_{cluster.name}.output" \
          --jobs 100
