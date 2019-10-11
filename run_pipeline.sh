#!/bin/bash

mkdir -p workup
mkdir -p workup/logs
mkdir -p workup/logs/cluster

snakemake \
--snakefile Snakefile \
--use-conda \
-j 32 \
--cluster-config cluster.yaml \
--cluster "sbatch -c {cluster.cpus} \
-t {cluster.time} -N {cluster.nodes} \
--mem {cluster.mem} \
--output {cluster.output} \
--error {cluster.error}"


#-t 24:00:00 --nodes 1 \
#--ntasks 1 --cpus-per-task 10 --mem 100G \
#--output slurms/job-%A_%a.out"

#--config bID="./config.txt" assembly="mm10" type="RNA-DNA" num_tags="5" \