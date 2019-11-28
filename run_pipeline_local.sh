#!/bin/bash

mkdir -p workup
mkdir -p workup/logs

snakemake \
--snakefile Snakefile \
--use-conda \
--cores 10
