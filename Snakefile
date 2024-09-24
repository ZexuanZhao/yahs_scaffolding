#!/usr/bin/env python

import os

# Get config file
configfile: "config.yaml"

# Dependencies
os.environ['PATH'] += ':' + os.path.abspath("./scripts")
os.environ['PATH'] += ':' + os.path.abspath("./bin")

## Make all scripts executable

def make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)

for filename in os.listdir("./scripts"):
    f = os.path.join("./scripts", filename)
    # checking if it is a file
    if os.path.isfile(f):
        make_executable(f)

# Get parameters
assembly_prefix = os.path.basename(config["assembly"]).replace(".fasta", "").replace(".fa","")
hic_R1_prefix = os.path.basename(config["hiC_read1"]).replace(".gz", "").replace(".fastq", "").replace(".fq", "")
hic_R2_prefix = os.path.basename(config["hiC_read2"]).replace(".gz", "").replace(".fastq", "").replace(".fq", "")
hic_reads_prefixs = [hic_R1_prefix, hic_R2_prefix]

# Opts
threads = config["threads"]
out_dir = config["outdir"]

# Include rule files
include: "rules/arima_hic_mapping.smk"
include: "rules/yahs.smk"

# Main
rule all:
    input:
        os.path.join(out_dir, "assembly", "assembly_after_scaffolding_scaffolds_final.fa"),
        os.path.join(out_dir, "assembly", "assembly_after_scaffolding_scaffolds_final.hic")

