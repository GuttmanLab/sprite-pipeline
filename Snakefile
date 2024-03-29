'''
Author: Peter Chovanec
Aim: A Snakemake workflow to process DNA-DNA SPRITE-seq data
'''

import os 
import sys
import datetime
from pathlib import Path
from shutil import which
import yaml

#Location of scripts
barcode_id_jar = "scripts/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "scripts/python/get_ligation_efficiency.py"
split_fq = "scripts/python/get_full_barcodes.py"
add_chr = "scripts/python/ensembl2ucsc.py"
get_clusters = "scripts/python/get_clusters.py"
get_cluster_size = "scripts/r/get_cluster_size_distribution.r"
hicorrector = "scripts/HiCorrector_1.2/bin/ic"
clusters_heatmap = "scripts/python/get_sprite_contacts.py"


# Check if file is executable
assert which(hicorrector) is not None, 'HiCorrector ic needs to be executable'

#Copy config file into logs
v = datetime.datetime.now()
run_date = v.strftime('%Y.%m.%d.')

#Load config.yaml file
INFO = ''

try:
    email = config['email']
except:
    INFO += "Won't send email on error\n"
    email = None

try:
    out_dir = config['output_dir']
    INFO += 'All data will be written to: ' + out_dir + '\n'
except:
    out_dir = ''
    INFO += 'Defaulting to working directory as output directory\n'

try:
    bid_config = config['bID']
    INFO += 'Using BarcodeID config ' + bid_config + '\n'
except:
    bid_config = 'workup/config.txt'
    INFO += 'Config "bID" not specified, looking for config at: ' + bid_config + '\n'

try:
    num_tags = config['num_tags']
    INFO += 'Using ' + num_tags + ' tags \n'
except:
    num_tags = "5"
    INFO += 'Config "num_tags" not specified, using:' + num_tags + '\n'

#Make pipeline compatible for multiple assemblies
try:
    assembly = config['assembly']
    assert assembly in ['mm10', 'hg38'], 'Only "mm10" or "hg38" currently supported'
    INFO += 'Using ' + assembly + '\n'
except:
    INFO += 'Config "assembly" not specified, defaulting to "mm10"\n'
    assembly = 'mm10'

try:
    bowtie2_index = config['bowtie2_index'][config['assembly']]
except:
    sys.exit('Bowtie2 index not specified in config.yaml') #no default, exit


try:
    mask = config['mask'][config['assembly']]
except:
    sys.exit('Mask path not specified in config.yaml') #no default, exit

try:
    chromosome = config['chromosome']
except:
    chromosome = 'genome'
    INFO += 'Defaulting to "genome"\n'
try:
    downweighting = config['downweighting']
except:
    downweighting = 'two_over_n'
    INFO += 'Defaulting to "two_over_n"\n'

try:
    ice_iterations = config['ice_iterations']
    resolution = config['resolution']
    min_cluster_size = config['min_cluster_size']
    max_cluster_size = config['max_cluster_size']
    max_val = config['max_value']
except:
    ice_iterations = 100
    resolution = 1000000
    min_cluster_size = 2
    max_cluster_size = 1000
    max_val = 255

try:
    samples = config['samples']
    # print('Using samples file:', samples)
except:
    samples = './samples.json'
    # print('Defaulting to working directory for samples json file')

#make output directories (aren't created automatically on cluster)
Path(out_dir + "workup/logs/cluster").mkdir(parents=True, exist_ok=True)
out_created = os.path.exists(out_dir + "workup/logs/cluster")
INFO += 'Output logs path created: ' + str(out_created) + '\n'

#get all samples from fastq Directory using the fastq2json.py scripts, then just
#load the json file with the samples
FILES = json.load(open(samples))
ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

CONFIG = [out_dir + "workup/logs/config_" + run_date + "yaml"]

#Shared
TRIM = expand(out_dir + "workup/trimmed/{sample}_{read}.fq.gz", sample = ALL_SAMPLES, 
              read = ["R1_val_1", "R2_val_2"])
TRIM_LOG = expand(out_dir + "workup/trimmed/{sample}_{read}.fastq.gz_trimming_report.txt", 
                  sample = ALL_SAMPLES, read = ["R1", "R2"])
LE_LOG_ALL = [out_dir + "workup/ligation_efficiency.txt"]
TRIM_RD = expand(out_dir + "workup/trimmed/{sample}_R1.barcoded.RDtrim.fastq.gz", 
                  sample = ALL_SAMPLES)
MASKED = expand(out_dir + "workup/alignments/{sample}.DNA.chr.masked.bam", sample=ALL_SAMPLES)
MULTI_QC = [out_dir + "workup/qc/multiqc_report.html"]

CHR_DNA = expand(out_dir + "workup/alignments/{sample}.DNA.chr.bam", sample=ALL_SAMPLES)
BARCODE_FULL = expand(out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz", sample=ALL_SAMPLES)

Bt2_DNA_ALIGN = expand(out_dir + "workup/alignments/{sample}.DNA.bowtie2.mapq20.bam", 
                       sample=ALL_SAMPLES)
#DNA-DNA
BARCODEID = expand(out_dir + "workup/fastqs/{sample}_{read}.barcoded.fastq.gz", 
                       sample = ALL_SAMPLES, read = ["R1", "R2"])
CLUSTERS = expand(out_dir + "workup/clusters/{sample}.DNA.clusters", sample=ALL_SAMPLES)
CLUSTERS_PLOT = [out_dir + "workup/clusters/cluster_sizes.pdf", out_dir + "workup/clusters/cluster_sizes.png"]
CLUSTERS_HP = expand([out_dir + "workup/heatmap/{sample}.DNA.iced.txt",
        out_dir + "workup/heatmap/{sample}.DNA.bias.txt",
        out_dir + "workup/heatmap/{sample}.DNA.raw.txt",
        out_dir + "workup/heatmap/{sample}.DNA.final.txt"], sample=ALL_SAMPLES)

PLOT = expand([out_dir + "workup/heatmap/{sample}.DNA.final.pdf",
        out_dir + "workup/heatmap/{sample}.DNA.final.png"], sample=ALL_SAMPLES)


################################################################################
# Execute before workflow starts
################################################################################
onstart:
    print(INFO)

################################################################################
# Rule all
################################################################################

rule all:
    input: CONFIG + ALL_FASTQ + TRIM + TRIM_LOG + BARCODE_FULL + BARCODEID + LE_LOG_ALL +
           TRIM_RD + Bt2_DNA_ALIGN + CHR_DNA + MASKED + CLUSTERS + MULTI_QC + CLUSTERS_PLOT +
           CLUSTERS_HP + PLOT

################################################################################
# Log
################################################################################

rule log_config:
    '''Copy config.yaml and place in logs folder with the date run
    '''
    output:
        out_dir + "workup/logs/config_" + run_date + "yaml"
    run:
        with open(output[0], 'w') as out:
            yaml.dump(config, out, default_flow_style=False)


####################################################################################################
#Trimming and barcode identification
####################################################################################################


rule adaptor_trimming_pe:
    '''
    #Trim adaptors
    #multiple cores requires pigz to be installed on the system
    '''
    input:
        [lambda wildcards: FILES[wildcards.sample]['R1'],
        lambda wildcards: FILES[wildcards.sample]['R2']]
    output:
         out_dir + "workup/trimmed/{sample}_R1_val_1.fq.gz",
         out_dir + "workup/trimmed/{sample}_R1.fastq.gz_trimming_report.txt",
         out_dir + "workup/trimmed/{sample}_R2_val_2.fq.gz",
         out_dir + "workup/trimmed/{sample}_R2.fastq.gz_trimming_report.txt"
    log:
        out_dir + "workup/logs/{sample}.trim_galore.logs"
    conda:
        "envs/sprite.yaml"
    shell:
        "trim_galore \
        --paired \
        --gzip \
        --quality 20 \
        --fastqc \
        -o {out_dir}workup/trimmed/ \
        {input} &> {log}"


#Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = out_dir + "workup/trimmed/{sample}_R1_val_1.fq.gz",
        r2 = out_dir + "workup/trimmed/{sample}_R2_val_2.fq.gz"
    output:
    #if statements have to be inline (each input is like a function)
        r1_barcoded = out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz",
        r2_barcoded = out_dir + "workup/fastqs/{sample}_R2.barcoded.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.bID.log"
    shell:
        "java -jar {barcode_id_jar} \
        --input1 {input.r1} --input2 {input.r2} \
        --output1 {output.r1_barcoded} --output2 {output.r2_barcoded} \
        --config {bid_config} &> {log}"


#Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        temp(out_dir + "workup/{sample}.ligation_efficiency.txt")
    conda:
        "envs/sprite.yaml"
    shell:
        "python {lig_eff} {input.r1} > {output}"


#Combine ligation efficiency from all samples into a single file
rule cat_ligation_efficiency:
    input:
        expand(out_dir + "workup/{sample}.ligation_efficiency.txt", sample=ALL_SAMPLES)
    output:
        out_dir + "workup/ligation_efficiency.txt"
    shell:
        "tail -n +1 {input} > {output}"



rule full_barcode:
    '''
    remove incomplete barcodes
    '''
    input:
        out_dir + "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.barcoded_short.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}_DPM.log"
    shell:
        "python {split_fq} --r1 {input} &> {log}"


rule cutadapt:
    '''
    Trim DPM if read through reads
    DPM from right GATCGGAAGAG
    DPM from left GGTGGTCTT ^ anchored (only appears at the start of read)
    '''
    input:
        out_dir + "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz"
    output:
        fastq=out_dir + "workup/trimmed/{sample}_R1.barcoded.RDtrim.fastq.gz",
        qc=out_dir + "workup/trimmed/{sample}_R1.barcoded.RDtrim.qc.txt"
    params:
        "-a GATCGGAAGAG -g file:dpm96.fasta"
    log:
        "logs/cutadapt/{sample}.log"
    threads: 
        10
    shell:
        '''
        cutadapt \
        {params} \
        -o {output.fastq} \
        -j {threads} \
        {input} \
        1> {output.qc} 2> {log}
        '''



############################################################################################
#DNA alignment
############################################################################################


rule bowtie2_align:
    '''
    MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
    -F: Do not output alignments with any bits set in INT present in the FLAG field
    '''
    input:
        fq=out_dir + "workup/trimmed/{sample}_R1.barcoded.RDtrim.fastq.gz"
    output:
        out_dir + "workup/alignments/{sample}.DNA.bowtie2.mapq20.bam"
    threads: 10
    log:
        out_dir + "workup/logs/{sample}.bowtie2.log"
    conda:
        "envs/sprite.yaml"
    shell:
        "(bowtie2 \
        -p 10 \
        -t \
        --phred33 \
        -x {bowtie2_index} \
        -U {input.fq} | \
        samtools view -bq 20 -F 4 -F 256 - > {output}) &> {log}"


rule add_chr:
    input:
        out_dir + "workup/alignments/{sample}.DNA.bowtie2.mapq20.bam"
    output:
        out_dir + "workup/alignments/{sample}.DNA.chr.bam"
    log:
        out_dir + "workup/logs/{sample}.DNA_chr.log"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        python {add_chr} -i {input} -o {output} --assembly {assembly} &> {log}
        '''

rule repeat_mask:
    input:
        out_dir + "workup/alignments/{sample}.DNA.chr.bam"
    output:
        out_dir + "workup/alignments/{sample}.DNA.chr.masked.bam"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        bedtools intersect -v -a {input} -b {mask} > {output}
        '''


rule make_clusters:
    input:
        out_dir + "workup/alignments/{sample}.DNA.chr.masked.bam",
    output:
        out_dir + "workup/clusters/{sample}.DNA.clusters"
    log:
        out_dir + "workup/logs/{sample}.make_clusters.log"
    conda:
        "envs/sprite.yaml"
    shell:
        "python {get_clusters} -i {input} -o {output} -n {num_tags} &> {log}"


rule multiqc:
    input:
        #needs to be the last file produced in the pipeline 
        expand(out_dir + "workup/clusters/{sample}.DNA.clusters", sample=ALL_SAMPLES)
    output:
        out_dir + "workup/qc/multiqc_report.html"
    log:
        out_dir + "workup/logs/multiqc.log"
    conda: 
        "envs/sprite.yaml"
    shell: 
        '''
        if [ ! -f {output} ]
        then
            multiqc {out_dir}workup -o {out_dir}workup/qc
        fi
        '''


rule plot_cluster_size:
    input:
        expand(out_dir + "workup/clusters/{sample}.DNA.clusters", sample=ALL_SAMPLES)
    output:
        out_dir + "workup/clusters/cluster_sizes.pdf",
        out_dir + "workup/clusters/cluster_sizes.png"
    log:
        out_dir + "workup/logs/cluster_sizes.log"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        Rscript scripts/r/get_cluster_size_distribution.r \
            {out_dir}workup/clusters/ \
            DNA.clusters
        '''

rule make_heatmap_matrix:
    input:
        out_dir + "workup/clusters/{sample}.DNA.clusters"
    output:
        out_dir + "workup/heatmap/{sample}.DNA.iced.txt",
        out_dir + "workup/heatmap/{sample}.DNA.bias.txt",
        out_dir + "workup/heatmap/{sample}.DNA.raw.txt",
        out_dir + "workup/heatmap/{sample}.DNA.final.txt",
    conda:
        "envs/sprite.yaml"
    log:
        out_dir + "workup/clusters/{sample}.DNA.matrix.log"
    shell:
        '''
        python {clusters_heatmap} \
        --clusters {out_dir}workup/clusters/{wildcards.sample}.DNA.clusters \
        --raw_contacts {out_dir}workup/heatmap/{wildcards.sample}.DNA.raw.txt \
        --biases {out_dir}workup/heatmap/{wildcards.sample}.DNA.bias.txt \
        --iced {out_dir}workup/heatmap/{wildcards.sample}.DNA.iced.txt \
        --output {out_dir}workup/heatmap/{wildcards.sample}.DNA.final.txt \
        --assembly {assembly} \
        --chromosome {chromosome} \
        --downweighting {downweighting} \
        --hicorrector {hicorrector} \
        --resolution {resolution} \
        --iterations {ice_iterations} \
        --max_cluster_size {max_cluster_size} \
        --min_cluster_size {min_cluster_size}
        '''

rule plot_heatmap:
    input:
        out_dir + "workup/heatmap/{sample}.DNA.final.txt"
    output:
        out_dir + "workup/heatmap/{sample}.DNA.final.pdf",
        out_dir + "workup/heatmap/{sample}.DNA.final.png"
    log:
        out_dir + "workup/logs/{sample}.heatmap.log"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        Rscript scripts/r/plot_heatmap.R \
            -i {input} \
            -m {max_val}
        '''


################################################################################
# If error occurs
################################################################################
#Send and email if an error occurs during execution
if email != None:
    onerror:
        shell('mail -s "an error occurred" ' + email + ' < {log}')
