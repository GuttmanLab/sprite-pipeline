'''
Author: Peter Chovanec
Aim: A Snakemake workflow to process DNA-DNA SPRITE-seq data
'''


import os 
import sys


#Location of scripts
barcode_id_jar = "scripts/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "scripts/python/get_ligation_efficiency.py"
split_fq = "scripts/python/get_full_barcodes.py"
add_chr = "scripts/python/ensembl2ucsc.py"
get_clusters = "scripts/python/get_clusters.py"

#Load config.yaml file

try:
    config_path = config["config_path"]
except:
    config_path = 'config.yaml'

configfile: config_path

try:
    email = config['email']
except:
    print("Won't send email on error")
    email = None

try:
    bid_config = config['bID']
    print('Using BarcodeID config', bid_config)
except:
    bid_config = 'workup/config.txt'
    print('Config "bID" not specified, looking for config at:', bid_config)

try:
    num_tags = config['num_tags']
    print('Using', num_tags, 'tags')
except:
    num_tags = "5"
    print('Config "num_tags" not specified, using:', num_tags)

#Make pipeline compatible for multiple assemblies
try:
    assembly = config['assembly']
    assert assembly in ['mm10', 'hg38'], 'Only "mm10" or "hg38" currently supported'
    print('Using', assembly)
except:
    print('Config "assembly" not specified, defaulting to "mm10"')
    assembly = 'mm10'

try:
    bowtie2_index = config['bowtie2_index'][config['assembly']]
except:
    print('Bowtie2 index not specified in config.yaml')
    sys.exit() #no default, exit


try:
    mask = config['mask'][config['assembly']]
except:
    print('Mask path not specified in config.yaml')
    sys.exit() #no default, exit


#get all samples from fastq Directory using the fastq2json.py scripts, then just
#load the json file with the samples
FILES = json.load(open("./samples.json"))
ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

#Shared
TRIM = expand("workup/trimmed/{sample}_{read}.fq.gz", sample = ALL_SAMPLES, 
              read = ["R1_val_1", "R2_val_2"])
TRIM_LOG = expand("workup/trimmed/{sample}_{read}.fastq.gz_trimming_report.txt", 
                  sample = ALL_SAMPLES, read = ["R1", "R2"])
LE_LOG_ALL = ["workup/ligation_efficiency.txt"]
MASKED = expand("workup/alignments/{sample}.DNA.chr.masked.bam", sample=ALL_SAMPLES)
MULTI_QC = ["workup/qc/multiqc_report.html"]

CHR_DNA = expand("workup/alignments/{sample}.DNA.chr.bam", sample=ALL_SAMPLES)

#Bowtie2 alignment
Bt2_DNA_ALIGN = expand("workup/alignments/{sample}.DNA.bowtie2.mapq20.bam", 
                       sample=ALL_SAMPLES)

#DNA-DNA
BARCODEID_DNA = expand("workup/fastqs/{sample}_{read}.barcoded.fastq.gz", 
                       sample = ALL_SAMPLES, read = ["R1", "R2"])
BCS_DNA = expand("workup/alignments/{sample}.DNA.chr.bam", sample=ALL_SAMPLES)
CLUSTERS_DNA = expand("workup/clusters/{sample}.DNA.clusters", sample=ALL_SAMPLES)



rule all:
    input: ALL_FASTQ + TRIM + TRIM_LOG + BARCODEID_DNA + LE_LOG_ALL + BARCODEID_DNA +
            Bt2_DNA_ALIGN + CHR_DNA + MASKED + CLUSTERS_DNA + MULTI_QC



#Send and email if an error occurs during execution
if email != None:
    onerror:
        shell('mail -s "an error occurred" ' + email + ' < {log}')


####################################################################################################
#Trimming and barcode identification
####################################################################################################

#Trim adaptors
#multiple cores requires pigz to be installed on the system
rule adaptor_trimming_pe:
    input:
        [lambda wildcards: FILES[wildcards.sample]['R1'],
        lambda wildcards: FILES[wildcards.sample]['R2']]
    output:
         "workup/trimmed/{sample}_R1_val_1.fq.gz",
         "workup/trimmed/{sample}_R1.fastq.gz_trimming_report.txt",
         "workup/trimmed/{sample}_R2_val_2.fq.gz",
         "workup/trimmed/{sample}_R2.fastq.gz_trimming_report.txt"
    log:
        "workup/logs/{sample}.trim_galore.logs"
    conda:
        "envs/trim_galore.yaml"
    shell:
        "trim_galore \
        --paired \
        --gzip \
        --quality 20 \
        --fastqc \
        -o workup/trimmed/ \
        {input} &> {log}"



#Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = "workup/trimmed/{sample}_R1_val_1.fq.gz",
        r2 = "workup/trimmed/{sample}_R2_val_2.fq.gz"
    output:
    #if statements have to be inline (each input is like a function)
        r1_barcoded = "workup/fastqs/{sample}_R1.barcoded.fastq.gz",
        r2_barcoded = "workup/fastqs/{sample}_R2.barcoded.fastq.gz"
    log:
        "workup/logs/{sample}.bID.log"
    shell:
        "java -jar {barcode_id_jar} \
        --input1 {input.r1} --input2 {input.r2} \
        --output1 {output.r1_barcoded} --output2 {output.r2_barcoded} \
        --config {bid_config} &> {log}"


#Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        temp("workup/{sample}.ligation_efficiency.txt")
    shell:
        "python {lig_eff} {input.r1} > {output}"


#Combine ligation efficiency from all samples into a single file
rule cat_ligation_efficiency:
    input:
        expand("workup/{sample}.ligation_efficiency.txt", sample=ALL_SAMPLES)
    output:
        "workup/ligation_efficiency.txt"
    shell:
        "tail -n +1 {input} > {output}"



rule full_barcode:
    '''
    remove incomplete barcodes
    '''
    input:
        "workup/fastqs/{sample}_R1.barcoded.fastq.gz"
    output:
        "workup/fastqs/{sample}_R1.barcoded_full.fastq.gz",
        "workup/fastqs/{sample}_R1.barcoded_short.fastq.gz"
    log:
        "workup/logs/{sample}_DPM.log"
    shell:
        "python {split_fq} --r1 {input} &> {log}"

############################################################################################
#DNA alignment
############################################################################################


rule bowtie2_align:
    '''
    MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
    -F: Do not output alignments with any bits set in INT present in the FLAG field
    '''
    input:
        fq="workup/fastqs/{sample}_R1.barcoded_full.fastq.gz"
    output:
        "workup/alignments/{sample}.DNA.bowtie2.mapq20.bam"
    threads: 10
    log:
        "workup/logs/{sample}.bowtie2.log"
    conda:
        "envs/bowtie2.yaml"
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
        "workup/alignments/{sample}.DNA.bowtie2.mapq20.bam"
    output:
        "workup/alignments/{sample}.DNA.chr.bam"
    log:
        "workup/logs/{sample}.DNA_bcs.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        '''
        python {add_chr} -i {input} -o {output} --assembly {assembly} &> {log}
        '''

rule repeat_mask:
    input:
        "workup/alignments/{sample}.DNA.chr.bam"
    output:
        "workup/alignments/{sample}.DNA.chr.masked.bam"
    conda:
        "envs/bedtools.yaml"
    shell:
        '''
        bedtools intersect -v -a {input} -b {mask} > {output}
        '''


rule make_clusters:
    input:
        "workup/alignments/{sample}.DNA.chr.masked.bam",
    output:
        "workup/clusters/{sample}.DNA.clusters"
    log:
        "workup/logs/{sample}.make_clusters.log"
    conda:
        "envs/python_dep.yaml"
    shell:
        "python {get_clusters} -i {input} -o {output} -n {num_tags} &> {log}"



rule multiqc:
    input:
        #needs to be the last file produced in the pipeline 
        expand("workup/clusters/{sample}.DNA.clusters", sample=ALL_SAMPLES)
    output:
        "workup/qc/multiqc_report.html"
    log:
        "workup/logs/multiqc.log"
    conda: 
        "envs/qc.yaml"
    shell: 
        "multiqc workup -o workup/qc"