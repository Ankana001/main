# Nextflow Process for Low Frequency Somatic Mutation Detection using GATK Mutect
### required version restrictions
Java = 11 or later
  <a href="https://github.com/broadinstitute/gatk/releases">GATK</a>

## reference sequence used: 
hg19

## directory structure
```bash
.
├── bwa_index
│   ├── hg19.fa
│   ├── hg19.fa.amb
│   ├── hg19.fa.ann
│   ├── hg19.fa.bwt
│   ├── hg19.fa.pac
│   └── hg19.fa.sa
├── fastq_files
│   ├── Accel-Amplicon-56G.tsv
│   ├── CNTRL-0000046-SFT-0000409_S4_MiSeq02-Run0488_L001_R1_001.fastq.gz
│   ├── CNTRL-0000046-SFT-0000409_S4_MiSeq02-Run0488_L001_R2_001.fastq.gz
│   ├── CNTRL-0000056-SFT-0000342_S3_MiSeq01-Run0591_L001_R1_001.fastq.gz
│   └── CNTRL-0000056-SFT-0000342_S3_MiSeq01-Run0591_L001_R2_001.fastq.gz
├── gatk_env.yaml
├── gatk_mutect2_dsl2.nf
└── gatk_ref_index
    ├── hg19.dict
    └── hg19.fa.fai
```
## conda environment creation and activation
To create environment from yaml file the following script should be run
```bash 
conda env create -f gatk-mutect2_env.yaml
conda activate gatk
```
## the index files for bwa was made previously using the following command and was stored in bwa_index folder
```bash
bwa index hg19.fa
```

## the index files for gatk was made previously using the following command and was stored in gatk_ref_index folder
```bash
gatk CreateSequenceDictionary -R hg19.fa
samtools faidx hg19.fa
```
## Running the pipeline
```bash
nextflow run gatk_mutect2_dsl2.nf
```
after doing any changes if you want to resume the flow from where you left you can add -resume paameter with this command

## Stages of the process
### 1. Indexing
This takes already present reference fasta and it's index files as input and creates a tuple channel which is used as reference input for alignment step.
  
### 2. Alignment
This takes the paired end read fastq files and reference fasta channel as input.
The script-
1. alignes the sample reads to the reference fatsa sequence and gets sam file as output
2. the sam file is converted into bam file using samtools

### 3. picard addReadGroups
picard AddOrReplaceReadGroups adds read group tag to the alignment outfile file.
  
#### usage example:
```bash 
java -jar picard.jar AddOrReplaceReadGroups \
      I=input.bam \
      O=output.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
```
This step is necesssary because for mutect2 to run the software needs the @RG tags.
