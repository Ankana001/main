## Nextflow Process for Low Frequency Somatic Mutation Detection using GATK Mutect
### required version restrictions
Java = 11 or later
\<a href="https://github.com/broadinstitute/gatk/releases">GATK</a>

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
### conda environment creation and activation
To create environment from yaml file the following script should be run
```bash 
conda env create -f gatk-mutect2_env.yaml
conda activate gatk
```
### the index files for bwa was made previously using the following command and was stored in bwa_index folder
```bash
bwa index hg19.fa
```

### the index files for gatk was made previously using the following command and was stored in gatk_ref_index folder
```bash
gatk CreateSequenceDictionary -R hg19.fa
samtools faidx hg19.fa
```
### Running the pipeline
```bash
nextflow run gatk_mutect2_dsl2.nf
```
after doing any changes if you want to resume the flow from where you left you can add -resume paameter with this command

