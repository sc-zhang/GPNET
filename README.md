## Introduction

This software is a toolkit for indentifying alleles associate with phenotypes, then create network for modeling allele
pyramiding with phenotypes.

## Dependencies

### Software

- Python >= 3.7
- GMAP >= 2021
- MAFFT
- rust

**Notice** GMAP and MAFFT are required by MATE module, rust is required by WPNET module for building the
`fast_ridge` library.

### Python Modules

Python modules used by GPNet are listed in `requirements.txt`.

- bioplotz>=0.1.0.dev16
- matplotlib
- networkx
- numpy
- outlier_utils
- pathos>=0.3.0
- pysam
- scikit-learn
- scipy

## Installation

- Download from release

[GPNET v1.1.0](https://github.com/sc-zhang/GPNET/releases/download/v1.1.0/GPNET-v1.1.0.zip)

- Build from source

```bash
cd /path/to/install
git clone https://github.com/sc-zhang/GPNET.git
cd GPNET
chmod +x gpnet.py

# build rust lib
bash ./build.sh

# install python packages
pip install -r requirements.txt
```

## Usage

### Main program

```bash
usage: gpnet.py [-h] [-v] {mate,wpnet,report} ...

positional arguments:
  {mate,wpnet,report}
    mate               MATE module of GPNet
    wpnet              WPNET module of GPNet
    report             Report module of GPNet

options:
  -h, --help           show this help message and exit
  -v, --version        show program's version number and exit
```

#### 1. First stage with MATE module

Details of MATE module could be found in [MATE](mate_module/README.md)

Creating allele matrix with MATE module

- Run with pan-genome data

```bash
gpnet.py mate -q genomes/ -r genes.fa -p pheno/ -o wrkdir/ -t 12
```

> **Notice**
> - **genomes** is a folder which contained all fasta files of samples, like: ABC.fa.
> - **genes.fa** is a fasta file contain all candidate CDS sequences, the id of genes in reference cds file must not
    contain invalid characters that cannot use in path, like '/', '\', '?', et al.
> - **pheno** is a folder which contained all phenotypes files, the phenotype file should be a tsv file named with
    target
> - **phenotype**, and with two columns, [sample, trait value], genome_name should be matched with the name of genome
    file without last suffix, like: ABC.

- Run with assembled CDS data

```bash
gpnet.py mate -q query_cds/ --query_type cds -r genes.fa -p pheno/ -o wrkdir/ -t 12
```

> **Notice**
> - **query_cds** is a folder which contained all cds files of samples, must with suffix .cds, like: ABC.cds.
> - **genes.fa** is a fasta file contain all candidate CDS sequences.
> - **pheno** is a folder which contained all phenotypes files, the phenotype file should be a tsv file named with
    target
> - **phenotype**, and with two columns, [sample_name, trait value], sample_name should be matched with the name of cds
    file without last suffix, like: ABC.

The mat files in wrkdir/0*.VariantMatrix/02.SignificantAlleles could be used for next stage analysis.

#### 2. Second stage with WPNET

Details of WPNET module could be found in [WPNET](wpnet_module/README.md)

##### Creating network

- Run for single phenotype

```bash
mkdir gpnet_wrkdir
gpnet.py wpnet predictor -i wrkdir/08.VariantMatrix/02.SignificantAlleles/Phenotype1.mat -o gpnet_wrkdir/Phenotype1.txt -m CWNET_FR --single -t 12 
```

> **Notice**
> - **08.VariantMatrix** is the output directory with MATE on pan-genome data, for CDS data, it should be
    **07.VariantMatrix**

- Run for all phenotypes

```bash
readlink -f wrkdir/08.VariantMatrix/02.SignificantAlleles/*.mat > mat.list
gpnet.py wpnet predictor -i mat.list -o gpnet_wrkdir -m CWNET_FR -t 12
```

> **Notice**
> - **08.VariantMatrix** is the output directory with MATE on pan-genome data, for CDS data, it should be
    **07.VariantMatrix**

The output directory contains several files with same format of MATE mat file, contain the allele information
and best pyramiding sample, weights of nodes and edges.

##### Selecting top _k_ pyramiding alleles

```bash
gpnet.py wpnet selector -i gpnet_wrkdir/Phenotype1.txt -m gpnet_wrkdir/Phenotype1.gpnet.bin -n 3 -o Phenotype1_select3.txt
```

> **Notice**
> The bin file would only be found with CWNET_FR method (the recommend and default method of WPNET)

The output file is a text file for which each line is a allele id.

#### 3. Final stage for getting result

- Get report of predicted best individual

```bash
gpnet.py report -i gpnet_wrkdir/Phenotype1.txt -m wrkdir/08.VariantMatrix/02.SignificantAlleles/Phenotype1.mat -o report
```

- Get report of selected top _k_ pyramiding alleles

```bash
gpnet.py report -i Phenotype1_select3.txt -m wrkdir/08.VariantMatrix/02.SignificantAlleles/Phenotype1.mat -o report
```

##### Result files

- **allele_source.tsv**: a tsv file contain two or three columns, for predict file, there are 3 columns
  [Allele, AlleleWeight, SourceSamples]; for select file, there are 2 columns [Allele, SourceSamples]. The source
  samples are the samples contain current Allele.
- **sample_score.tsv**: a tsv file contain three or four columns, for predict file, there are 4
  columns [Sample, ContainAlleleCount, Score, Alleles]; for select file, there are three
  columns [Sample, ContainAlleleCount, Alleles]. The ContainAlleleCount means the counts of alleles in best pyramiding
  or be selected could be found in current sample. The score is the sum of weights of alleles in best pyramiding or be
  selected could be found in current sample. The Alleles are the list of best pyramiding or selected alleles in current
  sample.
