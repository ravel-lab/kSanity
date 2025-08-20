# kSanity

## Overview
kSanity is a k-Mer based approach for detection and quantification of known bacterial strains in shotgun metagenomic data. The expected use-case is for the detection and quantification of bacterial live-biotherapuetic strains in clinical trial data, but the application could be used in any scenario where the goal is to track a known bacterial strain. Required inputs are 1) the closed genome of the strain(s) that are to be tracked, and 2) additional closed genomes of additional strains beloning to the same species (n>10).

First, a set of unique (present in only one strain) and shared (present in all strains) k-Mers are identified in the reference set. Next, each shotgun metagenomic data set is processed indvidually to identify and count the k-Mers observed in the sequence reads. Then detection and quantification is achieved by examining the fraction of unique k-Mers observed in the sequence reads, and quantification by the count of the unique compared to that of the shared k-Mers.

## Dependencies

kSanity requires the following software to be install and in the PATH

-Gzip
-Python3 >v3.8

and the following python packages

-Pandas
-Numpy
-Argparse
-Subprocess
-Json
-Sys
-Biopython
-Collections
-Gzip
-gc

## Installation

VIRGO2 can be installed by cloning the repository.

```
git clone https://github.com/ravel-lab/kSanity.git


```

## Usage

All kSanity operations are performed by the script `kSanity.py` that has the commands `kSanity.py build-ref` , `kSanity.py reads` , `kSanity.py DnQ` , and and `kSanity.py compile`. The `build-ref` command is run once on a dataset of closed genome sequences for conspecific strains and produces two reference json files. The `uniqueKmers.json` file contains a dictionary where the keys are the strain ids and the values are a list of the unique k-Mer sequences. The `genomeSignatures.json` files contains another dictionary where the keys are k-Mer sequences and the values are lists of strains whose genomes contain the k-Mer.

The `reads`command is performed individually for each metagenomic dataset and produces a json file containing k-Mers as keys and the number of times they were observed in the metagenome as values. 

The `DnQ` performs strain detection and quantification using the reference and read json files as input. This command should also be run on each metagenome individually. Strains are detected based on the percentage of their unique k-Mers are observed in a metagenome and are quantified as the ratio of the median number of observations of unique versus shared k-Mers. This operation reports a percent abundance that is iterpretted as the abundance of the strain relative to the species' total population (e.g. Strain 1 makes up 50% of the population of E. coli in the sample).

The `compile` command concatenates the data from multiple analyzed metagenomes.

```
>python3 kSanity.py -h
usage: kSanity.py [-h] {build-ref,reads,DnQ,compile} ...

kSanity is a tool to detect and quantify bacterial strains in shotgun metagenomic data

positional arguments:
  {build-ref,reads,DnQ,compile}

optional arguments:
  -h, --help            show this help message and exit
```

### `kSanity.py build-ref`

This module will build reference k-Mer dictionaries describing the unique and shared sequences present in each conspecific strains. All genomes provided are assumed to be in a single circular piece. Extra-chromomasal elements (e.g. plasmids) should not be provided. The user can provide a value for `k`, but the default is 55bp. Testing demonstrated kSanity is not sensitive to the value of `k`. However `k` should be an odd number less than the length of most of the sequence reads in the metagenomes. 

```
>python3 kSanity.py build-ref -h
usage: kSanity.py build-ref [-h] -r REFDIR -m MANIFEST [-e {fa,fna,fasta}] [-k KMER] -o OUTPUTPREFIX

optional arguments:
  -h, --help            show this help message and exit
  -r REFDIR, --refdir REFDIR
                        Path to directory containing reference genomes
  -m MANIFEST, --manifest MANIFEST
                        File containing list of strains, 1 entry per line
  -e {fa,fna,fasta}, --ending {fa,fna,fasta}
                        File ending of reference genomes default: fa
  -k KMER, --kmer KMER  kMer length used in analysis, odd numbers recommended between 33-77, default:55
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        Prefix used in the filename for the mapping output
```

### `kSanity.py reads`

This module should be performed on each metagenome individually, using the same value of `k` as used to build the reference, and produces a k-Mer dictionary describing the sequences observed in the sample. A single gzip'd fastq is expected, the reads can be of any length and can vary in length across the dataset. Although, all reads shorter than the length of `k` are discarded.

NOTE: If matching values of `k` are not used in the `build-ref` and `reads` steps, kSanity will not function.

```
>python3 kSanity.py reads -h
usage: kSanity.py reads [-h] -r READS [-k KMER] -o OUTPUTPREFIX

optional arguments:
  -h, --help            show this help message and exit
  -r READS, --reads READS
                        Path to the reads file, expected to be a gzipped fastq file
  -k KMER, --kmer KMER  kMer length used in analysis, odd numbers recommended between 33-77, default:55, must match that used in building the reference
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        Prefix used in the filename for the mapping output
```

### `kSanity.py DnQ`

The `DnQ` command detects and quantifies the bacterial strains, and should be run idividually on each metagenome. Provided inputs include reference files and the sample.json file generated by the reads command. Detection is achieved based on the percentage of a strains' unique k-Mers observed in a metagenome and quantification by the number of observations of unique versus shared k-Mers. Two files are generated, one `strainDetect.csv` reports the percentage of unique k-Mers observed for each strain in the reference, and the median number of times their unique k-Mers were observed. The `strainRelAbund.csv` reports the abundance of each strain in the metagenome as a perctange of the total species population (e.g. a value of 0.5 indicates that strain makes of 50% of the species population). 

```
>python3 kSanity.py DnQ -h
usage: kSanity.py DnQ [-h] -r REF -p PROFILE -o OUTPUTPREFIX [-t THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  -r REF, --ref REF     Full path to the reference files generated by the build command with prefix but without ending e.g. /path/to/PREFIX where /path/to/ contains PREFIX.
  -p PROFILE, --profile PROFILE
                        Full path to the sample kMer profile json file produced by the reads command
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        Prefix used in the filename for the mapping output
  -t THRESHOLD, --threshold THRESHOLD
                        Percent of unique kMers needed to be identified for a strain to be considered detected, default=0.5
```

### `kSanity.py compile`

The `compile` command concatenates all of the `strainRelAbund.csv` files in a directory to produce a single output csv file describing the relative abundance of each strain.

```
>python3 kSanity.py compile -h
usage: kSanity.py compile [-h] -i INPUTDIR -o OUTPUTPREFIX

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTDIR, --inputDir INPUTDIR
                        Full path to the directory containing the per sample kSanity output files
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        Prefix used in the filename for the compiled output
```
