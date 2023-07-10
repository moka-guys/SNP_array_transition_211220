# Overlap Map Bed Generator

This script is used to generate the bed file for the overlap map array project which
contains the positions of all probes to be masked within the CHAS array sofware.

The bed file is used with the 'overlap map' feature in CHAS which hides calls in regions
specified in a BED file. This is required because the array team needs a way to hide
'problematic' regions that have been called numerous times but contain no protein coding
genes - this reduces the number of calls needing interpretation so speeds up the analysis and checking process.


## Protocol

The script takes several inputs:

```bash
usage: __main__.py [-h] -N NUM_PROBES -O OUTDIR -G GENES_AED -P PROBES_BED

Generate the BED file for the overlap map array project which contains the positions of all probes to be masked within the CHAS array sofware

options:
  -h, --help            show this help message and exit
  -N NUM_PROBES, --num_probes NUM_PROBES
                        Number of probes away from the nearest coding region that we do not want
                        to be masked. E.g. if a value of 51 is input, the BED output will contain
                        all probes fom the 52nd probe and higher away from the nearest coding
                        region
  -O OUTDIR, --outdir OUTDIR
                        Directory to output created files to
  -G GENES_AED, --genes_aed GENES_AED
                        Genes.aed file containing a curated list of genes and their coordinates
                        which specifies whether the genes are coding or non-coding. This is used
                        to identify which probes require masking
  -P PROBES_BED, --probes_bed PROBES_BED
                        BED file containing all probe names and coordinates
```


The script performs the following actions:

1. Generates a .genome file from the packaged [genome.fa.fai](data/genome.fa.fai) (this file originates from the 001_Tools project (GRCh38.noalt.tar.gz - file-G9k9f600jy1g2X9j37K5FGQ3))
2. Reads the $GENES_AED file supplied on the command line, and manipulates to produce a BED file containing only coding regions
3. Runs bedtools complement to find the regions in the genome file that are not represented in the coding regions BED file (i.e. the non coding regions outside protein coding genes)
4. Extracts the probes that occur within these non-coding regions outside protein coding genes
5. Removes any probes that are within $NUM_PROBES (command-line supplied) distance of a coding region
6. Cleans up any intermediate files produced during the above steps


## Usage

The script is configured to be run on the command line as a module (i.e. from the top-level directory using the -m flag). It is recommended that a virtual environment is set up.

Dependencies can be installed as follows:

```bash
pip3 install -r requirements.txt
```

The tool can be run as follows (with example of NUM_PROBES = 51):

```bash
 python3 -m overlap_map_bed_generator -N 51 -O $PATH_TO_OUTDIR -G $PATH_TO_GENES_AED -P $PATH_TO_PROBES_BED
```
