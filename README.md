# Copy legacy ISL genome resources to Refgenie assets

This workflow is designed to copy all of the base genome resources used by the gene expression group to new Refgenie assets. 

## Prerequisites

 - Conda
 - Nextflow

## Usage

```
nextflow run isl_refs_to_refgenie/main.nf \
    --irapConfigDir=<path to directory with species-wise configuration files as used by IRAP> \
    --irapDataDir=<path to top level IRAP data directory, with 'references' subdirectory> \
    --refgenieDir=<path to refgenie top directory> \
    --islGenomes=<path togenome_references.conf>
```

## What's done

### Input file

Input files are discovered for each species:

 1. Current reference genome file from species configs
 2. Current cDNA reference file from species configs
 3. Current GTF annotation file from species configs
 4. Newest cDNA FASTA matching pattern in genome_references.conf
 5. Newest GTF matching pattern in genome references.com

### Building Refgenie assets

From the previous section:

 - 1. is indexed as a top-level reference item, and 3. and 5 added as assets thereof
 - 2. and 4. added as top-level reference items

Further each genome and transcriptome reference is added again with each of the available spike sets (currently just ERCC).

#### Transcriptome indices

Each of the 4 transcriptome assets for each species (current and newest, with/without spikes) is indexed by both Kallisto and Salmon, with versions picked to match our single-cell pipelines.
