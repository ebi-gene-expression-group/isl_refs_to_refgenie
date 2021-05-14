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
 6. Bacterial, fungal and viral sequence files, to be used in contamination checks. 

### Building Refgenie assets

From the previous section:

 - 1. is indexed as a top-level reference item, and 3. and 5 added as assets thereof
 - 2., 4. and 6. added as top-level reference items

Further each genome and transcriptome reference is added again with each of the available spike sets (currently just ERCC).

### Genome indices

All genome resources are indexed by hisat2. Contamination indices are also indexed by Bowtie2.

#### Transcriptome indices

Each of the 4 transcriptome assets for each species (current and newest, with/without spikes) is indexed by both Kallisto and Salmon, with versions picked to match our single-cell pipelines.

## Result

Top-level resources are named like `species-assembly[-spikename]', with the spike name being optional. For release-specific resources (principally cDNAs), the assembly component is appended with the release. The resulting content of refgenie is like:

```
> refgenie list 
                                                             Local refgenie assets                                                              
                                              Server subscriptions: http://refgenomes.databio.org                                               
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ genome                                                                                                 ┃ assets                              ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ escherichia_coli-fungi--str_k_12_substr_mg1655_ASM584v2-refseq_fungi_53--ercc                          │ fasta, hisat2_index                 │
│ escherichia_coli-fungi--str_k_12_substr_mg1655_ASM584v2-refseq_fungi_53                                │ fasta, hisat2_index, bowtie2_index  │
│ escherichia_coli-fungi-viruses--str_k_12_substr_mg1655_ASM584v2-refseq_fungi_53-refseq_viral_2_1--ercc │ fasta, hisat2_index                 │
│ escherichia_coli-fungi-viruses--str_k_12_substr_mg1655_ASM584v2-refseq_fungi_53-refseq_viral_2_1       │ fasta, hisat2_index, bowtie2_index  │
│ escherichia_coli--str_k_12_substr_mg1655_ASM584v2--ercc                                                │ fasta, hisat2_index                 │
│ escherichia_coli--str_k_12_substr_mg1655_ASM584v2                                                      │ fasta, hisat2_index                 │
│ escherichia_coli-viruses--str_k_12_substr_mg1655_ASM584v2-refseq_viral_2_1--ercc                       │ fasta, hisat2_index                 │
│ escherichia_coli-viruses--str_k_12_substr_mg1655_ASM584v2-refseq_viral_2_1                             │ fasta, hisat2_index, bowtie2_index  │
│ fungi--refseq_fungi_53--ercc                                                                           │ fasta, hisat2_index                 │
│ fungi--refseq_fungi_53                                                                                 │ fasta, hisat2_index                 │
│ fungi-viruses--refseq_fungi_53-refseq_viral_2_1--ercc                                                  │ fasta, hisat2_index, bowtie2_index  │
│ fungi-viruses--refseq_fungi_53-refseq_viral_2_1                                                        │ fasta, hisat2_index                 │
│ homo_sapiens--GRCh38_cdna_e95--ercc                                                                    │ fasta, salmon_index, kallisto_index │
│ homo_sapiens--GRCh38_cdna_e95                                                                          │ fasta, salmon_index, kallisto_index │
│ homo_sapiens--GRCh38_cdna_e99--ercc                                                                    │ fasta, salmon_index, kallisto_index │
│ homo_sapiens--GRCh38_cdna_e99                                                                          │ fasta, salmon_index, kallisto_index │
│ homo_sapiens--GRCh38--ercc                                                                             │ fasta, ensembl_gtf, hisat2_index    │
│ homo_sapiens--GRCh38                                                                                   │ fasta, ensembl_gtf, hisat2_index    │                                   
│ viruses--refseq_viral_2_1--ercc                                                                        │ fasta, hisat2_index                 │
│ viruses--refseq_viral_2_1                                                                              │ fasta, hisat2_index                 │
└────────────────────────────────────────────────────────────────────────────────────────────────────────┴─────────────────────────────────────┘
                                              use refgenie list -g <genome> for more detailed view  
```
