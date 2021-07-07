#!/usr/bin/env nextflow

process find_references {

    conda "${baseDir}/envs/refgenie.yml"

    output:
        file('reference.csv') into REFERENCES

    """
    refgenie list | awk -F'│' '{print $2}' | sed 's/ //g' | awk -F',' '{print $1}' | grep '^[A-Za-z]' > reference.csv
    """
}

process find_cdnas {
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        val(reference) from REFERENCES

    output:
        tuple val(reference), stdout into CDNAS

    """
    refgenie list -g $reference | grep fasta_txome | awk -F'│' '{print $4}' | sed 's/ //g'
    """
}

CDNAS
    .map{tuple(it[0], it[1].split(,))}
    .transpose()
    .view()
