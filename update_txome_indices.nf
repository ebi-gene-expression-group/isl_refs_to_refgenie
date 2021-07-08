#!/usr/bin/env nextflow

process find_references {

    conda "${baseDir}/envs/refgenie.yml"

    output:
        file('reference.csv') into REFERENCES

    """
    env COLUMNS=500 refgenie list | awk -F'│' '{print \$2}' | sed 's/ //g' | awk -F',' '{print \$1}' | grep '^[A-Za-z]' > reference.csv
    """
}

process find_cdnas {
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        val(reference) from REFERENCES.splitText().map{it.trim()}

    output:
        tuple val(reference), stdout into CDNAS

    """
    env COLUMNS=500 refgenie list -g \$(echo $reference | tr -d '\\n') | grep fasta_txome | awk -F'│' '{print \$4}' | sed 's/ //g'
    """
}

CDNAS
    .map{tuple(it[0], tuple(it[1].trim().split(',')))}
    .transpose()
    .map{tuple(it[0], it[1], tuple('kallisto', 'salmon'))}
    .transpose()
    .view()
    .set{
        CDNAS_FOR_INDEXING
    }

process build_index {
 
    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry': 'ignore' }
    maxRetries 3

    input:
        tuple val(reference), val(cdnaName), val(indexer) from CDNAS_FOR_INDEXING

    output:
        tuple val(reference), val(cdnaName), val(indexer) into INDEX_DONE

    """
    indexer_version=\$(cat ${baseDir}/envs/refgenie.yml | grep $indexer | awk -F'=' '{print \$2}')
    cdna_asset="fasta=${reference}/fasta_txome:${cdnaName}"
    tag=${cdnaName}--${indexer}_\${indexer_version} 
    
    build_asset.sh \
        -a ${reference} \
        -r ${indexer}_index \
        -t \$tag \
        -m yes \
        -s \$cdna_asset
    """
}

process reduce {

    cache false

    conda "${baseDir}/envs/refgenie.yml"
    
    memory { 20.GB * task.attempt }
    
    maxForks 1
    
    input:
        val(collected) from INDEX_DONE.collect().map { r -> 'collected' }

    output:
        val('reduced') into REDUCED

    """
    refgenie build --reduce
    """ 
}
