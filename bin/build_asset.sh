#!/usr/bin/env bash

assembly=$1
recipe=$2
fileType=$3
filePath=$4
refgenieDir=$5
tag=${6:-'default'}

if [ "$recipe" = 'fasta' ]; then
    successMarker="${refgenieDir}/alias/$assembly/fasta/$tag/${assembly}.chrom.sizes"
elif [ "$recipe" = 'ensembl_gtf' ]; then
    successMarker="${refgenieDir}/alias/$assembly/ensembl_gtf/$tag/${assembly}_ensembl_gene_body.bed"
fi

#if [ -e "$successMarker" ]; then
#    echo "$successMarker already present, $recipe build cancelled" 1>&2
#else
    filePart=''
    if [ -n "$fileType" ]; then
        filePart="--files $fileType=$filePath "
    fi
    refgenieCommand="refgenie build $assembly/${recipe}:${tag} ${filePart}-c ${refgenieDir}/genome_config.yaml -R"
    echo $refgenieCommand
    eval $refgenieCommand
   
    if [ -e "$successMarker" ]; then
        echo "Refgenie build process successful"
    else
        echo "Refgenie build process failed, '$successMarker' not present" 1>&2
        exit 1
    fi
#fi

touch .done
