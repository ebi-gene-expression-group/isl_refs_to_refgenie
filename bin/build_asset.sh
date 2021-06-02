#!/usr/bin/env bash

assembly=$1
recipe=$2
fileType=$3
filePath=$4
refgenieDir=$5
tags=${6:-'default'}
tagPrefix=${7:-''}
assets=${8:-''}

# Assemblies or tags with '.' cause refgenie errors

if [ -e "$completedFlag" ]; then
    echo "$completedFlag already present, $recipe build cancelled" 1>&2
else
    filePart=''
    assetsPart=''
    if [ -n "$fileType" ]; then
        filePart="--files $fileType=$filePath "
    elif [ -n "$assets" ]; then
        assetsPart="--assets $assets "
    fi

    built=0
    firsttag=

    for tag in $(echo "$tags" | tr -d '\n' | sed 's/,/ /g'); do
        completedFlag="${refgenieDir}/alias/$assembly/$recipe/${tagPrefix}$tag/_refgenie_build/refgenie_completed.flag"

        if [ "$built" -eq 0 ]; then
            firsttag=$tag        
            refgenieCommand="refgenie build $assembly/${recipe}:${tagPrefix}${firsttag} ${filePart}${assetsPart}-c ${refgenieDir}/genome_config.yaml"
        else
            # See https://github.com/refgenie/refgenie/issues/252
            #refgenieCommand="refgenie tag $assembly/${recipe}:${firsttag} --tag $tag -c ${refgenieDir}/genome_config.yaml"
            refgenieCommand="refgenie build $assembly/${recipe}:${tagPrefix}${tag} ${filePart}${assetsPart}-c ${refgenieDir}/genome_config.yaml"
        fi            

        echo $refgenieCommand
        eval $refgenieCommand
       
        if [ -e "$completedFlag" ]; then
            echo "Refgenie build process successful"
            built=1
        else
            echo "Refgenie build process failed, '$completedFlag' not present" 1>&2
            exit 1
        fi
    done

fi

touch .done
