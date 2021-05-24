#!/usr/bin/env bash

assembly=$1
recipe=$2
fileType=$3
filePath=$4
refgenieDir=$5
tags=${6:-'default'}

# Assemblies or tags with '.' cause refgenie errors

if [ -e "$completedFlag" ]; then
    echo "$completedFlag already present, $recipe build cancelled" 1>&2
else
    filePart=''
    if [ -n "$fileType" ]; then
        filePart="--files $fileType=$filePath "
    fi

    built=0
    firsttag=

    for tag in $(echo "$tags" | tr -d '\n' | sed 's/,/ /g'); do
        completedFlag="${refgenieDir}/alias/$assembly/$recipe/$tag/_refgenie_build/refgenie_completed.flag"

        if [ "$built" -eq 0 ]; then
            firsttag=$tag        
            refgenieCommand="refgenie build $assembly/${recipe}:${firsttag} ${filePart}-c ${refgenieDir}/genome_config.yaml -R"
        else
            refgenieCommand="refgenie tag $assembly/${recipe}:${firsttag} --tag $tag"
        fi            

        echo $refgenieCommand
        eval $refgenieCommand
       
        if [ -e "$completedFlag" ]; then
            echo "Refgenie build process successful"
        else
            echo "Refgenie build process failed, '$completedFlag' not present" 1>&2
            exit 1
        fi
    done

fi

touch .done
