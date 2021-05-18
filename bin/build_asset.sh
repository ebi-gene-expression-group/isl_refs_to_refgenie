#!/usr/bin/env bash

assembly=$1
recipe=$2
fileType=$3
filePath=$4
refgenieDir=$5
tags=${6:-'default'}

completedFlag="${refgenieDir}/alias/$assembly/$recipe/$tag/_refgenie_build/refgenie_completed.flag"

# Assemblies or tags with '.' cause refgenie errors

assembly=$(echo -e "$assebly" | sed 's/\./_/g')
tags=$(echo -e "$tags" | sed 's/\./_/g' | sed 's/,/ /g')

if [ -e "$completedFlag" ]; then
    echo "$completedFlag already present, $recipe build cancelled" 1>&2
else
    filePart=''
    if [ -n "$fileType" ]; then
        filePart="--files $fileType=$filePath "
    fi

    built=0
    firsttag=

    for tag in $tags; do
        if [ "$built" -eq 0 ]; then
            firsttag=$tag        
            refgenieCommand="refgenie build $assembly/${recipe}:${firsttag} ${filePart}-c ${refgenieDir}/genome_config.yaml -R"
            echo $refgenieCommand
            eval $refgenieCommand
           
            if [ -e "$completedFlag" ]; then
                echo "Refgenie build process successful"
            else
                echo "Refgenie build process failed, '$completedFlag' not present" 1>&2
                exit 1
            fi

        else
            refgenie tag $assembly/${recipe}:${firsttag}--tag $tag
        fi
    done

fi

touch .done
