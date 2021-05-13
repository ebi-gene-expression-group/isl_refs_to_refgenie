#!/usr/bin/env bash

assembly=$1
recipe=$2
fileType=$3
filePath=$4
refgenieDir=$5
tag=${6:-'default'}

completedFlag="${refgenieDir}/alias/$assembly/$recipe/$tag/_refgenie_build/refgenie_completed.flag"

if [ -e "$completedFlag" ]; then
    echo "$completedFlag already present, $recipe build cancelled" 1>&2
else
    filePart=''
    if [ -n "$fileType" ]; then
        filePart="--files $fileType=$filePath "
    fi
    refgenieCommand="refgenie build $assembly/${recipe}:${tag} ${filePart}-c ${refgenieDir}/genome_config.yaml -R"
    echo $refgenieCommand
    eval $refgenieCommand
   
    if [ -e "$completedFlag" ]; then
        echo "Refgenie build process successful"
    else
        echo "Refgenie build process failed, '$completedFlag' not present" 1>&2
        exit 1
    fi
fi

touch .done
