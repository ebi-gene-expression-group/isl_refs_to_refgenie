#!/usr/bin/env bash

usage() { echo "Usage: $0 [ -a <assembly> ] [ -r <recipe> ] [ -f <file type> ] [ -p <file path> ] [ -d <refgenie directory> ] [ -t <comma-separated tags> ] [ -x <tag prefix> ] [ -s <assets> ] [ -b <force rebuild?> ]" 1>&2; }

assembly=
recipe=
fileType=
filePath=
refgenieDir=
tags=default
tagPrefix=
assets=
forceRebuild=

while getopts ":a:r:f:p:d:t:x:s:b:" o; do
    case "${o}" in
        a)
            assembly=${OPTARG}
            ;;
        r)
            recipe=${OPTARG}
            ;;
        f)
            fileType=${OPTARG}
            ;;
        p)
            filePath=${OPTARG}
            ;;
        d)
            refgenieDir=${OPTARG}
            ;;
        t)
            tags=${OPTARG}
            ;;
        x)
            tagPrefix=${OPTARG}
            ;;
        s)
            assets=${OPTARG}
            ;;
        b)
            forceRebuild=${OPTARG}
            ;;
        *)
            usage
            exit 0
            ;;
    esac
done
shift $((OPTIND-1))

filePart=''
assetsPart=''
rebuildPart=''
if [ -n "$fileType" ]; then
    filePart="--files $fileType=$filePath "
elif [ -n "$assets" ]; then
    assetsPart="--assets $assets "
fi
if [ -n "$forceRebuild" ]; then
    rebuildPart=' -R '
fi

built=0
firsttag=

for tag in $(echo "$tags" | tr -d '\n' | sed 's/,/ /g'); do
    completedFlag="${refgenieDir}/alias/$assembly/$recipe/${tagPrefix}$tag/_refgenie_build/refgenie_completed.flag"

    refgenie seek ${assembly}/${recipe}:${tag} > /dev/null 2>&1

    if [ $? -eq 0 ]; then
        echo "Asset ${assembly}/${recipe}:${tag} already present"
        continue
    fi

    if [ "$built" -eq 0 ]; then
        firsttag=$tag        
        refgenieCommand="refgenie build $assembly/${recipe}:${tagPrefix}${firsttag} ${filePart}${assetsPart}-c ${refgenieDir}/genome_config.yaml${rebuildPart}"
    else
        # See https://github.com/refgenie/refgenie/issues/252
        #refgenieCommand="refgenie tag $assembly/${recipe}:${firsttag} --tag $tag -c ${refgenieDir}/genome_config.yaml"
        refgenieCommand="refgenie build $assembly/${recipe}:${tagPrefix}${tag} ${filePart}${assetsPart}-c ${refgenieDir}/genome_config.yaml${rebuildPart}"
    fi            

    echo $refgenieCommand
    eval $refgenieCommand
  
    if [ $? -ne 0 ]; then
	echo "Refgenie build returned non-zero exit status" 1>&2
	exit 1
    else
        if [ -e "$completedFlag" ]; then
            echo "Refgenie build process successful"
            built=1
        else
            echo "Refgenie build process failed, '$completedFlag' not present" 1>&2
            exit 1
        fi
    fi

done

touch .done
