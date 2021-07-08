#!/usr/bin/env bash

usage() { echo "Usage: $0 [ -a <assembly> ] [ -r <recipe> ] [ -f <file type> ] [ -p <file path> ] [ -d <refgenie directory> ] [ -t <comma-separated tags> ] [ -x <tag prefix> ] [ -s <assets> ] [ -l <aliases>] [ -b <force rebuild?> ] [ -m <map only> ]" 1>&2; }

assembly=
recipe=
fileType=
filePath=
refgenieDir=
tags=default
tagPrefix=
assets=
forceRebuild=
aliases=
mapOnly=

while getopts ":a:r:f:p:d:t:x:s:l:b:m:" o; do
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
        l)
            aliases=${OPTARG}
            ;;
        b)
            forceRebuild=${OPTARG}
            ;;
        m)
            mapOnly=" --map"
            ;;
        *)
            usage
            exit 1
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

if [ -n "$refenieDir" ]; then
    export REFGENIE=${refgenieDir}/genome_config.yaml 
elif [ -n "$REFGENIE" ]; then
    refgenieDir=$(dirname $REFGENIE)
else
    echo "Need to set env var REFGENIE or supply Refgenie directory via -d" 1>&2
    exit 1
fi

for tag in $(echo "$tags" | tr -d '\n' | sed 's/,/ /g'); do
    completedFlag="${refgenieDir}/alias/$assembly/$recipe/${tagPrefix}$tag/_refgenie_build/refgenie_completed.flag"

    refgenie seek ${assembly}/${recipe}:${tag} > /dev/null 2>&1

    if [ $? -eq 0 ] && [ -z  "$forceRebuild" ]; then
        echo "Asset ${assembly}/${recipe}:${tag} already present"
        continue
    fi

    if [ "$built" -eq 0 ]; then
        firsttag=$tag        
        refgenieCommand="refgenie --verbosity 5 build $assembly/${recipe}:${tagPrefix}${firsttag} ${filePart}${assetsPart}-c ${refgenieDir}/genome_config.yaml${rebuildPart}${mapOnly}"
    else
        # See https://github.com/refgenie/refgenie/issues/252
        #refgenieCommand="refgenie tag $assembly/${recipe}:${firsttag} --tag $tag -c ${refgenieDir}/genome_config.yaml"
        refgenieCommand="refgenie --verbosity 5 build $assembly/${recipe}:${tagPrefix}${tag} ${filePart}${assetsPart}-c ${refgenieDir}/genome_config.yaml${rebuildPart}${mapOnly}"
    fi            

    echo $refgenieCommand

    # Check for errors in case the return code isn't useful

    eval $refgenieCommand > cmd.out
    statusCode=$?
    grep "Changed status from running to failed" cmd.out > /dev/null
    errorPresent=$? 
    cat cmd.out && rm cmd.out

    if [ $? -ne 0 ] || [ $errorPresent -eq 0 ]; then
	    echo "Refgenie build returned non-zero exit status" 1>&2
	    exit 1
    else
        echo "Refgenie build process successful"

        # If this is a genome build, set any provided aliases

        if [ "$built" -eq 0 ] && [ "$recipe" = 'fasta' ] && [ -n "$aliases" ]; then
            digest=$(refgenie alias get -a $assembly)
            refgenie alias set --aliases $(echo -e "$aliases" | sed 's/,/ /g') --digest $digest
            if [ $? -ne 0 ]; then
                echo "Aliasing $assembly to $aliases failed" 1>&2
                exit 1
            fi
        fi
        built=1
    fi

done

touch .done
