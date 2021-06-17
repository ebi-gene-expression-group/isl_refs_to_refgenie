#!/usr/bin/env bash

configFile=$1
islGenomesFile=$2

config=$(cat $configFile)

species=$(echo -e "$config" | grep "species=" | awk -F'=' '{print $2}' | tr -d '\n')
sourceResource=$(grep "^$species " $islGenomesFile | awk '{print $3}' | tr -d '\n')

if [ "$sourceResource" == 'wbps' ]; then
    version=$(echo -e "$config" | grep -Eo -m 1 "WBPS[0-9]+" | grep -Eo "[0-9]+")
else 
    gtf_file=$(echo -e "$config" | grep gtf_file | awk -F '=' '{print $2}')
    echo -e "$gtf_file" | grep "\.gz" > /dev/null
    if [ $? -eq 0 ]; then
        version=$(echo -e "$gtf_file" | awk -F'.' '{print $(NF-2)}')
    else
        version=$(echo -e "$gtf_file" | awk -F'.' '{print $(NF-1)}')
    fi
fi

if [ -z "$version" ]; then
    echo "Failed to detect current version" 1>&2
    exit 1
fi

echo -n ${sourceResource}${version}
