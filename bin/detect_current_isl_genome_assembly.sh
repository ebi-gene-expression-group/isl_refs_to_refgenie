#!/usr/bin/env bash

configFile=$1
islGenomesFile=$2

config=$(cat $configFile)

# In theory, the species-wise config file might be using a different reference.
# Check for explicit references, assume same as ISL config otherwise.

grep "#assembly" $configFile > /dev/null
if [ $? -eq 0 ]; then
    echo -e "$config" | grep '#assembly' | awk -F'=' '{print $2}' | tr -d '\n'
else
    species=$(echo -e "$config" | grep "species=" | awk -F'=' '{print $2}' | tr -d '\n')
    grep "^$species " $islGenomesFile | awk '{print $7}' | tr -d '\n'
fi 
