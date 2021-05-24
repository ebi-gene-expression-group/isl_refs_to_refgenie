#!/usr/bin/env bash

species=$1
data_dir=$2
gtf_pattern=$3

find_file_pattern=$(basename "$gtf_pattern" | sed 's/RELNO/*/')
extract_release_pattern=$(basename "$gtf_pattern" | sed 's/RELNO/\(\.\*)/')

newest_file=$(ls $data_dir/reference/$species/$find_file_pattern | sort -rV | head -n 1)

if [[ "$(basename $newest_file)" =~ $extract_release_pattern ]]; then
    release=${BASH_REMATCH[1]}
fi

echo -n "$release"

