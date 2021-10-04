#!/usr/bin/env bash

# This script determines if supplied reference files for a species are already
# in use in Refgenie with the same tags etc, so that they can be excluded from
# re-addition

species=$1
assembly=$2
release=$3
reference=$4
cdna_file=$5
gtf_file=$6
tag=$7

function find_orig_refgenie_asset_name() {     
    asset_path=$1     
    basename $(grep "cp " $(dirname $asset_path)/_refgenie_build/refgenie_commands.sh | head -n 1 | awk '{print $2}')
}

function is_asset_new() {
    local is_new=1

    local genome=$1
    local asset_type=$2
    local tag=$3
    local asset=$4

    filename=$(refgenie seek ${genome}/${asset_type}:${tag} 2>/dev/null)
    if [ $? -eq 0 ]; then
        orig_asset_name=$(find_orig_refgenie_asset_name $filename)

        if [ "$(basename $asset)" = "$orig_asset_name" ]; then
            is_new=0
        fi
    fi
    
    return $is_new
}

# We'll be checking status of all our reference files

genome_new=1
gtf_new=0
cdna_new=0

# 1. Check if reference already present and tagged the same

digest=$(refgenie alias get -a ${species}--${tag} 2>/dev/null)

if [ $? -eq 0 ]; then
    current_aliases=$(refgenie alias get | grep $digest | awk -F' │' '{print $2}' | sed 's/^ //' | tr ',' '\n' | sed 's/^ //')
    echo -e "$current_aliases" | grep "^${species}--${tag}" > /dev/null
    if [ $? -eq 0 ]; then
        is_asset_new "${species}--${tag}" 'fasta' 'genome' "$reference"
        genome_new=$?
    fi
fi

# If the genome reference is new we don't need to check anything else

if [ "$genome_new" = 1 ]; then
    echo "REFERENCE: ${species}--${tag}/fasta:genome is not new or diferent to Refgenie"
    exit 1
else
    echo "REFERENCE: ${species}--${tag}/fasta:genome is not new or different to Refgenie"
fi

# 2. Check if cDNA is present and at the relevant tags

for t in cdna_$tag cdna_$release; do
    is_asset_new "${species}--${tag}" 'fasta_txome' "$t" "$cdna_file"
    if [ $? -eq 1 ]; then
        cdna_new=1
        break
    fi    
done

if [ "$cdna_new" = 1 ]; then
    echo "CDNA: One of ${species}--${tag}/fasta_txome entries (cdna_${tag}, cdna_${release}) is new to Refgenie"
    exit 1
else
    echo "CDNA: All of ${species}--${tag}/fasta_txome entries (cdna_${tag}, cdna_${release}) are known to Refgenie"
fi

# 3. Check if GTF is present and at the relevant tags

for t in $tag $release; do
    is_asset_new "${species}--${tag}" 'ensembl_gtf' "$t" "$gtf_file"
    if [ $? -eq 1 ]; then
        gtf_new=1
        break
    fi    
done

if [ "$gtf_new" = 1 ]; then
    echo "GTF: One of ${species}--${tag}/enembl_gtf entries (${tag}, ${release}) is new to Refgenie"
    exit 1
else
    echo "GTF: All of ${species}--${tag}/ensembl_gtf entries (${tag}, ${release}) are known to Refgenie"
fi

exit 0