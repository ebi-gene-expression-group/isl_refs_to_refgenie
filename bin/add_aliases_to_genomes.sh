#!/usr/bin/env bash

# Alias genomes so they can be identified without knowing assemblies

list=$(refgenie list)
species=$(echo -e "$list" | sed 's/^│ //g' | awk -F'--' '{print $1}'| grep '^[a-z]' | sort | uniq)
assemblies_list=$( echo -e "$list" | sed 's/^│ //g' | grep '^[a-z]' | awk -F'│' '{print $1}' | sed 's/ +$//g' | awk -F',' '{print $1}' )

for spec in $species; do
    for spikeset in '' '--spikes_ercc'; do
        echo -e "${spec}-${spikeset}"
        needed_aliases=''
        for alias in newest current; do
            echo -e "$list" | grep "${spec}--${alias}${spikeset}" > /dev/null
            if [ $? -eq 0 ]; then
                echo "$spec already has $alias alias (spikes: $spikeset)"
            else
                echo "$spec needs $alias alias (spikes: $spikeset)"
                needed_aliases="$needed_aliases $alias"
            fi
        done

        if [ -n "$needed_aliases" ]; then
            if [ -z "$spikeset" ]; then
                candidate_genomes=$(echo -e "$assemblies_list" | grep "^$spec" | grep -v 'spikes')
            else
                spikesearch=$(echo -e "$spikeset" | sed 's/\-/\\-/g')
                candidate_genomes=$(echo -e "$assemblies_list" | grep "^$spec" | grep "$spikesearch")
            fi
            n_candidate_genomes=$(echo -e "$candidate_genomes" | wc -l)

            # Determine from the GTF if the build was used for 'newest',
            # 'current' or both

            for can_gen in $candidate_genomes; do
                gtf_line=$(refgenie list -g $can_gen | grep "ensembl_gtf")
                digest=$(refgenie alias get -a ${can_gen})

                for na in $needed_aliases; do
                    echo $gtf_line | grep $na > /dev/null
                    if [ $? -eq 0 ]; then
                        echo "Aliasing $can_gen as ${spec}--${na}${spikeset} using digest $digest"
                        refgenie alias set --aliases ${spec}--${na}${spikeset} --digest $digest
                    fi
                done
            done
            
        fi
    done
done
