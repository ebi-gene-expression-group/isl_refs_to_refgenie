#!/usr/bin/env bash

alias=$1
alias_table=$2

digest=$(grep -P " $alias[^\w-]" alias_table.txt | awk -F'│' '{print $2}' | sed 's/ //g')
if [ $? -eq 0 ] && [ -n "$digest" ]; then
    echo -n "$digest"
else
    echo "ERROR: Can't find digest for alias $alias" 1>&2
fi
