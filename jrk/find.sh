#!/bin/bash

value='*2456312.24339*'
obs=$(<'obsfits.txt')
echo $obs |wc -w
for i in "${!obs[@]}"; do
    if [[ "${obs[$i]}" = "${value}" ]]; then 
	echo "${i}"
    fi
done
