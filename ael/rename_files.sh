#!/bin/bash

key=$1
txt=$2
rep=$3

echo $key
echo $txt
echo $rep

find . -name $key"*" -exec bash -c 'mv "$1" "${1/'$txt'/'$rep'}"' -- {} \;
