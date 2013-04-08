#!/bin/bash

FILE_NAME=flows
OUTPUT=flows_out
sed -e "s/n7C/c7/g" \
    -e "s/n/o/g" \
    -e "s/h1/c /g" \
    -e "s/h2/co/g" \
    -e "s/he2/c2 /g" \
    -e "s/li3/c3 /g" \
    -e "s/be4/c4 /g" \
    -e "s/b5/c5/g" \
    -e "s/c6C/c6/g" \
    -e "s/o8/c8/g" \
    -e "s/oo/on/g" <${FILE_NAME} >${OUTPUT}

rm ${FILE_NAME}
