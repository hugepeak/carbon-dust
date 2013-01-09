#!/bin/bash

FILE_NAME=flows
OUTPUT=flows_out
sed "s/n/o/g" <${FILE_NAME} >${OUTPUT}
sed "s/h1/c/g" <${FILE_NAME} >${OUTPUT}
sed "s/h2/co/g" <${FILE_NAME} >${OUTPUT}
sed "s/he2/c2/g" <${FILE_NAME} >${OUTPUT}
sed "s/li3/c3/g" <${FILE_NAME} >${OUTPUT}
sed "s/be4/c4/g" <${FILE_NAME} >${OUTPUT}
sed "s/b5/c5/g" <${FILE_NAME} >${OUTPUT}
sed "s/c6/c6/g" <${FILE_NAME} >${OUTPUT}
sed "s/n7/c7/g" <${FILE_NAME} >${OUTPUT}
sed "s/o8/c8/g" <${FILE_NAME} >${OUTPUT}
sed "s/oo/on/g" <${FILE_NAME} >${OUTPUT}

rm ${FILE_NAME}
