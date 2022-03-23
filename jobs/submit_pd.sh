#!/bin/bash

TARFILE=$1
LIST_OF_INPUT=$2
OUTPUT_ROOTFILE=$3

INPUT_PATH=/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/Workspaces/bootstrap_output_v2/
echo "======================="
echo "UNTAR WITH ${LIST_OF_INPUT}"
echo "======================="

pushd ${INPUT_PATH}
tar -zxf  "${TARFILE}"  $( cat ${LIST_OF_INPUT} | egrep -o "[0-9]+\.[0-9]+\.root" )
popd

BASE_DIR=/scratchfs/atlas/bowenzhang/bbtautau-hists/
cd ${BASE_DIR} && source ${BASE_DIR}/setup.sh

echo "======================="
echo "EXECUTING PART ${PART}"
echo "======================="

./run-pd ${LIST_OF_INPUT} ${OUTPUT_ROOTFILE}

echo "======================="
echo "REMOVING TMP WORKSPACES"
echo "======================="

while read -r line
do
    rm -f $line
done < ${LIST_OF_INPUT}

echo "======================="
echo "Finito !"
echo ${OUTPUT_ROOTFILE}
echo "======================="
