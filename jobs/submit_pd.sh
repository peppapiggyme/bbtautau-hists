#!/bin/bash

INPUT_PATH=/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/bootstrap_output/
OUTPUT_PATH=/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/ntuple/bootstrap_ntuple/
TRUNK_IDX=$1
MASS=$2

pushd $INPUT_PATH

tar -zxf trunk-$TRUNK_IDX-$MASS.tar.gz

mkdir -p /tmp/bowenzhang/
fileTmpListPrefix=/tmp/bowenzhang/ListOfPseudoData-${TRUNK_IDX}

rm -f $fileTmpListPrefix-$MASS.txt
touch $fileTmpListPrefix-$MASS.txt
list=$( tar -tf trunk-${TRUNK_IDX}-$MASS.tar.gz )
for elem in $list
do 
    echo "$INPUT_PATH/$elem" >> $fileTmpListPrefix-$MASS.txt
done

split $fileTmpListPrefix-$MASS.txt $fileTmpListPrefix-$MASS.txt. -d -l 100

ls /tmp/bowenzhang/

popd

BASE_DIR=${PWD}
cd ${BASE_DIR} && source ${BASE_DIR}/setup.sh

cd ${BASE_DIR}/../build-bbtautau-hists/


for PART in 00 01 02 03 04 05 06 07 08 09
do
    echo "======================="
    echo "EXECUTING PART ${PART}"
    echo "======================="

    ./run-pd $fileTmpListPrefix-$MASS.txt.${PART} $OUTPUT_PATH/params_${TRUNK_IDX}_${MASS}.${PART}.root

    echo "======================="
    echo "REMOVING TMP WORKSPACES"
    echo "======================="

    while read -r line
    do
        rm -f $line
    done < $fileTmpListPrefix-$MASS.txt.${PART}

    echo "======================="
    echo "FINISHING PART ${PART}"
    echo "======================="
done

echo "======================="
echo "Finito !"
echo "======================="
