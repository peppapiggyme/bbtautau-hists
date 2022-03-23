#!/bin/bash

MYSCRATHFS=/scratchfs/atlas/bowenzhang/

# masses="251 260 280 300 325 350 375 400 450 500 550 600 700 800 900 1000 1100 1200 1400 1600"
masses="251 260 280 300 325 350 375 400 450 500 550 600 700 800 900 1000 1100 1200 1400 1600"
# masses="500 800 900 1000 1400 1600"

if [ -z "$1" ]; then
    echo "NO IDX"
    return
fi

INPUT_PATH=/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/Workspaces/bootstrap_output_v2/
OUTPUT_PATH=/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/ntuple/bootstrap_ntuple/v2p5/
TRUNK_IDX=$1

# v2 -> no glob
# v2p1 -> all globs -> strange -25 thing
# v2p2 -> only gamma -> still strange -25 thing
# v2p3 -> only alpha
# v2p4 -> fix fcc fit, all globs

for MASS in $masses
do
    pushd $INPUT_PATH

    mkdir -p ${OUTPUT_PATH}/tmp
    fileTmpListPrefix="${OUTPUT_PATH}/tmp/ListOfPseudoData-${TRUNK_IDX}"

    TARFILE="trunk-${TRUNK_IDX}-${MASS}.tar.gz"

    rm -f ${fileTmpListPrefix}-${MASS}.txt
    touch ${fileTmpListPrefix}-${MASS}.txt
    list=$( tar -tf ${TARFILE} )
    for elem in $list
    do 
        echo "${INPUT_PATH}/$elem" >> ${fileTmpListPrefix}-${MASS}.txt
    done

    split ${fileTmpListPrefix}-${MASS}.txt ${fileTmpListPrefix}-${MASS}.txt. -d -l 50

    ls ${OUTPUT_PATH}/tmp

    popd

    PIECES="00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19"
    # PIECES="00 03 04 05 10 12 15 17 19"
    echo "PIECES are ${PIECES}"
    for PART in $PIECES
    do
        hep_sub jobs/submit_pd.sh \
        -g atlas -os CentOS7 \
        -o ${MYSCRATHFS}/bbtautau-hists/output/logs/job.pd.${TRUNK_IDX}.${MASS}.${PART}.out \
        -e ${MYSCRATHFS}/bbtautau-hists/output/logs/job.pd.${TRUNK_IDX}.${MASS}.${PART}.err \
        -wt mid \
        -argu "${TARFILE}" "${fileTmpListPrefix}-${MASS}.txt.${PART}" "${OUTPUT_PATH}/params_${TRUNK_IDX}_${MASS}.${PART}.root"
    done
done
