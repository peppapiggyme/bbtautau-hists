#!/bin/bash

INPUT_PATH=$1
OUTPUT_TAG=$2

BASE_DIR=${PWD}
cd ${BASE_DIR} && source ${BASE_DIR}/setup.sh

cd ${BASE_DIR}/../build-bbtautau-hists/
./run-test ${INPUT_PATH} ${OUTPUT_TAG}
