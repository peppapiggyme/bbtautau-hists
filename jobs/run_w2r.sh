W2R_VERSION=$1
W2R_LEVEL=Postfit

# INPUT_PATH=~/public_store/WSMakerOutput/All_Comb_v${W2R_VERSION}.master_HH_13TeV_master_Systs_hadhad_*
INPUT_PATH=/scratchfs/atlas/bowenzhang/ResolvedStatAna/WSMaker_Draft2/output/All_Comb_v${W2R_VERSION}.master_HH_13TeV_master_Systs_hadhad_*
LIST_OF_INPUTS=$( ls -d ${INPUT_PATH} )

# Comment out other functions!
W2R_EXE=/scratchfs/atlas/bowenzhang/build-bbtautau-hists/run-test

OUTPUT_PATH=/scratchfs/atlas/bowenzhang/bbtautau-hists/output/W2R/MVAOutput-${W2R_LEVEL}/

for ELEM in ${LIST_OF_INPUTS}
do 
    echo "PROCESSING ${ELEM}"
    OUTPUT_FILE_NAME=$( echo "${ELEM}" | egrep -o "Systs(.*)" ).root
    ${W2R_EXE} ${ELEM}/workspaces/combined/*root ${ELEM}/fccs/FitCrossChecks.root ${OUTPUT_PATH}/${OUTPUT_FILE_NAME}
    echo "DONE, STORED TO ${OUTPUT_FLIE_NAME}"
done