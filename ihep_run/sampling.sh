# INVER=All_Comb_v15
# OUTVER=rqcd_smooth
# MASS=1600

MYSCRATHFS=/scratchfs/atlas/bowenzhang/
OUTPUTPATH=/scratchfs/atlas/bowenzhang/ResolvedStatAna/WSMaker_Draft2/output/

# WORKSPACEPATH=${OUTPUTPATH}/${INVER}.${OUTVER}_HH_13TeV_${OUTVER}_Systs_hadhad_2HDM_MVA_${MASS}/workspaces/combined/${MASS}.root
# OUTPUT_TAG=${MASS}.mu0p2.ntoy10k.try2

# echo $WORKSPACEPATH
# echo $OUTPUT_TAG

# hep_sub jobs/submit_sampling.sh \
#     -g atlas -os CentOS7 \
#     -mem 8000 -wt long \
#     -o ${MYSCRATHFS}/bbtautau-hists/output/job.sampling.${INVER}-${OUTVER}-${OUTPUT_TAG}.out \
#     -e ${MYSCRATHFS}/bbtautau-hists/output/job.sampling.${INVER}-${OUTVER}-${OUTPUT_TAG}.err \
#     -argu ${WORKSPACEPATH} ${OUTPUT_TAG}

for MASS in X1200 X1400 X2000 X3000
do
    hep_sub jobs/submit_sampling.sh \
        -g atlas -os CentOS7 \
        -mem 16000 -wt mid \
        -o ${MYSCRATHFS}/bbtautau-hists/output/job.sampling.Boosted-${MASS}.out \
        -e ${MYSCRATHFS}/bbtautau-hists/output/job.sampling.Boosted-${MASS}.err \
        -argu /publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/Boosted/BoostedStatAna/run_full_origin_xx/${MASS}/results/bbtautau${MASS}/SPlusB_combined_NormalMeasurement_model.root Boosted${MASS}.100k.mu4
done
