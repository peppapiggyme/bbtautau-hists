INVER=Bowen_HadHadWSI_v3p1
OUTVER=bin_fine
MYSCRATHFS=/scratchfs/atlas/bowenzhang/

for i in 260 300 350 400 500 600 700 800 900 1000 1200 1600
do
    hep_sub jobs/submit_pulls.sh \
        -g atlas -os CentOS7 \
        -o ${MYSCRATHFS}/bbtautau-hists/output/job.Pulls.${INVER}-${OUTVER}.$i.out \
        -e ${MYSCRATHFS}/bbtautau-hists/output/job.Pulls.${INVER}-${OUTVER}.$i.err \
        -argu \
        ${MYSCRATHFS}/ResolvedStatAna/WSMaker_Group/output/${INVER}.${OUTVER}_HH_13TeV_${OUTVER}_Systs_hadhad_2HDM_MVA_$i/workspaces/combined/$i.root \
        ${MYSCRATHFS}/bbtautau-hists/output/Pulls.${INVER}-${OUTVER}.$i
done

hep_sub jobs/submit_pulls.sh \
    -g atlas -os CentOS7 \
    -o ${MYSCRATHFS}/bbtautau-hists/output/job.Pulls.${INVER}-${OUTVER}.SMHH.out \
    -e ${MYSCRATHFS}/bbtautau-hists/output/job.Pulls.${INVER}-${OUTVER}.SMHH.err \
    -argu \
    ${MYSCRATHFS}/ResolvedStatAna/WSMaker_Group/output/${INVER}.${OUTVER}_HH_13TeV_${OUTVER}_Systs_hadhad_SM_MVA_0/workspaces/combined/0.root \
    ${MYSCRATHFS}/bbtautau-hists/output/Pulls.${INVER}-${OUTVER}.SMHH