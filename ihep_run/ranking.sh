INVER=All_Comb_v4p3
OUTVER=master
MYSCRATHFS=/scratchfs/atlas/bowenzhang/
CHANNEL="lephad"

for i in 260 280 300 325 350 400 450 500 550 600 700 800 900 1000 1100 1200 1400 1600
do
    hep_sub jobs/submit_ranking.sh \
        -g atlas -os CentOS7 \
        -o ${MYSCRATHFS}/bbtautau-hists/output/job.Ranking.${INVER}-${OUTVER}.$i.out \
        -e ${MYSCRATHFS}/bbtautau-hists/output/job.Ranking.${INVER}-${OUTVER}.$i.err \
        -argu \
        ${MYSCRATHFS}/ResolvedStatAna/WSMaker_Group/output/${INVER}.${OUTVER}_HH_13TeV_${OUTVER}_Systs_${CHANNEL}_2HDM_MVA_$i/workspaces/combined/$i.root \
        ${MYSCRATHFS}/bbtautau-hists/output/Ranking.${CHANNEL}.${INVER}-${OUTVER}.$i
        # -mem 10000 \
        # -wt mid
done

# hep_sub jobs/submit_ranking.sh \
#     -g atlas -os CentOS7 \
#     -o ${MYSCRATHFS}/bbtautau-hists/output/job.Ranking.${INVER}-${OUTVER}.SMHH.out \
#     -e ${MYSCRATHFS}/bbtautau-hists/output/job.Ranking.${INVER}-${OUTVER}.SMHH.err \
#     -argu \
#     ${MYSCRATHFS}/ResolvedStatAna/WSMaker_Group/output/${INVER}.${OUTVER}_HH_13TeV_${OUTVER}_Systs_${CHANNEL}_SM_MVA_0/workspaces/combined/0.root \
#     ${MYSCRATHFS}/bbtautau-hists/output/Ranking.${CHANNEL}.${INVER}-${OUTVER}.SMHH \
#     -mem 10000 \
#     -wt mid

# hep_sub jobs/submit_ranking.sh \
#     -g atlas -os CentOS7 \
#     -o ${MYSCRATHFS}/bbtautau-hists/output/job.Huirun.Ranking.out \
#     -e ${MYSCRATHFS}/bbtautau-hists/output/job.Huirun.Ranking.err \
#     -argu \
#     /scratchfs/atlas/bowenzhang/HuirunWS/vbf_cp_m00.root \
#     ${MYSCRATHFS}/bbtautau-hists/output/Huirun.Ranking \
#     -mem 6000 \
#     -wt mid