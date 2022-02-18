W2R_EXE=/scratchfs/atlas/bowenzhang/build-bbtautau-hists/run-test
OUTPUT_PATH=/scratchfs/atlas/bowenzhang/bbtautau-hists/output/W2R/MVAOutput-Postfit/
INPUT_PATH=/scratchfs/atlas/bowenzhang/ResolvedStatAna/WSMaker_Draft2/output/

# ${W2R_EXE} ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_lephad_SMggFVBF_MVA_0/workspaces/combined/*root ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_tautau_SMggFVBF_MVA_0/fccs/FitCrossChecks.root ${OUTPUT_PATH}/Systs_lephad_SMggFVBF_MVA_0.root
# ${W2R_EXE} ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_lephad_2HDM_MVA_500/workspaces/combined/*root ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_tautau_2HDM_MVA_500/fccs/FitCrossChecks.root ${OUTPUT_PATH}/Systs_lephad_2HDM_MVA_500.root
# ${W2R_EXE} ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_lephad_2HDM_MVA_1000/workspaces/combined/*root ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_tautau_2HDM_MVA_1000/fccs/FitCrossChecks.root ${OUTPUT_PATH}/Systs_lephad_2HDM_MVA_1000.root

${W2R_EXE} ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_hadhad_SMggFVBF_MVA_0/workspaces/combined/*root ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_tautau_SMggFVBF_MVA_0/fccs/FitCrossChecks.root ${OUTPUT_PATH}/Systs_hadhad_SMggFVBF_MVA_0.root
${W2R_EXE} ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_hadhad_2HDM_MVA_500/workspaces/combined/*root ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_tautau_2HDM_MVA_500/fccs/FitCrossChecks.root ${OUTPUT_PATH}/Systs_hadhad_2HDM_MVA_500.root
${W2R_EXE} ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_hadhad_2HDM_MVA_1000/workspaces/combined/*root ${INPUT_PATH}/All_Comb_v15.rqcd_smooth_HH_13TeV_rqcd_smooth_Systs_tautau_2HDM_MVA_1000/fccs/FitCrossChecks.root ${OUTPUT_PATH}/Systs_hadhad_2HDM_MVA_1000.root
