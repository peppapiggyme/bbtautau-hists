log="/scratchfs/atlas/bowenzhang/ResolvedStatAna/WSMaker_Group/output/Bowen_HadHadWSI_v3p1.corr_HH_13TeV_corr_Systs_hadhad_2HDM_MVA_${1}/logs/build_workspace.log"

list=$( egrep -o "Bin\s[0-9]+\s[0-9]+" $log | awk '{print $3}' )

echo $list