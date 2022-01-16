whoami

echo ${PWD}

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_97a x86_64-centos7-gcc8-opt" cmake git --quiet

alias bbtt-make="cd ${PWD}/../build-bbtautau-hists; make -j"
alias bbtt-output="cd ${PWD}/output"

mkdir -p ../build-bbtautau-hists
pushd ../build-bbtautau-hists
