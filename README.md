# bbtautau-hists

[![Build Status](https://travis-ci.com/peppapiggyme/bbtautau-hists.svg?branch=master)](https://travis-ci.com/peppapiggyme/bbtautau-hists)

> The pipeline won't work anymore

(not just) a plotting tool for bbtautau-ers

# dependency
c++17, ROOT, CMake

# setup
```bash
git clone git@github.com:peppapiggyme/bbtautau-hists.git
cd ../
mkdir build/
cd build/
cmake ../bbtautau-hists/
make
```

Two useful commands are provided
```bash
bbtt-make
bbtt-output
```

# structrue of source code
```
.
├── examples        # example usage       
├── exec            # main functions
├── ext             # external code
└── src             # core analysis tools
```

# Examples
Example_bbll.cpp -> bbll Zhf CR, study Z+hf modelling
Example_CheckSyst.cpp -> check systematic uncertainties using SystCompTool
Example.cpp -> hadhad SR, basic usage of DrawStackTool
Example_klambda.cpp -> validate klam reweighting, basic usage of CompTool
Example_NonRes.cpp -> compare new/old non-resonant HH sample
Example_PseudoData.cpp -> fit workspaces built by pseudo-data, get a tree of local significances
Example_Pulls.cpp -> fit a workspace, usage of PullsTool
Example_Ranking.cpp -> get ranking of impact of NPs on mu, usage of RankingTool
Example_SamplingDist.cpp -> get toy sampling of a profiled likelihood test statistic
Example_ttbarGen.cpp -> draw ttbar uncertainty, SystCompTool 
Example_ttbarKine.cpp -> draw kinematic comparison of ttbar from different channels
Example_ttbarOTF.cpp -> draw ttbar uncertainty, SystCompTool
Example_ttbarRew.cpp -> draw ttbar uncertainty, SystCompTool
Example_ttbarSyst.cpp -> draw ttbar uncertainty, LepHadTTBarCompTool..
Example_W2R.cpp -> get objects for plotting from workspace, usage of WorkspaceToRoot 
Example_WS.cpp -> basic fit of workspace, usage of WorkSpace class
Example_WtDS.cpp -> draw single top Wt-chan uncertainty, SystCompTool and AutoBinningTool
Example_WtGen.cpp -> draw single top Wt-chan uncertainty
Example_WtOTF.cpp -> draw single top Wt-chan uncertainty
Example_yield.cpp -> an easy way to get full event yield
Example_ZtautauMG.cpp -> draw Z+hf uncertainty, SystCompTool and AutoBinningTool
Example_ZtautauRew.cpp -> extract reweighting of Z+hf, SystCompTool
Example_ZtautauUnc1.cpp -> draw Z+hf uncertainty
Example_ZtautauUnc2.cpp -> draw Z+hf uncertainty
Example_ZtautauUnc3.cpp -> draw Z+hf uncertainty
HadHadSR_v15_binning.cpp -> usage of AutoBinningTool to get optimal MVA output binning
HadHadSR_v15.cpp -> hadhad SR, basic usage of DrawStackTool
HadHadSR_v15_yield.cpp -> an easy way to get full event yield
thesis/* -> for figures in thesis

# run an example
```bash
# in `build/`, do
./run-main
# or
./bbtautau-hists
```

# run a test
```bash
# in `build/`, do
./run-test
```

# run pulls
```bash
# in `build/`, do
./run-pulls <path/to/workspace> <path/to/output/tag>
```

# run ranking
```bash
# in `build/`, do
./run-ranking <path/to/workspace> <path/to/output/tag>
```

# more info and ToDo items
see issues
