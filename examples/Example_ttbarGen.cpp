#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "CompTool.h"
#include "AutoBinningTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;
using BU = BinningUtils;

void test_hadhad_ttbarSys_Gen(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_SS",     "2 b-tag, 2 loose #tau, SS",        eRegionType::SR);
    // rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);

    Variables* vs_presel = new Variables();
    vs_presel->add("mBB",                  "m_{BB} [GeV]",                                      10);
    vs_presel->add("mMMC",                 "m_{#tau#tau} (MMC) [GeV]",                          24);
    vs_presel->add("mHH",                  "m_{HH} [GeV]",                                      24);
    vs_presel->add("dRBB",                 "#Delta R(B,B)",                                     10);
    vs_presel->add("dRTauTau",             "#Delta R(#tau_{had},#tau_{had})",                   10);
    vs_presel->add("MET",                  "E_{T}^{miss} [GeV]",                                10);
    vs_presel->add("nJets",                "Number of Jets",                                    1 );
    vs_presel->add("Yield",                "Event Yield",                                       1 );

    Variables* vs_pnn = new Variables();    
    auto binningFromFile = [](const std::string& fn)
    { 
        return BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>(fn));
    };

    vector<double> binning2HDM260   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM260.txt");
    vector<double> binning2HDM280   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM280.txt");
    vector<double> binning2HDM300   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM300.txt");
    vector<double> binning2HDM325   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM325.txt");
    vector<double> binning2HDM350   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM350.txt");
    vector<double> binning2HDM400   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM400.txt");
    vector<double> binning2HDM450   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM450.txt");
    vector<double> binning2HDM500   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM500.txt");
    vector<double> binning2HDM550   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM550.txt");
    vector<double> binning2HDM600   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM600.txt");
    vector<double> binning2HDM700   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM700.txt");
    vector<double> binning2HDM800   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM800.txt");
    vector<double> binning2HDM900   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM900.txt");
    vector<double> binning2HDM1000  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1000.txt");
    vector<double> binning2HDM1100  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1100.txt");
    vector<double> binning2HDM1200  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1200.txt");
    vector<double> binning2HDM1400  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1400.txt");
    vector<double> binning2HDM1600  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1600.txt");
    vector<double> binningSMBDT     =   BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_SMBDT.txt"), 1000, true);

    vs_pnn->add("PNN260",       "PNN260",       1,      &binning2HDM260[0],     binning2HDM260.size()-1);
    vs_pnn->add("PNN280",       "PNN280",       1,      &binning2HDM280[0],     binning2HDM280.size()-1);
    vs_pnn->add("PNN300",       "PNN300",       1,      &binning2HDM300[0],     binning2HDM300.size()-1);
    vs_pnn->add("PNN325",       "PNN325",       1,      &binning2HDM325[0],     binning2HDM325.size()-1);
    vs_pnn->add("PNN350",       "PNN350",       1,      &binning2HDM350[0],     binning2HDM350.size()-1);
    vs_pnn->add("PNN400",       "PNN400",       1,      &binning2HDM400[0],     binning2HDM400.size()-1);
    vs_pnn->add("PNN450",       "PNN450",       1,      &binning2HDM450[0],     binning2HDM450.size()-1);
    vs_pnn->add("PNN500",       "PNN500",       1,      &binning2HDM500[0],     binning2HDM500.size()-1);
    vs_pnn->add("PNN550",       "PNN550",       1,      &binning2HDM550[0],     binning2HDM550.size()-1);
    vs_pnn->add("PNN600",       "PNN600",       1,      &binning2HDM600[0],     binning2HDM600.size()-1);
    vs_pnn->add("PNN700",       "PNN700",       1,      &binning2HDM700[0],     binning2HDM700.size()-1);
    vs_pnn->add("PNN800",       "PNN800",       1,      &binning2HDM800[0],     binning2HDM800.size()-1);
    vs_pnn->add("PNN900",       "PNN900",       1,      &binning2HDM900[0],     binning2HDM900.size()-1);
    vs_pnn->add("PNN1000",      "PNN1000",      1,      &binning2HDM1000[0],    binning2HDM1000.size()-1);
    vs_pnn->add("PNN1100",      "PNN1100",      1,      &binning2HDM1100[0],    binning2HDM1100.size()-1);
    vs_pnn->add("PNN1200",      "PNN1200",      1,      &binning2HDM1200[0],    binning2HDM1200.size()-1);
    vs_pnn->add("PNN1400",      "PNN1400",      1,      &binning2HDM1400[0],    binning2HDM1400.size()-1);
    vs_pnn->add("PNN1600",      "PNN1600",      1,      &binning2HDM1600[0],    binning2HDM1600.size()-1);

    Variables* vs_bdt = new Variables();
    vs_bdt->add("SMBDT",        "SM BDT",       1,      &binningSMBDT[0],        binningSMBDT.size()-1);

    CompInfo* info = new CompInfo();
    info->ratio_high = 1.78;
    info->ratio_low = 0.22;
    info->shape_only = true;

    for (VariableInfo* v : *(vs_presel->content()))
    {
        Processes* ps = new Processes();
        ps->add("ttbarAF2",       "true-#tau_{h} t#bar{t} Nominal",  eProcessType::BKG,  eProcess::TTBAR_NOM,      "true-#tau_{h} t#bar{t} Wt Nominal",  kBlue+1);
        ps->add("ttbarAMCP8",     "true-#tau_{h} t#bar{t} aMC@NLO",  eProcessType::BKG,  eProcess::TTBAR_AMC,  "true-#tau_{h} t#bar{t} Wt aMC@NLO",  kRed+1);
        ps->add("ttbarPH7",       "true-#tau_{h} t#bar{t} Herwig7",  eProcessType::BKG,  eProcess::TTBAR_HW7,  "true-#tau_{h} t#bar{t} Wt Herwig7",  kViolet+1);
        ps->add("ttbarPP8RadHi",  "true-#tau_{h} t#bar{t} Hdamp",    eProcessType::BKG,  eProcess::TTBAR_HDAMP,  "true-#tau_{h} t#bar{t} Hdamp",       kGreen+1);

        Config* c = new Config(b, ps, rs, vs_presel);
        c->load(filename, "Preselection");
        info->parameter = "ttbar_Gen_Presel";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/TTbar-SS";
        if (ct->check(c))
        {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::Array);
            ct->makeYield(c, info->parameter);
            ct->paint(c);
            ct->run(c);
        }
        else 
        {
            clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
        }

        delete ps;
        delete ct;
        delete c;
    }

    info->logy = true;

    for (VariableInfo* v : *(vs_pnn->content()))
    {
        Processes* ps = new Processes();
        ps->add("ttbarAF2",       "true-#tau_{h} t#bar{t} Nominal",  eProcessType::BKG,  eProcess::TTBAR_NOM,      "true-#tau_{h} t#bar{t} Wt Nominal",  kBlue+1);
        ps->add("ttbarAMCP8",     "true-#tau_{h} t#bar{t} aMC@NLO",  eProcessType::BKG,  eProcess::TTBAR_AMC,  "true-#tau_{h} t#bar{t} Wt aMC@NLO",  kRed+1);
        ps->add("ttbarPH7",       "true-#tau_{h} t#bar{t} Herwig7",  eProcessType::BKG,  eProcess::TTBAR_HW7,  "true-#tau_{h} t#bar{t} Wt Herwig7",  kViolet+1);
        ps->add("ttbarPP8RadHi",  "true-#tau_{h} t#bar{t} Hdamp",    eProcessType::BKG,  eProcess::TTBAR_HDAMP,  "true-#tau_{h} t#bar{t} Hdamp",       kGreen+1);

        Config* c = new Config(b, ps, rs, vs_pnn);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "ttbar_Gen_Presel";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/TTbar-SS";
        if (ct->check(c))
        {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::Array, "using array and transformed", true);
            ct->makeYield(c, info->parameter);
            ct->paint(c);
            ct->run(c);
        }
        else 
        {
            clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
        }

        delete ps;
        delete ct;
        delete c;
    }

        for (VariableInfo* v : *(vs_bdt->content()))
    {
        Processes* ps = new Processes();
        ps->add("ttbarAF2",       "true-#tau_{h} t#bar{t} Nominal",  eProcessType::BKG,  eProcess::TTBAR_NOM,      "true-#tau_{h} t#bar{t} Wt Nominal",  kBlue+1);
        ps->add("ttbarAMCP8",     "true-#tau_{h} t#bar{t} aMC@NLO",  eProcessType::BKG,  eProcess::TTBAR_AMC,  "true-#tau_{h} t#bar{t} Wt aMC@NLO",  kRed+1);
        ps->add("ttbarPH7",       "true-#tau_{h} t#bar{t} Herwig7",  eProcessType::BKG,  eProcess::TTBAR_HW7,  "true-#tau_{h} t#bar{t} Wt Herwig7",  kViolet+1);
        ps->add("ttbarPP8RadHi",  "true-#tau_{h} t#bar{t} Hdamp",    eProcessType::BKG,  eProcess::TTBAR_HDAMP,  "true-#tau_{h} t#bar{t} Hdamp",       kGreen+1);

        Config* c = new Config(b, ps, rs, vs_bdt);
        c->load(filename, "BDTScorePreselection");
        info->parameter = "ttbar_Gen_Presel";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/TTbar-SS";
        if (ct->check(c))
        {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::Array, "using array and transformed", true);
            ct->makeYield(c, info->parameter);
            ct->paint(c);
            ct->run(c);
        }
        else 
        {
            clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
        }

        delete ps;
        delete ct;
        delete c;
    }

    delete b;
    delete rs;
    delete vs_presel;
    delete vs_pnn;
    delete vs_bdt;
    delete info;

}
