#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "CompTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;
using BU = BinningUtils;

void test_hadhad_ZtautauMG(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);

    Variables* vs_presel = new Variables();
    double binning[7] = {0, 50, 100, 150, 200, 300, 500};
    vs_presel->add("mBB",                  "m_{BB} [GeV]",                                      20, binning, 6);
    // vs_presel->add("mBB",                  "m_{BB} [GeV]",                                      10);
    vs_presel->add("mMMC",                 "m_{#tau#tau} (MMC) [GeV]",                          24);
    vs_presel->add("mHH",                  "m_{HH} [GeV]",                                      24);
    vs_presel->add("mHHScaled",            "m_{HH} (Scaled) [GeV]",                             24);
    vs_presel->add("dRBB",                 "#Delta R(B,B)",                                     10);
    vs_presel->add("pTBB",                 "p_{T}^{b,b}",                                       20);
    vs_presel->add("dRTauTau",             "#Delta R(#tau_{had},#tau_{had})",                   10);
    vs_presel->add("pTTauTau",             "p_{T}^{#tau,#tau}",                                 20);
    vs_presel->add("dPhiBBTauTau",         "#Delta #phi{bb,#tau#tau}",                          8 );
    vs_presel->add("dRBBTauTau",           "#Delta R{bb,#tau#tau}",                             5 );
    vs_presel->add("pTBalance",            "p_{T}^{b,b} / p_{T}^{#tau,#tau}",                   5 );
    vs_presel->add("MET",                  "E_{T}^{miss} [GeV]",                                10);
    vs_presel->add("nJets",                "Number of Jets",                                    1 );
    vs_presel->add("Yield",                "Event Yield",                                       1 );
    vs_presel->add("metSig",               "MET significance",                                  2 );
    vs_presel->add("metSigPU",             "MET significance (pile-up)",                        2 );
    vs_presel->add("Jet0Pt",               "Leading b-jet p_{T} (before corr.) [GeV]",          20);
    vs_presel->add("Jet1Pt",               "Sub-leading b-jet p_{T} (before corr.) [GeV]",      20);
    vs_presel->add("pTB0",                 "Leading b-jet p_{T} [GeV]",                         20);
    vs_presel->add("pTB1",                 "Sub-leading b-jet p_{T} [GeV]",                     20);
    vs_presel->add("LeadJetPt",            "Leading jet p_{T} [GeV]",                           20);
    vs_presel->add("SubleadJetPt",         "Sub-leading jet p_{T} [GeV]",                       20);
    vs_presel->add("LeadJetEta",           "Leading jet #eta",                                  2 );
    vs_presel->add("SubleadJetEta",        "Sub-leading jet #eta",                              2 );
    vs_presel->add("Jet0Eta",              "Leading b-jet #eta (before corr.)",                 2 );
    vs_presel->add("Jet1Eta",              "Sub-leading b-jet #eta (before corr.)",             2 );
    vs_presel->add("Tau0Pt",               "Leading #tau_{had} p_{T} [GeV]",                    10);
    vs_presel->add("Tau1Pt",               "Sub-leading #tau_{had} p_{T} [GeV]",                10);
    vs_presel->add("Tau0Eta",              "Leading #tau_{had} #eta",                           2 );
    vs_presel->add("Tau1Eta",              "Sub-leading #tau_{had} #eta",                       2 );
    vs_presel->add("Tau0Ntrk",             "Leading #tau_{had} N_{tracks}",                     1 );
    vs_presel->add("Tau1Ntrk",             "Sub-leading #tau_{had} N_{tracks}",                 1 );
    vs_presel->add("dPhiTauTau",           "#Delta #phi (#tau_{had},#tau_{had})",               4 );

    auto binningFromFile = [](const std::string& fn)
    { 
        return BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>(fn));
    };

    vector<double> binning2HDM260   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM260.txt");
    vector<double> binning2HDM300   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM300.txt");
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
    vector<double> binning2HDM1200  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1200.txt");
    vector<double> binning2HDM1400  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1400.txt");
    vector<double> binning2HDM1600  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1600.txt");
    vector<double> binningSMHH      =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_SMBDT.txt");

    Tools::printVector(binning2HDM300);

    Variables* vs_pnn = new Variables();
    vs_pnn->add("PNN260",                  "PNN260",                       100,   &binning2HDM260[0], binning2HDM260.size()-1);
    vs_pnn->add("PNN300",                  "PNN300",                       100,   &binning2HDM300[0], binning2HDM300.size()-1);
    vs_pnn->add("PNN350",                  "PNN350",                       100,   &binning2HDM350[0], binning2HDM350.size()-1);
    vs_pnn->add("PNN400",                  "PNN400",                       100,   &binning2HDM400[0], binning2HDM400.size()-1);
    vs_pnn->add("PNN450",                  "PNN450",                       100,   &binning2HDM450[0], binning2HDM450.size()-1);
    vs_pnn->add("PNN500",                  "PNN500",                       100,   &binning2HDM500[0], binning2HDM500.size()-1);
    vs_pnn->add("PNN550",                  "PNN550",                       100,   &binning2HDM550[0], binning2HDM550.size()-1);
    vs_pnn->add("PNN600",                  "PNN600",                       100,   &binning2HDM600[0], binning2HDM600.size()-1);
    vs_pnn->add("PNN700",                  "PNN700",                       100,   &binning2HDM700[0], binning2HDM700.size()-1);
    vs_pnn->add("PNN800",                  "PNN800",                       100,   &binning2HDM800[0], binning2HDM800.size()-1);
    vs_pnn->add("PNN900",                  "PNN900",                       100,   &binning2HDM900[0], binning2HDM900.size()-1);
    vs_pnn->add("PNN1000",                 "PNN1000",                      100,   &binning2HDM1000[0], binning2HDM1000.size()-1);
    vs_pnn->add("PNN1200",                 "PNN1200",                      100,   &binning2HDM1200[0], binning2HDM1200.size()-1);
    vs_pnn->add("PNN1400",                 "PNN1400",                      100,   &binning2HDM1400[0], binning2HDM1400.size()-1);
    vs_pnn->add("PNN1600",                 "PNN1600",                      100,   &binning2HDM1600[0], binning2HDM1600.size()-1);

    Variables* vs_bdt = new Variables();
    vs_bdt->add("SMBDT",                   "SM BDT",                       100,   &binningSMHH[0], binningSMHH.size()-1);

    Systematics* ss = new Systematics();
    ss->add("ZJETS_GENERATOR", "ZJETS_GENERATOR", eSystematicType::TwoSide);

    CompInfo* info = new CompInfo();
    info->ratio_high = 1.48;
    info->ratio_low = 0.52;
    info->shape_only = true;

    for (VariableInfo* v : *(vs_presel->content()))
    {
        Processes* ps = new Processes();
        ps->add("Zttbb",    "Zhf Nominal",  eProcessType::BKG,  eProcess::ZllHF,    "Zhf Nominal",  kBlue+1);
        ps->add("Zttbc",    "Zhf Nominal",  eProcessType::BKG,  eProcess::ZllHF,    "Zhf Nominal",  kBlue+1);
        ps->add("Zttcc",    "Zhf Nominal",  eProcessType::BKG,  eProcess::ZllHF,    "Zhf Nominal",  kBlue+1);

        ps->add("MadZttbb",   "Zhf MG5+Py8",  eProcessType::BKG,  eProcess::ZllHFMG,  "Zhf MG5+Py8",  kRed+1);
        ps->add("MadZttbc",   "Zhf MG5+Py8",  eProcessType::BKG,  eProcess::ZllHFMG,  "Zhf MG5+Py8",  kRed+1);
        ps->add("MadZttcc",   "Zhf MG5+Py8",  eProcessType::BKG,  eProcess::ZllHFMG,  "Zhf MG5+Py8",  kRed+1);

        Config* c = new Config(b, ps, rs, vs_presel, ss);
        c->load(filename, "Preselection");
        info->parameter = "Presel";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/tmp/bowenzhang/ZtautauMG";
        if (ct->check(c))
        {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::Array);
            ct->makeYield(c);
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
        ps->add("Zttbb",    "Zhf Nominal",  eProcessType::BKG,  eProcess::ZllHF,    "Zhf Nominal",  kBlue+1);
        ps->add("Zttbc",    "Zhf Nominal",  eProcessType::BKG,  eProcess::ZllHF,    "Zhf Nominal",  kBlue+1);
        ps->add("Zttcc",    "Zhf Nominal",  eProcessType::BKG,  eProcess::ZllHF,    "Zhf Nominal",  kBlue+1);

        ps->add("MadZttbb",   "Zhf MG5+Py8",  eProcessType::BKG,  eProcess::ZllHFMG,  "Zhf MG5+Py8",  kRed+1);
        ps->add("MadZttbc",   "Zhf MG5+Py8",  eProcessType::BKG,  eProcess::ZllHFMG,  "Zhf MG5+Py8",  kRed+1);
        ps->add("MadZttcc",   "Zhf MG5+Py8",  eProcessType::BKG,  eProcess::ZllHFMG,  "Zhf MG5+Py8",  kRed+1);

        Config* c = new Config(b, ps, rs, vs_pnn, ss);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "PNN";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/tmp/bowenzhang/ZtautauMG";
        if (ct->check(c))
        {
            ct->manipulate(c);
            // ct->rebin(c, eRebinOption::N_Rebin);
            ct->rebin(c, eRebinOption::Array, "using array and transformed", true);
            ct->makeYield(c);
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
        ps->add("Zttbb",    "Zhf Nominal",  eProcessType::BKG,  eProcess::ZllHF,    "Zhf Nominal",  kBlue+1);
        ps->add("Zttbc",    "Zhf Nominal",  eProcessType::BKG,  eProcess::ZllHF,    "Zhf Nominal",  kBlue+1);
        ps->add("Zttcc",    "Zhf Nominal",  eProcessType::BKG,  eProcess::ZllHF,    "Zhf Nominal",  kBlue+1);

        ps->add("MadZttbb",   "Zhf MG5+Py8",  eProcessType::BKG,  eProcess::ZllHFMG,  "Zhf MG5+Py8",  kRed+1);
        ps->add("MadZttbc",   "Zhf MG5+Py8",  eProcessType::BKG,  eProcess::ZllHFMG,  "Zhf MG5+Py8",  kRed+1);
        ps->add("MadZttcc",   "Zhf MG5+Py8",  eProcessType::BKG,  eProcess::ZllHFMG,  "Zhf MG5+Py8",  kRed+1);

        Config* c = new Config(b, ps, rs, vs_bdt, ss);
        c->load(filename, "BDTScorePreselection");
        info->parameter = "BDT";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/tmp/bowenzhang/ZtautauMG";
        if (ct->check(c))
        {
            ct->manipulate(c);
            // ct->rebin(c, eRebinOption::N_Rebin);
            ct->rebin(c, eRebinOption::Array, "using array and transformed", true);
            ct->makeYield(c);
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