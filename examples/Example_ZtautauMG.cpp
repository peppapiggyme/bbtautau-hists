#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "CompTool.h"
#include "SystCompTool.h"
#include "AutoBinningTool.h"

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
    //rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);
    rs->add("2tag2pjet_0ptv_LL_SS",     "2 b-tag, 2 loose #tau, SS",        eRegionType::SR);

    Variables* vs_presel = new Variables();
    double binning_mbb[15] = {0, 40, 60, 90, 110, 125, 140, 150, 170, 190, 220, 250, 300, 350, 400};
    double binning_mtt[14] = {60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 130, 150, 300};
    double binning_mhh[20] = {200, 240, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 520, 580, 640, 700, 800, 900, 1200};
    double binning_dR[12] = {0.4, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 2.8, 3.2};
    double binning_pTV[7] = {0., 70., 140., 280., 500., 1000., 2000.};
    double binning_pTxx[14] = {0., 40., 60, 80, 100., 120, 140, 160, 180, 200, 220, 250, 300, 500};

    vs_presel->add("mBB",                  "m_{BB} [GeV]",                                      1);
    vs_presel->add("mMMC",                 "m_{#tau#tau} (MMC) [GeV]",                          1);
    vs_presel->add("mHH",                  "m_{HH} [GeV]",                                      1);
    vs_presel->add("dRBB",                 "#Delta R(B,B)",                                     1);
    vs_presel->add("pTBB",                 "p_{T}^{b,b}",                                       1);
    vs_presel->add("dRTauTau",             "#Delta R(#tau_{had},#tau_{had})",                   1);
    vs_presel->add("pTTauTau",             "p_{T}^{#tau,#tau}",                                 1);

    // vs_presel->add("mBB",                  "m_{BB} [GeV]",                                      10, binning_mbb, 14);
    // vs_presel->add("mMMC",                 "m_{#tau#tau} (MMC) [GeV]",                          12, binning_mtt, 13);
    // vs_presel->add("mHH",                  "m_{HH} [GeV]",                                      12, binning_mhh, 19);
    // // vs_presel->add("mHHScaled",            "m_{HH} (Scaled) [GeV]",                             12);
    // vs_presel->add("dRBB",                 "#Delta R(B,B)",                                     4 , binning_dR, 11);
    // vs_presel->add("pTBB",                 "p_{T}^{b,b}",                                       10, binning_pTxx, 13);
    // vs_presel->add("dRTauTau",             "#Delta R(#tau_{had},#tau_{had})",                   4 , binning_dR, 11);
    // vs_presel->add("pTTauTau",             "p_{T}^{#tau,#tau}",                                 10, binning_pTxx, 13);
    // vs_presel->add("dPhiBBTauTau",         "#Delta #phi{bb,#tau#tau}",                          4 );
    // vs_presel->add("dRBBTauTau",           "#Delta R{bb,#tau#tau}",                             5 );
    // vs_presel->add("pTBalance",            "p_{T}^{b,b} / p_{T}^{#tau,#tau}",                   5 );
    // vs_presel->add("MET",                  "E_{T}^{miss} [GeV]",                                5 );
    // vs_presel->add("nJets",                "Number of Jets",                                    1 );
    // vs_presel->add("Yield",                "Event Yield",                                       1 );
    // vs_presel->add("metSig",               "MET significance",                                  2 );
    // vs_presel->add("metSigPU",             "MET significance (pile-up)",                        2 );
    // vs_presel->add("Jet0Pt",               "Leading b-jet p_{T} (before corr.) [GeV]",          10);
    // vs_presel->add("Jet1Pt",               "Sub-leading b-jet p_{T} (before corr.) [GeV]",      10);
    // vs_presel->add("pTB0",                 "Leading b-jet p_{T} [GeV]",                         10);
    // vs_presel->add("pTB1",                 "Sub-leading b-jet p_{T} [GeV]",                     10);
    // vs_presel->add("LeadJetPt",            "Leading jet p_{T} [GeV]",                           10);
    // vs_presel->add("SubleadJetPt",         "Sub-leading jet p_{T} [GeV]",                       10);
    // vs_presel->add("LeadJetEta",           "Leading jet #eta",                                  2 );
    // vs_presel->add("SubleadJetEta",        "Sub-leading jet #eta",                              2 );
    // vs_presel->add("Jet0Eta",              "Leading b-jet #eta (before corr.)",                 2 );
    // vs_presel->add("Jet1Eta",              "Sub-leading b-jet #eta (before corr.)",             2 );
    // vs_presel->add("Tau0Pt",               "Leading #tau_{had} p_{T} [GeV]",                    10);
    // vs_presel->add("Tau1Pt",               "Sub-leading #tau_{had} p_{T} [GeV]",                10);
    // vs_presel->add("Tau0Eta",              "Leading #tau_{had} #eta",                           2 );
    // vs_presel->add("Tau1Eta",              "Sub-leading #tau_{had} #eta",                       2 );
    // vs_presel->add("Tau0Ntrk",             "Leading #tau_{had} N_{tracks}",                     1 );
    // vs_presel->add("Tau1Ntrk",             "Sub-leading #tau_{had} N_{tracks}",                 1 );
    // vs_presel->add("dPhiTauTau",           "#Delta #phi (#tau_{had},#tau_{had})",               4 );
    // vs_presel->add("SherpapTV",            "Sherpa p_{T}^{Z} [GeV]",                            1 , binning_pTV, 6);

    auto binningFromFile = [](const std::string& fn)
    { 
        return BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>(fn));
    };

    vector<double> binning2HDM251   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM260.txt");
    vector<double> binning2HDM260   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM260.txt");
    vector<double> binning2HDM280   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM280.txt");
    vector<double> binning2HDM300   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM300.txt");
    vector<double> binning2HDM325   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM300.txt");
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

    // vector<double> binning2HDM900  =   {0.005, 0.01, 0.027, 0.085, 0.25, 0.5, 0.7, 0.878, 0.994, 1}; 
    // vector<double> binning2HDM800  =   {0.005, 0.01, 0.027, 0.085, 0.25, 0.5, 0.7, 0.878, 0.994, 1}; 
    // vector<double> binning2HDM1000  =   {0.005, 0.01, 0.027, 0.085, 0.25, 0.5, 0.7, 0.878, 0.994, 1}; 
    // vector<double> binning2HDM1100  =   {0.005, 0.01, 0.027, 0.085, 0.25, 0.5, 0.7, 0.878, 0.994, 1}; 
    // vector<double> binning2HDM1200  =   {0.005, 0.01, 0.027, 0.085, 0.25, 0.5, 0.7, 0.878, 0.994, 1}; 
    // vector<double> binning2HDM1400  =   {0.005, 0.01, 0.027, 0.085, 0.25, 0.5, 0.7, 0.878, 0.994, 1}; 
    // vector<double> binning2HDM1600  =   {0.005, 0.01, 0.027, 0.085, 0.25, 0.5, 0.7, 0.878, 0.994, 1}; 

    Tools::printVector(binningSMBDT); 

    Variables* vs_pnn = new Variables();
    vs_pnn->add("PNN251",                  "PNN251",                       100,   &binning2HDM251[0], binning2HDM251.size()-1);
    vs_pnn->add("PNN260",                  "PNN260",                       100,   &binning2HDM260[0], binning2HDM260.size()-1);
    vs_pnn->add("PNN280",                  "PNN280",                       100,   &binning2HDM280[0], binning2HDM280.size()-1);
    vs_pnn->add("PNN300",                  "PNN300",                       100,   &binning2HDM300[0], binning2HDM300.size()-1);
    vs_pnn->add("PNN325",                  "PNN325",                       100,   &binning2HDM325[0], binning2HDM325.size()-1);
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
    vs_pnn->add("PNN1100",                 "PNN1100",                      100,   &binning2HDM1100[0], binning2HDM1100.size()-1);
    vs_pnn->add("PNN1200",                 "PNN1200",                      100,   &binning2HDM1200[0], binning2HDM1200.size()-1);
    vs_pnn->add("PNN1400",                 "PNN1400",                      100,   &binning2HDM1400[0], binning2HDM1400.size()-1);
    vs_pnn->add("PNN1600",                 "PNN1600",                      100,   &binning2HDM1600[0], binning2HDM1600.size()-1);

    Variables* vs_bdt = new Variables();
    vs_bdt->add("SMBDT",                   "SM BDT",                       100,   &binningSMBDT[0], binningSMBDT.size()-1);

    Systematics* ss = new Systematics();
    // ss->add("ZJETS_GENERATOR", "ZJETS_GENERATOR", eSystematicType::TwoSide);
    // ss->add("ZJETS_ACC_GENERATOR_mBB", "ZJETS_ACC_GENERATOR_mBB", eSystematicType::OneSide);
    // ss->add("ZJETS_ACC_GENERATOR_pTBB", "ZJETS_ACC_GENERATOR_pTBB", eSystematicType::OneSide);
    // ss->add("ZJETS_ACC_GENERATOR_pTTauTau", "ZJETS_ACC_GENERATOR_pTTauTau", eSystematicType::OneSide);
    // ss->add("ZJETS_ACC_GENERATOR_dRBB", "ZJETS_ACC_GENERATOR_dRBB", eSystematicType::OneSide);
    // ss->add("ZJETS_ACC_GENERATOR_dRTauTau", "ZJETS_ACC_GENERATOR_dRTauTau", eSystematicType::OneSide);
    // // ss->add("ZJETS_ACC_GENERATOR", "ZJETS_GENERATOR", eSystematicType::OneSide);
    // ss->add("ZJETS_ACC_GENERATOR_MVA_SMBDT", "ZJETS_ACC_GENERATOR_MVA_SMBDT", eSystematicType::OneSide);
    // ss->add("ZJETS_ACC_GENERATOR_MVA_PNN300", "ZJETS_ACC_GENERATOR_MVA_PNN300", eSystematicType::OneSide);
    // ss->add("ZJETS_ACC_GENERATOR_MVA_PNN500", "ZJETS_ACC_GENERATOR_MVA_PNN500", eSystematicType::OneSide);
    // ss->add("ZJETS_ACC_GENERATOR_MVA_PNN800", "ZJETS_ACC_GENERATOR_MVA_PNN800", eSystematicType::OneSide);
    // ss->add("ZJETS_ACC_GENERATOR_MVA_PNN1000", "ZJETS_ACC_GENERATOR_MVA_PNN1000", eSystematicType::OneSide);
    // ss->add("ZJETS_ACC_GENERATOR_MVA_PNN1600", "ZJETS_ACC_GENERATOR_MVA_PNN1600", eSystematicType::OneSide);
    // ss->add("ZJETS_CR_DATA", "ZJETS_CR_DATA", eSystematicType::TwoSide);

    CompInfo* info = new CompInfo();
    info->ratio_high = 1.88;
    info->ratio_low = 0.12;
    info->shape_only = true;
    info->save_ratio = false;

    AutoBinningInfo* abi = new AutoBinningInfo();
    // abi->n_bins = 16;
    abi->n_bins = 5;
    AutoBinningTool* abt = new AutoBinningTool_v2(abi, eProcess::ZllHFMG);

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
        // SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/Ztautau-SS/MG";

        if (HistToolHelper::check(c)) {
            ct->manipulate(c);
            // ct->rebin(c, eRebinOption::Array);
            abt->run(c);
            abt->rebin(c);
            Tools::printVector(abt->binning());
        }

        if (ct->check(c))
        {
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
    info->save_ratio = false;
    abi->n_bins = 28;

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
        // SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/Ztautau-SS/MG";

        if (HistToolHelper::check(c)) {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::Array, "trans", true);
            // abt->run(c);
            // abt->rebin(c);
            Tools::printVector(abt->binning());
        }

        if (ct->check(c))
        {
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
        // SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/Ztautau-SS/MG";
        if (HistToolHelper::check(c)) {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::Array, "trans", true);
            // abt->run(c);
            // abt->rebin(c);
            Tools::printVector(abt->binning());
        }

        if (ct->check(c))
        {
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