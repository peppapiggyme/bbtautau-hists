#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "DrawStackTool.h"
#include "AutoBinningTool.h"


#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;
using BU = BinningUtils;

void test_bbll(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    // rs->add("2tag2pjet_0ptv_LL_OS",     "bbll Z+HF CR, 50 < mBB < 150 GeV",        eRegionType::SR);
    // rs->add("2tag2pjet_0ptv_LL_OS",     "bbll Z+HF CR, mBB sideband",        eRegionType::SR);
    // rs->add("2tag2pjet_0ptv_LL_OS",     "bbll Z+HF CR, 60 < mBB < 95 GeV",        eRegionType::SR);
    rs->add("2tag2pjet_0ptv_LL_OS",     "bbll Z+HF CR, 135 < mBB < 160 GeV",        eRegionType::SR);

    auto binningFromFile = [](const std::string& fn)
    { 
        return BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>(fn));
    };
    
    vector<double> binningMhh{200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1400, 1600};
    vector<double> binningMbb{0, 20, 30, 40, 50, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 150, 160, 180, 200, 240};
    vector<double> binningMtt{75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110};
    vector<double> binningdRBB{0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 5.0};
    vector<double> binningPtTauTau{0., 60., 120., 180., 240., 300., 400, 1000.};
    vector<double> binning2HDM600   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM600.txt");
    vector<double> binning2HDM700   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM700.txt");
    vector<double> binning2HDM800   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM800.txt");
    vector<double> binning2HDM900   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM900.txt");
    vector<double> binning2HDM1000  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1000.txt");
    vector<double> binning2HDM1200  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1200.txt");
    vector<double> binning2HDM1400  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1400.txt");
    vector<double> binning2HDM1600  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1600.txt");
    vector<double> binningSMBDT     =   BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_SMBDT.txt"), 1000, true);

    Variables* vs = new Variables();
    vs->add("mBB",                  "m_{BB} [GeV]",                                      1,         &binningMbb[0],         binningMbb.size()-1);
    vs->add("mMMC",                 "m_{ll} [GeV]",                                      1,         &binningMtt[0],         binningMtt.size()-1);
    vs->add("mHH",                  "m_{HH} [GeV]",                                      1,         &binningMhh[0],         binningMhh.size()-1);
    vs->add("pTBB",                 "p_{T}^{b, b}",                                      1,         &binningPtTauTau[0],    binningPtTauTau.size()-1);
    vs->add("pTTauTau",             "p_{T}^{l, l}",                                      1,         &binningPtTauTau[0],    binningPtTauTau.size()-1);
    // vs->add("dRBB",                 "#Delta R(B,B)",                                     1,      &binningdRBB[0],    binningdRBB.size()-1);
    vs->add("dRBB",                 "#Delta R(B,B)",                                     2);
    vs->add("dRTauTau",             "#Delta R(l,l)",                                     4 );
    // vs->add("dPhiBBTauTau",         "#Delta #phi{bb,#tau#tau}",                          4 );
    // vs->add("dRBBTauTau",           "#Delta R{bb,#tau#tau}",                             5 );
    // vs->add("pTBalance",            "p_{T}^{b,b} / p_{T}^{#tau,#tau}",                   5 );
    vs->add("MET",                  "E_{T}^{miss} [GeV]",                                5 );
    // vs->add("nJets",                "Number of Jets",                                    1 );
    // vs->add("Yield",                "Event Yield",                                       1 );
    // vs->add("metSig",               "MET significance",                                  2 );
    // vs->add("metSigPU",             "MET significance (pile-up)",                        2 );
    // vs->add("Jet0Pt",               "Leading b-jet p_{T} (before corr.) [GeV]",          10);
    // vs->add("Jet1Pt",               "Sub-leading b-jet p_{T} (before corr.) [GeV]",      10);
    vs->add("pTB0",                 "Leading b-jet p_{T} [GeV]",                         10);
    vs->add("pTB1",                 "Sub-leading b-jet p_{T} [GeV]",                     10);
    // vs->add("LeadJetPt",            "Leading jet p_{T} [GeV]",                           10);
    // vs->add("SubleadJetPt",         "Sub-leading jet p_{T} [GeV]",                       10);
    // vs->add("LeadJetEta",           "Leading jet #eta",                                  2 );
    // vs->add("SubleadJetEta",        "Sub-leading jet #eta",                              2 );
    // vs->add("Jet0Eta",              "Leading b-jet #eta (before corr.)",                 2 );
    // vs->add("Jet1Eta",              "Sub-leading b-jet #eta (before corr.)",             2 );
    vs->add("Tau0Pt",               "Leading lepton p_{T} [GeV]",                        5 );
    vs->add("Tau1Pt",               "Sub-leading lepton p_{T} [GeV]",                    5 );
    // vs->add("Tau0Eta",              "Leading #tau_{had} #eta",                           2 );
    // vs->add("Tau1Eta",              "Sub-leading #tau_{had} #eta",                       2 );
    // vs->add("Tau0Ntrk",             "Leading #tau_{had} N_{tracks}",                     1 );
    // vs->add("Tau1Ntrk",             "Sub-leading #tau_{had} N_{tracks}",                 1 );
    // vs->add("dPhiTauTau",           "#Delta #phi (#tau_{had},#tau_{had})",               4 );
    vs->add("PNN600",       "PNN600",       1,      &binning2HDM600[0],     binning2HDM600.size()-1);
    vs->add("PNN700",       "PNN700",       1,      &binning2HDM700[0],     binning2HDM700.size()-1);
    vs->add("PNN800",       "PNN800",       1,      &binning2HDM800[0],     binning2HDM800.size()-1);
    vs->add("PNN900",       "PNN900",       1,      &binning2HDM900[0],     binning2HDM900.size()-1);
    vs->add("PNN1000",      "PNN1000",      1,      &binning2HDM1000[0],    binning2HDM1000.size()-1);
    vs->add("PNN1200",      "PNN1200",      1,      &binning2HDM1200[0],    binning2HDM1200.size()-1);
    vs->add("PNN1400",      "PNN1400",      1,      &binning2HDM1400[0],    binning2HDM1400.size()-1);
    vs->add("PNN1600",      "PNN1600",      1,      &binning2HDM1600[0],    binning2HDM1600.size()-1);
    vs->add("SMBDT",        "SM BDT",       1,      &binningSMBDT[0],       binningSMBDT.size()-1);

    DrawStackInfo* info = new DrawStackInfo();
    info->blind = false;
    info->ratio_high = 1.68;
    info->ratio_low = 0.32;
    info->signal_scale = 1;

    info->logy = false;

    for (VariableInfo* v : *(vs->content()))
    {
        Processes* ps = new Processes();
        ps->add("data",                 "data",                         eProcessType::DATA,     eProcess::DATA,         "data",                         kBlack);
        ps->add("Ztt",                  "Z#tau#tau",                    eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-8);
        ps->add("Zeebb",                "Zll + bb",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("ZeeLowbb",             "Zll + bb",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zeecc",                "Zll + cc",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zeell",                "Zll + ll",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("ZeeLowllcc",           "Zll + ll",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("ZeeEW",                "Zll + ll",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zmmbb",                "Zll + bb",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("ZmmLowbb",             "Zll + bb",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zmmcc",                "Zll + cc",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zmmll",                "Zll + ll",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("ZmmLowllcc",           "Zll + ll",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("ZmmEW",                "Zll + ll",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("ttbar",                "t#bar{t}",                     eProcessType::BKG,      eProcess::TOP,          "top",                          kWhite);
        ps->add("fakes_top_ttbar",      "t#bar{t}",                     eProcessType::BKG,      eProcess::TOP,          "top",                          kWhite);
        ps->add("stop",                 "single top",                   eProcessType::BKG,      eProcess::TOP,          "top",                          kOrange+2);
        ps->add("fakes_top_stop",       "single top",                   eProcessType::BKG,      eProcess::TOP,          "top",                          kOrange+2);
        ps->add("diboson",              "diboson",                      eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-8);
        ps->add("fakes_Wenu",           "W+jets",                       eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-8);
        ps->add("fakes_Wmunu",          "W+jets",                       eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-8);
        ps->add("fakes_Wtaunu",         "W+jets",                       eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-8);
        ps->add("fakes_diboson",        "W+jets",                       eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-8);
        ps->add("ttV",                  "ttV",                          eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-8);
        ps->add("ttH",                  "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WH",                   "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ZH",                   "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ZHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggFHbb",               "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("VHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("VBFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("VHWW",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggFHWW",               "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("VBFHWW",               "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggFHZZ",               "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("VBFHZZ",               "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        // ps->add("hhwwbb",               "X1000",                        eProcessType::SIG,      eProcess::X1000,        "X1000",                        kRed);
        // ps->add("hhzzbb",               "X1000",                        eProcessType::SIG,      eProcess::X1000,        "X1000",                        kRed);
        // ps->add("hhttbb",               "X1000",                        eProcessType::SIG,      eProcess::X1000,        "X1000",                        kRed);

        Systematics* ss = new Systematics();
        // WARN: you must pass OneSide for DrawStackTool, otherwise it will double count!!!
        ss->add("ZJETS_GEN", "ZJETS_GEN", eSystematicType::OneSide, kBlue);
        ss->add("ZJETS_NORM", "ZJETS_NORM", eSystematicType::OneSide, kBlue);
        ss->add("ZJETS_XSEC", "ZJETS_XSEC", eSystematicType::OneSide, kBlue);

        Config* c = new Config(b, ps, rs, vs, ss);
        c->load(filename, "Preselection");
        info->parameter = "0716-Right";

        c->updateHistogramPtr(rs->content()->front(), v);
        DrawStackTool* ds = new DrawStackTool(info);
        ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/bbllCR-0716/";
        if (ds->check(c))
        {
            ds->manipulate(c);
            if (v->name.find("PNN") == std::string::npos)
            {
                info->logy = false;
                ds->rebin(c, eRebinOption::Array);
            }
            else
            {
                info->logy = true;
                ds->rebin(c, eRebinOption::Array, "", true);
            }
            ds->paint(c);
            ds->makeYield(c);
            ds->run(c);
        }
        else 
        {
            clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
        }

        delete ps;
        // delete ss;
        delete ds;
        delete c;
    }

    delete b;
    delete rs;
    delete vs;
    delete info;
}