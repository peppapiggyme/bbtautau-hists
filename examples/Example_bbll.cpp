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
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, Zll+HF CR",        eRegionType::SR);

    vector<double> binningMhh{200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1400, 1600};
    vector<double> binningMbb{0, 20, 30, 40, 50, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 150, 160, 180, 200, 240};
    vector<double> binningMtt{75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110};
    vector<double> binningPtTauTau{0., 60., 120., 180., 240., 300., 400, 1000.};
    Variables* vs = new Variables();
    vs->add("mBB",                  "m_{BB} [GeV]",                                      1,         &binningMbb[0],         binningMbb.size()-1);
    vs->add("mMMC",                 "m_{ll} [GeV]",                                      1,         &binningMtt[0],         binningMtt.size()-1);
    vs->add("dRBB",                 "#Delta R(B,B)",                                     4 );
    vs->add("pTBB",                 "p_{T}^{b,b}",                                       1,         &binningPtTauTau[0],    binningPtTauTau.size()-1);
    vs->add("dRTauTau",             "#Delta R(l,l)",                                     4 );
    vs->add("pTTauTau",             "p_{T}^{l, l}",                                      1,         &binningPtTauTau[0],    binningPtTauTau.size()-1);
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

        Config* c = new Config(b, ps, rs, vs, ss);
        c->load(filename, "Preselection");
        info->parameter = "Presel";

        c->updateHistogramPtr(rs->content()->front(), v);
        DrawStackTool* ds = new DrawStackTool(info);
        ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/bbllCR/";
        if (ds->check(c))
        {
            ds->manipulate(c);
            ds->rebin(c, eRebinOption::Array);
            ds->paint(c);
            ds->makeYield(c);
            ds->run(c);
        }
        else 
        {
            clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
        }

        delete ps;
        delete ss;
        delete ds;
        delete c;
    }

    delete b;
    delete rs;
    delete vs;
    delete info;
}