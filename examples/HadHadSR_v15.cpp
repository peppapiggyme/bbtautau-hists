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

void hadhadsr_v15(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    // rs->add("2tag2pjet_0ptv_LL_OS",     "0.046 < PNN900 < 0.994, PNN1000 < 0.994",         eRegionType::SR);
    // rs->add("2tag2pjet_0ptv_LL_OS",     "0.027 < PNN1000 < 0.994",         eRegionType::SR);
    // rs->add("2tag2pjet_0ptv_LL_OS",     "0.878 < PNN1000 < 0.994",         eRegionType::SR);
    // rs->add("2tag2pjet_0ptv_LL_OS",     "0.017 < PNN1100 < 0.994, PNN1000 < 0.994",        eRegionType::SR);
    // rs->add("2tag2pjet_0ptv_LL_OS",     "0.012 < PNN1200 < 0.992, PNN1000 < 0.994",        eRegionType::SR);
    // rs->add("2tag2pjet_0ptv_LL_OS",     "0.005 < PNN1400 < 0.987, PNN1000 < 0.994",        eRegionType::SR);
    // rs->add("2tag2pjet_0ptv_LL_OS",     "0.002 < PNN1600 < 0.981, PNN1000 < 0.994",        eRegionType::SR);
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::CR);
    rs->add("2tag2pjet_0ptv_LL_SS",     "2 b-tag, 2 loose #tau, SS",        eRegionType::CR);
    // rs->add("1tag2pjet_0ptv_LL_OS",     "1 b-tag, 2 loose #tau, OS",        eRegionType::CR);
    // rs->add("1tag2pjet_0ptv_LL_SS",     "1 b-tag, 2 loose #tau, SS",        eRegionType::CR);
    // rs->add("2tag2pjet_0ptv_LL_OS",     "ZVR Tight",         eRegionType::SR);

    vector<double> binningMhh{400, 450, 500, 550, 600, 630, 660, 690, 720, 750, 780, 810, 840, 870, 900, 930, 960, 1000, 1050, 1100, 1150, 1200, 1300, 1400};
    vector<double> binningMbb{0, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 200, 220, 240, 280, 320, 400};
    vector<double> binningMtt{70, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 110};
    vector<double> binningdRBB{0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.4, 4.2, 5.0};
    Variables* vs = new Variables();
    // vs->add("mBB",                  "m_{BB} [GeV]",                                      1,      &binningMbb[0],     binningMbb.size()-1);
    // vs->add("mMMC",                 "m_{#tau#tau} (MMC) [GeV]",                          1,      &binningMtt[0],     binningMtt.size()-1);
    // vs->add("mHH",                  "m_{HH} [GeV]",                                      1,      &binningMhh[0],     binningMhh.size()-1);
    // vs->add("mHHScaled",            "m_{HH} (Scaled) [GeV]",                             1,      &binningMhh[0],     binningMhh.size()-1);
    // vs->add("dRBB",                 "#Delta R(B,B)",                                     1,      &binningdRBB[0],    binningdRBB.size()-1);
    // vs->add("dRTauTau",             "#Delta R(#tau_{had},#tau_{had})",                   1,      &binningdRBB[0],    binningdRBB.size()-1);
    // vs->add("dRBB_rapidity",        "#Delta R(B,B) using Rapidity",                      4 );
    // vs->add("dRBB_uncorrected",     "#Delta R(B,B) (before corr.)",                      4 );
    // vs->add("pTBB",                 "p_{T}^{b,b}",                                       10);
    // vs->add("pTTauTau",             "p_{T}^{#tau,#tau}",                                 10);
    // vs->add("dPhiBBTauTau",         "#Delta #phi{bb,#tau#tau}",                          4 );
    // vs->add("dRBBTauTau",           "#Delta R{bb,#tau#tau}",                             5 );
    // vs->add("pTBalance",            "p_{T}^{b,b} / p_{T}^{#tau,#tau}",                   5 );
    // vs->add("MET",                  "E_{T}^{miss} [GeV]",                                5 );
    // vs->add("nJets",                "Number of Jets",                                    1 );
    // vs->add("Yield",                "Event Yield",                                       1 );
    // vs->add("metSig",               "MET significance",                                  2 );
    // vs->add("metSigPU",             "MET significance (pile-up)",                        2 );
    // vs->add("Jet0Pt",               "Leading b-jet p_{T} (before corr.) [GeV]",          10);
    // vs->add("Jet1Pt",               "Sub-leading b-jet p_{T} (before corr.) [GeV]",      10);
    // vs->add("pTB0",                 "Leading b-jet p_{T} [GeV]",                         10);
    // vs->add("pTB1",                 "Sub-leading b-jet p_{T} [GeV]",                     10);
    // vs->add("LeadJetPt",            "Leading jet p_{T} [GeV]",                           10);
    // vs->add("SubleadJetPt",         "Sub-leading jet p_{T} [GeV]",                       10);
    // vs->add("LeadJetEta",           "Leading jet #eta",                                  2 );
    // vs->add("SubleadJetEta",        "Sub-leading jet #eta",                              2 );
    // vs->add("Jet0Eta",              "Leading b-jet #eta (before corr.)",                 2 );
    // vs->add("Jet1Eta",              "Sub-leading b-jet #eta (before corr.)",             2 );
    // vs->add("Tau0Pt",               "Leading #tau_{had} p_{T} [GeV]",                    10);
    // vs->add("Tau1Pt",               "Sub-leading #tau_{had} p_{T} [GeV]",                10);
    // vs->add("Tau0Eta",              "Leading #tau_{had} #eta",                           2 );
    // vs->add("Tau1Eta",              "Sub-leading #tau_{had} #eta",                       2 );
    // vs->add("Tau0Ntrk",             "Leading #tau_{had} N_{tracks}",                     1 );
    // vs->add("Tau1Ntrk",             "Sub-leading #tau_{had} N_{tracks}",                 1 );
    // vs->add("dPhiTauTau",           "#Delta #phi (#tau_{had},#tau_{had})",               4 );

    Variables* vs_pnn = new Variables();
    auto binningFromFile = [](const std::string& fn)
    { 
        return BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>(fn));
    };

    vector<double> binning2HDM251   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM251.txt");
    vector<double> binning2HDM260   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM260.txt");
    vector<double> binning2HDM280   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM280.txt");
    vector<double> binning2HDM300   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM300.txt");
    vector<double> binning2HDM325   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM325.txt");
    vector<double> binning2HDM350   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM350.txt");
    vector<double> binning2HDM375   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM375.txt");
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
    // vector<double> binningSMBDT     =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_SMBDT.txt");
    vector<double> binningSMBDT      =   BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_SMBDT.txt"), 1000, true);

    vs_pnn->add("PNN251",       "PNN251",       1,      &binning2HDM251[0],     binning2HDM251.size()-1);
    vs_pnn->add("PNN260",       "PNN260",       1,      &binning2HDM260[0],     binning2HDM260.size()-1);
    vs_pnn->add("PNN280",       "PNN280",       1,      &binning2HDM280[0],     binning2HDM280.size()-1);
    vs_pnn->add("PNN300",       "PNN300",       1,      &binning2HDM300[0],     binning2HDM300.size()-1);
    vs_pnn->add("PNN325",       "PNN325",       1,      &binning2HDM325[0],     binning2HDM325.size()-1);
    vs_pnn->add("PNN350",       "PNN350",       1,      &binning2HDM350[0],     binning2HDM350.size()-1);
    vs_pnn->add("PNN375",       "PNN375",       1,      &binning2HDM375[0],     binning2HDM375.size()-1);
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

    Systematics* ss = new Systematics();
    // WARN: you must pass OneSide for DrawStackTool, otherwise it will double count!!!
    // ss->add("ZJETS_ACC_GENERATOR_MBB", "ZJETS_ACC_GENERATOR_MBB", eSystematicType::OneSide, kBlue);
    // ss->add("ZJETS_ACC_SCALE_MVA_PNN1000", "ZJETS_ACC_SCALE_MVA_PNN1000", eSystematicType::OneSide, kBlue);
    // ss->add("ZJETS_NORM_CR", "ZJETS_NORM_CR", eSystematicType::OneSide, kBlue);
    // ss->add("ZJETS_XSEC", "ZJETS_XSEC", eSystematicType::OneSide, kBlue);
    // ss->add("TAUS_TRUEHADTAU_SME_TES_INSITUEXP", "TAUS_TRUEHADTAU_SME_TES_INSITUEXP", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_SME_TES_INSITUFIT", "TAUS_TRUEHADTAU_SME_TES_INSITUFIT", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_SME_TES_MODEL_CLOSURE", "TAUS_TRUEHADTAU_SME_TES_MODEL_CLOSURE", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_SME_TES_PHYSICSLIST", "TAUS_TRUEHADTAU_SME_TES_PHYSICSLIST", eSystematicType::OneSide, kBlack);
    // ss->add("PRW_DATASF", "PRW_DATASF", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEELECTRON_EFF_ELEBDT_STAT", "TAUS_TRUEELECTRON_EFF_ELEBDT_STAT", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEELECTRON_EFF_ELEBDT_SYST", "TAUS_TRUEELECTRON_EFF_ELEBDT_SYST", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL", "TAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RECO_TOTAL", "TAUS_TRUEHADTAU_EFF_RECO_TOTAL", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPT2025", "TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPT2025", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPT2530", "TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPT2530", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPT3040", "TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPT3040", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPTGE40", "TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPTGE40", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RNNID_3PRONGSTATSYSTPT2025", "TAUS_TRUEHADTAU_EFF_RNNID_3PRONGSTATSYSTPT2025", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RNNID_3PRONGSTATSYSTPT2530", "TAUS_TRUEHADTAU_EFF_RNNID_3PRONGSTATSYSTPT2530", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RNNID_3PRONGSTATSYSTPT3040", "TAUS_TRUEHADTAU_EFF_RNNID_3PRONGSTATSYSTPT3040", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RNNID_3PRONGSTATSYSTPTGE40", "TAUS_TRUEHADTAU_EFF_RNNID_3PRONGSTATSYSTPTGE40", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RNNID_HIGHPT", "TAUS_TRUEHADTAU_EFF_RNNID_HIGHPT", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_RNNID_SYST", "TAUS_TRUEHADTAU_EFF_RNNID_SYST", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA161718", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA161718", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA1718", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA1718", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2016", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2016", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018AFTTS1", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018AFTTS1", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC161718", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC161718", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC1718", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC1718", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2016", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2016", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2018", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2018", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2018AFTTS1", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2018AFTTS1", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYST161718", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYST161718", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYST1718", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYST1718", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYST2016", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYST2016", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYST2018", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYST2018", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYST2018AFTTS1", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYST2018AFTTS1", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU161718", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU161718", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU1718", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU1718", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU2016", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU2016", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU2018", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU2018", eSystematicType::OneSide, kBlack);
    // ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU2018AFTTS1", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU2018AFTTS1", eSystematicType::OneSide, kBlack);
    // ss->add("FT_EFF_Eigen_B_0", "FT_EFF_Eigen_B_0", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_Eigen_B_1", "FT_EFF_Eigen_B_1", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_Eigen_B_2", "FT_EFF_Eigen_B_2", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_Eigen_C_0", "FT_EFF_Eigen_C_0", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_Eigen_C_1", "FT_EFF_Eigen_C_1", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_Eigen_C_2", "FT_EFF_Eigen_C_2", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_Eigen_C_3", "FT_EFF_Eigen_C_3", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_Eigen_Light_0", "FT_EFF_Eigen_Light_0", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_Eigen_Light_1", "FT_EFF_Eigen_Light_1", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_Eigen_Light_2", "FT_EFF_Eigen_Light_2", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_Eigen_Light_3", "FT_EFF_Eigen_Light_3", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_extrapolation", "FT_EFF_extrapolation", eSystematicType::OneSide, kBlue);
    // ss->add("FT_EFF_extrapolation_from_charm", "FT_EFF_extrapolation_from_charm", eSystematicType::OneSide, kBlue);

    // // // ss->add("ZJETS_PTV_DEPENDENCY_CR", "ZJETS_PTV_DEPENDENCY_CR", eSystematicType::OneSide, kBlue);

    DrawStackInfo* info = new DrawStackInfo();
    info->blind = false;
    info->ratio_high = 1.68;
    info->ratio_low = 0.32;
    info->signal_scale = 1;
    // info->signal_scale = 0.04;

    info->logy = false;

    for (VariableInfo* v : *(vs->content()))
    {
        Processes* ps = new Processes();
        ps->add("data",                 "data",                         eProcessType::DATA,     eProcess::DATA,         "data",                         kBlack);
        ps->add("Zttbb",                "Z#tau#tau + bb",               eProcessType::BKG   ,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zttbc",                "Z#tau#tau + bc",               eProcessType::BKG   ,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zttbl",                "Z#tau#tau + bl",               eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zttcc",                "Z#tau#tau + cc",               eProcessType::BKG   ,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zttcl",                "Z#tau#tau + cl",               eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zttl",                 "Z#tau#tau + l",                eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zbb",                  "Zll + bb",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zbc",                  "Zll + bc",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zbl",                  "Zll + bl",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zcc",                  "Zll + cc",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zcl",                  "Zll + cl",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zl",                   "Zll + l",                      eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("ttbar",                "t#bar{t}",                     eProcessType::BKG,      eProcess::TTBAR,        "t#bar{t}",                     kWhite);
        ps->add("ttbarSFTF",            "t#bar{t} TF",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
        ps->add("ttbarSFFT",            "t#bar{t} FT",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
        ps->add("ttbarSFFF",            "t#bar{t} FF",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
        ps->add("stops",                "single top",                   eProcessType::BKG,      eProcess::STOP,         "single top",                   kOrange+2);
        ps->add("stopt",                "single top",                   eProcessType::BKG,      eProcess::STOP,         "single top",                   kOrange+2);
        ps->add("stopWt",               "single top",                   eProcessType::BKG,      eProcess::STOP,         "single top",                   kOrange+2);
        ps->add("WW",                   "WW",                           eProcessType::BKG,      eProcess::DIBOSON,      "Diboson",                      kGreen-8);
        ps->add("WZ",                   "WZ",                           eProcessType::BKG,      eProcess::DIBOSON,      "Diboson",                      kGreen-8);
        ps->add("ZZ",                   "ZZ",                           eProcessType::BKG,      eProcess::DIBOSON,      "Diboson",                      kGreen-8);
        ps->add("Wtt",                  "W+jets",                       eProcessType::BKG,      eProcess::WJETS,        "W+jets",                       kGreen-10);
        ps->add("W",                    "W+jets",                       eProcessType::BKG,      eProcess::WJETS,        "W+jets",                       kGreen-10);
        ps->add("Fake",                 "Multi-jet",                    eProcessType::BKG,      eProcess::FAKE,         "Multi-jet",                    kMagenta-10);
        ps->add("ttW",                  "ttV",                          eProcessType::BKG,      eProcess::ttV,          "ttV",                          kYellow+3);
        ps->add("ttZ",                  "ttV",                          eProcessType::BKG,      eProcess::ttV,          "ttV",                          kYellow+3);
        ps->add("ttH",                  "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggZHbb",               "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("qqZHbb",               "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggZHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("qqZHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("VBFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        // ps->add("Hhhbbtautau900",       "X900",                         eProcessType::SIG,      eProcess::X900,        "X900",                         kRed);
        // ps->add("Hhhbbtautau1000",      "X1000",                        eProcessType::SIG,      eProcess::X1000,        "X1000",                        kMagenta);
        // ps->add("hhttbb",               "Non-res HH",                   eProcessType::SIG,      eProcess::SMHH,         "Non-res HH",                   kMagenta+2);

        Config* c = new Config(b, ps, rs, vs, ss);
        c->load(filename, "Preselection");
        info->parameter = "Preselection";

        c->updateHistogramPtr(rs->content()->back(), v);
        DrawStackTool* ds = new DrawStackTool(info);
        ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/SR-OSSS-0822/";
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
        delete ds;
        delete c;
    }

    info->logy = true;
    info->ratio_high = 1.88;
    info->ratio_low = 0.12;

    for (VariableInfo* v : *(vs_pnn->content()))
    {
        Processes* ps = new Processes();
        ps->add("data",                 "data",                         eProcessType::DATA,     eProcess::DATA,         "data",                         kBlack);
        ps->add("Zttbb",                "Z#tau#tau + bb",               eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zttbc",                "Z#tau#tau + bc",               eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zttbl",                "Z#tau#tau + bl",               eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zttcc",                "Z#tau#tau + cc",               eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zttcl",                "Z#tau#tau + cl",               eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zttl",                 "Z#tau#tau + l",                eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zbb",                  "Zll + bb",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zbc",                  "Zll + bc",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zbl",                  "Zll + bl",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zcc",                  "Zll + cc",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zcl",                  "Zll + cl",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zl",                   "Zll + l",                      eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("ttbar",                "t#bar{t}",                     eProcessType::BKG,      eProcess::TTBAR,        "t#bar{t}",                     kWhite);
        ps->add("ttbarSFTF",            "t#bar{t} TF",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
        ps->add("ttbarSFFT",            "t#bar{t} FT",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
        ps->add("ttbarSFFF",            "t#bar{t} FF",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
        ps->add("stops",                "single top",                   eProcessType::BKG,      eProcess::STOP,         "single top",                   kOrange+2);
        ps->add("stopt",                "single top",                   eProcessType::BKG,      eProcess::STOP,         "single top",                   kOrange+2);
        ps->add("stopWt",               "single top",                   eProcessType::BKG,      eProcess::STOP,         "single top",                   kOrange+2);
        ps->add("WW",                   "WW",                           eProcessType::BKG,      eProcess::DIBOSON,      "Diboson",                      kGreen-8);
        ps->add("WZ",                   "WZ",                           eProcessType::BKG,      eProcess::DIBOSON,      "Diboson",                      kGreen-8);
        ps->add("ZZ",                   "ZZ",                           eProcessType::BKG,      eProcess::DIBOSON,      "Diboson",                      kGreen-8);
        ps->add("Wtt",                  "W+jets",                       eProcessType::BKG,      eProcess::WJETS,        "W+jets",                       kGreen-10);
        ps->add("W",                    "W+jets",                       eProcessType::BKG,      eProcess::WJETS,        "W+jets",                       kGreen-10);
        ps->add("Fake",                 "Multi-jet",                    eProcessType::BKG,      eProcess::FAKE,         "Multi-jet",                    kMagenta-10);
        ps->add("ttW",                  "ttV",                          eProcessType::BKG,      eProcess::ttV,          "ttV",                          kYellow+3);
        ps->add("ttZ",                  "ttV",                          eProcessType::BKG,      eProcess::ttV,          "ttV",                          kYellow+3);
        ps->add("ttH",                  "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ZHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ZHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("VBFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        // ps->add("Hhhbbtautau900",       "X900",                         eProcessType::SIG,      eProcess::X900,        "X900",                         kRed);
        // ps->add("Hhhbbtautau1000",      "X1000",                        eProcessType::SIG,      eProcess::XtoHH,        "X1000",                        kMagenta+2);
        // ps->add("hhttbb",               "Non-res HH",                   eProcessType::SIG,      eProcess::SMHH,         "Non-res HH",                   kMagenta+2);

        Config* c = new Config(b, ps, rs, vs, ss);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "PNNScorePreselection";

        c->updateHistogramPtr(rs->content()->back(), v);
        DrawStackTool* ds = new DrawStackTool(info);
        ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/SR-OSSS-0822";
        if (ds->check(c))
        {
            ds->manipulate(c);
            // ds->rebin(c, eRebinOption::Array);
            ds->rebin(c, eRebinOption::Array, "using array and transformed", true);
            ds->paint(c);
            ds->makeYield(c);
            ds->run(c);
        }
        else 
        {
            clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
        }

        delete ps;
        delete ds;
        delete c;
    }

    for (VariableInfo* v : *(vs_bdt->content()))
    {
        Processes* ps = new Processes();
        ps->add("data",                 "data",                         eProcessType::DATA,     eProcess::DATA,         "data",                         kBlack);
        ps->add("Zttbb",                "Z#tau#tau + bb",               eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zttbc",                "Z#tau#tau + bc",               eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zttbl",                "Z#tau#tau + bl",               eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zttcc",                "Z#tau#tau + cc",               eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zttcl",                "Z#tau#tau + cl",               eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zttl",                 "Z#tau#tau + l",                eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zbb",                  "Zll + bb",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zbc",                  "Zll + bc",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zbl",                  "Zll + bl",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zcc",                  "Zll + cc",                     eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10);
        ps->add("Zcl",                  "Zll + cl",                     eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("Zl",                   "Zll + l",                      eProcessType::BKG,      eProcess::ZllLF,        "Z+lf",                         kBlue-5);
        ps->add("ttbar",                "t#bar{t}",                     eProcessType::BKG,      eProcess::TTBAR,        "t#bar{t}",                     kWhite);
        ps->add("ttbarSFTF",            "t#bar{t} TF",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
        ps->add("ttbarSFFT",            "t#bar{t} FT",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
        ps->add("ttbarSFFF",            "t#bar{t} FF",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
        ps->add("stops",                "single top",                   eProcessType::BKG,      eProcess::STOP,         "single top",                   kOrange+2);
        ps->add("stopt",                "single top",                   eProcessType::BKG,      eProcess::STOP,         "single top",                   kOrange+2);
        ps->add("stopWt",               "single top",                   eProcessType::BKG,      eProcess::STOP,         "single top",                   kOrange+2);
        ps->add("WW",                   "WW",                           eProcessType::BKG,      eProcess::DIBOSON,      "Diboson",                      kGreen-8);
        ps->add("WZ",                   "WZ",                           eProcessType::BKG,      eProcess::DIBOSON,      "Diboson",                      kGreen-8);
        ps->add("ZZ",                   "ZZ",                           eProcessType::BKG,      eProcess::DIBOSON,      "Diboson",                      kGreen-8);
        ps->add("Wtt",                  "W+jets",                       eProcessType::BKG,      eProcess::WJETS,        "W+jets",                       kGreen-10);
        ps->add("W",                    "W+jets",                       eProcessType::BKG,      eProcess::WJETS,        "W+jets",                       kGreen-10);
        ps->add("Fake",                 "Multi-jet",                    eProcessType::BKG,      eProcess::FAKE,         "Multi-jet",                    kMagenta-10);
        ps->add("ttW",                  "ttV",                          eProcessType::BKG,      eProcess::ttV,          "ttV",                          kYellow+3);
        ps->add("ttZ",                  "ttV",                          eProcessType::BKG,      eProcess::ttV,          "ttV",                          kYellow+3);
        ps->add("ttH",                  "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ZHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ZHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("VBFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        // ps->add("Hhhbbtautau900",       "X900",                         eProcessType::SIG,      eProcess::X900,        "X900",                         kRed);
        // ps->add("Hhhbbtautau1000",      "X1000",                        eProcessType::SIG,      eProcess::XtoHH,        "X1000",                        kMagenta+2);
        // ps->add("hhttbb",               "Non-res HH",                   eProcessType::SIG,      eProcess::SMHH,         "Non-res HH",                   kMagenta+2);

        Config* c = new Config(b, ps, rs, vs, ss);
        c->load(filename, "BDTScorePreselection");
        info->parameter = "BDTScorePreselection";

        c->updateHistogramPtr(rs->content()->back(), v);
        DrawStackTool* ds = new DrawStackTool(info);
        ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/SR-OSSS-0822";
        if (ds->check(c))
        {
            ds->manipulate(c);
            // ds->rebin(c, eRebinOption::Array);
            ds->rebin(c, eRebinOption::Array, "using array and transformed", true);
            ds->paint(c);
            ds->makeYield(c);
            ds->run(c);
        }
        else 
        {
            clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
        }

        delete ps;
        delete ds;
        delete c;
    }

    delete b;
    delete rs;
    delete vs;
    delete ss;
    delete info;

}
