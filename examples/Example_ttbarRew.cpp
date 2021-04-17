/**
 * This is for lephad channel
 */

#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "SystCompTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;

void test_lephad_ttbarRew(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_OS",     "2 b-tag, loose #tau, OS",        eRegionType::SR);

    Variables* vs_presel = new Variables();
    vs_presel->add("mBB",                  "m_{BB} [GeV]",                                      20);
    vs_presel->add("mMMC",                 "m_{#tau#tau} (MMC) [GeV]",                          20);
    vs_presel->add("mHH",                  "m_{HH} [GeV]",                                      50);
    vs_presel->add("dRBB",                 "#Delta R(B,B)",                                     5 );
    vs_presel->add("pTBB",                 "p_{T}^{b,b}",                                       20);
    vs_presel->add("dRTauTau",             "#Delta R(#tau_{had},#tau_{had})",                   5 );
    vs_presel->add("pTTauTau",             "p_{T}^{#tau,#tau}",                                 20);
    vs_presel->add("MET",                  "E_{T}^{miss} [GeV]",                                20);
    vs_presel->add("nJets",                "Number of Jets",                                    1 );
    vs_presel->add("Yield",                "Event Yield",                                       1 );
    vs_presel->add("pTB0",                 "Leading b-jet p_{T} [GeV]",                         20);
    vs_presel->add("pTB1",                 "Sub-leading b-jet p_{T} [GeV]",                     20);
    vs_presel->add("pTTau",                "#tau_{had} p_{T} [GeV]",                            20);
    vs_presel->add("pTLep",                "Lepton p_{T} [GeV]",                                20);
    vs_presel->add("pTTauLep",             "p_{T}^{#tau,lep}",                                  20);
    vs_presel->add("dRTauLep",             "#Delta R(#tau_{had},lep)",                          5 );

    Systematics* ss = new Systematics();
    // OK the name is misleading. It is ttbar with the nominal ttbar reweighting applied
    ss->add("TTBarNominal", "t#bar{t} ad-hoc reweight", eSystematicType::TwoSide);
    ss->add("TTBAR_REWEIGHT_NNLO", "t#bar{t} NNLO reweight", eSystematicType::TwoSide);

    CompInfo* info = new CompInfo();
    info->ratio_high = 1.38;
    info->ratio_low = 0.62;
    info->shape_only = false;

    for (VariableInfo* v : *(vs_presel->content()))
    {
        Processes* ps = new Processes();
        ps->add("ttbar",            "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);

        Config* c = new Config(b, ps, rs, vs_presel, ss);
        c->load(filename, "Preselection");
        info->parameter = "Presel";
        c->updateHistogramPtr(rs->content()->front(), v);
        SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/ttbarRew/";
        if (ct->check(c))
        {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::N_Rebin);
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
    delete info;

}

void test_hadhad_ttbarRew(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);

    Variables* vs_presel = new Variables();
    vs_presel->add("mBB",                  "m_{BB} [GeV]",                                      20);
    vs_presel->add("mMMC",                 "m_{#tau#tau} (MMC) [GeV]",                          20);
    vs_presel->add("mHH",                  "m_{HH} [GeV]",                                      12);
    vs_presel->add("mHHScaled",            "m_{HH} (Scaled) [GeV]",                             12);
    vs_presel->add("dRBB",                 "#Delta R(B,B)",                                     5 );
    vs_presel->add("pTBB",                 "p_{T}^{b,b}",                                       10);
    vs_presel->add("dRTauTau",             "#Delta R(#tau_{had},#tau_{had})",                   5 );
    vs_presel->add("pTTauTau",             "p_{T}^{#tau,#tau}",                                 10);
    vs_presel->add("dPhiBBTauTau",         "#Delta #phi{bb,#tau#tau}",                          8 );
    vs_presel->add("dRBBTauTau",           "#Delta R{bb,#tau#tau}",                             5 );
    vs_presel->add("pTBalance",            "p_{T}^{b,b} / p_{T}^{#tau,#tau}",                   5 );
    vs_presel->add("MET",                  "E_{T}^{miss} [GeV]",                                5 );
    vs_presel->add("nJets",                "Number of Jets",                                    1 );
    vs_presel->add("Yield",                "Event Yield",                                       1 );
    vs_presel->add("metSig",               "MET significance",                                  2 );
    vs_presel->add("metSigPU",             "MET significance (pile-up)",                        2 );
    vs_presel->add("Jet0Pt",               "Leading b-jet p_{T} (before corr.) [GeV]",          10);
    vs_presel->add("Jet1Pt",               "Sub-leading b-jet p_{T} (before corr.) [GeV]",      10);
    vs_presel->add("pTB0",                 "Leading b-jet p_{T} [GeV]",                         10);
    vs_presel->add("pTB1",                 "Sub-leading b-jet p_{T} [GeV]",                     10);
    vs_presel->add("LeadJetPt",            "Leading jet p_{T} [GeV]",                           10);
    vs_presel->add("SubleadJetPt",         "Sub-leading jet p_{T} [GeV]",                       10);
    vs_presel->add("LeadJetEta",           "Leading jet #eta",                                  2 );
    vs_presel->add("SubleadJetEta",        "Sub-leading jet #eta",                              2 );
    vs_presel->add("Jet0Eta",              "Leading b-jet #eta (before corr.)",                 2 );
    vs_presel->add("Jet1Eta",              "Sub-leading b-jet #eta (before corr.)",             2 );
    vs_presel->add("Tau0Pt",               "Leading #tau_{had} p_{T} [GeV]",                    5 );
    vs_presel->add("Tau1Pt",               "Sub-leading #tau_{had} p_{T} [GeV]",                5 );
    vs_presel->add("Tau0Eta",              "Leading #tau_{had} #eta",                           2 );
    vs_presel->add("Tau1Eta",              "Sub-leading #tau_{had} #eta",                       2 );
    vs_presel->add("Tau0Ntrk",             "Leading #tau_{had} N_{tracks}",                     1 );
    vs_presel->add("Tau1Ntrk",             "Sub-leading #tau_{had} N_{tracks}",                 1 );
    vs_presel->add("dPhiTauTau",           "#Delta #phi (#tau_{had},#tau_{had})",               4 );

    Variables* vs_pnn = new Variables();
    vs_pnn->add("PNN260",                  "PNN260",                                            100);
    vs_pnn->add("PNN300",                  "PNN300",                                            100);
    vs_pnn->add("PNN350",                  "PNN350",                                            100);
    vs_pnn->add("PNN400",                  "PNN400",                                            100);
    vs_pnn->add("PNN450",                  "PNN450",                                            100);
    vs_pnn->add("PNN500",                  "PNN500",                                            100);
    vs_pnn->add("PNN550",                  "PNN550",                                            100);
    vs_pnn->add("PNN600",                  "PNN600",                                            100);
    vs_pnn->add("PNN700",                  "PNN700",                                            100);
    vs_pnn->add("PNN800",                  "PNN800",                                            100);
    vs_pnn->add("PNN900",                  "PNN900",                                            100);
    vs_pnn->add("PNN1000",                 "PNN1000",                                           100);
    vs_pnn->add("PNN1200",                 "PNN1200",                                           100);
    vs_pnn->add("PNN1400",                 "PNN1400",                                           100);
    vs_pnn->add("PNN1600",                 "PNN1600",                                           100);

    Variables* vs_bdt = new Variables();
    vs_bdt->add("SMBDT",                   "SM BDT",                                            100);

    Systematics* ss = new Systematics();
    // OK the name is misleading. It is ttbar with the nominal ttbar reweighting applied
    ss->add("TTBAR_REWEIGHT_ADHOC", "t#bar{t} ad-hoc reweight", eSystematicType::TwoSide);
    ss->add("TTBAR_REWEIGHT_NNLO", "t#bar{t} NNLO reweight", eSystematicType::TwoSide);

    CompInfo* info = new CompInfo();
    info->ratio_high = 1.38;
    info->ratio_low = 0.62;
    info->shape_only = true;

    for (VariableInfo* v : *(vs_presel->content()))
    {
        Processes* ps = new Processes();
        ps->add("ttbar",            "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);
        ps->add("ttbarTF",          "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);
        ps->add("ttbarFT",          "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);
        ps->add("ttbarFF",          "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);

        Config* c = new Config(b, ps, rs, vs_presel, ss);
        c->load(filename, "Preselection");
        info->parameter = "Presel_HadHad";
        c->updateHistogramPtr(rs->content()->front(), v);
        SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/ttbarRew/";
        
        if (HistToolHelper::check(c)) {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::N_Rebin);
        }

        if (ct->check(c))
        {
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
        ps->add("ttbar",            "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);
        ps->add("ttbarTF",          "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);
        ps->add("ttbarFT",          "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);
        ps->add("ttbarFF",          "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);

        Config* c = new Config(b, ps, rs, vs_presel, ss);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "PNN_HadHad";
        c->updateHistogramPtr(rs->content()->front(), v);
        SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/ttbarRew/";
        
        if (HistToolHelper::check(c)) {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::N_Rebin);
        }

        if (ct->check(c))
        {
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
        ps->add("ttbar",            "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);
        ps->add("ttbarTF",          "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);
        ps->add("ttbarFT",          "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);
        ps->add("ttbarFF",          "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlack);

        Config* c = new Config(b, ps, rs, vs_presel, ss);
        c->load(filename, "BDTScorePreselection");
        info->parameter = "BDT_HadHad";
        c->updateHistogramPtr(rs->content()->front(), v);
        SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/ttbarRew/";
        
        if (HistToolHelper::check(c)) {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::N_Rebin);
        }

        if (ct->check(c))
        {
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
    delete info;

}