#include "Config.h"
#include "Utils.h"
#include "CompTool.h"

#include "ExamplesInclude.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;

void test_hadhad_NonRes(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);

    Variables* vs_presel = new Variables();
    double binning_mass[17] = {0, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300};
    vs_presel->add("mBB",                  "m_{BB} [GeV]",                                      1, binning_mass, 16);
    vs_presel->add("mMMC",                 "m_{#tau#tau} (MMC) [GeV]",                          1, binning_mass, 16);
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

    Variables* vs_bdt = new Variables();
    double binning_bdt[14] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0};
    vs_bdt->add("SMBDT",                   "SM BDT Score",                                      50, binning_bdt, 13);

    CompInfo* info = new CompInfo();
    info->ratio_high = 1.48;
    info->ratio_low = 0.52;
    info->shape_only = true;
    info->ratio_tex = "Herwig7 / Pythia8";

    for (VariableInfo* v : *(vs_presel->content()))
    {
        Processes* ps = new Processes();
        // ps->add("hhttbb",       "Pythia8",  eProcessType::SIG,  eProcess::HH,    "Pythia8",  kBlue+1);
        // ps->add("hhttbbPhH7",   "Herwig7",  eProcessType::SIG,  eProcess::HH,    "Herwig7",  kRed+1);
        ps->add("hhttbbL10",       "Pythia8",  eProcessType::SIG,  eProcess::HH,    "Pythia8",  kBlue+1);
        ps->add("hhttbbL10PhH7",   "Herwig7",  eProcessType::SIG,  eProcess::HH,    "Herwig7",  kRed+1);

        Config* c = new Config(b, ps, rs, vs_presel);
        c->load(filename, "Preselection");
        // info->parameter = "Presel_ShapeOnly";
        info->parameter = "Presel_l10_ShapeOnly";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/NonRes";
        if (ct->check(c))
        {
            ct->makeYield(c, info->parameter);
            ct->rebin(c, eRebinOption::Array);
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
    // info->auto_binning = true;

    for (VariableInfo* v : *(vs_bdt->content()))
    {
        Processes* ps = new Processes();
        // ps->add("hhttbb",       "Pythia8",  eProcessType::SIG,  eProcess::HH,    "Pythia8",  kBlue+1);
        // ps->add("hhttbbPhH7",   "Herwig7",  eProcessType::SIG,  eProcess::HH,    "Herwig7",  kRed+1);
        ps->add("hhttbbL10",       "Pythia8",  eProcessType::SIG,  eProcess::HH,    "Pythia8",  kBlue+1);
        ps->add("hhttbbL10PhH7",   "Herwig7",  eProcessType::SIG,  eProcess::HH,    "Herwig7",  kRed+1);

        Config* c = new Config(b, ps, rs, vs_bdt);
        c->load(filename, "BDTScorePreselection");
        // info->parameter = "BDT_ShapeOnly";
        info->parameter = "BDT_l10_ShapeOnly";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/NonRes";
        if (ct->check(c))
        {
            ct->makeYield(c, info->parameter);
            ct->rebin(c, eRebinOption::Array);
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
    delete vs_bdt;
    delete info;

}
