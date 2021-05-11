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

void test_ttbarKine(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);

    Variables* vs_presel = new Variables();
    vs_presel->add("mBB",                  "m_{BB} [GeV]",                                      10);
    vs_presel->add("mMMC",                 "m_{#tau#tau} (MMC) [GeV]",                          10);
    vs_presel->add("dRBB",                 "#Delta R(B,B)",                                     5);
    vs_presel->add("dRTauTau",             "#Delta R(#tau_{had},#tau_{had})",                   5);
    vs_presel->add("pTB0",                 "Leading b-jet p_{T} [GeV]",                         5);
    vs_presel->add("pTB1",                 "Sub-leading b-jet p_{T} [GeV]",                     5);
    vs_presel->add("Tau0Pt",               "Leading #tau_{had} p_{T} [GeV]",                    5);
    vs_presel->add("Tau1Pt",               "Sub-leading #tau_{had} p_{T} [GeV]",                5);
    vs_presel->add("HT",                   "Sum of pT (jet) [GeV]",                             2);
    vs_presel->add("sT",                   "Sum of pT (all) [GeV]",                             2);

    CompInfo* info = new CompInfo();
    info->ratio_high = 2;
    info->ratio_low = 0;
    info->shape_only = true;
    info->save_ratio = true;

    // AutoBinningInfo* abi = new AutoBinningInfo();
    // abi->n_bins = 20;
    // AutoBinningTool* abt = new AutoBinningTool_v2(abi, eProcess::TTBARHH);

    for (VariableInfo* v : *(vs_presel->content()))
    {
        Processes* ps = new Processes();
        ps->add("ttbarFromSLT",      "ttbar (lephad SLT)",   eProcessType::BKG,   eProcess::TTBARLH,     "ttbar (lephad)",      kRed+1);
        ps->add("ttbarFromLTT",      "ttbar (lephad LTT)",   eProcessType::BKG,   eProcess::TTBARLH,     "ttbar (lephad)",      kOrange+1);
        ps->add("ttbar",            "ttbar (hadhad)",   eProcessType::BKG,   eProcess::TTBARHH,     "ttbar (hadhad)",      kBlue);
        ps->add("ttbarFromZCR",     "ttbar (ZCR)",      eProcessType::BKG,   eProcess::TTBARZCR,    "ttbar (ZCR)",         kBlack);

        Config* c = new Config(b, ps, rs, vs_presel);
        c->load(filename, "Preselection");
        info->parameter = "ttbar";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/ttbar";
        if (ct->check(c))
        {
            // abt->run(c);
            // abt->rebin(c);
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
}