#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "LephadTTBarCompTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;

void test_lephad_ttbarSys(const std::string& filename)
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
    ss->add("TTBarNominal", "TTBarNominal", eSystematicType::TwoSide);
    // ttbar reweighting variations
    ss->add("TTBarReweight", "TTBarReweight", eSystematicType::TwoSide);
    ss->add("TTBarReweight_Closure_TauPt", "TTBarReweight_Closure_TauPt", eSystematicType::TwoSide);
    ss->add("TTBarReweight_Stat", "TTBarReweight_Stat", eSystematicType::TwoSide);
    // ttbar internal weights, some are missing? e.g. PDFs?
    ss->add("TTBAR_ACC_MUF", "TTBAR_ACC_MUF", eSystematicType::TwoSide);
    ss->add("TTBAR_ACC_MUR", "TTBAR_ACC_MUR", eSystematicType::TwoSide);
    ss->add("TTBAR_ACC_ISR_ALPHAS", "TTBAR_ACC_ISR_ALPHAS", eSystematicType::TwoSide);
    ss->add("TTBAR_ACC_FSR", "TTBAR_ACC_FSR", eSystematicType::TwoSide);

    CompInfo* info = new CompInfo();
    info->ratio_high = 1.38;
    info->ratio_low = 0.62;
    info->shape_only = true;

    for (VariableInfo* v : *(vs_presel->content()))
    {
        Processes* ps = new Processes();
        ps->add("ttbar",            "ttbar Nominal",        eProcessType::BKG,      eProcess::TTBAR,            "ttbar Nominal",        kBlue+1);
        ps->add("ttbarPH7",         "ttbar Herwig",         eProcessType::BKG,      eProcess::TTBARPH7,         "ttbar Herwig",         kRed-2);
        ps->add("ttbarAMCP8",       "ttbar aMC@NLO",        eProcessType::BKG,      eProcess::TTBARAMCP8,       "ttbar aMC@NLO",        kRed);
        ps->add("ttbarPP8RadHi",    "ttbar hdamp=3mtop",    eProcessType::BKG,      eProcess::TTBARPP8RADHI,    "ttbar hdamp=3mtop",    kRed+2);

        Config* c = new Config(b, ps, rs, vs_presel, ss);
        c->load(filename, "Preselection");
        info->parameter = "Presel";
        c->updateHistogramPtr(rs->content()->front(), v);
        LephadTTBarCompTool* ct = new LephadTTBarCompTool(info);
        ct->output_path = "/tmp/bowenzhang/ttbarSys";
        if (ct->check(c))
        {
            ct->manipulate(c);
            ct->combineReweightSysts(c);
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