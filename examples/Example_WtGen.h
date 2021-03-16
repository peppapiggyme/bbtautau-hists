#include "Config.h"
#include "Utils.h"
#include "CompTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;

void test_hadhad_WtGen(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);

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

    CompInfo* info = new CompInfo();
    info->ratio_high = 1.78;
    info->ratio_low = 0.22;
    info->shape_only = false;

    for (VariableInfo* v : *(vs_presel->content()))
    {
        Processes* ps = new Processes();
        ps->add("stopWtAF2",       "s-top Wt Nominal",  eProcessType::BKG,  eProcess::STOPWT,      "s-top Wt Nominal",  kBlue+1);
        ps->add("stopWtaMCAF2",    "s-top Wt aMC@NLO",  eProcessType::BKG,  eProcess::STOPWT_AMC,  "s-top Wt aMC@NLO",  kRed+1);
        ps->add("stopWtHw7AF2",    "s-top Wt Herwig7",  eProcessType::BKG,  eProcess::STOPWT_HW7,  "s-top Wt Herwig7",  kViolet+1);

        Config* c = new Config(b, ps, rs, vs_presel);
        c->load(filename, "Preselection");
        info->parameter = "Wt_Gen_Presel";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/tmp/zhangbw/Stop";
        if (ct->check(c))
        {
            ct->manipulate(c);
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
        ps->add("stopWtAF2",       "s-top Wt Nominal",  eProcessType::BKG,  eProcess::STOPWT,      "s-top Wt Nominal",  kBlue+1);
        ps->add("stopWtaMCAF2",    "s-top Wt aMC@NLO",  eProcessType::BKG,  eProcess::STOPWT_AMC,  "s-top Wt aMC@NLO",  kRed+1);
        ps->add("stopWtHw7AF2",    "s-top Wt Herwig7",  eProcessType::BKG,  eProcess::STOPWT_HW7,  "s-top Wt Herwig7",  kViolet+1);

        Config* c = new Config(b, ps, rs, vs_pnn);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "Wt_Gen_PNN";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/tmp/zhangbw/Stop";
        if (ct->check(c))
        {
            ct->manipulate(c);
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
        ps->add("stopWtAF2",       "s-top Wt Nominal",  eProcessType::BKG,  eProcess::STOPWT,      "s-top Wt Nominal",  kBlue+1);
        ps->add("stopWtaMCAF2",    "s-top Wt aMC@NLO",  eProcessType::BKG,  eProcess::STOPWT_AMC,  "s-top Wt aMC@NLO",  kRed+1);
        ps->add("stopWtHw7AF2",    "s-top Wt Herwig7",  eProcessType::BKG,  eProcess::STOPWT_HW7,  "s-top Wt Herwig7",  kViolet+1);

        Config* c = new Config(b, ps, rs, vs_bdt);
        c->load(filename, "BDTScorePreselection");
        info->parameter = "Wt_Gen_BDT";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/tmp/zhangbw/Stop";
        if (ct->check(c))
        {
            ct->manipulate(c);
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