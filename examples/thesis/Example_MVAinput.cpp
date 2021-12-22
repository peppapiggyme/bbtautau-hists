#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "DistributionTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;
using BU = BinningUtils;

void thesis_mvainput(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("SR", "");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "",        eRegionType::SR);

    Variables* vs_presel = new Variables();
    vs_presel->add("mBB",                  "m_{bb} [GeV]",                           10);
    vs_presel->add("mMMC",                 "m_{#tau#tau}^{MMC} [GeV]",               10);
    vs_presel->add("mHH",                  "m_{HH} [GeV]",                           20);
    vs_presel->add("dRBB",                 "#Delta R(b, b)",                         4);
    vs_presel->add("dRTauTau",             "#Delta R(#tau, #tau)",                   4);

    DistributionInfo* info = new DistributionInfo();
    info->atlas = false;
    info->shape_only = true;

    for (VariableInfo* v : *(vs_presel->content()))
    {
        Processes* ps = new Processes();
        ps->add("hhttbb",      "SM HH",             eProcessType::SIG,  eProcess::P1,    "SM HH",                TColor::GetColor(0, 20, 30));
        ps->add("hhttbbVBFSM", "SM HH",             eProcessType::SIG,  eProcess::P1,    "SM HH",                TColor::GetColor(0, 20, 30));
        ps->add("Hhhbbtautau300",  "X(300)",        eProcessType::SIG,  eProcess::P2,    "X(300)",               TColor::GetColor(31, 95, 102));
        ps->add("Hhhbbtautau1000", "X(1000)",       eProcessType::SIG,  eProcess::P3,    "X(1000)",              TColor::GetColor(93, 147, 153));
        ps->add("Zttbb",    "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zttbc",    "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zttcc",    "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zttbl",    "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zttcl",    "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zttl",     "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zbb",      "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zbc",      "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zbl",      "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zcc",      "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zcl",      "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("Zl",       "Z + jets",             eProcessType::BKG,  eProcess::P6,    "Z + jets",             TColor::GetColor(105, 103, 255));
        ps->add("ttbar",    "Top-quark",            eProcessType::BKG,  eProcess::P8,    "Top-quark",            TColor::GetColor(218, 87, 235));
        ps->add("ttbarTF",  "Top-quark",            eProcessType::BKG,  eProcess::P8,    "Top-quark",            TColor::GetColor(218, 87, 235));
        ps->add("ttbarFT",  "Top-quark",            eProcessType::BKG,  eProcess::P8,    "Top-quark",            TColor::GetColor(218, 87, 235));
        ps->add("ttbarFF",  "Top-quark",            eProcessType::BKG,  eProcess::P8,    "Top-quark",            TColor::GetColor(218, 87, 235));
        ps->add("stops",    "Top-quark",            eProcessType::BKG,  eProcess::P8,    "Top-quark",            TColor::GetColor(218, 87, 235));
        ps->add("stopt",    "Top-quark",            eProcessType::BKG,  eProcess::P8,    "Top-quark",            TColor::GetColor(218, 87, 235));
        ps->add("stopWt",   "Top-quark",            eProcessType::BKG,  eProcess::P8,    "Top-quark",            TColor::GetColor(218, 87, 235));
        ps->add("ttH",        "SM Higgs",           eProcessType::BKG,  eProcess::P10,   "Single Higgs",         TColor::GetColor(250, 146, 132));
        ps->add("WHbb",       "SM Higgs",           eProcessType::BKG,  eProcess::P10,   "Single Higgs",         TColor::GetColor(250, 146, 132));
        ps->add("ggZHbb",     "SM Higgs",           eProcessType::BKG,  eProcess::P10,   "Single Higgs",         TColor::GetColor(250, 146, 132));
        ps->add("qqZHbb",     "SM Higgs",           eProcessType::BKG,  eProcess::P10,   "Single Higgs",         TColor::GetColor(250, 146, 132));
        ps->add("WHtautau",   "SM Higgs",           eProcessType::BKG,  eProcess::P10,   "Single Higgs",         TColor::GetColor(250, 146, 132));
        ps->add("ggZHtautau", "SM Higgs",           eProcessType::BKG,  eProcess::P10,   "Single Higgs",         TColor::GetColor(250, 146, 132));
        ps->add("qqZHtautau", "SM Higgs",           eProcessType::BKG,  eProcess::P10,   "Single Higgs",         TColor::GetColor(250, 146, 132));
        ps->add("ggFHtautau", "SM Higgs",           eProcessType::BKG,  eProcess::P10,   "Single Higgs",         TColor::GetColor(250, 146, 132));
        ps->add("VBFHtautau", "SM Higgs",           eProcessType::BKG,  eProcess::P10,   "Single Higgs",         TColor::GetColor(250, 146, 132));
        ps->add("Fake",     "Multi-jet",            eProcessType::BKG,  eProcess::P12,   "Multi-jet",            TColor::GetColor(170, 209, 250));

        Config* c = new Config(b, ps, rs, vs_presel);
        c->load(filename, "Preselection");
        info->parameter = "Presel";
        c->updateHistogramPtr(rs->content()->front(), v);
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/MVAinput";

        if (HistToolHelper::check(c)) {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::N_Rebin);
        }

        if (ct->check(c))
        {
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