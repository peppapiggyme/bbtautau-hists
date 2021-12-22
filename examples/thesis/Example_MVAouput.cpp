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

void add_bkg_mvaoutput(const Processes* ps)
{
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
}

void thesis_mvaoutput(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("SR", "");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "",        eRegionType::SR);

    auto binningFromFile = [](const std::string& fn, bool BDT=false)
    { 
        return BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>(fn), 1000, BDT);
    };

    vector<double> binning2HDM300   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM300.txt");
    vector<double> binning2HDM500   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM500.txt");
    vector<double> binning2HDM1000  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1000.txt");
    vector<double> binning2HDM1600  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1600.txt");
    vector<double> binningSMBDT     =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_SMBDT.txt", true);

    Tools::printVector(binning2HDM300);
    Tools::printVector(binning2HDM500);
    Tools::printVector(binning2HDM1000);
    Tools::printVector(binning2HDM1600);
    Tools::printVector(binningSMBDT);

    Variables* vs_pnn300 = new Variables();
    Variables* vs_pnn500 = new Variables();
    Variables* vs_pnn1000 = new Variables();
    Variables* vs_pnn1600 = new Variables();

    int nRebin = 5;

    vs_pnn300->add("PNN300",        "PNN300 score",       nRebin); //,      &binning2HDM300[0],     binning2HDM300.size()-1);
    vs_pnn500->add("PNN500",        "PNN500 score",       nRebin); //,      &binning2HDM500[0],     binning2HDM500.size()-1);
    vs_pnn1000->add("PNN1000",      "PNN1000 score",      nRebin); //,      &binning2HDM1000[0],    binning2HDM1000.size()-1);
    vs_pnn1600->add("PNN1600",      "PNN1600 score",      nRebin); //,      &binning2HDM1600[0],    binning2HDM1600.size()-1);

    Variables* vs_bdt = new Variables();
    vs_bdt->add("SMBDT",            "BDT score",           nRebin); //,      &binningSMBDT[0],        binningSMBDT.size()-1);

    DistributionInfo* info = new DistributionInfo();
    info->atlas = false;
    info->shape_only = true;
    info->logy = true;

    for (VariableInfo* v : *(vs_pnn300->content()))
    {
        Processes* ps = new Processes();
        ps->add("Hhhbbtautau300",  "X(300)",        eProcessType::SIG,  eProcess::P1,    "X(300)",               TColor::GetColor(0, 20, 30));
        add_bkg_mvaoutput(ps);

        Config* c = new Config(b, ps, rs, vs_pnn300);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "PNN300";
        c->updateHistogramPtr(rs->content()->front(), v);
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/MVAoutput";

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


    for (VariableInfo* v : *(vs_pnn500->content()))
    {
        Processes* ps = new Processes();
        ps->add("Hhhbbtautau500",  "X(500)",        eProcessType::SIG,  eProcess::P1,    "X(500)",               TColor::GetColor(0, 20, 30));
        add_bkg_mvaoutput(ps);
        
        Config* c = new Config(b, ps, rs, vs_pnn500);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "PNN500";
        c->updateHistogramPtr(rs->content()->front(), v);
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/MVAoutput";

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


    for (VariableInfo* v : *(vs_pnn1000->content()))
    {
        Processes* ps = new Processes();
        ps->add("Hhhbbtautau1000", "X(1000)",       eProcessType::SIG,  eProcess::P1,    "X(1000)",              TColor::GetColor(0, 20, 30));
        add_bkg_mvaoutput(ps);

        Config* c = new Config(b, ps, rs, vs_pnn1000);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "PNN1000";
        c->updateHistogramPtr(rs->content()->front(), v);
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/MVAoutput";

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


    for (VariableInfo* v : *(vs_pnn1600->content()))
    {
        Processes* ps = new Processes();
        ps->add("Hhhbbtautau1600", "X(1600)",       eProcessType::SIG,  eProcess::P1,    "X(1600)",              TColor::GetColor(0, 20, 30));
        add_bkg_mvaoutput(ps);

        Config* c = new Config(b, ps, rs, vs_pnn1600);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "PNN1600";
        c->updateHistogramPtr(rs->content()->front(), v);
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/MVAoutput";

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


    for (VariableInfo* v : *(vs_bdt->content()))
    {
        Processes* ps = new Processes();
        ps->add("hhttbb",      "SM HH",             eProcessType::SIG,  eProcess::P1,    "SM HH",                TColor::GetColor(0, 20, 30));
        ps->add("hhttbbVBFSM", "SM HH",             eProcessType::SIG,  eProcess::P1,    "SM HH",                TColor::GetColor(0, 20, 30));
        add_bkg_mvaoutput(ps);

        Config* c = new Config(b, ps, rs, vs_bdt);
        c->load(filename, "BDTScorePreselection");
        info->parameter = "BDT";
        c->updateHistogramPtr(rs->content()->front(), v);
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/MVAoutput";

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
    delete vs_pnn300;
    delete vs_pnn500;
    delete vs_pnn1000;
    delete vs_pnn1600;
    delete vs_bdt;
    delete info;
}