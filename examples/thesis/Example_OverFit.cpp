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


void add_bkg_overfit(const Processes* ps, int color)
{
    ps->add("Zttbb",                "Z#tau#tau + bb",               eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zttbc",                "Z#tau#tau + bc",               eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zttbl",                "Z#tau#tau + bl",               eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zttcc",                "Z#tau#tau + cc",               eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zttcl",                "Z#tau#tau + cl",               eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zttl",                 "Z#tau#tau + l",                eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zbb",                  "Zll + bb",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zbc",                  "Zll + bc",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zbl",                  "Zll + bl",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zcc",                  "Zll + cc",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zcl",                  "Zll + cl",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Zl",                   "Zll + l",                      eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ttbar",                "t#bar{t}",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ttbarSFTF",            "t#bar{t} TF",                  eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ttbarSFFT",            "t#bar{t} FT",                  eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ttbarSFFF",            "t#bar{t} FF",                  eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("stops",                "single top",                   eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("stopt",                "single top",                   eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("stopWt",               "single top",                   eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("WW",                   "WW",                           eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("WZ",                   "WZ",                           eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ZZ",                   "ZZ",                           eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Wtt",                  "W+jets",                       eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("W",                    "W+jets",                       eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("Fake",                 "Multi-jet",                    eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ttW",                  "ttV",                          eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ttZ",                  "ttV",                          eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ttH",                  "SM Higgs",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("WHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ggZHbb",               "SM Higgs",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("qqZHbb",               "SM Higgs",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("WHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ggZHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("qqZHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("ggFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
    ps->add("VBFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::P10,      "Total Bkg",       color);
}

void thesis_overfit(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("SR", "");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "",        eRegionType::SR);

    auto binningFromFile = [](const std::string& fn, bool BDT=false)
    { 
        return BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>(fn), 1000, BDT);
    };

    Variables* vs_pnn300 = new Variables();
    Variables* vs_pnn500 = new Variables();
    Variables* vs_pnn1000 = new Variables();
    Variables* vs_pnn1600 = new Variables();

    int nRebin = 5;

    vs_pnn300->add("PNN300Train",        "PNN300 score",       nRebin);
    vs_pnn300->add("PNN300Test",        "PNN300 score",       nRebin);
    vs_pnn500->add("PNN500Train",        "PNN500 score",       nRebin);
    vs_pnn500->add("PNN500Test",        "PNN500 score",       nRebin);
    vs_pnn1000->add("PNN1000Train",      "PNN1000 score",      nRebin);
    vs_pnn1000->add("PNN1000Test",      "PNN1000 score",      nRebin);
    vs_pnn1600->add("PNN1600Train",      "PNN1600 score",      nRebin);
    vs_pnn1600->add("PNN1600Test",      "PNN1600 score",      nRebin);

    Variables* vs_bdt = new Variables();
    vs_bdt->add("SMBDTTrain",            "BDT score",           nRebin);
    vs_bdt->add("SMBDTTest",            "BDT score",           nRebin);

    DistributionInfo* info = new DistributionInfo();
    info->atlas = false;
    info->shape_only = true;
    info->logy = true;

    {
        Processes* ps_train = new Processes();
        Processes* ps_test = new Processes();
        ps_train->add("Hhhbbtautau300",  "X(300)",        eProcessType::SIG,  eProcess::P1,    "X(300)",         TColor::GetColor(196, 110, 57));
        ps_test ->add("Hhhbbtautau300",  "X(300)",        eProcessType::SIG,  eProcess::P1,    "X(300)",         TColor::GetColor(120, 46, 0));
        add_bkg_overfit(ps_train, TColor::GetColor(73, 163, 166));
        add_bkg_overfit(ps_test,  TColor::GetColor(19, 83, 84));

        Config* c_train = new Config(b, ps_train, rs, vs_pnn300);
        Config* c_test = new Config(b, ps_test, rs, vs_pnn300);
        c_train->load(filename, "PNNScoreOvertrainingPreselection");
        c_test->load(filename, "PNNScoreOvertrainingPreselection");
        info->parameter = "PNN300";
        c_train->updateHistogramPtr(rs->content()->front(), vs_pnn300->content()->at(0)); // 0->Train
        c_test->updateHistogramPtr(rs->content()->front(), vs_pnn300->content()->at(1)); // 1->Test
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/OverFit";

        if (HistToolHelper::check(c_train) && HistToolHelper::check(c_test)) {
            ct->manipulate(c_train);
            ct->manipulate(c_test);
            ct->rebin(c_train, eRebinOption::N_Rebin);
            ct->rebin(c_test, eRebinOption::N_Rebin);
        }

        if (ct->check(c_train) && ct->check(c_test))
        {
            ct->paint(c_train);
            ct->paint(c_test);
            ct->run_overfit(c_train, c_test);
        }
        else 
        {
            clog << "Can not draw " << c_train->current_region->name << " " << c_train->current_variable->name << '\n';
        }

        delete ps_train;
        delete ps_test;
        delete ct;
        delete c_train;
        delete c_test;
    }

    {
        Processes* ps_train = new Processes();
        Processes* ps_test = new Processes();
        ps_train->add("Hhhbbtautau500",  "X(500)",        eProcessType::SIG,  eProcess::P1,    "X(500)",         TColor::GetColor(196, 110, 57));
        ps_test ->add("Hhhbbtautau500",  "X(500)",        eProcessType::SIG,  eProcess::P1,    "X(500)",         TColor::GetColor(120, 46, 0));
        add_bkg_overfit(ps_train, TColor::GetColor(73, 163, 166));
        add_bkg_overfit(ps_test,  TColor::GetColor(19, 83, 84));

        Config* c_train = new Config(b, ps_train, rs, vs_pnn500);
        Config* c_test = new Config(b, ps_test, rs, vs_pnn500);
        c_train->load(filename, "PNNScoreOvertrainingPreselection");
        c_test->load(filename, "PNNScoreOvertrainingPreselection");
        info->parameter = "PNN500";
        c_train->updateHistogramPtr(rs->content()->front(), vs_pnn500->content()->at(0)); // 0->Train
        c_test->updateHistogramPtr(rs->content()->front(), vs_pnn500->content()->at(1)); // 1->Test
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/OverFit";

        if (HistToolHelper::check(c_train) && HistToolHelper::check(c_test)) {
            ct->manipulate(c_train);
            ct->manipulate(c_test);
            ct->rebin(c_train, eRebinOption::N_Rebin);
            ct->rebin(c_test, eRebinOption::N_Rebin);
        }

        if (ct->check(c_train) && ct->check(c_test))
        {
            ct->paint(c_train);
            ct->paint(c_test);
            ct->run_overfit(c_train, c_test);
        }
        else 
        {
            clog << "Can not draw " << c_train->current_region->name << " " << c_train->current_variable->name << '\n';
        }

        delete ps_train;
        delete ps_test;
        delete ct;
        delete c_train;
        delete c_test;
    }

    {
        Processes* ps_train = new Processes();
        Processes* ps_test = new Processes();
        ps_train->add("Hhhbbtautau1000",  "X(1000)",        eProcessType::SIG,  eProcess::P1,    "X(1000)",         TColor::GetColor(196, 110, 57));
        ps_test ->add("Hhhbbtautau1000",  "X(1000)",        eProcessType::SIG,  eProcess::P1,    "X(1000)",         TColor::GetColor(120, 46, 0));
        add_bkg_overfit(ps_train, TColor::GetColor(73, 163, 166));
        add_bkg_overfit(ps_test,  TColor::GetColor(19, 83, 84));

        Config* c_train = new Config(b, ps_train, rs, vs_pnn1000);
        Config* c_test = new Config(b, ps_test, rs, vs_pnn1000);
        c_train->load(filename, "PNNScoreOvertrainingPreselection");
        c_test->load(filename, "PNNScoreOvertrainingPreselection");
        info->parameter = "PNN1000";
        c_train->updateHistogramPtr(rs->content()->front(), vs_pnn1000->content()->at(0)); // 0->Train
        c_test->updateHistogramPtr(rs->content()->front(), vs_pnn1000->content()->at(1)); // 1->Test
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/OverFit";

        if (HistToolHelper::check(c_train) && HistToolHelper::check(c_test)) {
            ct->manipulate(c_train);
            ct->manipulate(c_test);
            ct->rebin(c_train, eRebinOption::N_Rebin);
            ct->rebin(c_test, eRebinOption::N_Rebin);
        }

        if (ct->check(c_train) && ct->check(c_test))
        {
            ct->paint(c_train);
            ct->paint(c_test);
            ct->run_overfit(c_train, c_test);
        }
        else 
        {
            clog << "Can not draw " << c_train->current_region->name << " " << c_train->current_variable->name << '\n';
        }

        delete ps_train;
        delete ps_test;
        delete ct;
        delete c_train;
        delete c_test;
    }

    {
        Processes* ps_train = new Processes();
        Processes* ps_test = new Processes();
        ps_train->add("Hhhbbtautau1600",  "X(1600)",        eProcessType::SIG,  eProcess::P1,    "X(1600)",         TColor::GetColor(196, 110, 57));
        ps_test ->add("Hhhbbtautau1600",  "X(1600)",        eProcessType::SIG,  eProcess::P1,    "X(1600)",         TColor::GetColor(120, 46, 0));
        add_bkg_overfit(ps_train, TColor::GetColor(73, 163, 166));
        add_bkg_overfit(ps_test,  TColor::GetColor(19, 83, 84));

        Config* c_train = new Config(b, ps_train, rs, vs_pnn1600);
        Config* c_test = new Config(b, ps_test, rs, vs_pnn1600);
        c_train->load(filename, "PNNScoreOvertrainingPreselection");
        c_test->load(filename, "PNNScoreOvertrainingPreselection");
        info->parameter = "PNN1600";
        c_train->updateHistogramPtr(rs->content()->front(), vs_pnn1600->content()->at(0)); // 0->Train
        c_test->updateHistogramPtr(rs->content()->front(), vs_pnn1600->content()->at(1)); // 1->Test
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/OverFit";

        if (HistToolHelper::check(c_train) && HistToolHelper::check(c_test)) {
            ct->manipulate(c_train);
            ct->manipulate(c_test);
            ct->rebin(c_train, eRebinOption::N_Rebin);
            ct->rebin(c_test, eRebinOption::N_Rebin);
        }

        if (ct->check(c_train) && ct->check(c_test))
        {
            ct->paint(c_train);
            ct->paint(c_test);
            ct->run_overfit(c_train, c_test);
        }
        else 
        {
            clog << "Can not draw " << c_train->current_region->name << " " << c_train->current_variable->name << '\n';
        }

        delete ps_train;
        delete ps_test;
        delete ct;
        delete c_train;
        delete c_test;
    }

    {
        Processes* ps_train = new Processes();
        Processes* ps_test = new Processes();
        ps_train->add("hhttbb",       "SM HH",        eProcessType::SIG,  eProcess::P1,    "SM HH",         TColor::GetColor(196, 110, 57));
        ps_train->add("hhttbbSMVBF",  "SM HH",        eProcessType::SIG,  eProcess::P1,    "SM HH",         TColor::GetColor(196, 110, 57));
        ps_test ->add("hhttbb",       "SM HH",        eProcessType::SIG,  eProcess::P1,    "SM HH",         TColor::GetColor(120, 46, 0));
        ps_test ->add("hhttbbSMVBF",  "SM HH",        eProcessType::SIG,  eProcess::P1,    "SM HH",         TColor::GetColor(120, 46, 0));
        add_bkg_overfit(ps_train, TColor::GetColor(73, 163, 166));
        add_bkg_overfit(ps_test,  TColor::GetColor(19, 83, 84));

        Config* c_train = new Config(b, ps_train, rs, vs_bdt);
        Config* c_test = new Config(b, ps_test, rs, vs_bdt);
        c_train->load(filename, "BDTScoreOvertrainingPreselection");
        c_test->load(filename, "BDTScoreOvertrainingPreselection");
        info->parameter = "SMBDT";
        c_train->updateHistogramPtr(rs->content()->front(), vs_bdt->content()->at(0)); // 0->Train
        c_test->updateHistogramPtr(rs->content()->front(), vs_bdt->content()->at(1)); // 1->Test
        DistributionTool* ct = new DistributionTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/OverFit";

        if (HistToolHelper::check(c_train) && HistToolHelper::check(c_test)) {
            ct->manipulate(c_train);
            ct->manipulate(c_test);
            ct->rebin(c_train, eRebinOption::N_Rebin);
            ct->rebin(c_test, eRebinOption::N_Rebin);
        }

        if (ct->check(c_train) && ct->check(c_test))
        {
            ct->paint(c_train);
            ct->paint(c_test);
            ct->run_overfit(c_train, c_test);
        }
        else 
        {
            clog << "Can not draw " << c_train->current_region->name << " " << c_train->current_variable->name << '\n';
        }

        delete ps_train;
        delete ps_test;
        delete ct;
        delete c_train;
        delete c_test;
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