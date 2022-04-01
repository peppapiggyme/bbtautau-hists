/**
 * This program was used for 
 * estimate the variiable ranking for MVA
 * 
 * It calculates the optimal binning for each 
 * MVA with one variable missing
 * 
 * Therefore test the importance of the
 * missing variable
 */

#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "AutoBinningTool.h"
#include "DrawStackTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;

static void add_bkg_binning(const Processes* ps)
{
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
    ps->add("ttbarTF",            "t#bar{t} TF",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
    ps->add("ttbarFT",            "t#bar{t} FT",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
    ps->add("ttbarFF",            "t#bar{t} FF",                  eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9);
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
}

void test_bukin_binning(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");
    b->name_convention = Utils::NameConvention::WSMaker;

    Regions* rs = new Regions();
    rs->add("sdfsdf",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);


    Variables* vs_bukin = new Variables();        
    vs_bukin->add("Bukin",              "Likelihood Ratio",         1);

    AutoBinningInfo* info = new AutoBinningInfo();
    info->min_mcstats = 0.2; 

    DrawStackInfo* ds_info = new DrawStackInfo();
    ds_info->atlas = false;
    ds_info->blind = false;
    ds_info->ratio_high = 1.65;
    ds_info->ratio_low = 0.35;
    ds_info->signal_scale = 1;
    ds_info->signal_linestyle = 2;
    ds_info->signal_linewidth = 3;
    ds_info->show_scaling = true;
    ds_info->legend_ncolumns = 1;
    ds_info->legend_scaling_horizontal = 1.7;
    ds_info->logy = false;
    ds_info->output_format = "png";
    ds_info->draw_overflow = false;
    ds_info->ymax_ratio_nolog = 2.4;

    for (VariableInfo* v : *(vs_bukin->content()))
    { 
        info->parameter = "semiboost"; 
    
        Processes* ps = new Processes(); 
        ps->add("data",      "data",      eProcessType::DATA, eProcess::DATA, "data", kBlack); 
        ps->add("hhttbbVBFl1cvv3cv1",      "VBF kvv=3.0",      eProcessType::SIG, eProcess::HH, "VBF kvv=3.0", kMagenta+2); 
        add_bkg_binning(ps);

        Config* c = new Config(b, ps, rs, vs_bukin); 
        c->load(filename, ""); 
    
        c->updateHistogramPtr(rs->content()->front(), v); 

        AutoBinningTool* ab = new AutoBinningTool_v1(info, BinningCriteria::CaseThree); 
        ab->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/"; 
        DrawStackTool* ds = new DrawStackTool(ds_info);
        ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/";
        if (ds->check(c) && ab->check(c))
        {
            ds->manipulate(c);
            // ds->rebin(c, eRebinOption::Array);
            ab->run(c); 
            ab->rebin(c);
            ds->paint(c);
            ds->makeYield(c);
            ds->run(c);
        }
        else 
        {
            clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
        }
    
        delete ps; 
        delete ab; 
        delete ds;
        delete c; 
    }

    delete b;
    delete rs;
    delete vs_bukin;
    delete info;
    delete ds_info;
}





