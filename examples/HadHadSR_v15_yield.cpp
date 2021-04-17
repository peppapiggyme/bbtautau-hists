#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "DrawStackTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;

void hadhadsr_v15_yield(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);

    Variables* vs = new Variables();
    // vs->add("Yield",                "Event Yield",                                       1 );
    vs->add("SMBDT",                "SM BDT",                                           1 );

    DrawStackInfo* info = new DrawStackInfo();
    info->blind = false;
    info->ratio_high = 1.38;
    info->ratio_low = 0.62;
    info->signal_scale = 50;

    for (RegionInfo* r : *(rs->content()))
    {
        Processes* ps = new Processes();
        // Use normalisation factor ?
        ps->assign_norm_factors = true;

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
        ps->add("ttV",                  "ttV",                          eProcessType::BKG,      eProcess::ttV,          "ttV",                          kYellow+3);
        ps->add("ttH",                  "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ZHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ZHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("VBFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("hhttbb",               "Non-res HH",                   eProcessType::SIG,      eProcess::SMHH,         "Non-res HH",                   kMagenta+2);

        Config* c = new Config(b, ps, rs, vs);
        c->load(filename, "BDTScorePreselection");

        c->updateHistogramPtr(r, vs->content()->front());
        DrawStackTool* ds = new DrawStackTool(info);
        ds->output_path = "/tmp/bowenzhang/";
        if (ds->check(c))
        {
            ds->makeYield(c);
        }
        else 
        {
            clog << "Can not make yield of " << c->current_region->name << " " << c->current_variable->name << '\n';
        }

        delete ps;
        delete ds;
        delete c;
    }

    delete b;
    delete rs;
    delete vs;
    delete info;

}