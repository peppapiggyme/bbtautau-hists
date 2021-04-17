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

void test_hadhad_yield(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);
    rs->add("1tag2pjet_0ptv_LL_OS",     "1 b-tag, 2 loose #tau, OS",        eRegionType::CR);
    rs->add("2tag2pjet_0ptv_LL_SS",     "2 b-tag, 2 loose #tau, SS",        eRegionType::CR);
    rs->add("1tag2pjet_0ptv_LL_SS",     "1 b-tag, 2 loose #tau, SS",        eRegionType::CR);

    Variables* vs = new Variables();
    vs->add("Yield",                "Event Yield",                                       1 );

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
        ps->add("ttH",                  "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ZHbb",                 "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("WHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ZHtautau",             "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("ggFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("VBFHtautau",           "SM Higgs",                     eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);
        ps->add("hhttbbKLn20p0from1p0",     "#kappa_{#lambda}=-20 (from #kappa_{#lambda}=1)",       eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=-20 (from #kappa_{#lambda}=1)",       kRed);
        ps->add("hhttbbKLn15p0from1p0",     "#kappa_{#lambda}=-15 (from #kappa_{#lambda}=1)",       eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=-15 (from #kappa_{#lambda}=1)",       kRed);
        ps->add("hhttbbKLn10p0from1p0",     "#kappa_{#lambda}=-10 (from #kappa_{#lambda}=1)",       eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=-10 (from #kappa_{#lambda}=1)",       kRed);
        ps->add("hhttbbKLn5p0from1p0",      "#kappa_{#lambda}=-5 (from #kappa_{#lambda}=1)",        eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=-5 (from #kappa_{#lambda}=1)",        kRed);
        ps->add("hhttbbKLn4p0from1p0",      "#kappa_{#lambda}=-4 (from #kappa_{#lambda}=1)",        eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=-4 (from #kappa_{#lambda}=1)",        kRed);
        ps->add("hhttbbKLn3p0from1p0",      "#kappa_{#lambda}=-3 (from #kappa_{#lambda}=1)",        eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=-3 (from #kappa_{#lambda}=1)",        kRed);
        ps->add("hhttbbKLn2p0from1p0",      "#kappa_{#lambda}=-2 (from #kappa_{#lambda}=1)",        eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=-2 (from #kappa_{#lambda}=1)",        kRed);
        ps->add("hhttbbKLn1p0from1p0",      "#kappa_{#lambda}=-1 (from #kappa_{#lambda}=1)",        eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=-1 (from #kappa_{#lambda}=1)",        kRed);
        ps->add("hhttbbKL0p0from1p0",       "#kappa_{#lambda}=0 (from #kappa_{#lambda}=1)",         eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=0 (from #kappa_{#lambda}=1)",         kRed);
        ps->add("hhttbbKL0p4from1p0",       "#kappa_{#lambda}=0.4 (from #kappa_{#lambda}=1)",       eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=0.4 (from #kappa_{#lambda}=1)",       kRed);
        ps->add("hhttbbKL0p8from1p0",       "#kappa_{#lambda}=0.8 (from #kappa_{#lambda}=1)",       eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=0.8 (from #kappa_{#lambda}=1)",       kRed);
        ps->add("hhttbbKL1p2from1p0",       "#kappa_{#lambda}=1.2 (from #kappa_{#lambda}=1)",       eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=1.2 (from #kappa_{#lambda}=1)",       kRed);
        ps->add("hhttbbKL1p6from1p0",       "#kappa_{#lambda}=1.6 (from #kappa_{#lambda}=1)",       eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=1.6 (from #kappa_{#lambda}=1)",       kRed);
        ps->add("hhttbbKL2p0from1p0",       "#kappa_{#lambda}=2 (from #kappa_{#lambda}=1)",         eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=2 (from #kappa_{#lambda}=1)",         kRed);
        ps->add("hhttbbKL2p4from1p0",       "#kappa_{#lambda}=2.4 (from #kappa_{#lambda}=1)",       eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=2.4 (from #kappa_{#lambda}=1)",       kRed);
        ps->add("hhttbbKL3p0from1p0",       "#kappa_{#lambda}=3 (from #kappa_{#lambda}=1)",         eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=3 (from #kappa_{#lambda}=1)",         kRed);
        ps->add("hhttbbKL4p0from1p0",       "#kappa_{#lambda}=4 (from #kappa_{#lambda}=1)",         eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=4 (from #kappa_{#lambda}=1)",         kRed);
        ps->add("hhttbbKL5p0from1p0",       "#kappa_{#lambda}=5 (from #kappa_{#lambda}=1)",         eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=5 (from #kappa_{#lambda}=1)",         kRed);
        ps->add("hhttbbKL10p0from1p0",      "#kappa_{#lambda}=10 (from #kappa_{#lambda}=1)",        eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=10 (from #kappa_{#lambda}=1)",        kRed);
        ps->add("hhttbbKL15p0from1p0",      "#kappa_{#lambda}=15 (from #kappa_{#lambda}=1)",        eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=15 (from #kappa_{#lambda}=1)",        kRed);
        ps->add("hhttbbKL20p0from1p0",      "#kappa_{#lambda}=20 (from #kappa_{#lambda}=1)",        eProcessType::SIG,      eProcess::HHKLXFROM1,       "#kappa_{#lambda}=20 (from #kappa_{#lambda}=1)",        kRed);
        ps->add("hhttbbKLn20p0from10p0",    "#kappa_{#lambda}=-20 (from #kappa_{#lambda}=10)",      eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=-20 (from #kappa_{#lambda}=10)",      kRed);
        ps->add("hhttbbKLn15p0from10p0",    "#kappa_{#lambda}=-15 (from #kappa_{#lambda}=10)",      eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=-15 (from #kappa_{#lambda}=10)",      kRed);
        ps->add("hhttbbKLn10p0from10p0",    "#kappa_{#lambda}=-10 (from #kappa_{#lambda}=10)",      eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=-10 (from #kappa_{#lambda}=10)",      kRed);
        ps->add("hhttbbKLn5p0from10p0",     "#kappa_{#lambda}=-5 (from #kappa_{#lambda}=10)",       eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=-5 (from #kappa_{#lambda}=10)",       kRed);
        ps->add("hhttbbKLn4p0from10p0",     "#kappa_{#lambda}=-4 (from #kappa_{#lambda}=10)",       eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=-4 (from #kappa_{#lambda}=10)",       kRed);
        ps->add("hhttbbKLn3p0from10p0",     "#kappa_{#lambda}=-3 (from #kappa_{#lambda}=10)",       eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=-3 (from #kappa_{#lambda}=10)",       kRed);
        ps->add("hhttbbKLn2p0from10p0",     "#kappa_{#lambda}=-2 (from #kappa_{#lambda}=10)",       eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=-2 (from #kappa_{#lambda}=10)",       kRed);
        ps->add("hhttbbKLn1p0from10p0",     "#kappa_{#lambda}=-1 (from #kappa_{#lambda}=10)",       eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=-1 (from #kappa_{#lambda}=10)",       kRed);
        ps->add("hhttbbKL0p0from10p0",      "#kappa_{#lambda}=0 (from #kappa_{#lambda}=10)",        eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=0 (from #kappa_{#lambda}=10)",        kRed);
        ps->add("hhttbbKL0p4from10p0",      "#kappa_{#lambda}=0.4 (from #kappa_{#lambda}=10)",      eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=0.4 (from #kappa_{#lambda}=10)",      kRed);
        ps->add("hhttbbKL0p8from10p0",      "#kappa_{#lambda}=0.8 (from #kappa_{#lambda}=10)",      eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=0.8 (from #kappa_{#lambda}=10)",      kRed);
        ps->add("hhttbbKL1p0from10p0",      "#kappa_{#lambda}=1 (from #kappa_{#lambda}=10)",        eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=1 (from #kappa_{#lambda}=10)",        kRed);
        ps->add("hhttbbKL1p2from10p0",      "#kappa_{#lambda}=1.2 (from #kappa_{#lambda}=10)",      eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=1.2 (from #kappa_{#lambda}=10)",      kRed);
        ps->add("hhttbbKL1p6from10p0",      "#kappa_{#lambda}=1.6 (from #kappa_{#lambda}=10)",      eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=1.6 (from #kappa_{#lambda}=10)",      kRed);
        ps->add("hhttbbKL2p0from10p0",      "#kappa_{#lambda}=2 (from #kappa_{#lambda}=10)",        eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=2 (from #kappa_{#lambda}=10)",        kRed);
        ps->add("hhttbbKL2p4from10p0",      "#kappa_{#lambda}=2.4 (from #kappa_{#lambda}=10)",      eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=2.4 (from #kappa_{#lambda}=10)",      kRed);
        ps->add("hhttbbKL3p0from10p0",      "#kappa_{#lambda}=3 (from #kappa_{#lambda}=10)",        eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=3 (from #kappa_{#lambda}=10)",        kRed);
        ps->add("hhttbbKL4p0from10p0",      "#kappa_{#lambda}=4 (from #kappa_{#lambda}=10)",        eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=4 (from #kappa_{#lambda}=10)",        kRed);
        ps->add("hhttbbKL5p0from10p0",      "#kappa_{#lambda}=5 (from #kappa_{#lambda}=10)",        eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=5 (from #kappa_{#lambda}=10)",        kRed);
        ps->add("hhttbbKL15p0from10p0",     "#kappa_{#lambda}=15 (from #kappa_{#lambda}=10)",       eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=15 (from #kappa_{#lambda}=10)",       kRed);
        ps->add("hhttbbKL20p0from10p0",     "#kappa_{#lambda}=20 (from #kappa_{#lambda}=10)",       eProcessType::SIG,      eProcess::HHKLXFROM10,      "#kappa_{#lambda}=20 (from #kappa_{#lambda}=10)",       kRed);
        ps->add("hhttbbKL1p0",              "#kappa_{#lambda}=1 ",          eProcessType::SIG,      eProcess::HHKL1,    "#kappa_{#lambda}=1 ",      kMagenta+2);
        ps->add("hhttbbKL10p0",             "#kappa_{#lambda}=10 ",         eProcessType::SIG,      eProcess::HHKL10,   "#kappa_{#lambda}=10 ",     kMagenta+2);

        Config* c = new Config(b, ps, rs, vs);
        c->load(filename, "Preselection");

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