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

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;

void add_bkg_binning(const Processes* ps)
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

void thesis_binning(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);

    Variables* vs_bdt = new Variables();        
    vs_bdt->add("SMBDT",              "BDT score",         1);
    vs_bdt->add("SMBDT0",             "BDT score",         1);
    vs_bdt->add("SMBDT1",             "BDT score",         1);
    vs_bdt->add("SMBDT2",             "BDT score",         1);
    vs_bdt->add("SMBDT3",             "BDT score",         1);
    vs_bdt->add("SMBDT4",             "BDT score",         1);

    AutoBinningInfo* info = new AutoBinningInfo();

    for (VariableInfo* v : *(vs_bdt->content()))
    { 
        info->parameter = "BDTScorePreselection"; 
        info->required_mcstats = 1.1; /* 0.2 to be more conservative */
        info->min_mcstats = 0.2; 
        info->min_bkg = 5; 
        info->n_bins = 100; 
        info->eff_factor = 1.0; /* only for CaseTwo */
        info->for_bdt = true;
    
        Processes* ps = new Processes(); 
        ps->add("hhttbb",      "hhttbb",      eProcessType::SIG, eProcess::HH, "SM HH", kMagenta+2); 
        ps->add("hhttbbSMVBF", "hhttbbSMVBF", eProcessType::SIG, eProcess::HH, "SM HH", kMagenta+2); 
        add_bkg_binning(ps);

        Config* c = new Config(b, ps, rs, vs_bdt); 
        c->load(filename, "BDTScorePreselection"); 
    
        c->updateHistogramPtr(rs->content()->front(), v); 
        AutoBinningTool* ab = new AutoBinningTool_v1(info, BinningCriteria::CaseTwo); 
        ab->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/var-impo-bdt/"; 
        if (ab->check(c)) 
        { 
            ab->manipulate(c); 
            ab->paint(c); 
            ab->makeYield(c, info->parameter); 
            ab->run(c); 
        } 
        else  
        { 
            clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << 'n'; 
        } 
    
        delete ps; 
        delete ab; 
        delete c; 
    }

    delete b;
    delete rs;
    delete vs_bdt;
    delete info;
}




