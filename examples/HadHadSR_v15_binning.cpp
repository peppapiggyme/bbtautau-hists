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

#define ADDBKGTOPS() \
ps->add("Zttbb",                "Z#tau#tau + bb",         eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10); \
ps->add("Zttbc",                "Z#tau#tau + bc",         eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10); \
ps->add("Zttbl",                "Z#tau#tau + bl",         eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kBlue-5); \
ps->add("Zttcc",                "Z#tau#tau + cc",         eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10); \
ps->add("Zttcl",                "Z#tau#tau + cl",         eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kBlue-5); \
ps->add("Zttl",                 "Z#tau#tau + l",          eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kBlue-5); \
ps->add("Zbb",                  "Zll + bb",               eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10); \
ps->add("Zbc",                  "Zll + bc",               eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10); \
ps->add("Zbl",                  "Zll + bl",               eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kBlue-5); \
ps->add("Zcc",                  "Zll + cc",               eProcessType::BKG,      eProcess::ZllHF,        "Z+hf",                         kBlue-10); \
ps->add("Zcl",                  "Zll + cl",               eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kBlue-5); \
ps->add("Zl",                   "Zll + l",                eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kBlue-5); \
ps->add("ttbar",                "t#bar{t}",               eProcessType::BKG,      eProcess::TOP,          "top",                          kWhite); \
ps->add("ttbarSFTF",            "t#bar{t} TF",            eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9); \
ps->add("ttbarSFFT",            "t#bar{t} FT",            eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9); \
ps->add("ttbarSFFF",            "t#bar{t} FF",            eProcessType::BKG,      eProcess::TTBARFAKE,    "fake-#tau t#bar{t}",           kYellow-9); \
ps->add("stops",                "single top",             eProcessType::BKG,      eProcess::TOP,          "top",                          kOrange+2); \
ps->add("stopt",                "single top",             eProcessType::BKG,      eProcess::TOP,          "top",                          kOrange+2); \
ps->add("stopWt",               "single top",             eProcessType::BKG,      eProcess::TOP,          "top",                          kOrange+2); \
ps->add("WW",                   "WW",                     eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-8); \
ps->add("WZ",                   "WZ",                     eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-8); \
ps->add("ZZ",                   "ZZ",                     eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-8); \
ps->add("Wtt",                  "Others",                 eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-10); \
ps->add("W",                    "Others",                 eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kGreen-10); \
ps->add("Fake",                 "Multi-jet",              eProcessType::BKG,      eProcess::FAKE,         "Multi-jet",                    kMagenta-10); \
ps->add("ttV",                  "ttV",                    eProcessType::BKG,      eProcess::OTHERS,       "Others",                       kYellow+3); \
ps->add("ttH",                  "SM Higgs",               eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5); \
ps->add("WHbb",                 "SM Higgs",               eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5); \
ps->add("ZHbb",                 "SM Higgs",               eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5); \
ps->add("WHtautau",             "SM Higgs",               eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5); \
ps->add("ZHtautau",             "SM Higgs",               eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5); \
ps->add("ggFHtautau",           "SM Higgs",               eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5); \
ps->add("VBFHtautau",           "SM Higgs",               eProcessType::BKG,      eProcess::H,            "SM Higgs",                     kYellow-5);

#define TESTAUTOBINNING(VARSET, ROOTDIR, SIGNALNAME, NBINS) \
for (VariableInfo* v : *(VARSET->content())) \
{ \
    info->parameter = ROOTDIR; \
    info->required_mcstats = 0.2; /* 0.2 to be more conservative */\
    info->min_mcstats = 0.2; \
    info->min_bkg = 5; \
    info->n_bins = NBINS; \
    info->eff_factor = 1.0; /* only for CaseTwo */\
 \
    Processes* ps = new Processes(); \
    ADDBKGTOPS(); \
    ps->add(SIGNALNAME,             SIGNALNAME,               eProcessType::SIG,      eProcess::HH,           SIGNALNAME,                     kMagenta+2); \
 \
    Config* c = new Config(b, ps, rs, VARSET); \
    c->load(filename, ROOTDIR); \
 \
    c->updateHistogramPtr(rs->content()->front(), v); \
    AutoBinningTool* ab = new AutoBinningTool_v1(info, BinningCriteria::CaseTwo); \
    ab->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/binning/CaseTwo"; \
    if (ab->check(c)) \
    { \
        ab->manipulate(c); \
        ab->paint(c); \
        ab->makeYield(c, info->parameter); \
        ab->run(c); \
    } \
    else  \
    { \
        clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n'; \
    } \
 \
    delete ps; \
    delete ab; \
    delete c; \
}


void hadhadsr_v15_binning(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);

    Variables* vs_bdt = new Variables();        vs_bdt->add("SMBDT",              "BDT Score",         1);
    Variables* vs_260 = new Variables();        vs_260->add("PNN260",             "PNN260",            1);
    Variables* vs_280 = new Variables();        vs_280->add("PNN280",             "PNN280",            1);
    Variables* vs_300 = new Variables();        vs_300->add("PNN300",             "PNN300",            1);
    Variables* vs_325 = new Variables();        vs_325->add("PNN325",             "PNN325",            1);
    Variables* vs_350 = new Variables();        vs_350->add("PNN350",             "PNN350",            1);
    Variables* vs_400 = new Variables();        vs_400->add("PNN400",             "PNN400",            1);
    Variables* vs_450 = new Variables();        vs_450->add("PNN450",             "PNN450",            1);
    Variables* vs_500 = new Variables();        vs_500->add("PNN500",             "PNN500",            1);
    Variables* vs_550 = new Variables();        vs_550->add("PNN550",             "PNN550",            1);
    Variables* vs_600 = new Variables();        vs_600->add("PNN600",             "PNN600",            1);
    Variables* vs_700 = new Variables();        vs_700->add("PNN700",             "PNN700",            1);
    Variables* vs_800 = new Variables();        vs_800->add("PNN800",             "PNN800",            1);
    Variables* vs_900 = new Variables();        vs_900->add("PNN900",             "PNN900",            1);
    Variables* vs_1000 = new Variables();       vs_1000->add("PNN1000",           "PNN1000",           1);
    Variables* vs_1100 = new Variables();       vs_1100->add("PNN1100",           "PNN1100",           1);
    Variables* vs_1200 = new Variables();       vs_1200->add("PNN1200",           "PNN1200",           1);
    Variables* vs_1400 = new Variables();       vs_1400->add("PNN1400",           "PNN1400",           1);
    Variables* vs_1600 = new Variables();       vs_1600->add("PNN1600",           "PNN1600",           1);

    AutoBinningInfo* info = new AutoBinningInfo();

    TESTAUTOBINNING(vs_bdt,     "BDTScorePreselection",     "hhttbb",                   20);
    TESTAUTOBINNING(vs_260,     "PNNScorePreselection",     "Hhhbbtautau260",           20);
    TESTAUTOBINNING(vs_280,     "PNNScorePreselection",     "Hhhbbtautau280",           20);
    TESTAUTOBINNING(vs_300,     "PNNScorePreselection",     "Hhhbbtautau300",           20);
    TESTAUTOBINNING(vs_325,     "PNNScorePreselection",     "Hhhbbtautau325",           20);
    TESTAUTOBINNING(vs_350,     "PNNScorePreselection",     "Hhhbbtautau350",           20);
    TESTAUTOBINNING(vs_400,     "PNNScorePreselection",     "Hhhbbtautau400",           20);
    TESTAUTOBINNING(vs_450,     "PNNScorePreselection",     "Hhhbbtautau450",           20);
    TESTAUTOBINNING(vs_500,     "PNNScorePreselection",     "Hhhbbtautau500",           20);
    TESTAUTOBINNING(vs_550,     "PNNScorePreselection",     "Hhhbbtautau550",           14);
    TESTAUTOBINNING(vs_600,     "PNNScorePreselection",     "Hhhbbtautau600",           14);
    TESTAUTOBINNING(vs_700,     "PNNScorePreselection",     "Hhhbbtautau700",           14);
    TESTAUTOBINNING(vs_800,     "PNNScorePreselection",     "Hhhbbtautau800",            8);
    TESTAUTOBINNING(vs_900,     "PNNScorePreselection",     "Hhhbbtautau900",            8);
    TESTAUTOBINNING(vs_1000,    "PNNScorePreselection",    "Hhhbbtautau1000",            8);
    TESTAUTOBINNING(vs_1100,    "PNNScorePreselection",    "Hhhbbtautau1100",            8);
    TESTAUTOBINNING(vs_1200,    "PNNScorePreselection",    "Hhhbbtautau1200",            8);
    TESTAUTOBINNING(vs_1400,    "PNNScorePreselection",    "Hhhbbtautau1400",            8);
    TESTAUTOBINNING(vs_1600,    "PNNScorePreselection",    "Hhhbbtautau1600",            8);

    delete b;
    delete rs;
    delete vs_bdt;
    delete vs_260;
    delete vs_280;
    delete vs_300;
    delete vs_325;
    delete vs_350;
    delete vs_400;
    delete vs_450;
    delete vs_500;
    delete vs_550;
    delete vs_600;
    delete vs_700;
    delete vs_800;
    delete vs_900;
    delete vs_1000;
    delete vs_1100;
    delete vs_1200;
    delete vs_1400;
    delete vs_1600;
    delete info;

}