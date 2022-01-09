#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "DrawStackTool.h"
#include "AutoBinningTool.h"
#include "Colors.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;
using BU = BinningUtils;

void thesis_ffcr(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");
    b->name_convention = Utils::NameConvention::CxAODReader;

    Regions* rs = new Regions();
    rs->add("1tag2pjet_0ptv_LL_SS",     "1-b-tagged SS ID Region",         eRegionType::SR);
    rs->add("1tag2pjet_0ptv_1A_SS",     "1-b-tagged SS Anti-ID Region",    eRegionType::SR);
    rs->add("1tag2pjet_0ptv_LL_OS",     "1-b-tagged OS ID Region",         eRegionType::SR);
    rs->add("1tag2pjet_0ptv_1A_OS",     "1-b-tagged OS Anti-ID Region",    eRegionType::SR);

    vector<double> binning_1p{30, 35, 40, 45, 50, 60, 75, 1000};
    vector<double> binning_3p{30, 40, 50, 65, 1000};

    Variables* vs = new Variables();
    vs->add("TauPt1P",              "1-prong #tau_{had} p_{T} [GeV]",                       1,          &binning_1p[0],         binning_1p.size()-1);
    vs->add("TauPt3P",              "3-prong #tau_{had} p_{T} [GeV]",                       1,          &binning_3p[0],         binning_3p.size()-1);

    DrawStackInfo* info = new DrawStackInfo();
    info->atlas = false;
    info->blind = false;
    info->show_scaling = true;
    info->legend_ncolumns = 1;
    info->legend_scaling_horizontal = 1.7;
    info->legend_scaling_vertical = 0.85;
    info->logx = true;
    info->logy = true;
    info->output_format = "pdf";
    info->draw_overflow = false;
    info->draw_ratio = false;
    info->xmax = 1000;

    for (RegionInfo* r : *(rs->content()))
    {
        for (VariableInfo* v : *(vs->content()))
        {
            Processes* ps = new Processes();
            ps->add("data",       "data",               eProcessType::DATA,     eProcess::DATA,         "Data",                             kBlack);
            ps->add("ttbar",      "t#bar{t}",           eProcessType::BKG,      eProcess::P2,           "True-#tau_{had} t#bar{t}",         bbtt_kGOLDEN);
            ps->add("ttbarTF",    "t#bar{t} TF",        eProcessType::BKG,      eProcess::P1,           "Fake-#tau_{had} t#bar{t}",         bbtt_kORANGE);
            ps->add("ttbarFT",    "t#bar{t} FT",        eProcessType::BKG,      eProcess::P1,           "Fake-#tau_{had} t#bar{t}",         bbtt_kORANGE);
            ps->add("ttbarFF",    "t#bar{t} FF",        eProcessType::BKG,      eProcess::P1,           "Fake-#tau_{had} t#bar{t}",         bbtt_kORANGE);
            ps->add("Zttbl",      "Z#tau#tau + bl",     eProcessType::BKG,      eProcess::P4,           "Z #rightarrow #tau#tau + lf",      bbtt_kBLUE_D);
            ps->add("Zttcl",      "Z#tau#tau + cl",     eProcessType::BKG,      eProcess::P4,           "Z #rightarrow #tau#tau + lf",      bbtt_kBLUE_D);
            ps->add("Zttl",       "Z#tau#tau + l",      eProcessType::BKG,      eProcess::P4,           "Z #rightarrow #tau#tau + lf",      bbtt_kBLUE_D);
            ps->add("Wtt",        "W+jets",             eProcessType::BKG,      eProcess::P6,           "W #rightarrow #tau#nu",            bbtt_kPURPLE);
            ps->add("stops",      "single top",         eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("stopt",      "single top",         eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("stopWt",     "single top",         eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zttbb",      "Z#tau#tau + bb",     eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zttbc",      "Z#tau#tau + bc",     eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zttcc",      "Z#tau#tau + cc",     eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zbb",        "Zll + bb",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zbc",        "Zll + bc",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zbl",        "Zll + bl",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zcc",        "Zll + cc",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zcl",        "Zll + cl",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zl",         "Zll + l",            eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("WW",         "WW",                 eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("WZ",         "WZ",                 eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ZZ",         "ZZ",                 eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("W",          "W+jets",             eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ttW",        "ttV",                eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ttZ",        "ttV",                eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ttH",        "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("WHbb",       "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ggZHbb",     "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("qqZHbb",     "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("WHtautau",   "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ggZHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("qqZHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ggFHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("VBFHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);

            Config* c = new Config(b, ps, rs, vs);
            c->load(filename, "FF");
            info->parameter = "FF";

            c->updateHistogramPtr(r, v);
            DrawStackTool* ds = new DrawStackTool(info);
            ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/FFCR/";
            if (ds->check(c))
            {
                ds->manipulate(c);
                ds->rebin(c, eRebinOption::Array);
                ds->paint(c);
                ds->makeYield(c);
                ds->run(c);
            }
            else 
            {
                clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
            }

            delete ps;
            delete ds;
            delete c;
        }
    }
}

void thesis_ffcr_2tag(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");
    b->name_convention = Utils::NameConvention::CxAODReader;

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_SS",     "2-b-tagged SS ID Region",         eRegionType::SR);
    rs->add("2tag2pjet_0ptv_1A_SS",     "2-b-tagged SS Anti-ID Region",    eRegionType::SR);
    rs->add("2tag2pjet_0ptv_LL_OS",     "2-b-tagged OS ID Region",         eRegionType::SR);
    rs->add("2tag2pjet_0ptv_1A_OS",     "2-b-tagged OS Anti-ID Region",    eRegionType::SR);

    vector<double> binning_1p{30, 35, 40, 45, 50, 60, 75, 1000};
    vector<double> binning_3p{30, 40, 50, 65, 1000};

    Variables* vs = new Variables();
    vs->add("TauPt1P",              "1-prong #tau_{had} p_{T} [GeV]",                       1,          &binning_1p[0],         binning_1p.size()-1);
    vs->add("TauPt3P",              "3-prong #tau_{had} p_{T} [GeV]",                       1,          &binning_3p[0],         binning_3p.size()-1);

    DrawStackInfo* info = new DrawStackInfo();
    info->atlas = false;
    info->blind = false;
    info->show_scaling = true;
    info->legend_ncolumns = 1;
    info->legend_scaling_horizontal = 1.7;
    info->legend_scaling_vertical = 0.95;
    info->logx = true;
    info->logy = false;
    info->output_format = "pdf";
    info->draw_overflow = false;
    info->draw_ratio = false;
    info->xmax = 1000;
    info->ymax_ratio_nolog = 1.8;

    for (RegionInfo* r : *(rs->content()))
    {
        for (VariableInfo* v : *(vs->content()))
        {
            Processes* ps = new Processes();
            ps->add("data",       "data",               eProcessType::DATA,     eProcess::DATA,         "Data",                             kBlack);
            ps->add("ttbar",      "t#bar{t}",           eProcessType::BKG,      eProcess::P2,           "True-#tau_{had} t#bar{t}",         bbtt_kGOLDEN);
            ps->add("ttbarTF",    "t#bar{t} TF",        eProcessType::BKG,      eProcess::P1,           "Fake-#tau_{had} t#bar{t}",         bbtt_kORANGE);
            ps->add("ttbarFT",    "t#bar{t} FT",        eProcessType::BKG,      eProcess::P1,           "Fake-#tau_{had} t#bar{t}",         bbtt_kORANGE);
            ps->add("ttbarFF",    "t#bar{t} FF",        eProcessType::BKG,      eProcess::P1,           "Fake-#tau_{had} t#bar{t}",         bbtt_kORANGE);
            ps->add("Zttbb",      "Z#tau#tau + bb",     eProcessType::BKG,      eProcess::P4,           "Z #rightarrow #tau#tau + hf",      bbtt_kBLUE_L);
            ps->add("Zttbc",      "Z#tau#tau + bc",     eProcessType::BKG,      eProcess::P4,           "Z #rightarrow #tau#tau + hf",      bbtt_kBLUE_L);
            ps->add("Zttcc",      "Z#tau#tau + cc",     eProcessType::BKG,      eProcess::P4,           "Z #rightarrow #tau#tau + hf",      bbtt_kBLUE_L);
            ps->add("Zttbl",      "Z#tau#tau + bl",     eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zttcl",      "Z#tau#tau + cl",     eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zttl",       "Z#tau#tau + l",      eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("stops",      "single top",         eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("stopt",      "single top",         eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("stopWt",     "single top",         eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zbb",        "Zll + bb",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zbc",        "Zll + bc",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zbl",        "Zll + bl",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zcc",        "Zll + cc",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zcl",        "Zll + cl",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Zl",         "Zll + l",            eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("WW",         "WW",                 eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("WZ",         "WZ",                 eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ZZ",         "ZZ",                 eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("Wtt",        "W+jets",             eProcessType::BKG,      eProcess::P10,           "Others",                           bbtt_kGREEN);
            ps->add("W",          "W+jets",             eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ttW",        "ttV",                eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ttZ",        "ttV",                eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ttH",        "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("WHbb",       "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ggZHbb",     "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("qqZHbb",     "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("WHtautau",   "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ggZHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("qqZHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("ggFHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
            ps->add("VBFHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);

            Config* c = new Config(b, ps, rs, vs);
            c->load(filename, "FF");
            info->parameter = "FF";

            c->updateHistogramPtr(r, v);
            DrawStackTool* ds = new DrawStackTool(info);
            ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/FFCR/";
            if (ds->check(c))
            {
                ds->manipulate(c);
                ds->rebin(c, eRebinOption::Array);
                ds->paint(c);
                ds->makeYield(c);
                ds->run(c);
            }
            else 
            {
                clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
            }

            delete ps;
            delete ds;
            delete c;
        }
    }
}