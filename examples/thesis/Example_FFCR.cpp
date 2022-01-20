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

#define EQUALWIDTH_BINNING(B, START, END, STEP) \
    for (int i = (START); i <= (END); i += (STEP)) (B).push_back(i);

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

void thesis_ffcr_closure(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");
    b->name_convention = Utils::NameConvention::CxAODReader;

    Regions* rs = new Regions();
    rs->add("1tag2pjet_0ptv_LL_SS",     "1-b-tagged SS ID Region",         eRegionType::SR);

    vector<double> bTau0Pt;     EQUALWIDTH_BINNING(bTau0Pt  , 40 , 240  , 5);
    vector<double> bTau1Pt;     EQUALWIDTH_BINNING(bTau1Pt  , 25 , 125  , 5);

    Variables* vs = new Variables();
    vs->add("Tau0Pt",              "Leading #tau_{had} p_{T} [GeV]",                    1,      &bTau0Pt[0],     bTau0Pt.size()-1);
    vs->add("Tau1Pt",              "Sub-leading #tau_{had} p_{T} [GeV]",                1,      &bTau1Pt[0],     bTau1Pt.size()-1);
    vs->add("Tau0Eta",             "Leading #tau_{had} #eta",                           2);
    vs->add("Tau1Eta",             "Sub-leading #tau_{had} #eta",                       2);

    Systematics* ss = new Systematics();

    ss->add("EG_RESOLUTION_ALL", "EG_RESOLUTION_ALL", eSystematicType::TwoSide, 0);
    ss->add("EG_SCALE_ALL", "EG_SCALE_ALL", eSystematicType::TwoSide, 0);
    ss->add("EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR", eSystematicType::TwoSide, 0);
    ss->add("EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR", eSystematicType::TwoSide, 0);
    ss->add("EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_Eigen_B_0", "FT_EFF_Eigen_B_0", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_Eigen_B_1", "FT_EFF_Eigen_B_1", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_Eigen_B_2", "FT_EFF_Eigen_B_2", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_Eigen_C_0", "FT_EFF_Eigen_C_0", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_Eigen_C_1", "FT_EFF_Eigen_C_1", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_Eigen_C_2", "FT_EFF_Eigen_C_2", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_Eigen_Light_0", "FT_EFF_Eigen_Light_0", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_Eigen_Light_1", "FT_EFF_Eigen_Light_1", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_Eigen_Light_2", "FT_EFF_Eigen_Light_2", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_Eigen_Light_3", "FT_EFF_Eigen_Light_3", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_extrapolation", "FT_EFF_extrapolation", eSystematicType::TwoSide, 0);
    ss->add("FT_EFF_extrapolation_from_charm", "FT_EFF_extrapolation_from_charm", eSystematicType::TwoSide, 0);
    ss->add("JET_BJES_Response", "JET_BJES_Response", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Detector1", "JET_EffectiveNP_Detector1", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Detector2", "JET_EffectiveNP_Detector2", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Mixed1", "JET_EffectiveNP_Mixed1", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Mixed2", "JET_EffectiveNP_Mixed2", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Mixed3", "JET_EffectiveNP_Mixed3", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Modelling1", "JET_EffectiveNP_Modelling1", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Modelling2", "JET_EffectiveNP_Modelling2", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Modelling3", "JET_EffectiveNP_Modelling3", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Modelling4", "JET_EffectiveNP_Modelling4", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Statistical1", "JET_EffectiveNP_Statistical1", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Statistical2", "JET_EffectiveNP_Statistical2", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Statistical3", "JET_EffectiveNP_Statistical3", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Statistical4", "JET_EffectiveNP_Statistical4", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Statistical5", "JET_EffectiveNP_Statistical5", eSystematicType::TwoSide, 0);
    ss->add("JET_EffectiveNP_Statistical6", "JET_EffectiveNP_Statistical6", eSystematicType::TwoSide, 0);
    ss->add("JET_EtaIntercalibration_Modelling", "JET_EtaIntercalibration_Modelling", eSystematicType::TwoSide, 0);
    ss->add("JET_EtaIntercalibration_NonClosure_2018data", "JET_EtaIntercalibration_NonClosure_2018data", eSystematicType::TwoSide, 0);
    ss->add("JET_EtaIntercalibration_NonClosure_highE", "JET_EtaIntercalibration_NonClosure_highE", eSystematicType::TwoSide, 0);
    ss->add("JET_EtaIntercalibration_NonClosure_posEta", "JET_EtaIntercalibration_NonClosure_posEta", eSystematicType::TwoSide, 0);
    ss->add("JET_EtaIntercalibration_TotalStat", "JET_EtaIntercalibration_TotalStat", eSystematicType::TwoSide, 0);
    ss->add("JET_Flavor_Composition", "JET_Flavor_Composition", eSystematicType::TwoSide, 0);
    ss->add("JET_Flavor_Response", "JET_Flavor_Response", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_DataVsMC_MC16", "JET_JER_DataVsMC_MC16", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_1", "JET_JER_EffectiveNP_1", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_10", "JET_JER_EffectiveNP_10", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_11", "JET_JER_EffectiveNP_11", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_12restTerm", "JET_JER_EffectiveNP_12restTerm", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_2", "JET_JER_EffectiveNP_2", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_3", "JET_JER_EffectiveNP_3", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_4", "JET_JER_EffectiveNP_4", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_5", "JET_JER_EffectiveNP_5", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_6", "JET_JER_EffectiveNP_6", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_7", "JET_JER_EffectiveNP_7", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_8", "JET_JER_EffectiveNP_8", eSystematicType::TwoSide, 0);
    ss->add("JET_JER_EffectiveNP_9", "JET_JER_EffectiveNP_9", eSystematicType::TwoSide, 0);
    ss->add("JET_JvtEfficiency", "JET_JvtEfficiency", eSystematicType::TwoSide, 0);
    ss->add("JET_Pileup_OffsetMu", "JET_Pileup_OffsetMu", eSystematicType::TwoSide, 0);
    ss->add("JET_Pileup_OffsetNPV", "JET_Pileup_OffsetNPV", eSystematicType::TwoSide, 0);
    ss->add("JET_Pileup_PtTerm", "JET_Pileup_PtTerm", eSystematicType::TwoSide, 0);
    ss->add("JET_Pileup_RhoTopology", "JET_Pileup_RhoTopology", eSystematicType::TwoSide, 0);
    ss->add("JET_PunchThrough_MC16", "JET_PunchThrough_MC16", eSystematicType::TwoSide, 0);
    ss->add("MET_SoftTrk_ResoPara", "MET_SoftTrk_ResoPara", eSystematicType::OneSide, 0);
    ss->add("MET_SoftTrk_ResoPerp", "MET_SoftTrk_ResoPerp", eSystematicType::OneSide, 0);
    ss->add("MET_SoftTrk_Scale", "MET_SoftTrk_Scale", eSystematicType::TwoSide, 0);
    ss->add("MUON_EFF_ISO_STAT", "MUON_EFF_ISO_STAT", eSystematicType::TwoSide, 0);
    ss->add("MUON_EFF_ISO_SYS", "MUON_EFF_ISO_SYS", eSystematicType::TwoSide, 0);
    ss->add("MUON_EFF_RECO_SYS_LOWPT", "MUON_EFF_RECO_SYS_LOWPT", eSystematicType::TwoSide, 0);
    ss->add("MUON_EFF_TrigStatUncertainty", "MUON_EFF_TrigStatUncertainty", eSystematicType::TwoSide, 0);
    ss->add("MUON_EFF_TrigSystUncertainty", "MUON_EFF_TrigSystUncertainty", eSystematicType::TwoSide, 0);
    ss->add("MUON_ID", "MUON_ID", eSystematicType::TwoSide, 0);
    ss->add("MUON_MS", "MUON_MS", eSystematicType::TwoSide, 0);
    ss->add("MUON_SCALE", "MUON_SCALE", eSystematicType::TwoSide, 0);
    ss->add("PRW_DATASF", "PRW_DATASF", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEELECTRON_EFF_ELEBDT_STAT", "TAUS_TRUEELECTRON_EFF_ELEBDT_STAT", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEELECTRON_EFF_ELEBDT_SYST", "TAUS_TRUEELECTRON_EFF_ELEBDT_SYST", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL", "TAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_RECO_TOTAL", "TAUS_TRUEHADTAU_EFF_RECO_TOTAL", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPTGE40", "TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPTGE40", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_RNNID_HIGHPT", "TAUS_TRUEHADTAU_EFF_RNNID_HIGHPT", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_RNNID_SYST", "TAUS_TRUEHADTAU_EFF_RNNID_SYST", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA161718", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA161718", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA1718", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA1718", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2016", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2016", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018AFTTS1", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018AFTTS1", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2018", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2018", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU161718", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU161718", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_SME_TES_INSITUEXP", "TAUS_TRUEHADTAU_SME_TES_INSITUEXP", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_SME_TES_INSITUFIT", "TAUS_TRUEHADTAU_SME_TES_INSITUFIT", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_SME_TES_MODEL_CLOSURE", "TAUS_TRUEHADTAU_SME_TES_MODEL_CLOSURE", eSystematicType::TwoSide, 0);
    ss->add("TAUS_TRUEHADTAU_SME_TES_PHYSICSLIST", "TAUS_TRUEHADTAU_SME_TES_PHYSICSLIST", eSystematicType::TwoSide, 0);

    DrawStackInfo* info = new DrawStackInfo();
    info->atlas = false;
    info->blind = false;
    info->show_scaling = true;
    info->legend_ncolumns = 1;
    info->legend_scaling_horizontal = 1.7;
    info->legend_scaling_vertical = 0.85;
    info->ratio_high = 1.35;
    info->ratio_low = 0.65;
    info->logx = false;
    info->logy = false;
    info->output_format = "pdf";
    info->draw_overflow = true;
    info->draw_ratio = true;

    for (RegionInfo* r : *(rs->content()))
    {
        for (VariableInfo* v : *(vs->content()))
        {
            info->draw_overflow = v->name.find("Eta") == std::string::npos;
            info->ymax_ratio_nolog = v->name.find("Eta") == std::string::npos ? 1.4 : 1.7;

            Processes* ps = new Processes();
            ps->add("data",       "data",               eProcessType::DATA,     eProcess::DATA,         "Data",                             kBlack);
            ps->add("Fake",       "Fake",               eProcessType::BKG,      eProcess::P1,           "Mutli-jet",                        bbtt_kPINK);
            ps->add("ttbar",      "t#bar{t}",           eProcessType::BKG,      eProcess::P2,           "True-#tau_{had} t#bar{t}",         bbtt_kGOLDEN);
            ps->add("ttbarSFTF",  "t#bar{t} TF",        eProcessType::BKG,      eProcess::P3,           "Fake-#tau_{had} t#bar{t}",         bbtt_kORANGE);
            ps->add("ttbarSFFT",  "t#bar{t} FT",        eProcessType::BKG,      eProcess::P3,           "Fake-#tau_{had} t#bar{t}",         bbtt_kORANGE);
            ps->add("ttbarSFFF",  "t#bar{t} FF",        eProcessType::BKG,      eProcess::P3,           "Fake-#tau_{had} t#bar{t}",         bbtt_kORANGE);
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

            Config* c = new Config(b, ps, rs, vs, ss);
            c->load(filename, "Preselection");
            info->parameter = "Closure";

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