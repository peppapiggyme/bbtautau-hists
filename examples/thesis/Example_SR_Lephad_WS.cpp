#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "DrawStackTool_WS.h"
#include "Colors.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;
using BU = BinningUtils;

void thesis_sr_lephad_ws_bdt(const std::string& filename, bool prefit, bool SLT)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");
    b->name_convention = Utils::NameConvention::WSMaker;

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",    "#tau_{lep}#tau_{had} " + (SLT ? string("SLT") : string("LTT")) + " SR",         eRegionType::SR);

    Variables* vs = new Variables();
    vs->add("SMBDT",                   "BDT score",                       1);

    DrawStackInfo* info = new DrawStackInfo();
    info->atlas = false;
    info->blind = false;
    info->ratio_high = 1.35;
    info->ratio_low = 0.65;
    info->signal_scale = 100;
    info->signal_linestyle = 2;
    info->signal_linewidth = 3;
    info->show_scaling = true;
    info->legend_ncolumns = 1;
    info->legend_scaling_horizontal = 1.5;
    info->legend_scaling_vertical = 1.2;
    info->logy = false;
    info->output_format = "pdf";
    info->logy = true;
    info->draw_overflow = false;
    info->ratio_high = 1.45;
    info->ratio_low = 0.55;
    info->prefit = prefit;

    for (VariableInfo* v : *(vs->content()))
    {
        Processes* ps = new Processes();
        ps->add("data",       "data",               eProcessType::DATA,     eProcess::DATA,         "Data",                             kBlack);
        ps->add("ttbar",      "t#bar{t}",           eProcessType::BKG,      eProcess::P1,           "True-#tau_{had} t#bar{t}",         bbtt_kGOLDEN);
        ps->add("Fake",       "Fake #tau_{had}",    eProcessType::BKG,      eProcess::P4,           "Fake #tau_{had}",                        bbtt_kPINK);
        ps->add("Ztthf",      "Z#tau#tau + hf",     eProcessType::BKG,      eProcess::P8,           "Z #rightarrow #tau#tau + hf",      bbtt_kBLUE_L);
        ps->add("Zttlf",      "Z#tau#tau + lf",     eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("stop",       "single top",         eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("Zhf",        "Zll + hf",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("Zlf",        "Zll + lf",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("diboson",     "Diboson",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("Wtt",        "W+jets",             eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("W",          "W+jets",             eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("ttW",        "ttV",                eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("ttZ",        "ttV",                eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("DY",         "DY",                 eProcessType::BKG,      eProcess::P10,          "Others",                     bbtt_kCYAN);
        ps->add("DYtt",       "DYtt",               eProcessType::BKG,      eProcess::P10,          "Others",                     bbtt_kCYAN);
        ps->add("ttH",        "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("WHbb",       "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("ggZHbb",     "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("qqZHbb",     "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("WHtautau",   "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("ggZHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("qqZHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("ggFHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("VBFHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("hhttbb",      "SM HH",             eProcessType::SIG,      eProcess::P101,         "SM HH",                            kRed);
        ps->add("hhttbbVBFSM", "SM HH",             eProcessType::SIG,      eProcess::P101,         "SM HH",                            kRed);

        Config* c = new Config(b, ps, rs, vs);
        string folder = SLT ? "Region_BMin0_incJet1_distNN_J2_DSM_T2_SpcTauLH_Y2015_LTT0_L1" : "Region_BMin0_incJet1_distNN_J2_DSM_T2_SpcTauLH_Y2015_LTT1_L1";
        c->load(filename, folder);
        info->parameter = "BDT_" + (SLT ? string("SLT") : string("LTT"));

        c->updateHistogramPtr(rs->content()->back(), v);
        DrawStackTool* ds = new DrawStackTool_WS(info);
        if (info->prefit)
            ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SR_Comb_WS/";
        else
            ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SR_Comb_WS_Post/";
        if (ds->check(c))
        {
            ds->manipulate(c);
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

void thesis_sr_lephad_ws_pnn(const std::string& filename, const std::string& mass, bool prefit, bool SLT)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");
    b->name_convention = Utils::NameConvention::WSMaker;

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "#tau_{lep}#tau_{had} " + (SLT ? string("SLT") : string("LTT")) + " SR",         eRegionType::SR);

    Variables* vs = new Variables();
    vs->add("PNN" + mass,               "PNN" + mass + " score",                1);

    DrawStackInfo* info = new DrawStackInfo();
    info->atlas = false;
    info->blind = false;
    info->ratio_high = 1.35;
    info->ratio_low = 0.65;
    info->signal_scale = 0.1;
    info->signal_linestyle = 2;
    info->signal_linewidth = 3;
    info->show_scaling = true;
    info->legend_ncolumns = 1;
    info->legend_scaling_horizontal = 1.7;
    info->logy = false;
    info->output_format = "pdf";
    info->logy = true;
    info->draw_overflow = false;
    info->ratio_high = 1.55;
    info->ratio_low = 0.45;
    info->auto_ratio = true;
    info->prefit = prefit;
    
    for (VariableInfo* v : *(vs->content()))
    {
        Processes* ps = new Processes();
        ps->add("data",       "data",               eProcessType::DATA,     eProcess::DATA,         "Data",                             kBlack);
        ps->add("ttbar",      "t#bar{t}",           eProcessType::BKG,      eProcess::P1,           "True-#tau_{had} t#bar{t}",         bbtt_kGOLDEN);
        ps->add("Fake",       "Fake #tau_{had}",    eProcessType::BKG,      eProcess::P4,           "Fake #tau_{had}",                        bbtt_kPINK);
        ps->add("Ztthf",      "Z#tau#tau + hf",     eProcessType::BKG,      eProcess::P8,           "Z #rightarrow #tau#tau + hf",      bbtt_kBLUE_L);
        ps->add("Zttlf",      "Z#tau#tau + lf",     eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("stop",       "single top",         eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("Zhf",        "Zll + hf",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("Zlf",        "Zll + lf",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("diboson",     "Diboson",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("Wtt",        "W+jets",             eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("W",          "W+jets",             eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("ttW",        "ttV",                eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("ttZ",        "ttV",                eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("DY",         "DY",                 eProcessType::BKG,      eProcess::P10,          "Others",                     bbtt_kCYAN);
        ps->add("DYtt",       "DYtt",               eProcessType::BKG,      eProcess::P10,          "Others",                     bbtt_kCYAN);
        ps->add("ttH",        "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("WHbb",       "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("ggZHbb",     "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("qqZHbb",     "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("WHtautau",   "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("ggZHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("qqZHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("ggFHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("VBFHtautau", "SM Higgs",           eProcessType::BKG,      eProcess::P12,          "Single Higgs",                     bbtt_kCYAN);
        ps->add("Hhhbbtautau"+mass, "X("+mass+")",  eProcessType::SIG,      eProcess::P101,         "X("+mass+")",                      kRed);

        Config* c = new Config(b, ps, rs, vs);
        string folder = SLT ? "Region_BMin0_incJet1_dist" + mass + "_J2_D2HDMPNN_T2_SpcTauLH_Y2015_LTT0_L1" : "Region_BMin0_incJet1_dist" + mass + "_J2_D2HDMPNN_T2_SpcTauLH_Y2015_LTT1_L1";
        c->load(filename, folder);
        info->parameter = "PNN_" + (SLT ? string("SLT") : string("LTT"));

        c->updateHistogramPtr(rs->content()->back(), v);
        DrawStackTool* ds = new DrawStackTool_WS(info);
        if (info->prefit)
            ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SR_Comb_WS/";
        else
            ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SR_Comb_WS_Post/";
        if (ds->check(c))
        {
            ds->manipulate(c);
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