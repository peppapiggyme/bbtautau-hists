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

static void thesis_sr_ws_input_helper(const std::string& filename,
    const std::string& varname, const std::string& varname_tex, const std::string& folder);

void thesis_sr_ws_input(const std::string& filename)
{
    thesis_sr_ws_input_helper(filename, "dRBB", "#Delta R(b, b)", "Region_BMin0_incJet1_distdRBB_J2_Y2015_DLLOS_T2_SpcTauHH_L0");
    thesis_sr_ws_input_helper(filename, "dRTauTau", "#Delta R(#tau, #tau)", "Region_BMin0_incJet1_distdRTauTau_J2_Y2015_DLLOS_T2_SpcTauHH_L0");
    thesis_sr_ws_input_helper(filename, "mBB", "m_{bb} [GeV]", "Region_BMin0_incJet1_distmBB_J2_Y2015_DLLOS_T2_SpcTauHH_L0");
    thesis_sr_ws_input_helper(filename, "mHH", "m_{HH} [GeV]", "Region_BMin0_incJet1_distmHH_J2_Y2015_DLLOS_T2_SpcTauHH_L0");
    thesis_sr_ws_input_helper(filename, "mMMC", "m_{#tau#tau}^{MMC} [GeV]", "Region_BMin0_incJet1_distmMMC_J2_Y2015_DLLOS_T2_SpcTauHH_L0");
}

void thesis_sr_ws_input_helper(const std::string& filename, 
    const std::string& varname, const std::string& varname_tex, const std::string& folder)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");
    b->name_convention = Utils::NameConvention::WSMaker;

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "SR",         eRegionType::SR);

    Variables* vs = new Variables();
    vs->add(varname,                    varname_tex,                       1);

    DrawStackInfo* info = new DrawStackInfo();
    info->atlas = false;
    info->blind = false;
    info->signal_scale = 100;
    info->signal_linestyle = 2;
    info->signal_linewidth = 3;
    info->show_scaling = true;
    info->legend_ncolumns = 1;
    info->legend_scaling_horizontal = 1.7;
    info->logy = false;
    info->output_format = "pdf";
    info->ymax_ratio_nolog = 2.5;
    info->draw_overflow = false;
    info->draw_ratio = false;
    info->ratio_high = 1.45;
    info->ratio_low = 0.55;
    info->prefit = false;

    for (VariableInfo* v : *(vs->content()))
    {
        // if (v->name == "mBB" || v->name == "mMMC")
        // {
        //     info->xmax = 600.;
        // }
        // else if (v->name == "mHH")
        // {
        //     info->xmax = 2000.;
        // }
        // else 
        // {
        //     info->xmax = 3.;
        // }

        Processes* ps = new Processes();
        ps->add("data",       "data",               eProcessType::DATA,     eProcess::DATA,         "Data",                             kBlack);
        ps->add("Ztthf",      "Z#tau#tau + hf",     eProcessType::BKG,      eProcess::P1,           "Z #rightarrow #tau#tau + hf",      bbtt_kBLUE_L);
        ps->add("ttbar",      "t#bar{t}",           eProcessType::BKG,      eProcess::P4,           "True-#tau_{had} t#bar{t}",         bbtt_kGOLDEN);
        ps->add("ttbarFake",  "t#bar{t} Fake",      eProcessType::BKG,      eProcess::P6,           "Fake-#tau_{had} t#bar{t}",         bbtt_kORANGE);
        ps->add("Fake",       "Multi-jet",          eProcessType::BKG,      eProcess::P8,           "Multi-jet",                        bbtt_kPINK);
        ps->add("Zttlf",      "Z#tau#tau + lf",     eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("stop",       "single top",         eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("Zhf",        "Zll + hf",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("Zlf",        "Zll + lf",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("diboson",     "Diboson",           eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("Wtt",        "W+jets",             eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("W",          "W+jets",             eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("ttW",        "ttV",                eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
        ps->add("ttZ",        "ttV",                eProcessType::BKG,      eProcess::P10,          "Others",                           bbtt_kGREEN);
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
        c->load(filename, folder);
        info->parameter = "Input_SMBDT";

        c->updateHistogramPtr(rs->content()->back(), v);
        DrawStackTool* ds = new DrawStackTool_WS(info);
        ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SR_MVAInput_WS_Post/";
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
