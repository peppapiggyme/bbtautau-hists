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

void thesis_zcr(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");
    b->name_convention = Utils::NameConvention::WSMaker;

    Regions* rs = new Regions();
    rs->add("ZhfCR",     "Z+hf CR",         eRegionType::SR);

    Variables* vs = new Variables();
    vs->add("mll",       "m_{ll} [GeV]",    1);

    Systematics* ss = new Systematics();

    // WARN: you must pass OneSide for DrawStackTool, otherwise it will double count!!!
    ss->add("EG_RESOLUTION_ALL", "EG_RESOLUTION_ALL", eSystematicType::OneSide, 0);
    ss->add("EG_SCALE_ALL", "EG_SCALE_ALL", eSystematicType::OneSide, 0);
    ss->add("EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR", eSystematicType::OneSide, 0);
    ss->add("EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR", eSystematicType::OneSide, 0);
    ss->add("EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_Eigen_B_0", "FT_EFF_Eigen_B_0", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_Eigen_B_1", "FT_EFF_Eigen_B_1", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_Eigen_B_2", "FT_EFF_Eigen_B_2", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_Eigen_C_0", "FT_EFF_Eigen_C_0", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_Eigen_C_1", "FT_EFF_Eigen_C_1", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_Eigen_C_2", "FT_EFF_Eigen_C_2", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_Eigen_Light_0", "FT_EFF_Eigen_Light_0", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_Eigen_Light_1", "FT_EFF_Eigen_Light_1", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_Eigen_Light_2", "FT_EFF_Eigen_Light_2", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_Eigen_Light_3", "FT_EFF_Eigen_Light_3", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_extrapolation", "FT_EFF_extrapolation", eSystematicType::OneSide, 0);
    ss->add("FT_EFF_extrapolation_from_charm", "FT_EFF_extrapolation_from_charm", eSystematicType::OneSide, 0);
    ss->add("JET_BJES_Response", "JET_BJES_Response", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Detector1", "JET_EffectiveNP_Detector1", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Detector2", "JET_EffectiveNP_Detector2", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Mixed1", "JET_EffectiveNP_Mixed1", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Mixed2", "JET_EffectiveNP_Mixed2", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Mixed3", "JET_EffectiveNP_Mixed3", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Modelling1", "JET_EffectiveNP_Modelling1", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Modelling2", "JET_EffectiveNP_Modelling2", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Modelling3", "JET_EffectiveNP_Modelling3", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Modelling4", "JET_EffectiveNP_Modelling4", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Statistical1", "JET_EffectiveNP_Statistical1", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Statistical2", "JET_EffectiveNP_Statistical2", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Statistical3", "JET_EffectiveNP_Statistical3", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Statistical4", "JET_EffectiveNP_Statistical4", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Statistical5", "JET_EffectiveNP_Statistical5", eSystematicType::OneSide, 0);
    ss->add("JET_EffectiveNP_Statistical6", "JET_EffectiveNP_Statistical6", eSystematicType::OneSide, 0);
    ss->add("JET_EtaIntercalibration_Modelling", "JET_EtaIntercalibration_Modelling", eSystematicType::OneSide, 0);
    ss->add("JET_EtaIntercalibration_NonClosure_2018data", "JET_EtaIntercalibration_NonClosure_2018data", eSystematicType::OneSide, 0);
    ss->add("JET_EtaIntercalibration_NonClosure_highE", "JET_EtaIntercalibration_NonClosure_highE", eSystematicType::OneSide, 0);
    ss->add("JET_EtaIntercalibration_NonClosure_posEta", "JET_EtaIntercalibration_NonClosure_posEta", eSystematicType::OneSide, 0);
    ss->add("JET_EtaIntercalibration_TotalStat", "JET_EtaIntercalibration_TotalStat", eSystematicType::OneSide, 0);
    ss->add("JET_Flavor_Composition", "JET_Flavor_Composition", eSystematicType::OneSide, 0);
    ss->add("JET_Flavor_Response", "JET_Flavor_Response", eSystematicType::OneSide, 0);
    ss->add("JET_JER_DataVsMC_MC16", "JET_JER_DataVsMC_MC16", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_1", "JET_JER_EffectiveNP_1", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_10", "JET_JER_EffectiveNP_10", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_11", "JET_JER_EffectiveNP_11", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_12restTerm", "JET_JER_EffectiveNP_12restTerm", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_2", "JET_JER_EffectiveNP_2", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_3", "JET_JER_EffectiveNP_3", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_4", "JET_JER_EffectiveNP_4", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_5", "JET_JER_EffectiveNP_5", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_6", "JET_JER_EffectiveNP_6", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_7", "JET_JER_EffectiveNP_7", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_8", "JET_JER_EffectiveNP_8", eSystematicType::OneSide, 0);
    ss->add("JET_JER_EffectiveNP_9", "JET_JER_EffectiveNP_9", eSystematicType::OneSide, 0);
    ss->add("JET_JvtEfficiency", "JET_JvtEfficiency", eSystematicType::OneSide, 0);
    ss->add("JET_Pileup_OffsetMu", "JET_Pileup_OffsetMu", eSystematicType::OneSide, 0);
    ss->add("JET_Pileup_OffsetNPV", "JET_Pileup_OffsetNPV", eSystematicType::OneSide, 0);
    ss->add("JET_Pileup_PtTerm", "JET_Pileup_PtTerm", eSystematicType::OneSide, 0);
    ss->add("JET_Pileup_RhoTopology", "JET_Pileup_RhoTopology", eSystematicType::OneSide, 0);
    ss->add("JET_PunchThrough_MC16", "JET_PunchThrough_MC16", eSystematicType::OneSide, 0);
    ss->add("MET_SoftTrk_ResoPara", "MET_SoftTrk_ResoPara", eSystematicType::OneSide, 0);
    ss->add("MET_SoftTrk_ResoPerp", "MET_SoftTrk_ResoPerp", eSystematicType::OneSide, 0);
    ss->add("MET_SoftTrk_Scale", "MET_SoftTrk_Scale", eSystematicType::OneSide, 0);
    ss->add("MUON_EFF_ISO_STAT", "MUON_EFF_ISO_STAT", eSystematicType::OneSide, 0);
    ss->add("MUON_EFF_ISO_SYS", "MUON_EFF_ISO_SYS", eSystematicType::OneSide, 0);
    ss->add("MUON_EFF_RECO_SYS_LOWPT", "MUON_EFF_RECO_SYS_LOWPT", eSystematicType::OneSide, 0);
    ss->add("MUON_EFF_TrigStatUncertainty", "MUON_EFF_TrigStatUncertainty", eSystematicType::OneSide, 0);
    ss->add("MUON_EFF_TrigSystUncertainty", "MUON_EFF_TrigSystUncertainty", eSystematicType::OneSide, 0);
    ss->add("MUON_ID", "MUON_ID", eSystematicType::OneSide, 0);
    ss->add("MUON_MS", "MUON_MS", eSystematicType::OneSide, 0);
    ss->add("MUON_SCALE", "MUON_SCALE", eSystematicType::OneSide, 0);
    ss->add("PRW_DATASF", "PRW_DATASF", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEELECTRON_EFF_ELEBDT_STAT", "TAUS_TRUEELECTRON_EFF_ELEBDT_STAT", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEELECTRON_EFF_ELEBDT_SYST", "TAUS_TRUEELECTRON_EFF_ELEBDT_SYST", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL", "TAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_RECO_TOTAL", "TAUS_TRUEHADTAU_EFF_RECO_TOTAL", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPTGE40", "TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPTGE40", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_RNNID_HIGHPT", "TAUS_TRUEHADTAU_EFF_RNNID_HIGHPT", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_RNNID_SYST", "TAUS_TRUEHADTAU_EFF_RNNID_SYST", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA161718", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA161718", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA1718", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA1718", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2016", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2016", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018AFTTS1", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018AFTTS1", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2018", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2018", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU161718", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU161718", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_SME_TES_INSITUEXP", "TAUS_TRUEHADTAU_SME_TES_INSITUEXP", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_SME_TES_INSITUFIT", "TAUS_TRUEHADTAU_SME_TES_INSITUFIT", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_SME_TES_MODEL_CLOSURE", "TAUS_TRUEHADTAU_SME_TES_MODEL_CLOSURE", eSystematicType::OneSide, 0);
    ss->add("TAUS_TRUEHADTAU_SME_TES_PHYSICSLIST", "TAUS_TRUEHADTAU_SME_TES_PHYSICSLIST", eSystematicType::OneSide, 0);

    DrawStackInfo* info = new DrawStackInfo();
    info->atlas = false;
    info->blind = false;
    info->ratio_high = 1.22;
    info->ratio_low = 0.78;
    info->signal_scale = 1;
    info->signal_linestyle = 2;
    info->signal_linewidth = 3;
    info->show_scaling = false;
    info->legend_ncolumns = 1;
    info->legend_scaling_horizontal = 1.7;
    info->logy = false;
    info->output_format = "pdf";
    info->draw_overflow = false;

    for (VariableInfo* v : *(vs->content()))
    {
        Processes* ps = new Processes();
        map<eProcess, double> scaleFactors = {
            {eProcess::P1, 1.38},
            {eProcess::P4, 0.96},
        };
        ps->setScaleFactorsMap(scaleFactors);

        ps->add("data",       "data",       eProcessType::DATA,  eProcess::DATA,    "Data",                    kBlack);
        ps->add("Zbb",        "Zll + bb",   eProcessType::BKG,   eProcess::P1,      "Z #rightarrow ll + hf",       TColor::GetColor(55, 169, 222));
        ps->add("Zbc",        "Zll + bc",   eProcessType::BKG,   eProcess::P1,      "Z #rightarrow ll + hf",       TColor::GetColor(55, 169, 222));
        ps->add("Zcc",        "Zll + cc",   eProcessType::BKG,   eProcess::P1,      "Z #rightarrow ll + hf",       TColor::GetColor(55, 169, 222));
        ps->add("Zbl",        "Zll + bl",   eProcessType::BKG,   eProcess::P2,      "Z #rightarrow ll + lf",       TColor::GetColor(62, 88, 189));
        ps->add("Zcl",        "Zll + cl",   eProcessType::BKG,   eProcess::P2,      "Z #rightarrow ll + lf",       TColor::GetColor(62, 88, 189));
        ps->add("Zl",         "Zll + l",    eProcessType::BKG,   eProcess::P2,      "Z #rightarrow ll + lf",       TColor::GetColor(62, 88, 189));
        ps->add("ttbar",      "t#bar{t}",   eProcessType::BKG,   eProcess::P4,      "t#bar{t}",                TColor::GetColor(222, 173, 33));
        // ps->add("WW",         "WW",         eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        ps->add("stopWt",     "single top", eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        // ps->add("WZ",         "WZ",         eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        ps->add("ZZ",         "ZZ",         eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        ps->add("DY",         "DY",         eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        ps->add("DYtt",       "DYtt",       eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        ps->add("W",          "W+jets",     eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        ps->add("ttH",        "SM Higgs",   eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        ps->add("HtoVV",      "SM Higgs",   eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        ps->add("VHtautau",   "SM Higgs",   eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        ps->add("ggHtautau",  "SM Higgs",   eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));
        ps->add("VBFHtautau", "SM Higgs",   eProcessType::BKG,   eProcess::P8,      "Others",                  TColor::GetColor(35, 140, 121));

        ps->add("Zbb",        "Zll + bb",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("Zbc",        "Zll + bc",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("Zcc",        "Zll + cc",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("Zbl",        "Zll + bl",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("Zcl",        "Zll + cl",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("Zl",         "Zll + l",    eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("ttbar",      "t#bar{t}",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        // ps->add("WW",         "WW",         eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("stopWt",     "single top", eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        // ps->add("WZ",         "WZ",         eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("ZZ",         "ZZ",         eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("DY",         "DY",         eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("DYtt",       "DYtt",       eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("W",          "W+jets",     eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("ttH",        "SM Higgs",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("HtoVV",      "SM Higgs",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("VHtautau",   "SM Higgs",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("ggHtautau",  "SM Higgs",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        ps->add("VBFHtautau", "SM Higgs",   eProcessType::SIG,   eProcess::P16,     "Pre-norm. Bkg.",               TColor::GetColor(255, 0, 0));
        
        Config* c = new Config(b, ps, rs, vs, ss);
        c->load(filename, "");
        info->parameter = "ZVR";

        c->updateHistogramPtr(rs->content()->back(), v);
        DrawStackTool* ds = new DrawStackTool(info);
        ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/ZhfCR/";
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