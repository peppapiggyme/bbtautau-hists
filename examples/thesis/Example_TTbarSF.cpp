#include "ExamplesInclude.h"

#include "Utils.h"
#include "Colors.h"
#include "Config.h"
#include "DrawStackTool.h"
#include "AutoBinningTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;
using BU = BinningUtils;

void thesis_ttbarsf(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");
    b->name_convention = Utils::NameConvention::CxAODReader;

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_OS",     "t#bar{t}-ID-CR",         eRegionType::SR);

    Variables* vs = new Variables();
    vs->add("mTW",       "m_{T}^{W} [GeV]",           2);
    vs->add("TauPt",     "#tau_{had} p_{T} [GeV]",    4);

    Systematics* ss = new Systematics();

    ss->add("EG_RESOLUTION_ALL", "EG_RESOLUTION_ALL", eSystematicType::TwoSide, 0, false, true);
    ss->add("EG_SCALE_ALL", "EG_SCALE_ALL", eSystematicType::TwoSide, 0, false, true);
    ss->add("EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR", eSystematicType::TwoSide, 0, false, true);
    ss->add("EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR", eSystematicType::TwoSide, 0, false, true);
    ss->add("EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR", eSystematicType::TwoSide, 0, false, true);
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
    ss->add("JET_BJES_Response", "JET_BJES_Response", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Detector1", "JET_EffectiveNP_Detector1", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Detector2", "JET_EffectiveNP_Detector2", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Mixed1", "JET_EffectiveNP_Mixed1", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Mixed2", "JET_EffectiveNP_Mixed2", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Mixed3", "JET_EffectiveNP_Mixed3", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Modelling1", "JET_EffectiveNP_Modelling1", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Modelling2", "JET_EffectiveNP_Modelling2", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Modelling3", "JET_EffectiveNP_Modelling3", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Modelling4", "JET_EffectiveNP_Modelling4", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Statistical1", "JET_EffectiveNP_Statistical1", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Statistical2", "JET_EffectiveNP_Statistical2", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Statistical3", "JET_EffectiveNP_Statistical3", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Statistical4", "JET_EffectiveNP_Statistical4", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Statistical5", "JET_EffectiveNP_Statistical5", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EffectiveNP_Statistical6", "JET_EffectiveNP_Statistical6", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EtaIntercalibration_Modelling", "JET_EtaIntercalibration_Modelling", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EtaIntercalibration_NonClosure_2018data", "JET_EtaIntercalibration_NonClosure_2018data", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EtaIntercalibration_NonClosure_highE", "JET_EtaIntercalibration_NonClosure_highE", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EtaIntercalibration_NonClosure_posEta", "JET_EtaIntercalibration_NonClosure_posEta", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_EtaIntercalibration_TotalStat", "JET_EtaIntercalibration_TotalStat", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_Flavor_Composition", "JET_Flavor_Composition", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_Flavor_Response", "JET_Flavor_Response", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_DataVsMC_MC16", "JET_JER_DataVsMC_MC16", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_1", "JET_JER_EffectiveNP_1", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_10", "JET_JER_EffectiveNP_10", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_11", "JET_JER_EffectiveNP_11", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_12restTerm", "JET_JER_EffectiveNP_12restTerm", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_2", "JET_JER_EffectiveNP_2", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_3", "JET_JER_EffectiveNP_3", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_4", "JET_JER_EffectiveNP_4", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_5", "JET_JER_EffectiveNP_5", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_6", "JET_JER_EffectiveNP_6", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_7", "JET_JER_EffectiveNP_7", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_8", "JET_JER_EffectiveNP_8", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JER_EffectiveNP_9", "JET_JER_EffectiveNP_9", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_JvtEfficiency", "JET_JvtEfficiency", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_Pileup_OffsetMu", "JET_Pileup_OffsetMu", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_Pileup_OffsetNPV", "JET_Pileup_OffsetNPV", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_Pileup_PtTerm", "JET_Pileup_PtTerm", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_Pileup_RhoTopology", "JET_Pileup_RhoTopology", eSystematicType::TwoSide, 0, false, true);
    ss->add("JET_PunchThrough_MC16", "JET_PunchThrough_MC16", eSystematicType::TwoSide, 0, false, true);
    ss->add("MET_SoftTrk_ResoPara", "MET_SoftTrk_ResoPara", eSystematicType::OneSide, 0, false, true);
    ss->add("MET_SoftTrk_ResoPerp", "MET_SoftTrk_ResoPerp", eSystematicType::OneSide, 0, false, true);
    ss->add("MET_SoftTrk_Scale", "MET_SoftTrk_Scale", eSystematicType::TwoSide, 0, false, true);
    ss->add("MUON_EFF_ISO_STAT", "MUON_EFF_ISO_STAT", eSystematicType::TwoSide, 0, false, true);
    ss->add("MUON_EFF_ISO_SYS", "MUON_EFF_ISO_SYS", eSystematicType::TwoSide, 0, false, true);
    ss->add("MUON_EFF_RECO_SYS_LOWPT", "MUON_EFF_RECO_SYS_LOWPT", eSystematicType::TwoSide, 0, false, true);
    ss->add("MUON_EFF_TrigStatUncertainty", "MUON_EFF_TrigStatUncertainty", eSystematicType::TwoSide, 0, false, true);
    ss->add("MUON_EFF_TrigSystUncertainty", "MUON_EFF_TrigSystUncertainty", eSystematicType::TwoSide, 0, false, true);
    ss->add("MUON_ID", "MUON_ID", eSystematicType::TwoSide, 0, false, true);
    ss->add("MUON_MS", "MUON_MS", eSystematicType::TwoSide, 0, false, true);
    ss->add("MUON_SCALE", "MUON_SCALE", eSystematicType::TwoSide, 0, false, true);
    ss->add("PRW_DATASF", "PRW_DATASF", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEELECTRON_EFF_ELEBDT_STAT", "TAUS_TRUEELECTRON_EFF_ELEBDT_STAT", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEELECTRON_EFF_ELEBDT_SYST", "TAUS_TRUEELECTRON_EFF_ELEBDT_SYST", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL", "TAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_RECO_TOTAL", "TAUS_TRUEHADTAU_EFF_RECO_TOTAL", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPTGE40", "TAUS_TRUEHADTAU_EFF_RNNID_1PRONGSTATSYSTPTGE40", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_RNNID_HIGHPT", "TAUS_TRUEHADTAU_EFF_RNNID_HIGHPT", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_RNNID_SYST", "TAUS_TRUEHADTAU_EFF_RNNID_SYST", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA161718", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA161718", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA1718", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA1718", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2016", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2016", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018AFTTS1", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2018AFTTS1", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2018", "TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2018", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU161718", "TAUS_TRUEHADTAU_EFF_TRIGGER_SYSTMU161718", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_SME_TES_INSITUEXP", "TAUS_TRUEHADTAU_SME_TES_INSITUEXP", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_SME_TES_INSITUFIT", "TAUS_TRUEHADTAU_SME_TES_INSITUFIT", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_SME_TES_MODEL_CLOSURE", "TAUS_TRUEHADTAU_SME_TES_MODEL_CLOSURE", eSystematicType::TwoSide, 0, false, true);
    ss->add("TAUS_TRUEHADTAU_SME_TES_PHYSICSLIST", "TAUS_TRUEHADTAU_SME_TES_PHYSICSLIST", eSystematicType::TwoSide, 0, false, true);
    ss->add("TTBarReweight", "TTBarReweight", eSystematicType::OneSide, 0);
    ss->add("TTBarReweight_Stat", "TTBarReweight_Stat", eSystematicType::OneSide, 0);
    ss->add("TTBarReweight_Closure_TauPt", "TTBarReweight_Closure_TauPt", eSystematicType::OneSide, 0);

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
    info->legend_scaling_horizontal = 1.9;
    info->logy = false;
    info->output_format = "pdf";
    info->draw_overflow = true;

    for (VariableInfo* v : *(vs->content()))
    {
        Processes* ps = new Processes();
        map<eProcess, double> scaleFactors = {
            {eProcess::P1, 0.93},
            {eProcess::P3, 0.93},
        };
        ps->setScaleFactorsMap(scaleFactors);

        ps->add("data",       "data",                       eProcessType::DATA,  eProcess::DATA,    "Data",                         kBlack);
        ps->add("ttbar",      "true-#tau_{had} t#bar{t}",   eProcessType::BKG,   eProcess::P1,      "True-#tau_{had} t#bar{t}",     bbtt_kGOLDEN);
        ps->add("ttbarFake",  "fake-#tau_{had} t#bar{t}",   eProcessType::BKG,   eProcess::P3,      "Fake-#tau_{had} t#bar{t}",     bbtt_kORANGE);
        ps->add("stopWt",     "single top",                 eProcessType::BKG,   eProcess::P5,      "Single top",             bbtt_kBROWN);
        ps->add("stops",      "single top",                 eProcessType::BKG,   eProcess::P5,      "Single top",             bbtt_kBROWN);
        ps->add("stopt",      "single top",                 eProcessType::BKG,   eProcess::P5,      "Single top",             bbtt_kBROWN);
        ps->add("Z",          "Z+jets",                     eProcessType::BKG,   eProcess::P6,      "Z+jets",                       bbtt_kBLUE_L);
        ps->add("W",          "W+jets",                     eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("WW",         "WW",                         eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("WZ",         "WZ",                         eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("ZZ",         "ZZ",                         eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("DY",         "DY",                         eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("DYtt",       "DYtt",                       eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("ttW",        "ttW",                        eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("ttZ",        "ttZ",                        eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("ttH",        "SM Higgs",                   eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("WHbb",       "SM Higgs",                   eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("WHtautau",   "SM Higgs",                   eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("ggZHtautau", "SM Higgs",                   eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("qqZHtautau", "SM Higgs",                   eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("ggZHbb",     "SM Higgs",                   eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("qqZHbb",     "SM Higgs",                   eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("ggFHtautau", "SM Higgs",                   eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);
        ps->add("VBFHtautau", "SM Higgs",                   eProcessType::BKG,   eProcess::P8,      "Others",                       bbtt_kGREEN);

        Config* c = new Config(b, ps, rs, vs, ss);
        c->load(filename, "Preselection_Mbb150");
        info->parameter = "TTbarSF";

        c->updateHistogramPtr(rs->content()->back(), v);
        DrawStackTool* ds = new DrawStackTool(info);
        ds->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/TtbarSF/";
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