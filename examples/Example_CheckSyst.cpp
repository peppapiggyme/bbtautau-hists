#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "CompTool.h"
#include "SystCompTool.h"
#include "AutoBinningTool.h"
#include "Colors.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;
using BU = BinningUtils;

void test_check_syst(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("SR", "#sqrt{s} = 13 TeV");
    b->name_convention = Utils::NameConvention::WSMaker;

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "LepHad SLT",        eRegionType::SR);

    auto binningFromFile = [](const std::string& fn)
    { 
        return BU::intToDoubleBinEdgesForMVAInverse_LepHad1090Version(BU::readBinningFromFile<int>(fn));
    };

    // vector<double> binningSMNN      =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_LepHad_SLT_SM_NN.txt");
    // vector<double> binning2HDM500   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_LepHad_SLT_2HDM500.txt");
    // vector<double> binning2HDM1000  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_LepHad_SLT_2HDM1000.txt");
    vector<double> binning2HDM1100  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_LepHad_SLT_2HDM1100.txt");

    // vector<double> binning2HDM1000  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_LepHad_LTT_2HDM1000.txt");
    // vector<double> binning2HDM1100  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_LepHad_LTT_2HDM1100.txt");

    Variables* vs_pnn = new Variables();
    // vs_pnn->add("SM_NN",                   "NN score",                           100,   &binningSMNN[0], binningSMNN.size()-1);
    // vs_pnn->add("PNN500",                  "PNN500 score",                       100,   &binning2HDM500[0], binning2HDM500.size()-1);
    // vs_pnn->add("PNN1000",                 "PNN1000 score",                      100,   &binning2HDM1000[0], binning2HDM1000.size()-1);
    vs_pnn->add("PNN1100",                 "PNN1100 score",                      100,   &binning2HDM1100[0], binning2HDM1100.size()-1);

    Systematics* ss = new Systematics();
    // ss->add("THEO_ACC_StopWt_FSR", "THEO_ACC_StopWt_FSR", eSystematicType::TwoSide, bbtt_kBLUE4COMP, true);
    ss->add("Subtraction_bkg", "Subtraction_bkg", eSystematicType::TwoSide, bbtt_kBLUE4COMP);
    ss->add("ttReweighting", "ttReweighting", eSystematicType::OneSide, bbtt_kRED4COMP);

    CompInfo* info = new CompInfo();
    info->atlas = false;

    // info->ratio_high = 1.55;
    // info->ratio_low = 0.45;

    info->ratio_high = 1.75;
    info->ratio_low = 0.25;

    info->shape_only = false;
    info->logy = false;
    info->save_ratio = false;
    info->output_format = "pdf";
    info->legend_scaling_vertical = 0.8;
    info->legend_scaling_horizontal = 1.2;

    for (VariableInfo* v : *(vs_pnn->content()))
    {
        Processes* ps = new Processes();
        // ps->add("WZ",       "WZ",  eProcessType::BKG,  eProcess::DIBOSON,      "Diboson",  kBlack);
        // ps->add("WW",       "WW",  eProcessType::BKG,  eProcess::DIBOSON,      "Diboson",  kBlack);
        // ps->add("ZZ",       "ZZ",  eProcessType::BKG,  eProcess::DIBOSON,      "Diboson",  kBlack);

        ps->add("Fake",       "Fake",  eProcessType::BKG,  eProcess::FAKE,      "Fake",  kBlack);

        // ps->add("stopWt",       "stopWt",  eProcessType::BKG,  eProcess::FAKE,      "Single top Wt",  kBlack);

        Config* c = new Config(b, ps, rs, vs_pnn, ss);
        c->load(filename, "");
        info->parameter = "LepHad_LTT_PNN";
        c->updateHistogramPtr(rs->content()->front(), v);
        SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/CheckSyst/";

        if (HistToolHelper::check(c)) {
            ct->manipulate(c);
            ct->rebin(c, eRebinOption::Array, "", true);
        }

        if (ct->check(c))
        {
            ct->makeYield(c);
            ct->paint(c);
            ct->run(c);
        }
        else 
        {
            clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
        }

        delete ps;
        delete ct;
        delete c;
    }

    delete b;
    delete rs;
    delete vs_pnn;
    delete info;

}