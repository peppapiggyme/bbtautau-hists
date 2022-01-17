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

static string channel = "SLT";
static string version = "v13";

void test_check_syst(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("SR", "#sqrt{s} = 13 TeV");
    b->name_convention = Utils::NameConvention::WSMaker;

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "LepHad " + channel + " " + version,        eRegionType::SR);

    auto binningFromFile = [](const std::string& fn)
    { 
        return BU::intToDoubleBinEdgesForMVAInverse_LepHad1090Version(BU::readBinningFromFile<int>(fn));
    };

    vector<double> binning2HDM1000 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM1000." + version + ".txt");
    vector<double> binning2HDM1100 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM1100." + version + ".txt");

    Tools::printVector(binning2HDM1000);
    Tools::printVector(binning2HDM1100);

    Variables* vs_pnn = new Variables();
    vs_pnn->add("PNN1000",                 "PNN1000 score",                      1,   &binning2HDM1000[0], binning2HDM1000.size()-1);
    // vs_pnn->add("PNN1100",                 "PNN1100 score",                      1,   &binning2HDM1100[0], binning2HDM1100.size()-1);

    Systematics* ss = new Systematics();
    // ss->add("THEO_ACC_StopWt_FSR", "THEO_ACC_StopWt_FSR", eSystematicType::TwoSide, bbtt_kBLUE4COMP, true);
    // ss->add("THEO_ACC_TTBAR_ME", "THEO_ACC_TTBAR_ME", eSystematicType::OneSide, bbtt_kBLUE4COMP);
    // ss->add("THEO_ACC_TTBAR_PS_SLT", "THEO_ACC_TTBAR_PS_SLT", eSystematicType::OneSide, bbtt_kRED4COMP);
    // ss->add("Subtraction_bkg", "Subtraction_bkg", eSystematicType::TwoSide, bbtt_kBLUE4COMP);
    ss->add("FFVarrQCD", "FFVarrQCD", eSystematicType::TwoSide, bbtt_kBLUE4COMP);
    ss->add("ttReweighting", "ttReweighting", eSystematicType::OneSide, bbtt_kRED4COMP);

    CompInfo* info = new CompInfo();
    info->atlas = false;

    // info->ratio_high = 1.55;
    // info->ratio_low = 0.45;

    info->ratio_high = 1.8;
    info->ratio_low = 0.2;

    info->shape_only = false;
    info->logy = true;
    info->save_ratio = false;
    info->output_format = "png";
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

        // ps->add("ttbar",       "ttbar",  eProcessType::BKG,  eProcess::TTBAR,      "ttbar",  kBlack);

        Config* c = new Config(b, ps, rs, vs_pnn, ss);
        c->load(filename, "");
        info->parameter = "LepHad_" + channel + "_PNN." + version + "_Input." + version + "_Binning";
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