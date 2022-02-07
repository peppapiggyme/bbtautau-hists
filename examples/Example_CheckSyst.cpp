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
#include <map>

using std::cout; 
using std::endl;
using std::clog;
using BU = BinningUtils;
using std::map;

static string channel = "SLT";
static string version = "v14";

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

    vector<double> binningSMNN = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_SM_NN." + version + ".txt");
    vector<double> binning2HDM251 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM251." + version + ".txt");
    vector<double> binning2HDM260 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM260." + version + ".txt");
    vector<double> binning2HDM280 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM280." + version + ".txt");
    vector<double> binning2HDM300 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM300." + version + ".txt");
    vector<double> binning2HDM325 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM325." + version + ".txt");
    vector<double> binning2HDM350 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM350." + version + ".txt");
    vector<double> binning2HDM375 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM375." + version + ".txt");
    vector<double> binning2HDM400 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM400." + version + ".txt");
    vector<double> binning2HDM450 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM450." + version + ".txt");
    vector<double> binning2HDM500 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM500." + version + ".txt");
    vector<double> binning2HDM550 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM550." + version + ".txt");
    vector<double> binning2HDM600 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM600." + version + ".txt");
    vector<double> binning2HDM700 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM700." + version + ".txt");
    vector<double> binning2HDM800 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM800." + version + ".txt");
    vector<double> binning2HDM900 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM900." + version + ".txt");
    vector<double> binning2HDM1000 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM1000." + version + ".txt");
    vector<double> binning2HDM1100 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM1100." + version + ".txt");
    vector<double> binning2HDM1200 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM1200." + version + ".txt");
    vector<double> binning2HDM1400 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM1400." + version + ".txt");
    vector<double> binning2HDM1600 = binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/" + version + "/Binning_LepHad_" + channel + "_2HDM1600." + version + ".txt");

    map<string, vector<double>> files = {
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_SM_NN.root", binningSMNN},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_251.root", binning2HDM251},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_260.root", binning2HDM260},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_280.root", binning2HDM280},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_300.root", binning2HDM300},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_325.root", binning2HDM325},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_350.root", binning2HDM350},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_375.root", binning2HDM375},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_400.root", binning2HDM400},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_450.root", binning2HDM450},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_500.root", binning2HDM500},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_550.root", binning2HDM550},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_600.root", binning2HDM600},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_700.root", binning2HDM700},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_800.root", binning2HDM800},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_900.root", binning2HDM900},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_1000.root", binning2HDM1000},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_1100.root", binning2HDM1100},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_1200.root", binning2HDM1200},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_1400.root", binning2HDM1400},
        {"/publicfs/atlas/atlasnew/higgs/hh2X/zhangbw/WSI/2022-01-22-Yanlin-v14/SLT211222/13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_1600.root", binning2HDM1600},
    };

    Tools::printVector(binningSMNN);
    Tools::printVector(binning2HDM1000);
    Tools::printVector(binning2HDM1100);

    Variables* vs_pnn = new Variables();
    string varname;
    if (filename.find("SM_NN") != string::npos) varname = "SM_NN";
    else varname = filename.substr(filename.find("PNN"), filename.find(".") - filename.find("PNN"));
    vs_pnn->add(varname,                 varname,                      1,   &(files[filename][0]), files[filename].size()-1);

    Systematics* ss = new Systematics();
    // ss->add("THEO_ACC_StopWt_FSR", "THEO_ACC_StopWt_FSR", eSystematicType::TwoSide, bbtt_kBLUE4COMP, true);
    // ss->add("THEO_ACC_TTBAR_ME", "THEO_ACC_TTBAR_ME", eSystematicType::OneSide, bbtt_kBLUE4COMP);
    // ss->add("THEO_ACC_TTBAR_PS_SLT", "THEO_ACC_TTBAR_PS_SLT", eSystematicType::OneSide, bbtt_kRED4COMP);
    // ss->add("Subtraction_bkg", "Subtraction_bkg", eSystematicType::TwoSide, bbtt_kGREEN);
    ss->add("FFVarrQCD", "FFVarrQCD", eSystematicType::TwoSide, kBlue);
    ss->add("FFVarrQCD", "FFVarrQCD_smoothed", eSystematicType::TwoSide, kRed);
    // ss->add("ttReweighting", "ttReweighting", eSystematicType::OneSide, bbtt_kRED4COMP);

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
        ps->add("Fake",       "Fake",  eProcessType::BKG,  eProcess::FAKE,      "Fake",  kBlack);

        Config* c = new Config(b, ps, rs, vs_pnn, ss);
        c->load(filename, "");
        // info->parameter = "rQCD_LepHad_" + channel + "_SMNN." + version + "_Input." + version + "_Binning";
        info->parameter = "rQCDSmoothing_" + channel + "_PNN." + version + "_Input." + version + "_Binning";
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