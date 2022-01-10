#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "CompTool.h"
#include "SystCompTool.h"
#include "AutoBinningTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;
using BU = BinningUtils;

void thesis_syst_zhf(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("SR", "#sqrt{s} = 13 TeV");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "",        eRegionType::SR);

    Variables* vs_presel = new Variables();

    auto binningFromFile = [](const std::string& fn, bool BDT=false)
    { 
        return BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>(fn), 1000, BDT);
    };

    vector<double> binning2HDM251   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM260.txt");
    vector<double> binning2HDM260   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM260.txt");
    vector<double> binning2HDM280   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM280.txt");
    vector<double> binning2HDM300   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM300.txt");
    vector<double> binning2HDM325   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM300.txt");
    vector<double> binning2HDM350   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM350.txt");
    vector<double> binning2HDM400   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM400.txt");
    vector<double> binning2HDM450   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM450.txt");
    vector<double> binning2HDM500   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM500.txt");
    vector<double> binning2HDM550   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM550.txt");
    vector<double> binning2HDM600   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM600.txt");
    vector<double> binning2HDM700   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM700.txt");
    vector<double> binning2HDM800   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM800.txt");
    vector<double> binning2HDM900   =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM900.txt");
    vector<double> binning2HDM1000  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1000.txt");
    vector<double> binning2HDM1100  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1100.txt");
    vector<double> binning2HDM1200  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1200.txt");
    vector<double> binning2HDM1400  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1400.txt");
    vector<double> binning2HDM1600  =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1600.txt");
    vector<double> binningSMBDT     =   binningFromFile("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_SMBDT.txt", true);

    Tools::printVector(binningSMBDT); 

    Variables* vs_pnn = new Variables();
    vs_pnn->add("PNN251",                  "PNN251",                       100,   &binning2HDM251[0], binning2HDM251.size()-1);
    vs_pnn->add("PNN260",                  "PNN260",                       100,   &binning2HDM260[0], binning2HDM260.size()-1);
    vs_pnn->add("PNN280",                  "PNN280",                       100,   &binning2HDM280[0], binning2HDM280.size()-1);
    vs_pnn->add("PNN300",                  "PNN300",                       100,   &binning2HDM300[0], binning2HDM300.size()-1);
    vs_pnn->add("PNN325",                  "PNN325",                       100,   &binning2HDM325[0], binning2HDM325.size()-1);
    vs_pnn->add("PNN350",                  "PNN350",                       100,   &binning2HDM350[0], binning2HDM350.size()-1);
    vs_pnn->add("PNN400",                  "PNN400",                       100,   &binning2HDM400[0], binning2HDM400.size()-1);
    vs_pnn->add("PNN450",                  "PNN450",                       100,   &binning2HDM450[0], binning2HDM450.size()-1);
    vs_pnn->add("PNN500",                  "PNN500",                       100,   &binning2HDM500[0], binning2HDM500.size()-1);
    vs_pnn->add("PNN550",                  "PNN550",                       100,   &binning2HDM550[0], binning2HDM550.size()-1);
    vs_pnn->add("PNN600",                  "PNN600",                       100,   &binning2HDM600[0], binning2HDM600.size()-1);
    vs_pnn->add("PNN700",                  "PNN700",                       100,   &binning2HDM700[0], binning2HDM700.size()-1);
    vs_pnn->add("PNN800",                  "PNN800",                       100,   &binning2HDM800[0], binning2HDM800.size()-1);
    vs_pnn->add("PNN900",                  "PNN900",                       100,   &binning2HDM900[0], binning2HDM900.size()-1);
    vs_pnn->add("PNN1000",                 "PNN1000",                      100,   &binning2HDM1000[0], binning2HDM1000.size()-1);
    vs_pnn->add("PNN1100",                 "PNN1100",                      100,   &binning2HDM1100[0], binning2HDM1100.size()-1);
    vs_pnn->add("PNN1200",                 "PNN1200",                      100,   &binning2HDM1200[0], binning2HDM1200.size()-1);
    vs_pnn->add("PNN1400",                 "PNN1400",                      100,   &binning2HDM1400[0], binning2HDM1400.size()-1);
    vs_pnn->add("PNN1600",                 "PNN1600",                      100,   &binning2HDM1600[0], binning2HDM1600.size()-1);

    Variables* vs_bdt = new Variables();
    vs_bdt->add("SMBDT",                   "SM BDT",                       100,   &binningSMBDT[0], binningSMBDT.size()-1);

    CompInfo* info = new CompInfo();
    info->atlas = false;
    info->ratio_high = 1.88;
    info->ratio_low = 0.12;
    info->shape_only = true;
    info->logy = true;
    info->save_ratio = false;

    for (VariableInfo* v : *(vs_bdt->content()))
    {
        Systematics* ss = new Systematics();
        ss->add("ZJETS_ACC_GENERATOR_MBB", "GENERATOR", eSystematicType::TwoSide, kMagenta);
        ss->add("ZJETS_ACC_SCALE_MVA_" + v->name, "SCALE", eSystematicType::TwoSide, kGreen+1);
        
        Processes* ps = new Processes();
        ps->add("Zttbb",    "Z+hf",  eProcessType::BKG,  eProcess::ZllHF,    "Baseline",  kBlue+1);
        ps->add("Zttbc",    "Z+hf",  eProcessType::BKG,  eProcess::ZllHF,    "Baseline",  kBlue+1);
        ps->add("Zttcc",    "Z+hf",  eProcessType::BKG,  eProcess::ZllHF,    "Baseline",  kBlue+1);

        Config* c = new Config(b, ps, rs, vs_presel, ss);
        c->load(filename, "BDTScorePreselection");
        info->parameter = "BDT";
        c->updateHistogramPtr(rs->content()->front(), v);
        SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SystZhf/";

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

    for (VariableInfo* v : *(vs_pnn->content()))
    {
        Systematics* ss = new Systematics();
        ss->add("ZJETS_ACC_GENERATOR_MBB", "GENERATOR", eSystematicType::TwoSide, kMagenta);
        ss->add("ZJETS_ACC_SCALE_MVA_" + v->name, "SCALE", eSystematicType::TwoSide, kGreen+1);
        
        Processes* ps = new Processes();
        ps->add("Zttbb",    "Z+hf",  eProcessType::BKG,  eProcess::ZllHF,    "Baseline",  kBlue+1);
        ps->add("Zttbc",    "Z+hf",  eProcessType::BKG,  eProcess::ZllHF,    "Baseline",  kBlue+1);
        ps->add("Zttcc",    "Z+hf",  eProcessType::BKG,  eProcess::ZllHF,    "Baseline",  kBlue+1);

        Config* c = new Config(b, ps, rs, vs_presel, ss);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "PNN";
        c->updateHistogramPtr(rs->content()->front(), v);
        SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SystZhf/";

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
    delete vs_presel;
    delete vs_pnn;
    delete vs_bdt;
    delete info;

}