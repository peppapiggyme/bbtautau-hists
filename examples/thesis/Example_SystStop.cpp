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

void thesis_syst_stop(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("SR", "#sqrt{s} = 13 TeV");
    b->make_systs_shape_only = true;

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "SR",        eRegionType::SR);

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
    vs_pnn->add("PNN251",                  "PNN251 score",                       100,   &binning2HDM251[0], binning2HDM251.size()-1);
    vs_pnn->add("PNN260",                  "PNN260 score",                       100,   &binning2HDM260[0], binning2HDM260.size()-1);
    vs_pnn->add("PNN280",                  "PNN280 score",                       100,   &binning2HDM280[0], binning2HDM280.size()-1);
    vs_pnn->add("PNN300",                  "PNN300 score",                       100,   &binning2HDM300[0], binning2HDM300.size()-1);
    vs_pnn->add("PNN325",                  "PNN325 score",                       100,   &binning2HDM325[0], binning2HDM325.size()-1);
    vs_pnn->add("PNN350",                  "PNN350 score",                       100,   &binning2HDM350[0], binning2HDM350.size()-1);
    vs_pnn->add("PNN400",                  "PNN400 score",                       100,   &binning2HDM400[0], binning2HDM400.size()-1);
    vs_pnn->add("PNN450",                  "PNN450 score",                       100,   &binning2HDM450[0], binning2HDM450.size()-1);
    vs_pnn->add("PNN500",                  "PNN500 score",                       100,   &binning2HDM500[0], binning2HDM500.size()-1);
    vs_pnn->add("PNN550",                  "PNN550 score",                       100,   &binning2HDM550[0], binning2HDM550.size()-1);
    vs_pnn->add("PNN600",                  "PNN600 score",                       100,   &binning2HDM600[0], binning2HDM600.size()-1);
    vs_pnn->add("PNN700",                  "PNN700 score",                       100,   &binning2HDM700[0], binning2HDM700.size()-1);
    vs_pnn->add("PNN800",                  "PNN800 score",                       100,   &binning2HDM800[0], binning2HDM800.size()-1);
    vs_pnn->add("PNN900",                  "PNN900 score",                       100,   &binning2HDM900[0], binning2HDM900.size()-1);
    vs_pnn->add("PNN1000",                 "PNN1000 score",                      100,   &binning2HDM1000[0], binning2HDM1000.size()-1);
    vs_pnn->add("PNN1100",                 "PNN1100 score",                      100,   &binning2HDM1100[0], binning2HDM1100.size()-1);
    vs_pnn->add("PNN1200",                 "PNN1200 score",                      100,   &binning2HDM1200[0], binning2HDM1200.size()-1);
    vs_pnn->add("PNN1400",                 "PNN1400 score",                      100,   &binning2HDM1400[0], binning2HDM1400.size()-1);
    vs_pnn->add("PNN1600",                 "PNN1600 score",                      100,   &binning2HDM1600[0], binning2HDM1600.size()-1);

    Variables* vs_bdt = new Variables();
    vs_bdt->add("SMBDT",                   "BDT score",                       100,   &binningSMBDT[0], binningSMBDT.size()-1);

    Systematics* ss = new Systematics();
    ss->add("SingleTop_Wtchan_TopInterference", "Top_Interference", eSystematicType::OneSide, bbtt_kRED4COMP);

    CompInfo* info = new CompInfo();
    info->atlas = false;
    info->ratio_high = 1.68;
    info->ratio_low = 0.32;
    info->shape_only = true;
    info->logy = true;
    info->save_ratio = false;
    info->output_format = "pdf";
    info->legend_scaling_vertical = 1.0;
    info->legend_scaling_horizontal = 1.2;

    for (VariableInfo* v : *(vs_bdt->content()))
    {
        Processes* ps = new Processes();
        ps->add("stopWt",       "s-top Wt Nominal",  eProcessType::BKG,  eProcess::STOPWT,      "Baseline",  kBlack);

        Config* c = new Config(b, ps, rs, vs_bdt, ss);
        c->load(filename, "BDTScorePreselection");
        info->parameter = "BDT";
        c->updateHistogramPtr(rs->content()->front(), v);
        SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SystStop/";

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
        Processes* ps = new Processes();
        ps->add("stopWt",       "s-top Wt Nominal",  eProcessType::BKG,  eProcess::STOPWT,      "Baseline",  kBlack);

        Config* c = new Config(b, ps, rs, vs_pnn, ss);
        c->load(filename, "PNNScorePreselection");
        info->parameter = "PNN";
        c->updateHistogramPtr(rs->content()->front(), v);
        SystCompTool* ct = new SystCompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SystStop/";

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
    delete vs_bdt;
    delete info;

}

void thesis_syst_stop_param_topint(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("SR", "#sqrt{s} = 13 TeV");
    b->make_systs_shape_only = true;

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "SR",        eRegionType::SR);

    Variables* vs = new Variables();

    vs->add("pTBB",                   "m_{bb}",                    1);

    Systematics* ss = new Systematics();

    CompInfo* info = new CompInfo();
    info->atlas = false;
    info->ratio_high = 1.88;
    info->ratio_low = 0.12;
    info->shape_only = true;
    info->logy = false;
    info->save_ratio = false;
    info->output_format = "pdf";
    info->legend_scaling_horizontal = 1.3;
    info->ymax_ratio_nolog = 2.0;

    AutoBinningInfo* abi = new AutoBinningInfo();
    abi->n_bins = 6;
    AutoBinningTool* abt = new AutoBinningTool_v2(abi, eProcess::ZllHFMG);
    
    for (VariableInfo* v : *(vs->content()))
    {
        Processes* ps = new Processes();
        ps->add("stopWt",       "s-top Wt Nominal",  eProcessType::BKG,  eProcess::STOPWT,      "Baseline (DR)",    kBlack);
        ps->add("stopWtDS",     "s-top Wt DS",       eProcessType::BKG,  eProcess::STOPWT_DS,   "Alternative (DS)", bbtt_kRED4COMP);

        Config* c = new Config(b, ps, rs, vs, ss);
        c->load(filename, "Preselection");
        info->parameter = "TopInt";
        c->updateHistogramPtr(rs->content()->front(), v);
        CompTool* ct = new CompTool(info);
        ct->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SystStop_Param/";

        if (HistToolHelper::check(c)) {
            ct->manipulate(c);
            abt->run(c);
            abt->rebin(c);
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
    delete vs;
    delete info;

}