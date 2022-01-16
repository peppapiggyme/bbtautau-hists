#include "Utils.h"

#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

using namespace std;
using namespace RooStats;
using namespace RooFit;
using namespace HistFactory;

int main(int argc, char** argv)
{
    if (argc != 3) 
    {
        cerr << "usage:\n\tul-roostat <index_calc> <workspace>\n";
    }

    string type{argv[1]};
    TFile* f = TFile::Open(argv[2]);
    RooWorkspace* ws = static_cast<RooWorkspace*>(f->Get("combined"));

    // Data
    // ----
    RooAbsData *data = ws->data("obsData");

    if (!data)
        throw std::runtime_error("No data!");

    // S+B model
    // ---------
    ModelConfig *sbModel = (ModelConfig*)ws->obj("ModelConfig");
    // no systematic uncertainties
    const RooArgSet *nuisPar = sbModel->GetNuisanceParameters();
    if (nuisPar && nuisPar->getSize() > 0)
    {
        cout << "Switching off nuisance parameters by setting them constant to their default values";
        SetAllConstant(*nuisPar);
    }

    // In this case, B model is S+B model with mu = 0
    // ----------------------------------------------
    ModelConfig *bModel = (ModelConfig*)sbModel->Clone();
    bModel->SetName(TString(sbModel->GetName()) + TString("_with_poi_0"));
    RooRealVar *var = dynamic_cast<RooRealVar *>(bModel->GetParametersOfInterest()->first());
    if (!var)
        throw std::logic_error("There must be at least 1 POI in the model!");
    /// @todo: Need to understand why SetSnapShot
    double oldval = var->getVal();
    var->setVal(0);
    bModel->SetSnapshot(RooArgSet(*var));
    var->setVal(oldval);

    // Define calculator for hypothesis test
    // Here we use FrequentistCalculator
    // One can choose from AsymptoticCalculator, HybridCalculator
    // ----------------------------------------------------------
    FrequentistCalculator fc(*data, *bModel, *sbModel);
    AsymptoticCalculator ac(*data, *bModel, *sbModel);

    // Define test statistic
    // Use toy MC to generate the distribution of profiled likelihood ratio test statistic
    // -----------------------------------------------------------------------------------
    ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
    profll.SetOneSided(true);
    profll.SetMinimizer("Minuit2");
    profll.SetStrategy(0);
    profll.SetTolerance(1e-7);

    if (type == "1")
    {
        ToyMCSampler *toymcs = (ToyMCSampler*)fc.GetTestStatSampler();
        toymcs->SetTestStatistic(&profll);
        toymcs->SetNToys(10000);
    }
    
    // Define hypothesis test inverter (for upper limits)
    // --------------------------------------------------
    HypoTestInverter* calc_fc = new HypoTestInverter(fc);
    HypoTestInverter* calc_ac = new HypoTestInverter(ac);
    HypoTestInverter* calc = type == "1" ? calc_fc : calc_ac;

    calc->SetConfidenceLevel(0.95);
    calc->UseCLs(true);

    // Configure the POI scanning
    // --------------------------
    int npoints = 21;
    double poimin = 0.0;
    double poimax = 2.0;

    calc->SetFixedScan(npoints, poimin, poimax);

    // Calculate interval (upper limits)
    // ---------------------------------
    HypoTestInverterResult *r = calc->GetInterval();

    // Extract results
    // ---------------
    double upperLimit = r->UpperLimit();
    double expectedLimit_m2 = r->GetExpectedUpperLimit(-2);
    double expectedLimit_m1 = r->GetExpectedUpperLimit(-1);
    double expectedLimit    = r->GetExpectedUpperLimit(0);
    double expectedLimit_p1 = r->GetExpectedUpperLimit(1);
    double expectedLimit_p2 = r->GetExpectedUpperLimit(2);

    Tools::println("Observed / Expected upper limits: % / [%, %, %, %, %]", 
        upperLimit, expectedLimit_m2, expectedLimit_m1, expectedLimit, expectedLimit_p1, expectedLimit_p2);

    HypoTestInverterPlot *plot = new HypoTestInverterPlot("hi", "", r);

    TCanvas *c = new TCanvas("hypoInv", "", 800, 600);
    plot->Draw();
    c->SaveAs("hypo_scan.pdf");
    //for (int i = 0; i < r->ArraySize(); i++)
    //{
        SamplingDistribution *limDist = calc->GetUpperLimitDistribution(true, 10000);

        SamplingDistPlot limPlot(100);
        limPlot.AddSamplingDistribution(limDist);
        limPlot.GetTH1F()->SetStats(true); // display statistics
        limPlot.SetLineColor(kBlue);
        new TCanvas("limPlot", "Upper Limit Distribution");
        limPlot.Draw();
 
        /// save result in a file
        limDist->SetName("RULDist");
        TFile *fileOut = new TFile("RULDist.root", "RECREATE");
        limDist->Write();
        fileOut->Close();

        //c->SaveAs(string("hypo_scan_" + to_string(i) + ".pdf").c_str());
    //}

    f->Close();
    f = nullptr;
    delete calc_ac;
    delete calc_fc;
    delete plot;
    delete c;

    return 0;
}
