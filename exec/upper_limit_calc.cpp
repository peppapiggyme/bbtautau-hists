#include "CommonInclude.h"
#include "Utils.h"

#include "TRandom3.h"

#define EXIT_HERE return;

using namespace std;

int main(int argc, char** argv)
{
    if (argc != 4) {
        fprintf(stderr, "usage:\n\tul <ndata> <nbkg> <nsig>\n");
        exit(0);
    }
    else {
        fprintf(stdout, "\n\nUpper limit on counting experiment using toys\n\n");
    }

    // plotting
    TCanvas * c1 = new TCanvas( "c1_ex4", "m_{4l}", 800, 600);

    // total number of events
    double nData{std::atoi(argv[1])};
    double nBkg{std::atoi(argv[2])};
    double nSig{std::atoi(argv[3])};

    double nBkgErr = 0;

    // Hypothesis test
    // H0 (null hypothesis) : mu = mu0
    // H1 (alternative hypothesis) : mu < mu0
    double mu = 0.0;
    double nScan = 15;
    double muMin = 0.0;
    double muMax = 1.4;
    int nCounter = 0;
    double x[(int)nScan];
    double yObs[(int)nScan];
    double yExp[(int)nScan];
    unsigned long nToys = 1000000; // 1e6 x 0.05 = 5e4

    while (nCounter < nScan)
        {
            Tools::println("Scanning mu = %", mu);
            // generate toy mc
            TRandom3 r(nCounter);
        
            // number of events is our test statistic here
            // generate N's distribution
            TH1D* hEvent = new TH1D("Ex4_h", "", 150, 0, 150);
            for (unsigned long i = 0; i < nToys; i++)
                {
                    double nuRandom = r.Gaus(0, 1);
                    double b = nBkg * (1. + nuRandom * nBkgErr);
                    double s = nSig;
                    double nRandom = r.Poisson(mu * s + b);
                    hEvent->Fill(nRandom);
                }
            TH1D* hEventBkgOnly = new TH1D("Ex4_h_bkgOnly", "", 150, 0, 150);
            for (unsigned long i = 0; i < nToys; i++)
                {
                    double nuRandom = r.Gaus(0, 1);
                    double b = nBkg * (1. + nuRandom * nBkgErr);
                    double nRandom = r.Poisson(b);
                    hEventBkgOnly->Fill(nRandom);
                }

            hEvent->SetTitle(";N;Event");
            hEvent->Draw();
            hEventBkgOnly->SetLineColor(4);
            hEventBkgOnly->Draw("same");
            // expected N
            double nExp{0}, q{0.5};
            hEventBkgOnly->GetQuantiles(1, &nExp, &q);
            TLine* lExp = new TLine(nExp, 0., nExp, hEvent->GetMaximum());
            lExp->SetLineColor(4);
            lExp->SetLineWidth(2);
            lExp->Draw("SAME");
            // observed N
            double nObs = nData;
            TLine* lObs = new TLine(nObs, 0., nObs, hEvent->GetMaximum());
            lObs->SetLineColor(2);
            lObs->SetLineWidth(2);
            lObs->Draw("SAME");

            TText tex(.5, .8, string("#mu = " + to_string(mu)).c_str());
            tex.SetNDC();
            tex.Draw();

            c1->SetLogy(true);
            c1->SaveAs(string("Ex4_N_" + to_string(nCounter) + ".png").c_str());

            auto pValue = [&hEvent](double n) {
                return hEvent->Integral(0, hEvent->FindBin(n)) / hEvent->Integral();
            };

            x[nCounter] = mu;
            yObs[nCounter] = pValue(nObs);
            yExp[nCounter] = pValue(nExp);

            Tools::println("Summary: mu=%, obs p-value=%, exp p-value=%", mu, yObs[nCounter], yExp[nCounter]);

            mu += (muMax - muMin) / (nScan - 1);
            nCounter++;

            delete hEvent;
            delete hEventBkgOnly;
            delete lExp;
            delete lObs;
        }

    TGraph* gObs = new TGraph((int)nScan, x, yObs);
    TGraph* gExp = new TGraph((int)nScan, x, yExp);
    TGraph* gObsForEval = new TGraph((int)nScan, yObs, x);
    TGraph* gExpForEval = new TGraph((int)nScan, yExp, x);
    gObs->SetMinimum(0.0);
    gObs->SetTitle(";#mu;p-value");

    // styles
    gObs->SetLineColor(kBlack);
    gObs->SetLineWidth(2);
    gObs->SetMarkerSize(0);
    gExp->SetLineColor(kBlue);
    gExp->SetLineWidth(2);
    gExp->SetMarkerSize(0);

    TLine* pValue_0p05 = new TLine(gObs->GetXaxis()->GetXmin(), 0.05, gObs->GetXaxis()->GetXmax(), 0.05);
    pValue_0p05->SetLineColor(kRed);
    pValue_0p05->SetLineWidth(2);

    TCanvas* c = new TCanvas("Ex4_pValue", "mu Scan", 800, 600);
    gObs->Draw();
    gExp->Draw("SAME");
    pValue_0p05->Draw("SAME");
    c->SaveAs("Ex4_MuScan.png");

    Tools::println("Upper limit = %, Expected upper limit = %", gObsForEval->Eval(0.05), gExpForEval->Eval(0.05));

    // clear up
    delete gObs;
    delete gExp;
    delete gObsForEval;
    delete gExpForEval;
    delete pValue_0p05;
    delete c1;
}
