
void SampDistHelper(const std::string& filename);

void SampDist()
{
    gROOT->SetStyle("ATLAS");
    gStyle->SetErrorX(0.5);

    vector<std::string> files = {
        "pll_test_stat_distribution.BoostedX1200.100k.mu4.pdf",
        "pll_test_stat_distribution.BoostedX1400.100k.mu4.pdf",
        "pll_test_stat_distribution.BoostedX2000.100k.mu4.pdf",
        "pll_test_stat_distribution.BoostedX3000.100k.mu4.pdf",
    };

    for (auto& f : files)
    {
        SampDistHelper(f);
    }
}

void SampDistHelper(const std::string& filename) 
{
    TFile* cRootFile = TFile::Open(("/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SamplingDist/" + filename + ".root").c_str(), "read");

    TCanvas *cCanvas = new TCanvas("Sampling", "", 1200, 900);
    cCanvas->SetLogy();
    TH1F* cHistPlot = (TH1F*)cRootFile->Get("plotting");
    cHistPlot->SetDirectory(0);

    TH1F* cHistDummy = new TH1F("dummy", "", 1, 0, 4.5);
    cHistDummy->GetXaxis()->SetLabelSize(0.04);
    cHistDummy->GetXaxis()->SetTitleSize(0.045);
    cHistDummy->GetXaxis()->SetTitleOffset(1.2);
    cHistDummy->GetXaxis()->SetTitle(Form("-log#lambda(#mu=%.2f)", 0.4));
    cHistDummy->GetXaxis()->SetRangeUser(0, 4.5);
    cHistDummy->GetYaxis()->SetLabelSize(0.04);
    cHistDummy->GetYaxis()->SetTitleSize(0.045);
    cHistDummy->GetYaxis()->SetTitle("f(-log#lambda(#mu)|#mu)");
    cHistDummy->Draw("AXIS");
    
    cHistPlot->SetLineWidth(2);
    cHistPlot->SetMarkerStyle(kFullCircle);
    cHistPlot->Draw("E0 SAME");

    TH1F* cHistHalfChi2 = (TH1F*)cRootFile->Get("half_chi2");

    double fSampleDist = cHistPlot->Integral();
    cHistHalfChi2->Scale(1. / (cHistHalfChi2->Integral() * cHistHalfChi2->GetBinWidth(1)));
    cHistHalfChi2->SetLineColor(kRed+1);
    cHistHalfChi2->SetMarkerColor(kRed+1);
    cHistHalfChi2->SetMarkerStyle(kOpenCircle);
    cHistHalfChi2->SetLineWidth(2);
    cHistHalfChi2->Draw("E0 SAME");

    TLegend* legend = new TLegend(0.40, 0.78, 0.90, 0.92);
    legend->SetNColumns(2);
    legend->SetTextFont(42);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.036);
    legend->SetTextAlign(12);
    legend->AddEntry(cHistPlot, "Pseudo experiment", "lep");
    legend->AddEntry(cHistHalfChi2, "Approximation", "lep");
    legend->AddEntry("", "Number of toys = 10^{5}", "");
    legend->Draw();

    cCanvas->SetLogy();
    cCanvas->Update();
    /// @todo fix hardcoded output path
    string sOutputPath = \
        "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SamplingDist/" + filename;
    cCanvas->SaveAs(sOutputPath.c_str());


    delete cCanvas, cHistDummy, legend;
}