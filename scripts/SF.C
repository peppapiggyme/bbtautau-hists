
///////////////////////////////////////////////////////////////////////////////
// Settings
///////////////////////////////////////////////////////////////////////////////

// #define MAX_VALUE_Y 1.4
#define MAX_VALUE_Y 1.6

using namespace std;

static vector<int> colors{1, TColor::GetColor(224, 61, 59), TColor::GetColor(137, 90, 145), TColor::GetColor(36, 146, 224)};
static string lumi;
static string outfilename;
static void Load();
static void PlotSF();
static void PlotSF_Systs();

void SF()
{
    gROOT->SetStyle("ATLAS");
    gStyle->SetErrorX(0.5);
    gStyle->SetEndErrorSize(5);
    
    // PlotSF();`
    PlotSF_Systs();
}

void PlotSF()
{
    lumi = "#sqrt{s}=13 TeV, L=139 fb^{-1}";
    outfilename = "ttbar_sf_1P.pdf";

    TFile* f1 = TFile::Open("/scratchfs/atlas/bowenzhang/ttbarscalefactors/fits/sf.root");

    TH1F* h1 = (TH1F*)(f1->Get("TTBAR_FAKESF_OFFL_1P")->Clone());
    TH1F* h2 = (TH1F*)(f1->Get("TTBAR_FAKESF_TAU25_1P")->Clone());
    TH1F* h3 = (TH1F*)(f1->Get("TTBAR_FAKESF_TAU25EF_1P")->Clone());
    TH1F* h4 = (TH1F*)(f1->Get("TTBAR_FAKESF_TAU25RNN_1P")->Clone());

    TCanvas* c = new TCanvas("c", "", 1200, 900);

    TLegend* legend = new TLegend(0.60, 0.66, 0.92, 0.92);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(16);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.040);

    h1->GetXaxis()->SetTitleOffset(1.4);
    h1->GetXaxis()->SetTitleSize(0.055);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetXaxis()->SetTitle("#tau_{had} p_{T}");
    h1->GetYaxis()->SetTitle("Fake-#tau_{had} t#bar{t} Scale Factor");
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->SetTitleSize(0.055);
    h1->GetYaxis()->SetLabelSize(0.045);
    h1->SetMinimum(0.0);
    h1->SetMaximum(MAX_VALUE_Y);

    h1->SetLineColor(colors[0]);
    h1->SetLineWidth(2);
    h1->SetMarkerColor(colors[0]);
    h1->SetMarkerStyle(kFullSquare);
    h1->SetMarkerSize(1.3);

    h2->SetLineColor(colors[1]);
    h2->SetLineWidth(2);
    h2->SetLineStyle(2);
    h2->SetMarkerColor(colors[1]);
    h2->SetMarkerStyle(kOpenSquare);
    h2->SetMarkerSize(1.3);

    h3->SetLineColor(colors[2]);
    h3->SetLineWidth(2);
    h3->SetLineStyle(4);
    h3->SetMarkerColor(colors[2]);
    h3->SetMarkerStyle(kFullCircle);
    h3->SetMarkerSize(1.3);

    h4->SetLineColor(colors[3]);
    h4->SetLineWidth(2);
    h4->SetLineStyle(7);
    h4->SetMarkerColor(colors[3]);
    h4->SetMarkerStyle(kOpenCircle);
    h4->SetMarkerSize(1.3);

    h1->Draw("E1");
    h2->Draw("E1 SAME");
    h3->Draw("E1 SAME");
    h4->Draw("E1 SAME");

    legend->AddEntry(h1, "SF(off)", "lep");
    legend->AddEntry(h2, "SF(tau25)", "lep");
    legend->AddEntry(h3, "SF(tau25 EF)", "lep");
    legend->AddEntry(h4, "SF(tau25 EForMVA)", "lep");
    legend->Draw();

    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.040);
    text->DrawLatex(0.20, 0.86, "t#bar{t}-ID-CR");
    text->DrawLatex(0.20, 0.96, lumi.c_str());

    stringstream output;

    output << "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SF/" << outfilename;

    c->Update();
    c->SaveAs(output.str().c_str());

    delete legend;
    delete text;
    delete c;
}

void PlotSF_Systs()
{
    lumi = "#sqrt{s}=13 TeV, L=139 fb^{-1}";
    outfilename = "ttbar_sf_3P_syst.pdf";

    TFile* f1 = TFile::Open("/scratchfs/atlas/bowenzhang/ttbarscalefactors/fits/sf.root");

    TH1F* h = (TH1F*)(f1->Get("TTBAR_FAKESF_TAU25_3P")->Clone());
    TH1F* herr_1_up = (TH1F*)(f1->Get("TTBAR_FAKESF_TAU25_3P_SysTTBAR_FAKESF_TAU25_EIGEN0__1up")->Clone());
    TH1F* herr_1_dn = (TH1F*)(f1->Get("TTBAR_FAKESF_TAU25_3P_SysTTBAR_FAKESF_TAU25_EIGEN0__1down")->Clone());
    TH1F* herr_2_up = (TH1F*)(f1->Get("TTBAR_FAKESF_TAU25_3P_SysTTBAR_FAKESF_TAU25_EIGEN1__1up")->Clone());
    TH1F* herr_2_dn = (TH1F*)(f1->Get("TTBAR_FAKESF_TAU25_3P_SysTTBAR_FAKESF_TAU25_EIGEN1__1down")->Clone());
    TH1F* herr_1 = (TH1F*)(h->Clone());
    TH1F* herr_2 = (TH1F*)(h->Clone());
    for (int i = 0; i < h->GetNbinsX()+2; i++)
    {
        herr_1->SetBinError(i, 0.5 * TMath::Abs(herr_1_up->GetBinContent(i) - herr_1_dn->GetBinContent(i)));
        herr_2->SetBinError(i, 0.5 * TMath::Abs(herr_2_up->GetBinContent(i) - herr_2_dn->GetBinContent(i)));
    }

    TCanvas* c = new TCanvas("c", "", 1200, 900);

    TLegend* legend = new TLegend(0.60, 0.66, 0.92, 0.92);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(16);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.040);

    herr_1->GetXaxis()->SetTitleOffset(1.4);
    herr_1->GetXaxis()->SetTitleSize(0.055);
    herr_1->GetXaxis()->SetLabelSize(0.045);
    herr_1->GetXaxis()->SetTitle("#tau_{had} p_{T}");
    herr_1->GetYaxis()->SetTitle("Fake-#tau_{had} t#bar{t} Scale Factor");
    herr_1->GetYaxis()->SetTitleOffset(1.4);
    herr_1->GetYaxis()->SetTitleSize(0.055);
    herr_1->GetYaxis()->SetLabelSize(0.045);
    herr_1->SetMinimum(0.0);
    herr_1->SetMaximum(MAX_VALUE_Y);

    h->SetLineColor(colors[0]);
    h->SetLineWidth(2);
    h->SetMarkerColor(colors[0]);
    h->SetMarkerStyle(kFullSquare);
    h->SetMarkerSize(1.3);

    herr_1->SetLineWidth(0);
    herr_2->SetLineWidth(0);
    herr_1->SetFillColorAlpha(colors[1], 0.5);
    herr_2->SetFillColorAlpha(colors[2], 0.5);

    herr_1->Draw("E2 SAME");
    herr_2->Draw("E2 SAME");
    h->Draw("E1 SAME");

    legend->AddEntry(h, "SF(tau25)", "lep");
    legend->AddEntry(herr_1, "Eigen-variation 0", "f");
    legend->AddEntry(herr_2, "Eigen-variation 1", "f");
    legend->Draw();

    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.040);
    text->DrawLatex(0.20, 0.86, "t#bar{t}-ID-CR");
    text->DrawLatex(0.20, 0.96, lumi.c_str());

    stringstream output;

    output << "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SF/" << outfilename;

    c->Update();
    c->SaveAs(output.str().c_str());

    delete legend;
    delete text;
    delete c;
}