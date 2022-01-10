
///////////////////////////////////////////////////////////////////////////////
// Settings
///////////////////////////////////////////////////////////////////////////////
using namespace std;

static map<string, string> filenames;
static map<string, string> lumis;
static const string format{".pdf"};
static vector<string> tags{"1"};
static vector<string> prongs{"1", "3"};
static vector<string> vars{"UnbinLeadTau", "UnbinSubleadTau"};

static void PlotFF_DTT(const pair<string, string>& p, const char* folder);
static void PlotFF_STT(const pair<string, string>& p, const char* folder);
static void PlotTF(const pair<string, string>& p, const char* folder);

///////////////////////////////////////////////////////////////////////////////
// For fun
///////////////////////////////////////////////////////////////////////////////
#define __BEGINLOOP_TAGS \
    for (auto& tag : tags) \
    {

#define __BEGINLOOP_PRONGS \
    for (auto& prong : prongs) \
    {

#define __BEGINLOOP_VARS \
    for (auto& var : vars) \
    {

#define __ENDLOOP \
    }


///////////////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////////////
void FF() 
{
    gROOT->SetStyle("ATLAS");
    gStyle->SetErrorX(0.5);

    const string datapath = \
        "/scratchfs/atlas/bowenzhang/CxAODReaderNonRes/source/CxAODReader_HH_bbtautau/data/";

    filenames = {
        {"2015-2016", datapath + "FF_1d_etadep_15_16.root"},
        {"2017", datapath + "FF_1d_etadep_17.root"},
        {"2018", datapath + "FF_1d_etadep_18.root"},
    };

    lumis = {
        {"2015-2016", "#sqrt{s}=13 TeV, L=36.2 fb^{-1}"},
        {"2017", "#sqrt{s}=13 TeV, L=44.3 fb^{-1}"},
        {"2018", "#sqrt{s}=13 TeV, L=58.5 fb^{-1}"},
    };

    for (const auto& p : filenames)
    {
        PlotFF_DTT(p, "TauFFHadHad1D_DTT");
        PlotFF_STT(p, "TauFFHadHad1D_STT");
        PlotTF(p, "TauFFHadHad1D");
    }
}

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////
static void PlotFF_DTT(const pair<string, string>& p, const char* folder)
{
    TFile* f = TFile::Open(p.second.c_str(), "READ");
    
    __BEGINLOOP_TAGS
    __BEGINLOOP_PRONGS

    TCanvas* c = new TCanvas("c", "", 1200, 900);
    c->SetLogx();

    stringstream ss_barrel;
    stringstream ss_endcap;
    ss_barrel << folder << '/' << "FF_" << tag << "tag_TauPtBarrel_" << prong << "P";
    ss_endcap << folder << '/' << "FF_" << tag << "tag_TauPtEndcap_" << prong << "P";
    TH1* h_barrel = (TH1*)f->Get(ss_barrel.str().c_str());        
    TH1* h_endcap = (TH1*)f->Get(ss_endcap.str().c_str());

    string xtitle = prong + "-prong #tau_{had} p_T";
    h_barrel->GetXaxis()->SetTitle(xtitle.c_str());
    h_barrel->GetXaxis()->SetTitleOffset(1.4);
    h_barrel->GetXaxis()->SetTitleSize(0.055);
    h_barrel->GetXaxis()->SetLabelSize(0.045);
    h_barrel->GetYaxis()->SetTitle("Fake Factor");
    h_barrel->GetYaxis()->SetTitleOffset(1.4);
    h_barrel->GetYaxis()->SetTitleSize(0.055);
    h_barrel->GetYaxis()->SetLabelSize(0.045);
    h_barrel->SetMinimum(0);
    h_barrel->SetMaximum((int)(150 * max(h_barrel->GetMaximum(), h_endcap->GetMaximum())) / 100.);

    auto color_1 = TColor::GetColor(83, 173, 176);
    auto color_2 = TColor::GetColor(10, 73, 94);

    h_barrel->SetLineColor(color_1);
    h_barrel->SetLineWidth(2);
    h_barrel->SetMarkerColor(color_1);
    h_barrel->SetMarkerStyle(kFullSquare);
    h_barrel->SetMarkerSize(1.3);
    h_endcap->SetLineColor(color_2);
    h_endcap->SetLineWidth(2);
    h_endcap->SetMarkerColor(color_2);
    h_endcap->SetMarkerStyle(kFullCircle);
    h_endcap->SetMarkerSize(1.3);
    
    h_barrel->Draw("E0 SAME");
    h_endcap->Draw("E0 SAME");

    TLegend* legend = new TLegend(0.4, 0.8, 0.92, 0.92);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(16);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.040);
    string str_barrel{"FF(DTT, 0<|#eta|<1.37, " + prong + "-prong)"};
    string str_endcap{"FF(DTT, 1.52<|#eta|<2.5, " + prong + "-prong)"};
    legend->AddEntry(h_barrel, str_barrel.c_str(), "lp");
    legend->AddEntry(h_endcap, str_endcap.c_str(), "lp");
    legend->Draw();

    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.040);
    text->DrawLatex(0.20, 0.86, "1-b-tagged SS");
    text->DrawLatex(0.20, 0.96, string("Data in " + p.first + ", " + lumis[p.first]).c_str());

    stringstream output;

    output << "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/FF/FF_DTT_"
            << tag << "tag_" << prong << "prong_" << p.first << format;

    c->SaveAs(output.str().c_str());

    delete legend;
    delete text;
    delete c;
    
    __ENDLOOP
    __ENDLOOP
}

static void PlotFF_STT(const pair<string, string>& p, const char* folder)
{
    TFile* f = TFile::Open(p.second.c_str(), "READ");

    __BEGINLOOP_TAGS

    TH1* h_out = new TH1F(p.first.c_str(), "", 4, 0, 4);
    TCanvas* c = new TCanvas("c", "", 1200, 900);

    // 2 x 2
    size_t iter = 1;
    __BEGINLOOP_PRONGS
    __BEGINLOOP_VARS

    stringstream ss_name;
    ss_name << folder << '/' << "FF_" << tag << "tag_" << var << "_" << prong << "P";
    TH1* h = (TH1*)f->Get(ss_name.str().c_str());

    string xlabel = \
        var == "UnbinLeadTau" ? string("Leading, ") + prong + "-prong" : string("Sub-leading, ") + prong + "-prong";

    h_out->SetBinContent(iter, h->GetBinContent(1));
    h_out->SetBinError(iter, h->GetBinError(1));
    h_out->GetXaxis()->SetBinLabel(iter, xlabel.c_str());

    iter++;

    __ENDLOOP
    __ENDLOOP

    h_out->GetXaxis()->SetTitleOffset(1.4);
    h_out->GetXaxis()->SetTitleSize(0.055);
    h_out->GetXaxis()->SetLabelSize(0.045);
    h_out->GetYaxis()->SetTitle("Fake Factor");
    h_out->GetYaxis()->SetTitleOffset(1.4);
    h_out->GetYaxis()->SetTitleSize(0.055);
    h_out->GetYaxis()->SetLabelSize(0.045);
    h_out->SetMinimum(0);
    h_out->SetMaximum((int)(150 * h_out->GetMaximum()) / 100.);

    auto color = TColor::GetColor(10, 73, 94);

    h_out->SetLineColor(color);
    h_out->SetLineWidth(2);
    h_out->SetMarkerColor(color);
    h_out->SetMarkerStyle(kFullSquare);
    h_out->SetMarkerSize(1.3);
    
    h_out->Draw("E0 SAME");

    TLegend* legend = new TLegend(0.66, 0.86, 0.92, 0.92);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(16);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.040);
    string str_barrel{"FF(STT)"};
    legend->AddEntry(h_out, str_barrel.c_str(), "lp");
    legend->Draw();

    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.040);
    text->DrawLatex(0.20, 0.86, "1-b-tagged SS");
    text->DrawLatex(0.20, 0.96, string("Data in " + p.first + ", " + lumis[p.first]).c_str());

    stringstream output;

    output << "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/FF/FF_STT_"
            << tag << "tag_" << p.first << format;

    c->SaveAs(output.str().c_str());

    delete h_out;
    delete legend;
    delete text;
    delete c;
    
    __ENDLOOP
}

static void PlotTF(const pair<string, string>& p, const char* folder)
{
    TFile* f = TFile::Open(p.second.c_str(), "READ");

    TH1* h_out = new TH1F(p.first.c_str(), "", 4, 0, 4);
    TCanvas* c = new TCanvas("c", "", 1200, 900);

    // 2 x 2
    size_t iter = 1;
    __BEGINLOOP_PRONGS
    __BEGINLOOP_VARS

    stringstream ss_name;
    ss_name << folder << '/' << "TF_1to2_" << var << "_" << prong << "P";
    TH1* h = (TH1*)f->Get(ss_name.str().c_str());

    string xlabel = \
        var == "UnbinLeadTau" ? string("Leading, ") + prong + "-prong" : string("Sub-leading, ") + prong + "-prong";

    h_out->SetBinContent(iter, h->GetBinContent(1));
    h_out->SetBinError(iter, h->GetBinError(1));
    h_out->GetXaxis()->SetBinLabel(iter, xlabel.c_str());

    iter++;

    __ENDLOOP
    __ENDLOOP

    h_out->GetXaxis()->SetTitleOffset(1.4);
    h_out->GetXaxis()->SetTitleSize(0.055);
    h_out->GetXaxis()->SetLabelSize(0.045);
    h_out->GetYaxis()->SetTitle("Transfer Factor");
    h_out->GetYaxis()->SetTitleOffset(1.4);
    h_out->GetYaxis()->SetTitleSize(0.055);
    h_out->GetYaxis()->SetLabelSize(0.045);
    h_out->SetMinimum(0.2);
    h_out->SetMaximum(1.8);

    TLine* line = new TLine(0, 1, 4, 1);
    line->SetLineColor(kRed+1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);

    auto color = TColor::GetColor(10, 73, 94);

    h_out->SetLineColor(color);
    h_out->SetLineWidth(2);
    h_out->SetMarkerColor(color);
    h_out->SetMarkerStyle(kFullSquare);
    h_out->SetMarkerSize(1.3);
    
    h_out->Draw("E0 SAME");
    line->Draw("SAME");
    h_out->Draw("E0 SAME");

    TLegend* legend = new TLegend(0.66, 0.86, 0.92, 0.92);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(16);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.040);
    string str_barrel{"FF(STT)"};
    legend->AddEntry(h_out, str_barrel.c_str(), "lp");
    legend->Draw();

    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.040);
    text->DrawLatex(0.20, 0.86, "1-b-tagged SS");
    text->DrawLatex(0.20, 0.96, string("Data in " + p.first + ", " + lumis[p.first]).c_str());

    stringstream output;

    output << "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/FF/TF_"
           << p.first << format;

    c->SaveAs(output.str().c_str());

    delete h_out;
    delete legend;
    delete text;
    delete line;
    delete c;
}

