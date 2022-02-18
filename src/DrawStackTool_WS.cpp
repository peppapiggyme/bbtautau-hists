

#include "DrawStackTool_WS.h"
#include "CommonInclude.h"

#include <sstream>
#include <algorithm>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

DrawStackTool_WS::DrawStackTool_WS(const DrawStackInfo* info)
    : DrawStackTool(info)
{


}

DrawStackTool_WS::~DrawStackTool_WS()
{

}

void DrawStackTool_WS::run(const Config* c) const
{
    if (!fs::exists(fs::path(output_path)))
    {
        throw std::runtime_error("Output path does not exist");
    }

    const unique_ptr<TFile>& fin = c->getInputTFile();
    const std::string& dirStr = c->getInputDirStr();

    RooHist* cRooHistData = (RooHist*)fin->Get((dirStr + "/RooData").c_str());
    for (int i = 0; i < cRooHistData->GetN(); i++)
    {
        cRooHistData->SetPointEXhigh(i, 0.);
        cRooHistData->SetPointEXlow(i, 0.);
        if (cRooHistData->GetPointY(i) < 1e-5) {
            cRooHistData->SetPointEYhigh(i, 0.);
            cRooHistData->SetPointEYlow(i, 0.);
        }
    }

    TH1F* cHistPostErr = (TH1F*)fin->Get((dirStr + "/error").c_str());

    TH1F* cHistPrefitBkg = (TH1F*)fin->Get((dirStr + "/PrefitBkg").c_str());
    TH1F* cHistChi2 = (TH1F*)fin->Get((dirStr + "/chi2").c_str());
    Tools::println("Chi2 = [%]", cHistChi2->GetBinContent(1));

    vector<ProcessInfo*>* ps = c->processes->content();

    gROOT->SetStyle("ATLAS");
    gStyle->SetErrorX(0.5);

    TCanvas* c1 = new TCanvas("c", "", m_info->draw_ratio ? 900 : 1200, 900);

    TPad* upper_pad = new TPad("upper_pad", "", 0, m_info->draw_ratio ? 0.35 : 0, 1, 1);
    upper_pad->SetBottomMargin(m_info->draw_ratio ? 0.03 : 0.15);
    upper_pad->SetTickx(false);
    upper_pad->SetTicky(false);
    upper_pad->SetLogx(m_info->logx);
    upper_pad->SetLogy(m_info->logy);
    upper_pad->Draw(); 

    TPad* lower_pad = new TPad("lower_pad", "", 0, 0, 1, m_info->draw_ratio ? 0.35 : 0);
    if (m_info->draw_ratio)
    {
        lower_pad->SetTopMargin(0);
        lower_pad->SetBottomMargin(0.4);
        lower_pad->SetGridy();
        lower_pad->SetLogx(m_info->logx);
        lower_pad->Draw();
    }

    upper_pad-> SetTicks();
    lower_pad-> SetTicks();

    upper_pad->cd();

    
    TH1* data = (*m_it_data)->histogram;

    THStack* stack = new THStack();
    map<string, THStack*> stacks;

    // magically sort the backgrounds by 
    // eProcess codes
    sort(m_it_bkg, m_it_sig, [](const ProcessInfo* left, const ProcessInfo* right) {
        return static_cast<unsigned>(left->process) > static_cast<unsigned>(right->process);
    });

    for_each(m_it_bkg, m_it_sig, [this, &stack, &stacks](const ProcessInfo* p) {
        p->histogram->Scale(p->norm_factor);
        if (m_info->draw_overflow) p->histogram->GetXaxis()->SetRange(1, p->histogram->GetNbinsX() + 1);
        stack->Add(p->histogram);
    });

    TH1* bkg = (TH1*)stack->GetStack()->Last()->Clone();
    Tools::println("bkg nominal = %", bkg->Integral());

    stack->Draw("HIST");
    if (m_info->draw_ratio)
    {
        stack->GetXaxis()->SetLabelSize(0);
        stack->GetXaxis()->SetTitleSize(0);
        stack->GetXaxis()->SetTitleOffset(1.3);
    }
    else
    {
        stack->GetXaxis()->SetTitle((*m_it_data)->current_variable->name_tex.c_str());
        stack->GetXaxis()->SetTitleOffset(1.3);
        stack->GetXaxis()->SetTitleSize(0.055);
        stack->GetXaxis()->SetLabelSize(0.045);
    }
    stack->GetYaxis()->SetTitle("Events");
    stack->GetYaxis()->SetLabelSize(0.045);
    stack->GetYaxis()->SetTitleSize(0.055);
    stack->SetMaximum(data->GetMaximum() * m_info->ymax_ratio_nolog);
    stack->SetMinimum(0);
    if (m_info->logy)
    {
        stack->SetMaximum(data->GetMaximum() * 400);
        stack->SetMinimum(std::max(data->GetMinimum() * 1e-1, 1.));
    }
    stack->GetYaxis()->ChangeLabel(1, -1, 0);
    if (m_info->draw_overflow) stack->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);
    
    // user defined maximum x must be defined without turning on the draw_overflow option!
    if (!m_info->draw_overflow && m_info->xmax < DBL_MAX) 
        stack->GetXaxis()->SetRangeUser(stack->GetXaxis()->GetXmin(), m_info->xmax);

    cHistPostErr->SetFillStyle(3254);
    cHistPostErr->SetFillColor(kGray+2);
    cHistPostErr->SetLineWidth(0);
    cHistPostErr->SetMarkerSize(0);
    cHistPostErr->SetName("Uncertainty");
    if (m_info->draw_overflow) 
        cHistPostErr->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);
    cHistPostErr->Draw("E2 SAME");

    cHistPrefitBkg->SetLineColor(kBlue);
    cHistPrefitBkg->SetLineStyle(7);
    cHistPrefitBkg->SetLineWidth(2);
    if (!m_info->prefit) 
        cHistPrefitBkg->Draw("HIST SAME");

    if (!m_info->blind) {
        if (m_info->draw_overflow) {
            data->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);
            cRooHistData->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);
        }
        cRooHistData->SetMarkerStyle(20);
        cRooHistData->SetMarkerSize(1.2);
        cRooHistData->SetMarkerColor(kBlack);
        cRooHistData->SetLineColor(kBlack);
        cRooHistData->Draw("P SAME");
    }

    for_each(m_it_sig, m_it_end, [this](const ProcessInfo* p) {
        p->histogram->Scale(m_info->signal_scale);
        p->histogram->Scale(p->norm_factor);
        p->histogram->SetLineStyle(m_info->signal_linestyle);
        p->histogram->SetLineWidth(m_info->signal_linewidth);
        if (m_info->draw_overflow) p->histogram->GetXaxis()->SetRange(1, p->histogram->GetNbinsX() + 1);
        p->histogram->Draw("HIST SAME"); 
    });

    double x = 0.92 - m_info->legend_scaling_horizontal * 0.15 * m_info->legend_ncolumns;
    double y = 0.92 - m_info->legend_scaling_vertical * (0.05 * (ps->size() / 2 + 1)) * 2 / m_info->legend_ncolumns;
    TLegend* legend = new TLegend(x, y, 0.92, 0.92);
    legend->SetNColumns(m_info->legend_ncolumns);
    legend->SetTextFont(42);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.045);
    legend->SetTextAlign(32);

    // reverse for legend
    sort(m_it_bkg, m_it_sig, [](const ProcessInfo* left, const ProcessInfo* right) {
        return static_cast<unsigned>(left->process) < static_cast<unsigned>(right->process);
    });

    legend->AddEntry(cRooHistData, "Data", "ep");
    for_each(m_it_bkg, m_it_sig, [&legend](const ProcessInfo* p) {
        legend->AddEntry(p->histogram, p->name_tex.c_str(), "f"); });
    if (!m_info->prefit) 
        legend->AddEntry(cHistPrefitBkg, "Pre-fit Bkg.", "l");
    for_each(m_it_sig, m_it_end, [&legend, this](const ProcessInfo* p) {
        string signal_name = m_info->show_scaling ? to_string((double)(m_info->signal_scale * p->norm_factor)).substr(0, 4) + " x " : string();
        signal_name += p->name_tex;
        legend->AddEntry(p->histogram, signal_name.c_str(), "l"); });
    legend->AddEntry(cHistPostErr, "Uncertainty", "f");

    legend->Draw("SAME");

    TLatex *text = new TLatex();
    text->SetNDC();
    float baseline = m_info->atlas ? 0.86 : 0.92;
    if (m_info->atlas) {
        text->SetTextFont(72);
        text->SetTextSize(0.055);
        text->DrawLatex(0.20, 0.86, "ATLAS");
        text->SetTextFont(42);
        text->DrawLatex(0.20 + 0.12, 0.86, m_info->atlas_label);
        text->SetTextSize(0.045);
    }
    ostringstream oss;
    oss << c->basic->ecm << ", " << c->basic->luminosity;
    text->DrawLatex(0.20, baseline - 0.06, oss.str().c_str());
    text->SetTextSize(0.045);
    text->DrawLatex(0.20, baseline - 0.12, (*m_it_data)->current_region->name_tex.c_str());

    if (m_info->draw_ratio) {
        lower_pad->cd();
        double resize = 0.65 / 0.35;

        TH1* err = (TH1*)cHistPostErr->Clone();
        for (int i = 0; i < err->GetNbinsX() + 2; ++i)
        {
            if (bkg->GetBinContent(i) < 1e-5) {
                err->SetBinError(i, 0.0);
                err->SetBinContent(i, 1.0);
            }
            else
            {
                err->SetBinError(i, err->GetBinError(i) / err->GetBinContent(i));
                err->SetBinContent(i, 1.0);
            }
        }
        err->SetFillStyle(3254);
        err->SetFillColor(kGray+2);
        err->SetLineWidth(0);    
        err->SetMarkerSize(0);
        err->SetName("Uncertainty");
        err->GetXaxis()->SetTitle((*m_it_data)->current_variable->name_tex.c_str());
        err->GetXaxis()->SetTitleOffset(0.8 * resize);
        err->GetXaxis()->SetTitleSize(0.055 * resize);
        err->GetXaxis()->SetLabelSize(0.045 * resize);
        err->GetYaxis()->SetTitle("Data / Pred.");
        err->GetYaxis()->SetTitleOffset(0.4 * resize);
        err->GetYaxis()->SetTitleSize(0.055 * resize);
        err->GetYaxis()->SetLabelSize(0.045 * resize);
        err->GetYaxis()->SetNdivisions(505);
        if (m_info->draw_overflow) err->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);

        double rat_max = 0;
        double rat_min = 2;
        RooHist* rat = (RooHist*)cRooHistData->Clone();
        for (int i = 0; i < rat->GetN(); i++)
        {
            rat->SetPointY(i, rat->GetPointY(i) / bkg->GetBinContent(i+1));
            rat->SetPointEYhigh(i, rat->GetErrorYhigh(i) / bkg->GetBinContent(i+1));
            rat->SetPointEYlow(i, rat->GetErrorYlow(i) / bkg->GetBinContent(i+1));
            rat_max = std::max(rat_max, rat->GetPointY(i) + 1.5 * rat->GetErrorYhigh(i));
            rat_min = std::min(rat_min, rat->GetPointY(i) - 1.5 * rat->GetErrorYlow(i));
        }
        rat->SetTitle("ratio");

        if (!m_info->auto_ratio)
        {
            err->SetMinimum(m_info->ratio_low);
            err->SetMaximum(m_info->ratio_high);
        }
        else
        {
            err->SetMinimum(rat_min);
            err->SetMaximum(rat_max);
        }

        // draw
        err->Draw("E2");
        if (!m_info->blind) {
            if (m_info->draw_overflow) rat->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);
            rat->Draw("P SAME");
        }
    }

    ostringstream oss_out;
    oss_out << output_path << "/" 
            << (*m_it_data)->current_region->name << "_"
            << (*m_it_data)->current_variable->name << "_"
            << m_info->parameter << "." 
            << m_info->output_format;
    c1->Update();
    c1->SaveAs(oss_out.str().c_str());

    delete upper_pad;
    delete lower_pad;
    delete stack;
    delete legend;
    delete text;
    delete c1;
}
