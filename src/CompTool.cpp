#include "CompTool.h"
#include "CommonInclude.h"

#include <sstream>
#include <algorithm>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

CompTool::CompTool(const CompInfo* info)
    : HistTool(), m_info(info)
{

}

CompTool::~CompTool()
{

}

bool CompTool::check(const Config* c) const
{
    if (!HistTool::check(c))
        return false;

    vector<ProcessInfo*>* ps = c->processes->content();

    return ps->size() > 1;
}

void CompTool::paint(const Config* c) const
{
    vector<ProcessInfo*>* ps = c->processes->content();
    for_each(ps->begin(), ps->end(), [&c](ProcessInfo* p) {
        p->histogram->SetLineWidth(2);
        p->histogram->SetLineStyle(1);
        p->histogram->SetMarkerSize(0);
        p->histogram->SetMarkerColor(p->color);
        p->histogram->SetLineColor(p->color);

        for (auto& pp : p->systematic_histograms)
        {
            pp.second->SetLineWidth(2);
            // pp.second->SetLineStyle(2);
            pp.second->SetMarkerSize(0);
        }
    });
}

void CompTool::run(const Config* c) const
{
    if (!fs::exists(fs::path(output_path)))
    {
        throw std::runtime_error("Output path does not exist");
    }
    
    vector<ProcessInfo*>* ps = c->processes->content();
    vector<SystematicInfo*>* ss = c->systematics->content();

    gROOT->SetStyle("ATLAS");
    gStyle->SetErrorX(0.5);

    double resize = 1.5;

    TCanvas* c1 = new TCanvas("c", "", 900, 900);
    c1->SetRightMargin(1.6 * c1->GetRightMargin());

    TPad* upper_pad = new TPad("upper_pad", "", 0, 0.5, 1, 1);
    upper_pad->SetRightMargin(1.6 * upper_pad->GetRightMargin());
    upper_pad->SetBottomMargin(0);
    upper_pad->SetTickx(false);
    upper_pad->SetTicky(false);
    upper_pad->SetLogx(m_info->logx);
    upper_pad->SetLogy(m_info->logy);
    upper_pad->Draw();

    TPad* lower_pad = new TPad("lower_pad", "", 0, 0, 1, 0.5);
    lower_pad->SetRightMargin(1.6 * lower_pad->GetRightMargin());
    lower_pad->SetTopMargin(0);
    lower_pad->SetBottomMargin(0.26);
    lower_pad->SetGridy();
    lower_pad->SetLogx(m_info->logx);
    lower_pad->Draw();

    upper_pad->cd();

    if (m_info->shape_only) {
        ps->front()->histogram->Scale(1.0 / ps->front()->histogram->Integral());
        for (auto &pp : ps->front()->systematic_histograms)
        {
            pp.second->Scale(1.0 / pp.second->Integral());
        }
    } else {
        ps->front()->histogram->Scale(ps->front()->norm_factor);
        for (auto &pp : ps->front()->systematic_histograms)
        {
            pp.second->Scale(ps->front()->norm_factor);
        }
    }
    TH1* base = ps->front()->histogram;
    base->Draw("HIST");
    base->GetXaxis()->SetLabelSize(0);
    base->GetXaxis()->SetTitleSize(0);
    base->GetYaxis()->SetTitle(m_info->shape_only ? "Arbitrary Unit" : "Events");
    base->GetYaxis()->SetLabelSize(0.045 * resize);
    base->GetYaxis()->SetTitleSize(0.055 * resize);
    base->GetYaxis()->SetTitleOffset(0.8);
    base->SetMaximum(base->GetMaximum() * m_info->ymax_ratio_nolog);
    base->SetMinimum(0.);
    if (m_info->logy)
    {
        base->SetMaximum(base->GetMaximum() * 10);
    }
    base->GetYaxis()->ChangeLabel(1, -1, 0);

    for (auto& pp : ps->front()->systematic_histograms)
    {
        // chi2 test
        TH1* against = (TH1*)(*(ps->begin()+1))->histogram->Clone();
        TH1* here = (TH1*)(pp.second)->Clone();
        against->GetXaxis()->SetRangeUser(m_info->xmin_for_test, m_info->xmax_for_test);
        here->GetXaxis()->SetRangeUser(m_info->xmin_for_test, m_info->xmax_for_test);
        double chi2 = 0;
        int ndf = 0;
        int igood = 0;
        double pvalue = against->Chi2TestX(here, chi2, ndf, igood, "WW");
        double pvalue_ks = against->KolmogorovTest(here, "X");
        
        Tools::println("name [%]\n - chi2/ndf [%]\n - p_chi2 [%]\n - p_ks [%]\n - good [%]", 
            pp.first, chi2/(double)ndf, pvalue, pvalue_ks, igood);

        pp.second->Draw("HIST SAME");

        delete against;
        delete here;
    }

    for_each(ps->begin()+1, ps->end(), [this, &base](const ProcessInfo* p) {
        if (m_info->shape_only) {
            p->histogram->Scale(1.0 / p->histogram->Integral());
        } else {
            p->histogram->Scale(p->norm_factor);
        }
        p->histogram->Draw("E0 SAME"); 
    });

    TH1* base_copy = (TH1*)base->Clone();
    base_copy->SetFillStyle(3254);
    base_copy->SetFillColor(kGray+2);
    base_copy->SetLineWidth(0);
    base_copy->SetMarkerSize(0);
    base_copy->SetName("Uncertainty");
    base_copy->Draw("E2 SAME");

    double x = 0.92 - m_info->legend_scaling_horizontal * 0.32 ;
    double y = 0.92 - m_info->legend_scaling_vertical * 0.06 * (ps->size() + ps->front()->systematic_histograms.size() + 1);
    TLegend* legend = new TLegend(x, y, 0.92, 0.92);
    legend->SetTextFont(42);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.040 * resize);
    legend->SetTextAlign(12);

    vector<TH1*> dummys;
    for (auto& s : *ss) {
        TH1* dummy = new TH1F();
        dummy->SetLineColor(s->color);
        dummy->SetMarkerColor(s->color);
        dummy->SetLineWidth(2);
        dummy->SetMarkerSize(0);
        dummys.push_back(dummy);
    }

    for_each(ps->begin(), ps->end(), [&legend, &ps, &ss, &base_copy, &dummys](const ProcessInfo* p) {
        legend->AddEntry(p->histogram, p->name_tex.c_str(), p == ps->front() ? "f" : "lep");
        if (p == ps->front()) 
            legend->AddEntry(base_copy, (p->name_tex + " Stat. Unc.").c_str(), "f");
        if (p == ps->front())
        for (size_t i = 0; i < ss->size(); i++)
        {
            legend->AddEntry(dummys[i], ss->at(i)->name_tex.c_str(), "f");
        }
    });
    legend->Draw("SAME");

    TLatex *text = new TLatex();
    text->SetNDC();
    if (m_info->atlas) {
        text->SetTextFont(72);
        text->SetTextSize(0.055 * resize);
        text->DrawLatex(0.20, 0.86, "ATLAS");
        text->SetTextFont(42);
        text->DrawLatex(0.20 + 0.12, 0.86, m_info->atlas_label);
        text->SetTextSize(0.045 * resize);
    }
    text->SetTextSize(0.045 * resize);
    text->DrawLatex(0.20, 0.86, ps->front()->current_region->name_tex.c_str());

    lower_pad->cd();

    TH1* err = (TH1*)base->Clone();
    TH1* base_scale = (TH1*)base->Clone();
    for (int i = 0; i < base_scale->GetNbinsX() + 2; ++i)
    {
        base_scale->SetBinError(i, 0.0);
    }
    err->Divide(base_scale);
    err->SetLineWidth(0);
    err->SetFillStyle(3254);
    err->SetFillColor(kGray+2);
    err->SetMarkerSize(0);
    err->SetName("Unc.");
    err->GetXaxis()->SetTitle(ps->front()->current_variable->name_tex.c_str());
    err->GetXaxis()->SetTitleOffset(1.2);
    err->GetXaxis()->SetTitleSize(0.055 * resize);
    err->GetXaxis()->SetLabelSize(0.045 * resize);
    err->GetYaxis()->SetTitle(m_info->ratio_tex.c_str());
    err->GetYaxis()->SetTitleOffset(0.8);
    err->GetYaxis()->SetTitleSize(0.055 * resize);
    err->GetYaxis()->SetLabelSize(0.045 * resize);
    err->GetYaxis()->SetNdivisions(505);
    err->SetMinimum(m_info->ratio_low);
    err->SetMaximum(m_info->ratio_high);
    err->Draw("E2");

    for_each(ps->begin()+1, ps->end(), [this, &base_scale, &err](const ProcessInfo* p) {
        TH1* rat = (TH1*)p->histogram->Clone();
        rat->Divide(base_scale);
        // rat->Fit("pol1", "", "", 0, 250);
        rat->Draw("E0 SAME"); 
        if (m_info->save_ratio)
        {
            ostringstream oss;
            oss << output_path << "/CompareTo_" << rat->GetName() << ".root";
            TFile f(oss.str().c_str(), "recreate");
            f.cd();
            rat->Write(p->current_variable->name.c_str());
            f.Close();
        }
    });
    
    for (auto& pp : ps->front()->systematic_histograms)
    {
        TH1* rat_pp = (TH1*)pp.second->Clone();
        rat_pp->Divide(base_scale);
        for (int i = 0; i < rat_pp->GetNbinsX() + 2; ++i)
        {
            if (base_scale->GetBinContent(i) < 1e-6)
            {
                rat_pp->SetBinContent(i, 1.);
            }
        }
        if (rat_pp->GetLineColor() == kGreen || rat_pp->GetLineColor() == kGreen+2 || rat_pp->GetLineColor() == kGreen+4)
        {
            rat_pp->SetLineWidth(4);
        }
        rat_pp->Draw("HIST SAME");
    }

    ostringstream oss_out;
    oss_out << output_path << "/" 
            << c->current_region->name << "_"
            << c->current_variable->name << "_"
            << m_info->parameter << "."
            << m_info->output_format;

    c1->Update();
    c1->SaveAs(oss_out.str().c_str());

    for (auto& d : dummys) delete d;
    delete upper_pad;
    delete lower_pad;
    delete legend;
    delete text;
    delete c1;
}
