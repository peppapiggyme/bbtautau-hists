#include "DistributionTool.h"
#include "CommonInclude.h"

#include <sstream>
#include <algorithm>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

DistributionTool::DistributionTool(const DistributionInfo* info)
    : CompTool(info), m_info(info)
{

}

DistributionTool::~DistributionTool()
{

}

void DistributionTool::run(const Config* c) const
{
    if (!fs::exists(fs::path(output_path)))
    {
        throw std::runtime_error("Output path does not exist");
    }
    
    vector<ProcessInfo*>* ps = c->processes->content();

    sort(ps->begin(), ps->end(), [](const ProcessInfo* left, const ProcessInfo* right) {
        return static_cast<unsigned>(left->process) < static_cast<unsigned>(right->process);
    });

    // pre-smoothing
    for_each(ps->begin(), ps->end(), [](const ProcessInfo* p) {
        for (int i = 0; i <= p->histogram->GetNbinsX()+1; i++) 
        {
            if (p->histogram->GetBinContent(i) < 0) {
                p->histogram->SetBinContent(i, 0);
            }
        }
    });

    gROOT->SetStyle("ATLAS");
    gStyle->SetErrorX(0.5);

    TCanvas* c1 = new TCanvas("c", "", 1200, 900);
    c1->SetRightMargin(1.6 * c1->GetRightMargin());
    c1->SetBottomMargin(0.14);
    c1->SetTickx(false);
    c1->SetTicky(false);
    c1->SetLogx(m_info->logx);
    c1->SetLogy(m_info->logy);
    c1->Draw();

    TH1* base = ps->front()->histogram;
    
    for_each(ps->begin(), ps->end(), [this, &base](const ProcessInfo* p) {
        if (m_info->shape_only) {
            p->histogram->Scale(1.0 / p->histogram->Integral());
        } else {
            p->histogram->Scale(p->norm_factor);
        }
        if (p->type == eProcessType::BKG)
        {
            p->histogram->SetLineWidth(1);
            p->histogram->SetFillColorAlpha(p->histogram->GetLineColor(), 0.2);
        }
    });
    
    base->Draw("HIST");
    base->GetXaxis()->SetLabelSize(0.04);
    base->GetXaxis()->SetTitleSize(0.045);
    base->GetXaxis()->SetTitleOffset(1.2);
    base->GetXaxis()->SetTitle(ps->front()->current_variable->name_tex.c_str());
    if (!m_info->shape_only)
    {
        base->GetYaxis()->SetTitle("Events");
    }
    else
    {
        base->GetYaxis()->SetTitle("Arbitrary Unit");
    }
    base->GetYaxis()->SetLabelSize(0.04);
    base->GetYaxis()->SetTitleSize(0.045);
    float maxHeight = -1.0;
    for_each(ps->begin(), ps->end(), [&maxHeight](const ProcessInfo* p) {
        maxHeight = max(maxHeight, (float)p->histogram->GetMaximum());
    });
    if (!m_info->logy)
    {
        base->SetMaximum(maxHeight * 1.3);
        base->SetMinimum(0);
    }
    else 
    {
        base->SetMaximum(maxHeight * 5);
        if (m_info->shape_only)
        {
            base->SetMinimum(8e-5);
        }
    }

    base->GetYaxis()->ChangeLabel(1, -1, 0);

    for_each(ps->begin()+1, ps->end(), [this, &base](const ProcessInfo* p) {
        p->histogram->Draw("HIST SAME"); 
    });

    // fine tune for mvainput
    for_each(ps->begin(), ps->end(), [this, &base](const ProcessInfo* p) {
        if (p->type == eProcessType::SIG) p->histogram->Draw("HIST SAME"); 
        if (p->process == eProcess::P2) p->histogram->SetLineStyle(2);
        if (p->process == eProcess::P3) p->histogram->SetLineStyle(3);
    });

    double y = 0.92 - 0.05 * (ps->size() + ps->front()->systematic_histograms.size() + 1);
    TLegend* legend = new TLegend(0.70, y, 0.90, 0.92);
    legend->SetTextFont(42);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetTextAlign(12);

    for_each(ps->begin(), ps->end(), [&legend, &ps](const ProcessInfo* p) {
        legend->AddEntry(p->histogram, p->name_tex.c_str());
    });

    legend->Draw("SAME");

    TLatex *text = new TLatex();
    text->SetNDC();
    if (m_info->atlas) {
        text->SetTextFont(72);
        text->SetTextSize(0.050);
        text->DrawLatex(0.20, 0.96, "ATLAS");
        text->SetTextFont(42);
        text->DrawLatex(0.20 + 0.12, 0.96, m_info->atlas_label);
        text->SetTextSize(0.045);
    }
    ostringstream oss{c->basic->ecm};
    text->DrawLatex(0.20, 0.86, oss.str().c_str());
    text->SetTextSize(0.040);
    
    ostringstream oss_out;
    oss_out << output_path << "/" 
            << c->current_region->name << "_"
            << c->current_variable->name << "_"
            << m_info->parameter << ".pdf";
    c1->Update();
    c1->SaveAs(oss_out.str().c_str());

    delete legend;
    delete text;
    delete c1;
}

void DistributionTool::run_overfit(const Config* c_train, const Config* c_test) const
{
    if (!fs::exists(fs::path(output_path)))
    {
        throw std::runtime_error("Output path does not exist");
    }
    
    vector<ProcessInfo*>* ps = c_train->processes->content(); // train
    vector<ProcessInfo*>* ps_test = c_test->processes->content(); // test

    sort(ps->begin(), ps->end(), [](const ProcessInfo* left, const ProcessInfo* right) {
        return static_cast<unsigned>(left->process) < static_cast<unsigned>(right->process);
    });

    sort(ps_test->begin(), ps_test->end(), [](const ProcessInfo* left, const ProcessInfo* right) {
        return static_cast<unsigned>(left->process) < static_cast<unsigned>(right->process);
    });

    // pre-smoothing
    for_each(ps->begin(), ps->end(), [](const ProcessInfo* p) {
        for (int i = 0; i <= p->histogram->GetNbinsX()+1; i++) 
        {
            if (p->histogram->GetBinContent(i) < 0) {
                p->histogram->SetBinContent(i, 0);
            }
        }
    });

    for_each(ps_test->begin(), ps_test->end(), [](const ProcessInfo* p) {
        for (int i = 0; i <= p->histogram->GetNbinsX()+1; i++) 
        {
            if (p->histogram->GetBinContent(i) < 0) {
                p->histogram->SetBinContent(i, 0);
            }
        }
    });

    gROOT->SetStyle("ATLAS");
    gStyle->SetErrorX(0.5);

    TCanvas* c1 = new TCanvas("c", "", 1200, 900);
    c1->SetRightMargin(1.6 * c1->GetRightMargin());
    c1->SetBottomMargin(0.14);
    c1->SetTickx(false);
    c1->SetTicky(false);
    c1->SetLogx(m_info->logx);
    c1->SetLogy(m_info->logy);
    c1->Draw();

    TH1* base = ps->front()->histogram;
    
    for_each(ps->begin(), ps->end(), [this, &base](const ProcessInfo* p) {
        if (m_info->shape_only) {
            p->histogram->Scale(1.0 / p->histogram->Integral());
        } else {
            p->histogram->Scale(p->norm_factor);
        }
        p->histogram->SetLineStyle(2);
    });

    for_each(ps_test->begin(), ps_test->end(), [this, &base](const ProcessInfo* p) {
        if (m_info->shape_only) {
            p->histogram->Scale(1.0 / p->histogram->Integral());
        } else {
            p->histogram->Scale(p->norm_factor);
        }
    });

    base->Draw("HIST");
    base->GetXaxis()->SetLabelSize(0.04);
    base->GetXaxis()->SetTitleSize(0.045);
    base->GetXaxis()->SetTitleOffset(1.2);
    base->GetXaxis()->SetTitle(ps->front()->current_variable->name_tex.c_str());
    if (!m_info->shape_only)
    {
        base->GetYaxis()->SetTitle("Events");
    }
    else
    {
        base->GetYaxis()->SetTitle("Arbitrary Unit");
    }
    base->GetYaxis()->SetLabelSize(0.04);
    base->GetYaxis()->SetTitleSize(0.045);

    float maxHeight = -1.0;
    for_each(ps->begin(), ps->end(), [&maxHeight](const ProcessInfo* p) {
        maxHeight = max(maxHeight, (float)p->histogram->GetMaximum());
    });
    for_each(ps_test->begin(), ps_test->end(), [&maxHeight](const ProcessInfo* p) {
        maxHeight = max(maxHeight, (float)p->histogram->GetMaximum());
    });

    if (!m_info->logy)
    {
        base->SetMaximum(maxHeight * 1.3);
        base->SetMinimum(0);
    }
    else 
    {
        base->SetMaximum(maxHeight * 5);
        if (m_info->shape_only)
        {
            base->SetMinimum(8e-5);
        }
    }

    base->GetYaxis()->ChangeLabel(1, -1, 0);

    for_each(ps_test->begin(), ps_test->end(), [this, &base](const ProcessInfo* p) {
        p->histogram->Draw("HIST SAME"); 
    });

    for_each(ps->begin(), ps->end(), [this, &base](const ProcessInfo* p) {
        p->histogram->Draw("HIST SAME"); 
    });

    // K-S tests
    // make sure p[0] => signal and p[1] => total bkg! see examples/thesis/Example_OverFit.cpp

    double psig = ps->at(0)->histogram->KolmogorovTest(ps_test->at(0)->histogram, "X");
    double pbkg = ps->at(1)->histogram->KolmogorovTest(ps_test->at(1)->histogram, "X");

    char psig_s[64];
    char pbkg_s[64];

    sprintf(psig_s, "p_{sig} = %.2f", psig);
    sprintf(pbkg_s, "p_{bkg} = %.2f", pbkg);

    double y = 0.92 - 0.05 * (ps->size()+2);
    TLegend* legend = new TLegend(0.40, y, 0.90, 0.92);
    legend->SetNColumns(2);
    legend->SetTextFont(42);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.036);
    legend->SetTextAlign(12);

    for_each(ps->begin(), ps->end(), [&legend, &ps](const ProcessInfo* p) {
        legend->AddEntry(p->histogram, (p->name_tex + string(" [Train]")).c_str());
    });

    for_each(ps_test->begin(), ps_test->end(), [&legend, &ps](const ProcessInfo* p) {
        legend->AddEntry(p->histogram, (p->name_tex + string(" [Test]")).c_str());
    });

    legend->AddEntry("", psig_s, "");
    legend->AddEntry("", pbkg_s, "");

    legend->Draw("SAME");

    TLatex *text = new TLatex();
    text->SetNDC();
    if (m_info->atlas) {
        text->SetTextFont(72);
        text->SetTextSize(0.050);
        text->DrawLatex(0.20, 0.96, "ATLAS");
        text->SetTextFont(42);
        text->DrawLatex(0.20 + 0.12, 0.96, m_info->atlas_label);
        text->SetTextSize(0.045);
    }
    ostringstream oss{c_train->basic->ecm};
    text->DrawLatex(0.20, 0.86, oss.str().c_str());
    text->SetTextSize(0.040);
    
    ostringstream oss_out;
    oss_out << output_path << "/" 
            << c_train->current_region->name << "_"
            << c_train->current_variable->name << "_"
            << m_info->parameter << ".pdf";
    c1->Update();
    c1->SaveAs(oss_out.str().c_str());

    delete legend;
    delete text;
    delete c1;
}
