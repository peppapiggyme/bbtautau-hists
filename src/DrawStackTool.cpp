#include "DrawStackTool.h"
#include "CommonInclude.h"

#include <sstream>
#include <algorithm>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

DrawStackTool::DrawStackTool(const DrawStackInfo* info)
    : HistTool(), m_info(info)
{

}

DrawStackTool::~DrawStackTool()
{

}

bool DrawStackTool::check(const Config* c) const
{
    if (!HistTool::check(c))
        return false;

    vector<ProcessInfo*>* ps = c->processes->content();
    int n_data_check = count_if(ps->begin(), ps->end(), [](const ProcessInfo* p) {
        return p->type == eProcessType::DATA; });

    if (n_data_check != 1) {
        cerr << "FAIL: must be one data!\n";
        return false;
    }

    int n_bkg_check = count_if(ps->begin(), ps->end(), [](const ProcessInfo* p) {
        return p->type == eProcessType::BKG; });

    if (n_bkg_check < 1) {
        cerr << "FAIL: no background!\n";
        return false;
    }

    return true;

}

void DrawStackTool::paint(const Config* c) const
{
    vector<ProcessInfo*>* ps = c->processes->content();
    for_each(ps->begin(), ps->end(), [&c](ProcessInfo* p) {
        switch (p->type)
        {
        case eProcessType::DATA:
            p->histogram->SetMarkerStyle(20);
            p->histogram->SetMarkerSize(1.2);
            p->histogram->SetMarkerColor(p->color);
            p->histogram->SetLineColor(p->color);
            break;
        case eProcessType::BKG:
            p->histogram->SetLineColor(kBlack);
            p->histogram->SetLineWidth(1);
            p->histogram->SetFillColor(p->color);
            break;
        case eProcessType::SIG:
            p->histogram->SetLineWidth(2);
            p->histogram->SetLineStyle(1);
            p->histogram->SetLineColor(p->color);
            break;
        default:
            throw runtime_error("this never happen");
            break;
        }
    });
}

void DrawStackTool::manipulate(Config* c)
{
    HistTool::manipulate(c); // always do this first

    vector<ProcessInfo*>* ps = c->processes->content();
    m_it_data = ps->begin();
    m_it_end = ps->end();

    m_it_bkg = partition(ps->begin(), ps->end(), [](const ProcessInfo* p) {
        return p->type == eProcessType::DATA; });

    m_it_sig = partition(m_it_bkg, ps->end(), [](const ProcessInfo* p) {
        return p->type == eProcessType::BKG; });

    sort(m_it_bkg, m_it_sig, [this](const ProcessInfo* p1, const ProcessInfo* p2) {
        return *p1 < *p2; });

    sort(m_it_sig, ps->end(), [this](const ProcessInfo* p1, const ProcessInfo* p2) {
        return *p1 < *p2; });
}

void DrawStackTool::run(const Config* c) const
{
    if (!fs::exists(fs::path(output_path)))
    {
        throw std::runtime_error("Output path does not exist");
    }
    
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

    upper_pad->cd();

    TH1* data = (*m_it_data)->histogram;
    if (m_info->use_poisson_data_error) data->SetBinErrorOption(TH1::kPoisson);
    THStack* stack = new THStack();
    map<string, THStack*> stacks;
    set<string> systs;

    // magically sort the backgrounds by 
    // eProcess codes
    sort(m_it_bkg, m_it_sig, [](const ProcessInfo* left, const ProcessInfo* right) {
        return static_cast<unsigned>(left->process) > static_cast<unsigned>(right->process);
    });

    for_each(m_it_bkg, m_it_sig, [&systs, &stacks](const ProcessInfo* p) {
        for (auto &pp : p->systematic_histograms)
        {
            systs.insert(pp.first);
            stacks[pp.first] = new THStack();
        }
    });
    for_each(m_it_bkg, m_it_sig, [this, &stack, &stacks](const ProcessInfo* p) {
        p->histogram->Scale(p->norm_factor);
        if (m_info->draw_overflow) p->histogram->GetXaxis()->SetRange(1, p->histogram->GetNbinsX() + 1);
        stack->Add(p->histogram); 
        set<string> systs_here;
        for (auto &pp : p->systematic_histograms)
        {
            systs_here.insert(pp.first);
            pp.second->Scale(p->norm_factor);
            if (m_info->draw_overflow) pp.second->GetXaxis()->SetRange(1, pp.second->GetNbinsX() + 1);
        }
        for (auto &pp : stacks)
        {
            if (systs_here.find(pp.first) == systs_here.end())
            {
                stacks[pp.first]->Add(p->histogram);
            }
            else
            {
                stacks[pp.first]->Add(p->systematic_histograms.at(pp.first));
            }
        }
    });

    TH1* bkg = (TH1*)stack->GetStack()->Last()->Clone();
    Tools::println("bkg nominal = %", bkg->Integral());
    map<string, TH1*> bkgs;
    for (auto &pp : stacks)
    {
        bkgs[pp.first] = (TH1*)pp.second->GetStack()->Last()->Clone();
        Tools::println("bkg % = %", pp.first, bkgs[pp.first]->Integral());
    }
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
    if (!m_info->draw_overflow && m_info->xmax > DBL_MAX) 
        stack->GetXaxis()->SetRangeUser(stack->GetXaxis()->GetXmin(), m_info->xmax);

    for (auto &pp : stacks)
    {
        if (m_info->draw_overflow && pp.second->GetXaxis()) pp.second->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);
    }

    for (auto& pp : bkgs)
    {
        pp.second->Add(bkg, -1);
    }

    vector<double> total_errors(bkg->GetNbinsX()+2, 0.);
    for (int i = 0; i < bkg->GetNbinsX()+2; i++)
    {
        double total_error_2 = bkg->GetBinError(i) * bkg->GetBinError(i);
        for (auto& pp : bkgs)
        {
            total_error_2 += pp.second->GetBinContent(i) * pp.second->GetBinContent(i);
        }
        total_errors[i] = TMath::Sqrt(total_error_2 / 2.0); // -> to average up and down
    }

    TH1* bkg_stat = (TH1*)bkg->Clone();
    TH1* bkg_copy = (TH1*)bkg->Clone();
    for (int i = 0; i < bkg->GetNbinsX()+2; i++)
    {
        bkg_copy->SetBinError(i, total_errors[i]);
    }
    
    bkg_copy->SetFillStyle(3254);
    bkg_copy->SetFillColor(kGray+2);
    bkg_copy->SetLineWidth(0);
    // bkg_copy->SetFillColorAlpha(kRed + 3, 0.2);
    bkg_copy->SetMarkerSize(0);
    bkg_copy->SetName("Uncertainty");
    if (m_info->draw_overflow) bkg_copy->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);
    bkg_copy->Draw("E2 SAME");

    if (!m_info->blind) {
        if (m_info->draw_overflow) data->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);
        data->Draw("E0 X0 SAME");
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

    legend->AddEntry(data, "Data", "ep");
    for_each(m_it_bkg, m_it_sig, [&legend](const ProcessInfo* p) {
        legend->AddEntry(p->histogram, p->name_tex.c_str(), "f"); });
    for_each(m_it_sig, m_it_end, [&legend, this](const ProcessInfo* p) {
        string signal_name = m_info->show_scaling ? to_string((double)(m_info->signal_scale * p->norm_factor)).substr(0, 4) + " x " : string();
        signal_name += p->name_tex;
        legend->AddEntry(p->histogram, signal_name.c_str(), "l"); });
    legend->AddEntry(bkg_copy, "Uncertainty", "f");

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

        TH1* err = (TH1*)bkg_copy->Clone();
        TH1* bkg_scale = (TH1*)bkg->Clone();
        for (int i = 0; i < bkg_scale->GetNbinsX() + 2; ++i)
        {
            bkg_scale->SetBinError(i, 0.0);
        }
        err->Divide(bkg_scale);
        err->SetFillStyle(3254);
        err->SetFillColor(kGray+2);
        err->SetLineWidth(0);    
        // err->SetFillColorAlpha(kRed + 3, 0.2);
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
        err->SetMinimum(m_info->ratio_low);
        err->SetMaximum(m_info->ratio_high);
        if (m_info->draw_overflow) err->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);
        err->Draw("E2");

        TH1* err_stat = (TH1*)bkg_stat->Clone();
        err_stat->Divide(bkg_scale);
        err_stat->SetLineWidth(0);
        // err_stat->SetFillStyle(3254);
        err_stat->SetFillColorAlpha(kGray + 3, 0.2);
        err_stat->SetMarkerSize(0);
        err_stat->SetName("Unc. Stat-Only");
        // err_stat->Draw("E2 SAME");

        /// @todo: other tool might also need this!
        // {
        //     TLegend* legend = new TLegend(0.62, 0.88, 0.90, 0.98);
        //     legend->SetNColumns(2);
        //     legend->SetTextFont(42);
        //     legend->SetFillStyle(0);
        //     legend->SetBorderSize(0);
        //     legend->SetTextSize(0.035 * resize);
        //     legend->SetTextAlign(12);
        //     // legend->AddEntry(err_stat, "Stat Unc.", "f");
        //     legend->AddEntry(err, "Uncertainty", "f");
        //     legend->Draw("SAME");
        // }

        TH1* rat = (TH1*)data->Clone();
        // if (m_info->use_poisson_data_error) bkg_scale->SetBinErrorOption(TH1::kPoisson);
        rat->Divide(bkg_scale);
        if (m_info->use_poisson_data_error) rat->SetBinErrorOption(TH1::kPoisson);
        rat->SetTitle("lower_pad");
        if (!m_info->blind) {
            if (m_info->draw_overflow) rat->GetXaxis()->SetRange(1, data->GetNbinsX() + 1);
            rat->Draw("E0 X0 SAME");
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

    // TCanvas* c2 = new TCanvas("c2", "", 900, 900);
    // TH1* bkg_clone = (TH1*)bkg->Clone();
    // TH1* datadriven = (TH1*)data->Clone();
    // datadriven->Add(bkg_clone, -1);
    // datadriven->Scale((*m_it_sig)->histogram->Integral() / datadriven->Integral());
    // datadriven->Divide((*m_it_sig)->histogram);
    // datadriven->SetName((*m_it_data)->current_variable->name.c_str());
    // datadriven->SetTitle((*m_it_data)->current_variable->name.c_str());
    // datadriven->GetXaxis()->SetTitle((*m_it_data)->current_variable->name_tex.c_str());
    // datadriven->GetXaxis()->SetTitleOffset(1.2);
    // datadriven->GetXaxis()->SetTitleSize(0.045);
    // datadriven->GetXaxis()->SetLabelSize(0.04);
    // datadriven->GetYaxis()->SetTitle("Scale Factor");
    // datadriven->GetYaxis()->SetTitleOffset(1.2);
    // datadriven->GetYaxis()->SetTitleSize(0.045);
    // datadriven->GetYaxis()->SetLabelSize(0.04);
    // datadriven->GetYaxis()->SetNdivisions(505);
    // datadriven->Draw("E1");
    // ostringstream oss_out2;
    // oss_out2 << output_path << "/" << "sf_"
    //         << (*m_it_data)->current_region->name << "_"
    //         << (*m_it_data)->current_variable->name << "_"
    //         << m_info->parameter << ".png";
    // c2->Update();
    // c2->SaveAs(oss_out2.str().c_str());

    // oss_out2 << ".root";
    // TFile *fout = TFile::Open(oss_out2.str().c_str(), "recreate");
    // datadriven->Write();
    // fout->Close();

    delete upper_pad;
    delete lower_pad;
    delete stack;
    for (auto& pp : stacks)
    {
        delete pp.second;
    }
    delete legend;
    delete text;
    delete c1;
    // delete c2;
}
