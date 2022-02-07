#include "Config.h"
#include "WSMakerBinning.h"

#include <iterator>
#include <iostream>
#include <exception>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

BasicInfo::BasicInfo(const string& ecm, const string& lumi) noexcept
    : ecm(ecm), luminosity(lumi)
{

}

Config::Config(const BasicInfo* b, const Processes* ps, 
               const Regions* rs, const Variables* vs, const Systematics* ss) noexcept
    : basic(b), processes(ps), regions(rs), variables(vs), systematics(ss)
{
    m_fin = nullptr;
    m_dir = "";
    m_loaded = false;
    m_manipulated = false;
}

Config::~Config() noexcept
{
    m_fin->Close();
}

void Config::load(const string& fn, const string& dir)
{
    if (!fs::exists(fs::path(fn)))
    {
        throw std::runtime_error("Input does not exist");
    }
    m_fin.reset(TFile::Open(fn.c_str(), "read"));
    if (!m_fin)
    {
        throw std::runtime_error("Input does not exist");
    }
    m_dir = dir;
    if (!m_fin->GetListOfKeys()->Contains(dir.c_str()) && dir != "")
    {
        m_fin->ls();
        throw std::runtime_error(dir + " does not exist in the input file");
    }
    m_loaded = true;
}


void Config::updateHistogramPtr(RegionInfo* r, VariableInfo* v)
{
    if (!m_loaded)
    {
        throw std::runtime_error("You must load first");
    }

    if (m_manipulated)
    {
        throw std::runtime_error("Update histogram pointer for manipulated histograms is not supported!");
    }

    current_region = r;
    current_variable = v;

    TDirectory* d = nullptr;
    d = m_fin->GetDirectory(m_dir.c_str());

    for (ProcessInfo* p : *(processes->content()))
    {
        const std::string& fullname = Utils::histString(p, r, v, basic->name_convention);
        if (d->GetListOfKeys()->Contains(fullname.c_str()))
        {
            TH1* h = (TH1*)d->Get(fullname.c_str());
            Utils::histAssign(h, p, r, v);
        }
        else
        {
            clog << "INFO: " << fullname << " is not in " << m_dir << " (skip it)\n";
        }
        // else that p.histogram will remain as nullptr
        // later when make plot this should checked
        // -> well handled by HistTool::check()

        if (systematics)
        {
            for (SystematicInfo* s : *(systematics->content()))
            {
                const std::string& fullname = Utils::histStringSyst(p, r, v, s, basic->name_convention);
                const std::string& fullnameWithSystUp = fullname + "__1up";
                const std::string& fullnameWithSystDown = fullname + "__1down";
                TDirectory* d_syst = nullptr;
                d_syst = (TDirectory*)d->Get("Systematics");
                // d_syst->ls();
                cout << fullnameWithSystUp << endl;
                if (s->type == eSystematicType::TwoSide &&
                    d_syst->GetListOfKeys()->Contains(fullnameWithSystUp.c_str()) &&
                    d_syst->GetListOfKeys()->Contains(fullnameWithSystDown.c_str()))
                {
                    TH1* hUp = (TH1*)d_syst->Get(fullnameWithSystUp.c_str());
                    TH1* hDown = (TH1*)d_syst->Get(fullnameWithSystDown.c_str());

                    if (basic->make_systs_shape_only) 
                    {
                        Tools::println(":( make_systs_shape_only does not work with two side yet");
                    }

                    if (s->smooth) {
                        WSMakerBinning::smoothSyst(p->histogram, hUp);
                        WSMakerBinning::smoothSyst(p->histogram, hDown);
                    }

                    if (s->average) {
                        for (int i = 0; i < p->histogram->GetNbinsX()+2; i++)
                        {
                            float mid = p->histogram->GetBinContent(i);
                            float rel = 0.5 * (hUp->GetBinContent(i) - hDown->GetBinContent(i));
                            hUp->SetBinContent(i, mid + rel);
                            hDown->SetBinContent(i, mid - rel);
                        }
                    }

                    // set styles here for simplicity
                    hUp->SetLineColor(s->color);
                    // hUp->SetFillColor(s->color);
                    hUp->SetMarkerColor(s->color);
                    hUp->SetLineStyle(1);
                    hDown->SetLineColor(s->color);
                    // hDown->SetFillColor(s->color);
                    hDown->SetMarkerColor(s->color);
                    hDown->SetLineStyle(2);

                    if (s->name == s->name_tex)
                    {
                        Utils::histAssignSyst(hUp, p, Utils::systString(s) + "__1up");
                        Utils::histAssignSyst(hDown, p, Utils::systString(s) + "__1down");
                    }
                    else
                    {
                        Utils::histAssignSyst(hUp, p, Utils::systString(s) + " (1up)");
                        Utils::histAssignSyst(hDown, p, Utils::systString(s) + " (1down)");
                    }
                }
                else if (s->type == eSystematicType::OneSide 
                    && (d_syst->GetListOfKeys()->Contains(fullnameWithSystUp.c_str()) 
                    || d_syst->GetListOfKeys()->Contains(fullname.c_str())))
                {
                    TH1* hUp = d_syst->GetListOfKeys()->Contains(fullname.c_str()) ? \
                        (TH1*)d_syst->Get(fullname.c_str()) : (TH1*)d_syst->Get(fullnameWithSystUp.c_str());

                    hUp->SetName(fullnameWithSystUp.c_str());
                    TH1* hDown = (TH1*)hUp->Clone(fullnameWithSystDown.c_str());
                    
                    if (basic->make_systs_shape_only)
                    {
                        Tools::println("! make sure the syst is symmetric");
                        hUp->Scale(p->histogram->Integral() / hUp->Integral());
                    }

                    for (int i = 0; i < hUp->GetNbinsX() + 2; i++)
                    {
                        hDown->SetBinContent(i, p->histogram->GetBinContent(i) * 2. - hUp->GetBinContent(i));
                        /// @todo no bin error implemented yet
                    }

                    if (s->smooth) {
                        WSMakerBinning::smoothSyst(p->histogram, hUp);
                        WSMakerBinning::smoothSyst(p->histogram, hDown);
                    }

                    // set styles here for simplicity
                    hUp->SetLineColor(s->color);
                    // hUp->SetFillColor(s->color);
                    hUp->SetMarkerColor(s->color);
                    hUp->SetLineStyle(1);
                    hDown->SetLineColor(s->color);
                    // hDown->SetFillColor(s->color);
                    hDown->SetMarkerColor(s->color);
                    hDown->SetLineStyle(2);

                    if (s->name == s->name_tex)
                    {
                        Utils::histAssignSyst(hUp, p, Utils::systString(s) + "__1up");
                        Utils::histAssignSyst(hDown, p, Utils::systString(s) + "__1down");
                    }
                    else
                    {
                        Utils::histAssignSyst(hUp, p, Utils::systString(s) + " (1up)");
                        Utils::histAssignSyst(hDown, p, Utils::systString(s) + " (1down)");
                    }
                }
                else
                {
                    clog << "INFO: Can not add systematic uncertainty [" << s->name << "] for [" << p->name << "] (skip it)\n";
                }
            }
        }
    }
    
}

