#ifndef PULLS_H
#define PULLS_H

#include "StatisticTool.h"
#include "WorkSpace.h"

#include <fstream>

using std::ofstream;

class PullsEngine : public StatisticTool
{
public:
    PullsEngine(const WorkspaceInfo* info, const string& name="PullsEngine") 
        :  StatisticTool(name), m_cInfo(info) 
    {
        Tools::println("Initializing [%]", this->sToolName);
        RooMsgService::instance().setGlobalKillBelow(ERROR);
    }
    virtual ~PullsEngine() override
    {
        delete m_cTool;
    }
    PullsEngine(PullsEngine& old) = delete;
    PullsEngine& operator=(PullsEngine& old) = delete;

public:
    virtual void Execute() override 
    {
        m_cTool = new WorkSpace(m_cInfo);
        Tools::println(">> All: ");
        m_cTool->FitWithFixedMu(0.);
        m_mapFittedNPs = m_cTool->GetFittedNPs();
    }

public:
    void WriteToTxt(const string& outname)
    {
        ofstream output(outname);
        for (auto& np : m_mapFittedNPs)
        {
            output << OUTPUT_TABLE_4(
                np.first,
                std::get<0>(np.second),
                std::get<1>(np.second),
                std::get<2>(np.second));
        }
    }

private:
    const WorkspaceInfo* m_cInfo;
    WorkSpace* m_cTool;
    map<string, tuple<double, double, double>> m_mapFittedNPs;
};

class PullsPlotter : public StatisticTool
{
public:
    PullsPlotter(const string& name="PullsPlotter") : StatisticTool(name)
    {

    }
    
    void LoadFromTxt(const std::string& sFileName)
    {
        m_fCounter = 0.5;
        constexpr double offset = 0.0; // maybe useful when making comparison plots
        ifstream input(sFileName, std::ios::in);
        if (input.is_open())
        {
            string sLine;
            while (getline(input, sLine))
            {
                istringstream iss(sLine);
                Parameter par;
                iss >> par.name >> par.val >> par.err_hi >> par.err_lo;
                if (par.name.find("gamma") == std::string::npos)
                {
                    vLabels.push_back(par.name);
                    vX.push_back(m_fCounter + offset);
                    vX_ErrHi.push_back(0.0);
                    vX_ErrLo.push_back(0.0);
                    vY.push_back(par.val);
                    vY_ErrHi.push_back(std::fabs(par.err_hi));
                    vY_ErrLo.push_back(std::fabs(par.err_lo));
                    m_fCounter += 1.0;
                }
            }
        }
    }

    void Draw(const string output)
    {
        gROOT->SetStyle("ATLAS");
        gStyle->SetErrorX(0.5);

        TCanvas* c1 = new TCanvas("c", "", 1800, 800);
        c1->SetBottomMargin(0.50);
        c1->SetLeftMargin(0.05);
        c1->SetRightMargin(0.025);

        TH1* h = new TH1D("h", "", vLabels.size(), 0, vLabels.size());
        TH1* h2 = new TH1D("h2", "", vLabels.size(), 0, vLabels.size());

        h2->GetYaxis()->SetTitle("Pull");
        h2->GetYaxis()->SetTitleOffset(0.5);
        h2->GetYaxis()->SetTitleSize(0.025);
        h2->GetYaxis()->SetLabelSize(0.025);
        h2->GetYaxis()->SetRangeUser(-4.2, 4.2);
        h2->GetXaxis()->SetRangeUser(0., m_fCounter - 0.5);
        h2->GetXaxis()->LabelsOption("v");
        h2->GetXaxis()->SetLabelSize(0.025);

        h->SetLineWidth(0);
        h->SetFillColor(kGreen+1);
        h->SetMarkerSize(0);

        h2->SetLineWidth(0);
        h2->SetFillColor(kYellow+1);
        h2->SetMarkerSize(0);

        for (std::size_t i = 1; i <= vLabels.size(); ++i)
        {
            cout << i << ",  " << vLabels.at(i-1).c_str() << endl;
            //auto ibin = h->GetXaxis()->FindBin((double)i - 0.5);
            h->SetBinContent(i, 0.0);
            h2->SetBinContent(i, 0.0);
            h->SetBinError(i, 1.0);
            h2->SetBinError(i, 2.0);
            h2->GetXaxis()->SetBinLabel(i, Utils::systStringShort(vLabels.at(i-1)).c_str());
        }

        h->LabelsDeflate("X");
        h->LabelsDeflate("Y");
        
        h2->GetXaxis()->Paint("R");
        h2->Draw("E2");
        h->Draw("E2 SAME");
        
        TGraphAsymmErrors* g_nps = new TGraphAsymmErrors(vLabels.size(), 
            &vX[0], &vY[0], &vX_ErrHi[0], &vX_ErrLo[0], &vY_ErrHi[0], &vY_ErrLo[0]);

        g_nps->SetMarkerSize(1.0);
        g_nps->SetMarkerColor(kBlack);
        g_nps->SetMarkerStyle(21);

        g_nps->Draw("P E0 SAME");

        TLegend* legend = new TLegend(0.75, 0.85, 0.90, 0.93);
        legend->SetTextFont(42);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.035);
        legend->SetTextAlign(12);
        legend->AddEntry(h, "1#sigma");
        legend->AddEntry(h2, "2#sigma");
        legend->Draw("same");

        TLatex *text = new TLatex();
        text->SetNDC();
        text->SetTextFont(72);
        text->SetTextSize(0.040);
        text->DrawLatex(0.10, 0.96, "ATLAS");
        text->SetTextFont(42);
        text->DrawLatex(0.10 + 0.08, 0.96, "Internal");
        text->SetTextSize(0.035);
        text->SetTextSize(0.030);
        text->DrawLatex(0.85, 0.96, "Fit to Data");

        c1->SaveAs(output.c_str());

        delete g_nps;
        delete h;
        delete h2;
        delete legend;
        delete text;
        delete c1;
    }

public:
    float m_fCounter;
    std::vector<std::string> vLabels;
    std::vector<double> vX;
    std::vector<double> vX_ErrHi;
    std::vector<double> vX_ErrLo;
    std::vector<double> vY;
    std::vector<double> vY_ErrHi;
    std::vector<double> vY_ErrLo;
};

#endif