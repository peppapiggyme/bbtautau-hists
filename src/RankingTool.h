/**
 * @brief: make this tool as a child of WorkSpace
 */

#ifndef RANKING_H
#define RANKING_H

#include "StatisticTool.h"
#include "WorkSpace.h"

#include <fstream>
using std::ofstream;


class RankingEngine : public StatisticTool
{
public:
    RankingEngine(const WorkspaceInfo* info, const string& name="RankingEngine") 
        :  StatisticTool(name), m_cInfo(info) 
    {
        Tools::println("Initializing [%]", this->sToolName);
        RooMsgService::instance().setGlobalKillBelow(ERROR);
    }
    virtual ~RankingEngine() override
    {
        for (auto& p : m_fits)
        {
            if (!p.second)
            {
                delete p.second;
            }
        }
        m_fits.clear();
    }
    RankingEngine(RankingEngine& old) = delete;
    RankingEngine& operator=(RankingEngine& old) = delete;

public:
    virtual void Execute() override
    {
        auto run = [this](const string& sNP, const double nMode) 
        {
            string sPostFix = nMode > 0 ? "_Hi" : "_Lo";
            m_fits[sNP + sPostFix] = new WorkSpace(m_cInfo);
            m_fits[sNP + sPostFix]->FitWithFixedNP(sNP, m_fits["base"]->GetFittedNPs(), nMode);
            // m_fits[sNP + sPostFix]->Check();
            m_mapAltPOIs[sNP + sPostFix] = m_fits[sNP + sPostFix]->GetCache(m_fits[sNP + sPostFix]->NameOfPOI());
            if (!m_bCache)
            {
                delete m_fits[sNP + sPostFix];
            }
        };

        m_fits["base"] = new WorkSpace(m_cInfo);
        Tools::println(">> All: ");
        m_fits["base"]->FitAll();
        m_nPOI = m_fits["base"]->GetCache(m_fits["base"]->NameOfPOI());
        m_fits["base"]->Check();
        Tools::println(">> Fixed NP: ");
        for (const auto& sNP: m_fits["base"]->GetNPs())
        {
            for (double nMode : {1.0f, -1.0f})
            {
                run(sNP, nMode);
            }
        }
        if (!m_bCache)
        {
            delete m_fits["base"];
        }
    }

    map<WorkSpace::ePOI, double> GetFittedPOI()
    {
        return m_nPOI;
    }

    map<string, map<WorkSpace::ePOI, double>> GetFittedPOIWithFixedNPs()
    {
        return m_mapAltPOIs;
    }

public:
    void WriteToTxt(const std::string& outname)
    {
        /// @todo fix the copying
        auto val = GetFittedPOI();
        auto res = GetFittedPOIWithFixedNPs();
        ofstream output(outname);
        output << OUTPUT_TABLE_4("FIXED NP", "POI_VALUE", "POI_ERROR_UP", "POI_ERROR_DOWN");

        Tools::println("% ->      VALUE %      ERROR_UP %      ERROR_DOWN %",
            "No fixed NP",
            val[WorkSpace::ePOI::VALUE],
            val[WorkSpace::ePOI::ERRORUP],
            val[WorkSpace::ePOI::ERRORDOWN]);
        output << OUTPUT_TABLE_4("None",
            val[WorkSpace::ePOI::VALUE],
            val[WorkSpace::ePOI::ERRORUP],
            val[WorkSpace::ePOI::ERRORDOWN])

        for (auto& p : res) {
            Tools::println("% Fixed ->      VALUE %      ERROR_UP %      ERROR_DOWN %",
                p.first,
                p.second[WorkSpace::ePOI::VALUE],
                p.second[WorkSpace::ePOI::ERRORUP],
                p.second[WorkSpace::ePOI::ERRORDOWN]);
            output << OUTPUT_TABLE_4(
                p.first,
                p.second[WorkSpace::ePOI::VALUE],
                p.second[WorkSpace::ePOI::ERRORUP],
                p.second[WorkSpace::ePOI::ERRORDOWN]);
        }
    }

public:
    void CacheWorkSpaces()
    {
        m_bCache = true;
    }

protected:
    bool m_bCache = false;
    const WorkspaceInfo* m_cInfo;
    map<string, WorkSpace*> m_fits;
    map<WorkSpace::ePOI, double> m_nPOI;
    map<string, map<WorkSpace::ePOI, double>> m_mapAltPOIs;
};

struct RankingData
{
    RankingData(const string& np, double hi, double lo)
        : fixedNP(np), deltaMuByHi(hi), deltaMuByLo(lo) {}
    string fixedNP;
    double deltaMuByHi;
    double deltaMuByLo;
};

class RankingPlotter : public StatisticTool
{
public:
    RankingPlotter(const string& name="RankingPlotter") : StatisticTool(name)
    {

    }
    
    void LoadFromTxt(const std::string& sFileName)
    {
        ifstream input(sFileName, std::ios::in);
        if (input.is_open())
        {
            string sLine;
            unsigned nRow = 0;
            while (getline(input, sLine))
            {
                istringstream iss(sLine);
                string tmp;
                if (nRow == 0) 
                {
                    Tools::println("Starting!");
                }
                else if (nRow == 1)
                {
                    iss >> tmp 
                        >> m_nPOI[WorkSpace::ePOI::VALUE] 
                        >> m_nPOI[WorkSpace::ePOI::ERRORUP] 
                        >> m_nPOI[WorkSpace::ePOI::ERRORDOWN]; // symmetric of up due to Hesse
                }
                else
                {
                    iss >> tmp 
                        >> m_mapAltPOIs[tmp][WorkSpace::ePOI::VALUE] 
                        >> m_mapAltPOIs[tmp][WorkSpace::ePOI::ERRORUP] 
                        >> m_mapAltPOIs[tmp][WorkSpace::ePOI::ERRORDOWN]; // symmetric of up due to Hesse
                }
                nRow++;
            }
        }

        for (auto& p : m_mapAltPOIs)
        {
            string sNameOfNP = p.first.substr(0, p.first.length()-3);
            string sDirection = p.first.substr(p.first.length()-2);
            if (m_mapLookUp.find(sNameOfNP) == m_mapLookUp.end()) 
            {
                m_mapLookUp[sNameOfNP] = m_vData.size();
                if (sDirection == "Hi")
                {
                    m_vData.emplace_back(sNameOfNP, p.second[WorkSpace::ePOI::VALUE]-m_nPOI[WorkSpace::ePOI::VALUE], -999);
                }
                else 
                {
                    m_vData.emplace_back(sNameOfNP, -999, p.second[WorkSpace::ePOI::VALUE]-m_nPOI[WorkSpace::ePOI::VALUE]);
                }
            }
            else 
            {
                if (sDirection == "Hi")
                {
                    m_vData[m_mapLookUp[sNameOfNP]].deltaMuByHi = p.second[WorkSpace::ePOI::VALUE]-m_nPOI[WorkSpace::ePOI::VALUE];
                }
                else
                {
                    m_vData[m_mapLookUp[sNameOfNP]].deltaMuByLo = p.second[WorkSpace::ePOI::VALUE]-m_nPOI[WorkSpace::ePOI::VALUE];
                }
            }
        }

        double fSumSQ = 0.;
        for (const RankingData& data : m_vData)
        {
            fSumSQ += 0.25 * (std::fabs(data.deltaMuByHi - data.deltaMuByLo) * std::fabs(data.deltaMuByHi - data.deltaMuByLo));
        }
        fSumSQ = std::sqrt(fSumSQ);

        std::sort(m_vData.begin(), m_vData.end(), [&fSumSQ](RankingData& a, RankingData& b) {
            return std::fabs(a.deltaMuByHi - a.deltaMuByLo) > std::fabs(b.deltaMuByHi - b.deltaMuByLo);
        });
        
        unsigned i = 1;
        for (auto& info : m_vData)
        {
            Tools::println("% % % %", i++, info.fixedNP, info.deltaMuByHi, info.deltaMuByLo);
        }
        Tools::println("Sum in quadrature = %", fSumSQ);
    }

    void Draw(const string output, uint8_t nNPs=15)
    {
        gROOT->SetStyle("ATLAS");
        gStyle->SetErrorX(0.5);

        TCanvas* c1 = new TCanvas("c", "", 900, 1200);
        c1->SetLeftMargin(0.20);
        c1->SetBottomMargin(0.50);
        TH1* h_one = new TH1D("one", "", nNPs, 0, nNPs);
        TH1* h_hi = new TH1D("hi", "", nNPs, 0, nNPs);
        TH1* h_lo = new TH1D("lo", "", nNPs, 0, nNPs);
        for (std::size_t i = 1; i <= nNPs; i++) {
            h_hi->GetXaxis()->SetBinLabel(i, Utils::systStringShort(m_vData[i-1].fixedNP).c_str());
            h_hi->SetBinContent(i, m_vData[i-1].deltaMuByHi);
            h_lo->SetBinContent(i, m_vData[i-1].deltaMuByLo);
            h_one->SetBinContent(i, 0);
        }

        h_hi->SetBarWidth(0.90);
        h_lo->SetBarWidth(0.90);

        h_hi->SetLineColor(kRed);
        h_hi->SetLineWidth(1);
        h_hi->SetFillColor(kRed);
        h_hi->SetFillStyle(3245);
        h_hi->SetMarkerSize(0);
        h_hi->GetXaxis()->SetLabelSize(0.036);
        h_hi->GetXaxis()->SetTitleSize(0.036);
        h_hi->GetXaxis()->LabelsOption("v");
        double max_of_hist = std::max(std::abs(h_hi->GetMaximum()), std::abs(h_lo->GetMaximum()));
        h_hi->GetYaxis()->SetRangeUser(-1.6*max_of_hist, 1.6*max_of_hist);
        // h_hi->GetYaxis()->SetRangeUser(-0.5, 0.5);
        h_hi->GetYaxis()->SetLabelSize(0.036);
        h_hi->GetYaxis()->SetTitleSize(0.036);
        h_hi->GetYaxis()->SetTitleOffset(2);
        h_hi->GetYaxis()->SetTitle("#Delta#mu");
        h_lo->SetLineColor(kBlue);
        h_lo->SetLineWidth(1);
        h_lo->SetFillColor(kBlue);
        h_lo->SetFillStyle(3254);
        h_lo->SetMarkerSize(0);

        h_one->SetLineColor(kRed+1);
        h_one->SetLineWidth(3);

        h_hi->Draw("hist b");
        h_lo->Draw("hist b same");
        h_one->Draw("hist same");

        TLegend* legend = new TLegend(0.50, 0.85, 0.80, 0.93);
        legend->SetTextFont(42);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.035);
        legend->SetTextAlign(12);
        legend->AddEntry(h_hi, "+1#sigma post-fit impact on #mu", "f");
        legend->AddEntry(h_lo, "-1#sigma post-fit impact on #mu", "f");
        legend->Draw("same");

        TLatex *text = new TLatex();
        text->SetNDC();
        text->SetTextFont(72);
        text->SetTextSize(0.050);
        text->DrawLatex(0.50, 0.52, "ATLAS");
        text->SetTextFont(42);
        text->DrawLatex(0.50 + 0.18, 0.52, "Internal");
        text->SetTextSize(0.045);
        text->SetTextSize(0.040);
        text->DrawLatex(0.60, 0.96, "Fit to Data");

        c1->SaveAs(output.c_str());

        delete h_hi;
        delete h_lo;
        delete legend;
        delete text;
        delete c1;
    }

private:
    map<WorkSpace::ePOI, double> m_nPOI;
    map<string, map<WorkSpace::ePOI, double>> m_mapAltPOIs;
    vector<RankingData> m_vData;
    map<string, std::size_t> m_mapLookUp;
};

#endif