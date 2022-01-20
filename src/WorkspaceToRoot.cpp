#include "WorkspaceToRoot.h"
#include "Utils.h"
#include "roofitUtils.h"
#include "RooExpandedFitResult.h"

#include <regex>

using namespace std;

WorkspaceToRoot::WorkspaceToRoot(const WorkspaceInfo* info, const W2RInfo* w2rInfo, const string& name) 
    : StatisticTool(name)
    , m_cInfo(info)
    , m_cW2RInfo(w2rInfo)
    , m_cWs(make_unique<WorkSpace>(info))
{
    Tools::println("Initializing [%]", this->sToolName);
    RooMsgService::instance().setGlobalKillBelow(ERROR);

    if (m_cW2RInfo->verbose) 
    {
        m_cWs->Check();
    }
}

WorkspaceToRoot::~WorkspaceToRoot()
{

}

void WorkspaceToRoot::Execute()
{
    FetchContents(); 
    GetFitResult();
    MakeYield();
}

void WorkspaceToRoot::WriteToRootfile(const std::string& sOutName)
{
    Tools::println("Objects will be stored in [%]", sOutName);
    TFile* cOutput = TFile::Open(sOutName.c_str(), "recreate");

    for (const auto& pp : m_mapContent)
    {
        Tools::println("Storing objects in %", pp.first);

        TDirectory* cDirCat = cOutput->mkdir(pp.first.c_str());
        cDirCat->cd();
        RooPlot* cFrame = m_mapObs[pp.first]->frame();

        // data
        Tools::println("  -> Data");
        ((RooAbsData*)m_mapContent[pp.first]["Data"])->plotOn(cFrame, RooFit::DataError(RooAbsData::Poisson));
        RooHist* cRooHistData = cFrame->getHist();
        cRooHistData->SetName("RooData");
        cRooHistData->Write();

        // chi2
        Tools::println("  -> Chi2");
        m_mapPdf[pp.first]->plotOn(cFrame, RooFit::Normalization(1, RooAbsReal::RelativeExpected));
        double fChi2 = cFrame->chiSquare();
        TH1F* cHistChi2 = new TH1F("chi2", "", 1, 0, 1);
        cHistChi2->SetBinContent(1, fChi2);
        cHistChi2->Write();
        delete cHistChi2;

        TH1F* cHistData = nullptr;
        TH1F* cHistBkg = nullptr;

        // components
        for (const auto& qq : pp.second)
        {
            if (qq.first == "Data")
            {
                continue;
            }
            if (qq.first == "Signal" && IsCR(pp.first))
            {
                continue;
            }
            Tools::println("  -> %", qq.first);
            TH1F* cHist = (TH1F*)((RooAbsReal*)qq.second)->createHistogram(qq.first.c_str(), *m_mapObs[pp.first]);
            if (cHist->Integral() != 0.)
            {
                cHist->Scale(m_mapYields[pp.first][qq.first].val / cHist->Integral());
                cHist->SetName(qq.first.c_str());
                cHist->Write();
            }
            if (!cHistData)
            {
                cHistData = (TH1F*)cHist->Clone();
            }
            if (!cHistBkg && qq.first == "Background")
            {
                cHistBkg = (TH1F*)cHist->Clone();
            }
        }

        // data histogram -> only for the sake of DST checkings
        Tools::println("Generating data histogram nbins [%]", cHistData->GetNbinsX());
        for (int i = 0; i < cRooHistData->GetN(); i++)
        {
            cHistData->SetBinContent(i + 1, cRooHistData->GetPointY(i));
        }
        cHistData->SetName("data");
        cHistData->Write();

        // error
        Tools::println("  -> Error");
        double fMu = ((RooRealVar*)m_cWs->GetRooPOIs()->first())->getVal();
        ((RooRealVar*)m_cWs->GetRooPOIs()->first())->setVal(0.);
        m_mapPdf[pp.first]->plotOn(cFrame, 
            RooFit::VisualizeError(*m_cPostFitResult, 1), 
            RooFit::Name("PostfitError"), 
            RooFit::Normalization(1, RooAbsReal::RelativeExpected));
        ((RooRealVar*)m_cWs->GetRooPOIs()->first())->setVal(fMu);
        RooCurve* cRooCurveError = cFrame->getCurve();
        cRooCurveError->Write();

        /// @todo investigate the error propagation, why it is symmetric?
        /// anyway put it into histogram to maintain underflow and overflow bins
        TH1F* cHistError = (TH1F*)cHistBkg->Clone();
        cHistError->SetName("error");
        const int nBins = cHistBkg->GetNbinsX();
        const int nPoints = cRooCurveError->GetN();
        if (nPoints != (((nBins+3)<<1)+1)<<1)
        {
            Tools::println("inconsistent nbins [%] and npoint [%]", nBins, nPoints);
            throw std::runtime_error("inconsistency in cHistError construction!");
        }
        
        Tools::println("Generating error histogram, nbins [%]", cHistError->GetNbinsX());
        for (int i = 1; i < nPoints - i - 1; i += 2)
        {
            double fErrorUp = cRooCurveError->GetPointY(i);
            double fErrorDown = cRooCurveError->GetPointY(nPoints-i-1);
            cHistError->SetBinError(i >> 1, 0.5 * (fErrorUp - fErrorDown));
            Tools::println(" -> Bin [%] % +- % [up-%, down-%]",
                i >> 1, cHistError->GetBinContent(i >> 1), cHistError->GetBinError(i >> 1),
                fErrorUp, fErrorDown
            );
        }
        cHistError->Write();
        
        // postfit curve
        Tools::println("  -> Pre-fit background");
        TH1F* cHistPrefit = m_mapHistPrefit[pp.first];
        cHistPrefit->SetName("PrefitBkg");
        cHistPrefit->Write();

        cOutput->cd();
    }

    cOutput->Close();
    cOutput = nullptr;
    Tools::println("Objects have been stored in [%]", sOutName);
}

void WorkspaceToRoot::FetchContents()
{
    RooSimultaneous* cSimPdf = (RooSimultaneous*)m_cWs->GetRooWorkspace()->pdf("simPdf");
    auto& cCategories = cSimPdf->indexCat();
    string sChannel{cSimPdf->indexCat().GetName()};
    
    Tools::println("Channel [%]", sChannel);
    for (auto cat : cCategories)
    {
        string sCategory{cat->GetName()};
        Tools::println("find category: %", sCategory);

        m_mapContent[sCategory];

        RooAbsData* cDataTmp = m_cWs->GetRooData()->reduce((sChannel+"=="+sChannel+"::"+sCategory).c_str());
        m_mapPdf[sCategory] = cSimPdf->getPdf(sCategory.c_str());
        m_mapObs[sCategory] = (RooRealVar*)m_mapPdf[sCategory]->getObservables(m_cWs->GetModelConfig()->GetObservables())->first();

        m_mapBinWidth[sCategory] = m_mapPdf[sCategory]->getVariables()->find(("binWidth_obs_x_"+sCategory+"_0").c_str());
        if (!m_mapBinWidth[sCategory])
        {
            throw std::runtime_error("BinWidth not found!");
        }

        const RooAbsBinning& cBinning = m_mapObs[sCategory]->getBinning();
        const double fBinSize = cBinning.averageBinWidth();
        const double fBinLow = cBinning.lowBound();
        const double fBinHigh = cBinning.highBound();
        Tools::println("Binning [Size=%, Low=%, High=%, nBins=%]",
            fBinSize, fBinLow, fBinHigh, (size_t)((fBinHigh - fBinLow) / fBinSize));

        double fBinEdge = fBinLow;
        while (fBinEdge < fBinHigh)
        {
            string sBinIndexStr{"bin"+to_string((int)(fBinEdge / fBinSize - fBinLow))};
            Tools::println("  -> Setting range for %", sBinIndexStr);
            m_mapObs[sCategory]->setRange(sBinIndexStr.c_str(), fBinEdge, fBinEdge+fBinSize);
            fBinEdge += fBinSize;
        }

        m_mapContent[sCategory]["Data"] = (RooAbsArg*)cDataTmp;

        RooArgList cListBkg, cListSig, cListTotal;
        string sModelName{sCategory+"_model"};
        RooRealSumPdf* cModel = (RooRealSumPdf*)m_mapPdf[sCategory]->getComponents()->find(sModelName.c_str());
        if (!cModel)
        {
            throw std::runtime_error("Model not found!");
        }

        Tools::println("Fetching components of the pdf");
        for (auto component : cModel->funcList())
        {
            string sComponent{component->GetName()};
            std::regex reComp("L_x_(.*)_Region");
            std::smatch reCompMatch;
            if (std::regex_search(sComponent, reCompMatch, reComp)) 
            {
                sComponent = reCompMatch[1]; // 1 is the sub match
            }
            else 
            {
                throw std::runtime_error("no match");
            }

            Tools::println("  -> Found component: %", sComponent.c_str());
            m_mapContent[sCategory][sComponent] = component;
            cListTotal.add(*component);
            if (IsSignal(sComponent)) 
            {
                Tools::println("    * This is signal!");
                cListSig.add(*component);
            }
            else 
            {
                cListBkg.add(*component);
            }
        }
        m_mapContent[sCategory]["Total"] = new RooAddition("Total", "tot_sum", cListTotal);
        m_mapContent[sCategory]["Signal"] = new RooAddition("Signal", "sig_sum", cListSig);
        m_mapContent[sCategory]["Background"] = new RooAddition("Background", "bkg_sum", cListBkg);
    }
}

void WorkspaceToRoot::MakeYield()
{
    set<string> setSpecialItems = {
        "Data", "Total", "Signal", "Background"
    };

    for (const auto& pp : m_mapContent)
    {
        Tools::println("Yield in Category [%]", pp.first);

        // data
        Yield yData{((RooAbsData*)pp.second.at("Data"))->sumEntries(), 0.};
        Tools::println("|- Data %", yData.val);

        // special        
        Yield yTotal = GetIntegralAndError(pp.first, pp.second.at("Total"));
        Yield ySignal = GetIntegralAndError(pp.first, pp.second.at("Signal"));
        Yield yBackground = GetIntegralAndError(pp.first, pp.second.at("Background"));
        Tools::println("|- Total % +- %", yTotal.val, yTotal.err);
        Tools::println("|- Signal % +- %", ySignal.val, ySignal.err);
        Tools::println("|- Background % +- %", yBackground.val, yBackground.err);

        // put values into a cache
        m_mapYields[pp.first];
        m_mapYields[pp.first]["Data"] = yData;
        m_mapYields[pp.first]["Total"] = yTotal;
        m_mapYields[pp.first]["Signal"] = ySignal;
        m_mapYields[pp.first]["Background"] = yBackground;

        // components
        for (const auto& qq : pp.second)
        {
            if (setSpecialItems.find(qq.first) != setSpecialItems.end())
            {
                continue;
            }
            Yield y = GetIntegralAndError(pp.first, qq.second);
            Tools::println("|- % % +- %", qq.first, y.val, y.err);
            m_mapYields[pp.first][qq.first] = y;
        }

        Tools::println("\n\n");
    }
}

Yield WorkspaceToRoot::GetIntegralAndError(const string& sCategory, RooAbsArg* cComponent) const
{
    Yield y;
    RooArgSet cObsSet(*m_mapObs.at(sCategory));
    auto cIntegral = ((RooAbsReal*)cComponent)->createIntegral(cObsSet);
    float fBinWidth = ((RooAbsReal*)m_mapBinWidth.at(sCategory))->getVal();
    y.val = cIntegral->getVal() * fBinWidth;
    y.err = RU::getPropagatedError(*cIntegral, *m_cPostFitResult.get()) * fBinWidth;

    return y;
}

void WorkspaceToRoot::GetFitResult()
{
    double fMu = m_cW2RInfo->mu;

    Tools::println("Forcing POI = %", fMu);
    ((RooRealVar*)m_cWs->GetRooPOIs()->first())->setVal(fMu);
    ((RooRealVar*)m_cWs->GetRooPOIs()->first())->setError(1e-6);
    ((RooRealVar*)m_cWs->GetRooPOIs()->first())->setConstant(true);

    for (auto& np : *m_cWs->GetRooNPs())
    {
        string sNameNP{np->GetName()};
        Tools::println("Init NP value and errors are [%] = [%] +- [%]",
            sNameNP, ((RooRealVar*)np)->getVal(), ((RooRealVar*)np)->getError());

        if (sNameNP.substr(0, 5) == "alpha")
        {
            ((RooRealVar*)np)->setVal(0);
            ((RooRealVar*)np)->setError(1);
        }
        else if (sNameNP.substr(0, 5) == "gamma")
        {
            RooAbsPdf* cGamma = (RooAbsPdf*)m_cWs->GetRooWorkspace()->pdf((sNameNP+"_constraint").c_str());
            RooRealVar* cTau = (RooRealVar*)m_cWs->GetRooWorkspace()->obj((sNameNP+"_tau").c_str());
            if (cGamma && cTau)
            {
                ((RooRealVar*)np)->setVal(1);
                ((RooRealVar*)np)->setError(TMath::Sqrt(1.0 / cTau->getVal()));
            }
            else 
            {
                throw runtime_error("gamma constrain not exists");
            }
        }
        else if (sNameNP == "ATLAS_norm_Zhf")
        {
            ((RooRealVar*)np)->setVal(1.);
            ((RooRealVar*)np)->setError(1e-6);
        }
        else if (sNameNP == "ATLAS_norm_ttbar")
        {
            ((RooRealVar*)np)->setVal(1.);
            ((RooRealVar*)np)->setError(1e-6);
        }        
    }
    RooArgList cParams(*m_cWs->GetRooPOIs());
    cParams.add(*m_cWs->GetRooNPs());

    m_cPreFitResult = make_unique<RooExpandedFitResult>(cParams);
    if (m_cW2RInfo->verbose)
    {
        Tools::print("\n========== PRE-FIT RESULT ==========");
        m_cPreFitResult->Print();
        Tools::print("========== PRE-FIT RESULT ==========\n");
    }

    /// @todo mutmp is same as m_W2RInfo->mu, no?
    double fMuTmp = ((RooRealVar*)m_cWs->GetRooPOIs()->first())->getVal();
    ((RooRealVar*)m_cWs->GetRooPOIs()->first())->setVal(0.);
    for (const auto& pp : m_mapObs)
    {
        RooArgSet cObsSet(*pp.second);
        double fPrefitYield = m_mapPdf[pp.first]->expectedEvents(cObsSet);
        TH1F* cHistPrefit = (TH1F*)m_mapPdf[pp.first]->createHistogram("PrefitBkg", *pp.second);
        cHistPrefit->Scale(fPrefitYield / cHistPrefit->Integral());
        m_mapHistPrefit[pp.first] = cHistPrefit;
    }
    ((RooRealVar*)m_cWs->GetRooPOIs()->first())->setVal(fMuTmp);

    if (m_cW2RInfo->prefit)
    {
        m_cPostFitResult.reset((RooFitResult*)m_cPreFitResult->Clone());
    }
    else
    {   
        if (!m_cW2RInfo->postfit_from_file)
        {
            throw std::runtime_error("adhoc postfit not supported yet");
        }
        TFile* cFilePostFit = TFile::Open(m_cW2RInfo->postfit_result_file.c_str());
        m_cPostFitResult.reset(
            new RooExpandedFitResult((RooFitResult*)cFilePostFit->Get(m_cW2RInfo->postfit_result_directory.c_str())->Clone(), RooArgList())
        );

        if (!m_cPostFitResult)
        {
            throw std::runtime_error("postfit result not found, check path and directory");
        }

        RooArgList* cListPars = (RooArgList*)m_cPostFitResult->floatParsFinal().Clone();
        RooArgList* cListConstPars = (RooArgList*)m_cPostFitResult->constPars().Clone();
        cListPars->add(*cListConstPars);

        if (m_cW2RInfo->verbose)
        {
            cListPars->Print();
        }

        for (auto& np : *m_cWs->GetRooNPs())
        {
            string sNameNP{np->GetName()};
            auto cTheirs = cListPars->find(*np);
            if (!cTheirs) continue;
            ((RooRealVar*)np)->setVal(((RooRealVar*)cTheirs)->getVal());
            ((RooRealVar*)np)->setError(((RooRealVar*)cTheirs)->getError());

            Tools::println("Post-fit NP value and errors are [%] = [%] +- [%]",
                sNameNP, ((RooRealVar*)np)->getVal(), ((RooRealVar*)np)->getError());
        }
    }

    if (m_cW2RInfo->verbose && !m_cW2RInfo->prefit)
    {
        Tools::print("\n========== POST-FIT RESULT ==========");
        m_cPostFitResult->Print();
        Tools::print("========== POST-FIT RESULT ==========\n");
    }
    // else
}

bool WorkspaceToRoot::IsSignal(const string& sCompName)
{
    const set<string> setSignalPatterns = { "hhttbb", "Hhhbbtautau" };

    for (const auto& p : setSignalPatterns)
    {
        if (sCompName.find(p) != string::npos)
        {
            return true;
        }
    }

    return false;
}

bool WorkspaceToRoot::IsCR(const string& sCateName)
{
    const set<string> setCRPatterns = {"DZllbbCR"};

    for (const auto& p : setCRPatterns)
    {
        if (sCateName.find(p) != string::npos)
        {
            return true;
        }
    }

    return false;
}
