#include "WorkSpace.h"

void WorkSpace::Check() {
    Tools::println("%", m_cWs);
    Tools::println("%", m_cSBModel->GetName());

    const RooArgSet* observables = static_cast<const RooArgSet*>(m_cSBModel->GetObservables());
    Tools::println("Observables:");
    PrintObservables(observables);
    Tools::println("POIs:");
    PrintParametersOfInterest();
    Tools::println("NPs:");
    PrintNuisanceParameters();
}

void WorkSpace::CheckParaMaps()
{
    PRINT_PARA_MAP(m_mapNPsInit);
    PRINT_PARA_MAP(m_mapNPsFitted);
    PRINT_PARA_MAP(m_mapPOIsFitted);
}

void WorkSpace::Load() {
    // use the singleton to load file and workspace
    // this makes sure that m_cWs is loaded again for new WorkSpace instance
    // while the root file is opened only once.
    WorkspaceLoad::instance().LoadRootFile(m_cInfo);
    m_cWs = WorkspaceLoad::instance().GetWorkspace();
    
    m_cSBModel = static_cast<ModelConfig*>(m_cWs->obj(m_cInfo->config_name.c_str()));
    m_cNPs = m_cSBModel->GetNuisanceParameters();
    m_cPOIs = m_cSBModel->GetParametersOfInterest();
    m_sPOIName = static_cast<RooRealVar*>(m_cPOIs->first())->GetName();
    m_cData = m_cWs->data(m_cInfo->data_name.c_str());

    if (m_cInfo->use_asimov) {
        Tools::println("Using Asimov data! with mu = %", m_cInfo->mu_asimov);
        static_cast<RooRealVar*>(m_cPOIs->first())->setVal(m_cInfo->mu_asimov);
        RooArgSet* allParams = m_cSBModel->GetPdf()->getParameters(*m_cData);
        RooArgSet globObs("globObs"); // ?
        m_cData = AsymptoticCalculator::MakeAsimovData(*m_cSBModel, *allParams, globObs);
    }

    GetSetOfNPs(m_cNPs);
}

RooFitResult* WorkSpace::OneLinerFit(RooArgSet& cConstrainParas, RooAbsData& cData)
{
    std::cout << "Fitting with: \n"
    "o------------------------------------------------------------o\n"
    "| [One liner version]                                        |\n"
    "o------------------------------------------------------------o\n"
    "| This is quicker than the customized version (*)            |\n"
    "| (*) need to make more tests                                |\n"
    "o------------------------------------------------------------o\n";

    RooFitResult* cRes = m_cSBModel->GetPdf()->fitTo(
        cData, InitialHesse(false), Minos(false), Minimizer("Minuit", "Migrad"),
        Strategy(1), PrintLevel(m_cInfo->logLevel), Constrain(cConstrainParas), Save(true),
        Offset(RooStats::IsNLLOffset()));

    return cRes;
}

RooFitResult* WorkSpace::CustomizedFit(RooArgSet& cConstrainParas, RooAbsData& cData)
{
    std::cout << "Fitting with: \n"
    "o------------------------------------------------------------o\n"
    "| [Customized version]                                       |\n"
    "o------------------------------------------------------------o\n"
    "| The speed depends on the tolarence                         |\n"
    "o------------------------------------------------------------o\n";

    RooMsgService::instance().setGlobalKillBelow(ERROR);

    // Define NLL
    m_cNLL = m_cSBModel->GetPdf()->createNLL(cData, 
        Constrain(cConstrainParas), 
        GlobalObservables(*(m_cSBModel->GetGlobalObservables())), 
        Offset(RooStats::IsNLLOffset()), NumCPU(4));

    // Staring NLL value
    Tools::println("Starting NLL value = %", m_cNLL->getVal());

    // Define Minimizer
    RooMinimizer cMinimizer(*m_cNLL);
    // Print level
    cMinimizer.setPrintLevel(m_cInfo->logLevel);
    // Strategy
    const int nStrategy = ROOT::Math::MinimizerOptions::DefaultStrategy();
    cMinimizer.setStrategy(nStrategy);
    // double fTolerance = ROOT::Math::MinimizerOptions::DefaultTolerance();
    double fTolerance = m_cInfo->use_relative_tolerance ? m_cInfo->tolerance * m_cNLL->getVal() : m_cInfo->tolerance;
    cMinimizer.setEps(fTolerance);
    int nStatus = -1;
    constexpr int nOptConstFlag = 2; // ?
    cMinimizer.optimizeConst(nOptConstFlag);
    const std::string sAlgorithm = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();
    // Max calls and iterations
    const int nMaxCalls = ROOT::Math::MinimizerOptions::DefaultMaxFunctionCalls();
    const int nMaxIterations = ROOT::Math::MinimizerOptions::DefaultMaxIterations();
    cMinimizer.setMaxFunctionCalls(nMaxCalls);
    cMinimizer.setMaxIterations(nMaxIterations);
    // Error engine
    constexpr bool bUseMinos = false;

    Tools::println("Tolerance = [%]", fTolerance);
    Tools::println("Max Calls = [%], Max Iterations = [%]", nMaxCalls, nMaxIterations);

    // Do minimization
    nStatus = cMinimizer.minimize("Minuit2", sAlgorithm.c_str());
    Tools::println("Minimize status = [%]", nStatus);

    // Save the fit result pre
    RooFitResult* cResPre = cMinimizer.save();

    // Covariant matrix
    const TMatrixDSym cCovMatrix = cResPre->covarianceMatrix();
    double fDeterminant = cCovMatrix.Determinant();
    if (m_cInfo->logLevel > 0)
    {
        Tools::println("Covariant matrix determinant = [%]", fDeterminant);
    }

    // Eigen value maker
    TMatrixDSymEigen cEigenValueMaker(cCovMatrix);
    TVectorT<double> cEigenValues = cEigenValueMaker.GetEigenValues();
    // TMatrixT<double> cEigenVectors = cEigenValueMaker.GetEigenVectors();
    if (m_cInfo->logLevel > 0)
    {
        for (int i = 0; i < cEigenValues.GetNrows(); ++i) 
        {
            Tools::println("Eigen values -> %", cEigenValues[i]);
        }
    }

    // Improve fit result
    if (nStatus % 100 == 0) 
    {
        cMinimizer.setMinimizerType("Minuit2");
        nStatus = cMinimizer.hesse();
        Tools::println("Hesse status = [%]", nStatus);
        if (bUseMinos)
        {
            cMinimizer.minos();
        }
    }
    else 
    {
        Tools::println("Minimize failed with status code [%]", nStatus);
    }

    // Final NLL value
    typedef std::numeric_limits<double> dbl;
    std::cout.precision(dbl::digits10);
    Tools::println("Final NLL value = %", m_cNLL->getVal());

    RooFitResult* cRes = cMinimizer.save();

    delete cResPre;

    return cRes;
}

void WorkSpace::Fit(RooAbsData* cData) {
    if (!cData)
    {
        cData = m_cData;
    }
    // Information from workspace
    RooRealVar* cPOI = static_cast<RooRealVar*>(m_cPOIs->first());
    auto sMinimizerType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
    Tools::println("POI [%] initial value is [%]", cPOI->GetName(), cPOI->getVal());
    Tools::println("Set [%] as minimizer", sMinimizerType);
    // Remove constant paras
    RooArgSet cConstrainParas;
    cConstrainParas.add(*m_cNPs);
    RooStats::RemoveConstantParameters(&cConstrainParas);
    auto timeStart = steady_clock::now();

    RooFitResult* cRes = nullptr;

    // >>> core of fitting <<< START
    if (m_cInfo->use_oneline_fit)
    {
        cRes = OneLinerFit(cConstrainParas, *cData);
    }
    else 
    {
        cRes = CustomizedFit(cConstrainParas, *cData);
    }
    // >>> core of fitting <<< END

    auto timeEnd = steady_clock::now();

    bFitted = true;
    if (m_cInfo->logLevel > 0)
        cRes->Print();
    Tools::println("Fit status [%], spent [%ms]",
            cRes->status(), duration_cast<milliseconds>(timeEnd-timeStart).count());
    Tools::println("POI [%] final value is [%  %  %]",
            cPOI->GetName(), cPOI->getVal(), cPOI->getErrorHi(), cPOI->getErrorLo());

    if (cRes)
    {
        delete cRes;
    }
}

void WorkSpace::FitAll() {
    // prefit paramters (not useful)
    UpdateMapNPsInit(m_cNPs);
    // Fitting
    Fit();
    // Store results in maps and take snapshot
    UpdateMapNPsFinal(m_cNPs);
    UpdateMapPOIsFitted(m_cPOIs);
}

void WorkSpace::FitWithFixedNP(const string& sPara, 
                        const map<string, tuple<double, double, double>>& mapNPsFromFitAll,
                        double nMode)
{
    const auto& tupleVal = mapNPsFromFitAll.at(sPara);
    Tools::println("Initial value: % | error_up % | error_down %",
            std::get<0>(tupleVal), std::get<1>(tupleVal), std::get<2>(tupleVal));
    double fFixedVal = std::get<0>(tupleVal);
    fFixedVal += nMode > 0 ? TMath::Abs(nMode) * std::get<1>(tupleVal) : TMath::Abs(nMode) * std::get<2>(tupleVal);
    RooRealVar* cNP = (RooRealVar*)m_cNPs->find(sPara.c_str());
    cNP->setVal(fFixedVal);
    cNP->setConstant(true);
    Fit();
    UpdateMapNPsFinal(m_cNPs);
    UpdateMapPOIsFitted(m_cPOIs);
}

void WorkSpace::FitWithAllNPFixed(const map<string, tuple<double, double, double>>& mapNPsFromFitAll,
                            double nDirection)
{
    for (const auto& pp : mapNPsFromFitAll)
    {
        const auto& tupleVal = pp.second;
        double fFixedVal = std::get<0>(tupleVal);
        double fPos = std::get<1>(tupleVal) > 0 ? std::get<1>(tupleVal) : std::get<2>(tupleVal);
        double fNeg = std::get<1>(tupleVal) < 0 ? std::get<1>(tupleVal) : std::get<2>(tupleVal);
        fFixedVal += nDirection > 0 ? fPos : fNeg;
        RooRealVar* cNP = (RooRealVar*)m_cNPs->find(pp.first.c_str());
        cNP->setVal(fFixedVal);
        cNP->setConstant(true);
    }
    Fit();
    UpdateMapNPsFinal(m_cNPs);
    UpdateMapPOIsFitted(m_cPOIs);
}

void WorkSpace::FitWithFixedMu(double mu_0)
{
    static_cast<RooRealVar*>(m_cPOIs->first())->setVal(mu_0);
    static_cast<RooRealVar*>(m_cPOIs->first())->setConstant(true);
    Fit();
    UpdateMapNPsFinal(m_cNPs);
    UpdateMapPOIsFitted(m_cPOIs);
}

void WorkSpace::SetStatOnly()
{
    const RooArgSet *cNPs = m_cSBModel->GetNuisanceParameters();
    if (cNPs && cNPs->getSize() > 0)
    {
        Tools::println("Switching off NPs by setting them constant to their default values");
        SetAllConstant(*cNPs);
    }
}

map<WorkSpace::ePOI, double> WorkSpace::GetCache(const string& nm)
{
    tuple<double, double, double> tupleVals;
    if (m_mapNPsFitted.find(nm) != m_mapNPsFitted.end())
        tupleVals = m_mapNPsFitted.at(nm);
    else if (m_mapPOIsFitted.find(nm) != m_mapPOIsFitted.end())
        tupleVals = m_mapPOIsFitted.at(nm);

    map<WorkSpace::ePOI, double> res = {
        {WorkSpace::ePOI::VALUE,   std::get<0>(tupleVals)},
        {WorkSpace::ePOI::ERRORUP, std::get<1>(tupleVals)},
        {WorkSpace::ePOI::ERRORDOWN, std::get<2>(tupleVals)}
    };

    return std::move(res);
}

void WorkSpace::DrawProfiledLikelihoodTestStatDist(double fMu, int nToys, int nWorkers) const
{
    gROOT->SetStyle("ATLAS");
    gStyle->SetErrorX(0.5);
    
    Tools::println("Will generate sampling distribution at % = %", m_cPOIs->first()->GetName(), fMu);

    ProfileLikelihoodTestStat* cPLLTestStat = new ProfileLikelihoodTestStat(*m_cSBModel->GetPdf());
    cPLLTestStat->SetOneSided(true); // ATLAS use one-side 
    cPLLTestStat->SetMinimizer("Minuit2");
    cPLLTestStat->SetStrategy(0);
    double fTolerance = m_cInfo->tolerance;
    cPLLTestStat->SetTolerance(fTolerance);

    static_cast<RooRealVar*>(m_cPOIs->first())->setVal(fMu); // set POI value for generation
    
    RooArgSet cPOIs_tmp;
    cPOIs_tmp.add(*m_cSBModel->GetParametersOfInterest());

    /// @todo the evaluate is -2loglambda ??
    double fTestStatData = cPLLTestStat->Evaluate(*m_cData, cPOIs_tmp) * 0.5; 
    Tools::println("Evaluate profiled likelihood on data: %", fTestStatData);

    ToyMCSampler* cToySampler = new ToyMCSampler(*cPLLTestStat, nToys);
    cToySampler->SetPdf(*m_cSBModel->GetPdf());
    cToySampler->SetObservables(*m_cSBModel->GetObservables());
    cToySampler->SetGlobalObservables(*m_cSBModel->GetGlobalObservables());
    cToySampler->SetParametersForTestStat(*m_cSBModel->GetParametersOfInterest()); // set POI value for evaluation

    ProofConfig cProofConfig(*m_cWs, nWorkers, "", false);
    cToySampler->SetProofConfig(&cProofConfig); // enable proof

    RooArgSet* cAllParams = new RooArgSet();
    cAllParams->add(*m_cSBModel->GetParametersOfInterest());
    cAllParams->add(*m_cSBModel->GetNuisanceParameters());

    SamplingDistribution *cSamplingDist = cToySampler->GetSamplingDistribution(*cAllParams);
    SamplingDistPlot cSamplingPlot;
    cSamplingPlot.AddSamplingDistribution(cSamplingDist);
    
    TCanvas *cCanvas = new TCanvas("Sampling", "", 1200, 900);
    cCanvas->SetLogy();
    TH1F* cHistPlot = (TH1F*)cSamplingPlot.GetTH1F(cSamplingDist)->Clone("plotting");

    double xmin = cHistPlot->GetXaxis()->GetXmin();
    double xmax = std::max(cHistPlot->GetXaxis()->GetXmax(), fTestStatData);
    TH1F* cHistDummy = new TH1F("dummy", "", 1, xmin, xmax * 1.2);
    cHistDummy->GetXaxis()->SetLabelSize(0.04);
    cHistDummy->GetXaxis()->SetTitleSize(0.045);
    cHistDummy->GetXaxis()->SetTitleOffset(1.2);
    cHistDummy->GetXaxis()->SetTitle(Form("-log#lambda(#mu=%.2f)", fMu));
    cHistDummy->GetXaxis()->SetRangeUser(xmin, xmax * 1.2);
    cHistDummy->GetYaxis()->SetLabelSize(0.04);
    cHistDummy->GetYaxis()->SetTitleSize(0.045);
    cHistDummy->GetYaxis()->SetTitle("f(-log#lambda(#mu)|#mu)");
    cHistDummy->Draw("AXIS");
    
    cHistPlot->SetLineWidth(1);
    cHistPlot->Draw("HIST E0 SAME");

    TF1 *cChi2 = new TF1("f", "ROOT::Math::chisquared_pdf(2 * x, 1)", xmin, xmax);
    cChi2->SetLineColor(kRed+1);
    cChi2->SetLineWidth(1);
    cChi2->Draw("SAME");

    TLine* cLineObsData = new TLine(fTestStatData, cHistPlot->GetMinimum(), fTestStatData, cHistPlot->GetMaximum());
    cLineObsData->SetLineColor(kBlue+1);
    cLineObsData->SetLineStyle(2);
    cLineObsData->SetLineWidth(2);
    cLineObsData->Draw("SAME");

    cCanvas->SetLogy();
    cCanvas->Update();
    string sOutputPath = \
        "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/pll_test_stat_distribution." + m_cInfo->output_tag + ".pdf";
    cCanvas->SaveAs(sOutputPath.c_str());

    delete cLineObsData;
    delete cChi2;
    delete cAllParams;
    delete cToySampler;
    delete cPLLTestStat;
    delete cCanvas;
}

void WorkSpace::DrawProfiledLogLikelihood()
{
    gROOT->SetStyle("ATLAS");
    gStyle->SetErrorX(0.5);

    // const RooArgSet* cGlobalObs = m_cSBModel->GetGlobalObservables();
    
    // run the fit
    FitAll();

    RooRealVar* cPOI = static_cast<RooRealVar*>(m_cPOIs->first());

    Tools::println("Best fit mu = %, error = [%, %]", cPOI->getVal(), cPOI->getErrorHi(), cPOI->getErrorLo());

    RooPlot* cFrame = cPOI->frame(RooFit::Bins(100), RooFit::Range(-5, 5), RooFit::Title(""));
    m_cNLL->plotOn(cFrame, RooFit::ShiftToZero(), RooFit::Name("deltaNLL"));

    cFrame->SetMinimum(0);
    cFrame->SetMaximum(3);
    cFrame->SetXTitle("#mu");
    cFrame->SetYTitle("-#Delta log(L)");

    Tools::println("Creating profiled likelihood");
    RooAbsReal* cPLL = m_cNLL->createProfile(*m_cPOIs);

    cPLL->plotOn(cFrame, RooFit::LineColor(kRed), RooFit::Precision(-1), RooFit::Name("PLL"));

    TCanvas* cCanvas = new TCanvas("c", "", 1200, 900);
    cFrame->Draw();

    cCanvas->Update();
    string sOutputPath = \
        "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/pll_wrt_mu." + m_cInfo->output_tag + ".pdf";
    cCanvas->SaveAs(sOutputPath.c_str());
}