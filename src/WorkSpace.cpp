#include "WorkSpace.h"
#include "Minimization.h"

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

    m_nStatus = cRes->status();
    m_fNLL = m_cNLL->getVal();
    
    return cRes;
}

RooFitResult* WorkSpace::CustomizedFit(RooArgSet& cConstrainParas, RooAbsData& cData)
{
    std::cout << "Fitting with: \n"
    "o------------------------------------------------------------o\n"
    "| [Customized version]                                       |\n"
    "o------------------------------------------------------------o\n"
    "| The speed depends on the tolarence                         |\n"
    "| This one is not stable                                     |\n"
    "o------------------------------------------------------------o\n";

    RooMsgService::instance().setGlobalKillBelow(ERROR);

    // Define NLL
    m_cNLL = m_cSBModel->GetPdf()->createNLL(cData, 
        Constrain(cConstrainParas), 
        GlobalObservables(*(m_cSBModel->GetGlobalObservables())), 
        Offset(true), NumCPU(4));

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

    Tools::println("Tolerance = [%]", fTolerance);
    Tools::println("Max Calls = [%], Max Iterations = [%]", nMaxCalls, nMaxIterations);

    // Do minimization
    nStatus = cMinimizer.minimize("Minuit2", sAlgorithm.c_str());
    Tools::println("Minimize status = [%]", nStatus);
    m_nStatus = nStatus;

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
        if (m_cInfo->use_minos)
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
    m_fNLL = m_cNLL->getVal();

    delete cResPre;

    return cRes;
}

RooFitResult* WorkSpace::CommonStatToolFit(RooArgSet& cConstrainParas, RooAbsData& cData)
{
    std::cout << "Fitting with: \n"
    "o------------------------------------------------------------o\n"
    "| [CommonStatTool version]                                   |\n"
    "o------------------------------------------------------------o\n"
    "| Stable                                                     |\n"
    "| Only MINOS confidence interval                             |\n"
    "o------------------------------------------------------------o\n";
    m_cNLL = EXOSTATS::createNLL(m_cSBModel->GetPdf(), &cData, &cConstrainParas, 4);

    RooFitResult *cRes;
    EXOSTATS::minimize(m_cNLL, 3, nullptr, "", "", 0, kTRUE, &cRes, m_cInfo->use_minos);

    m_nStatus = cRes->status();
    m_fNLL = m_cNLL->getVal();

    return cRes;
}

/// @todo need some clean up ..
RooFitResult* WorkSpace::FCCFit(RooArgSet& cConstrainParas, RooAbsData& cData)
{
    std::cout << "Fitting with: \n"
    "o------------------------------------------------------------o\n"
    "| [FitCrossCheck version]                                    |\n"
    "o------------------------------------------------------------o\n"
    "| Stable                                                     |\n"
    "| Default is HESSE, can choose MINOS                         |\n"
    "o------------------------------------------------------------o\n";
    
    bool fancy(false);
    bool retryOnHesseFailure(false); // a bit less fancy

    // RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

    Constrain(cConstrainParas);

    const RooArgSet* glbObs = m_cSBModel->GetGlobalObservables();

    // fix poi to zero
    RooRealVar * poi = (RooRealVar*) m_cSBModel->GetParametersOfInterest()->first();
    if(poi) {
        cout << "Constant POI ";
    if(poi->isConstant()) { cout << "YES" << endl; } else { cout << "NO" << endl; }
        cout << "Value of POI  " << poi->getVal() << endl;
    } else {
        cout << "No POI - assuming background only fit" << endl;
    }
    int nCPU = 2;
        const char* NCORE = getenv("NCORE");
    if (NCORE) {
        TString nCPU_str = NCORE;
        nCPU = nCPU_str.Atoi();
    }
    RooAbsReal * m_cNLL = m_cSBModel->GetPdf()->createNLL(cData, Constrain(cConstrainParas), GlobalObservables(*glbObs), Offset(1), Optimize(2) );
    double nllval = m_cNLL->getVal();
    cout << "Starting NLL value " << nllval << endl;

    double nllnom(0);
    double nllup(0);
    double nlldn(0);

    cout << "Change to the LHood" << endl;
    TIterator *nItr = m_cSBModel->GetNuisanceParameters()->createIterator();
    RooRealVar* arg = 0;
    while ((arg=(RooRealVar*)nItr->Next())) {
        if (!arg) continue;
        if(TString(arg->GetName()).Contains("ATLAS_norm")) { arg->setVal(1); continue; }
        if(TString(arg->GetName()).Contains("gamma_stat")) { arg->setVal(1); continue; }
        arg->setVal(5);
        nllup = m_cNLL->getVal();
        arg->setVal(-5);
        nlldn = m_cNLL->getVal();
        arg->setVal(0); // must re-set before fitting!!
        nllnom = m_cNLL->getVal();
        cout << "\t" << arg->GetName() << "\t" << nllup-nllnom << "\t" << nlldn-nllnom << endl;
    }

    std::cout << "initial parameters" << std::endl;
    cConstrainParas.Print("v");

    std::cout << "INITIAL NLL = " << m_cNLL->getVal() << std::endl;

    RooMinimizer minim(*m_cNLL);
    int strategy = ROOT::Math::MinimizerOptions::DefaultStrategy();
    minim.setStrategy( strategy );

    Tools::println(">>>>>>>> DEBUG <<<<<<<< Line [%] Strategy [%]", __LINE__, strategy);

    double tol =  ROOT::Math::MinimizerOptions::DefaultTolerance(); //currently 0.1
    tol = std::min(tol,0.1);
    minim.setEps( 0.01 );

    Tools::println(">>>>>>>> DEBUG <<<<<<<< Line [%] Tolerant [%]", __LINE__, tol);

    int status = -1;

    int optConstFlag = 2;
    minim.optimizeConst(optConstFlag);

    TStopwatch sw; sw.Start();

    TString minimizer = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
    minimizer = "Minuit2";
    TString algorithm = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();

    Tools::println(">>>>>>>> DEBUG <<<<<<<< Line [%] MINIMIZER Type [%]", __LINE__, minimizer);
    Tools::println(">>>>>>>> DEBUG <<<<<<<< Line [%] MINIMIZER Algo [%]", __LINE__, algorithm);

    int maxCalls = ROOT::Math::MinimizerOptions::DefaultMaxFunctionCalls();
    int maxIterations = ROOT::Math::MinimizerOptions::DefaultMaxIterations();

    Tools::println(">>>>>>>> DEBUG <<<<<<<< Line [%] max calls [%]", __LINE__, maxCalls);
    Tools::println(">>>>>>>> DEBUG <<<<<<<< Line [%] max iters [%]", __LINE__, maxIterations);

    minim.setMaxFunctionCalls(maxCalls);
    minim.setMaxIterations(maxIterations);

    cout << "MAX CALLS " << maxCalls  << endl;
    cout << "MAX Iterations " << maxIterations << endl;

    cout << "FitPDF" << endl;
    cout << "\t minimizer " << minimizer << endl;
    cout << "\t algorithm " << algorithm << endl;
    cout << "\t strategy  " << strategy << endl;
    cout << "\t tolerance " << tol << endl;

    // HistFitter
    if(fancy) {
        Tools::println(">>>>>>>> DEBUG <<<<<<<< Line [%] in FANCY", __LINE__);
        bool kickApplied(false);
        for (int tries = 1, maxtries = 4; tries <= maxtries; ++tries) {
            cout << "try " << tries << " " << kickApplied << endl;
            sw.Print();
            std::cout << "prefit values" << std::endl;
            cConstrainParas.Print("v");
            status = minim.minimize(minimizer, algorithm);
            cout << "minimizer status = " << status << endl;
            if (status%1000 == 0) {  // ignore erros from Improve
            break;
            } else {
                if (tries == 1) {
                    cout << "    ----> Doing a re-scan first" << endl;
                    status = minim.minimize(minimizer,"Scan");
                }
                if (tries == 2) {
                    if (ROOT::Math::MinimizerOptions::DefaultStrategy() == 1 ) {
                        cout << "    ----> trying with strategy = 2" << endl;
                        minim.setStrategy(2);
                    }
                    else { tries++; } // skip this trial if strategy
                }
                if (tries == 3) {
                    cout << "    ----> trying with improve" << endl;
                    minimizer = "Minuit2";
                    algorithm = "migradimproved";
                }
                if (tries == 4 && !kickApplied) {
                    cout << "    ----> trying fit with different starting values" << endl;
                    RooFitResult* tmpResult = minim.save();
                    const RooArgList& randList = tmpResult->randomizePars();
                    cConstrainParas = randList;
                    delete tmpResult;
                    tries=0;          // reset the fit cycle
                    kickApplied=true; // do kick only once
                    minim.setStrategy(ROOT::Math::MinimizerOptions::DefaultStrategy());
                }
            }
        }
    } else if (retryOnHesseFailure) {
        Tools::println(">>>>>>>> DEBUG <<<<<<<< Line [%] in RETRY", __LINE__);
        int maxtries = 3;
        int tries = 1;
        while (true) {
            cout << "Starting minization, try number " << tries << " of " << maxtries << endl;
            status = minim.minimize(minimizer, algorithm);
            cout << "Calling Hesse ..."  << endl;
            minim.setMinimizerType("Minuit2");
            int statusH = minim.hesse();
            cout << "Hesse Status : " << statusH << endl;
            if (statusH == 0 || tries >= maxtries) {
                break;
            }
            cout << "    ----> trying fit with different starting values" << endl;
            RooFitResult* tmpResult = minim.save();
            const RooArgList& randList = tmpResult->randomizePars();
            cConstrainParas = randList;
            delete tmpResult;
            tries++;
        }
    } else {
        Tools::println(">>>>>>>> DEBUG <<<<<<<< Line [%] in PLAIN", __LINE__);
        status = minim.minimize(minimizer, algorithm);
    }
    cout << "Minimize Status : " << status << endl;
    m_nStatus = status;

    //cout << endl;
    RooFitResult * tmpResult = minim.save();
    const TMatrixDSym covarMat = tmpResult->covarianceMatrix();
    Double_t det = covarMat.Determinant();
    cout << "Determinant " << det << endl;
    if(det < 0) { cout << "Determinant negative" << endl; }

    // get eigenvectors and eigenvalues
    TMatrixDSymEigen eigenValueMaker(covarMat);
    TVectorT<double> eigenValues   = eigenValueMaker.GetEigenValues();
    TMatrixT<double> eigenVectors  = eigenValueMaker.GetEigenVectors();
    cout << endl << "Eigenvalues  " << endl;
    for( int l=0; l<eigenValues.GetNrows(); l++ ) {
        cout << "\t" << l << "\t" << eigenValues[l] << endl;
    }
    cout << endl;

    if (status%100 == 0) { // ignore errors in Hesse or in Improve
        if (!retryOnHesseFailure) {
            cout << "Calling Hesse ..."  << endl;
            minim.setMinimizerType("Minuit2");
            status = minim.hesse();
            cout << "Hesse Status : " << status << endl;
        }
        if (m_cInfo->use_minos) minim.minos();
    } else {
        cout << "FIT FAILED !" << endl;
    }

    sw.Print();

    sw.Stop();

    RooFitResult * r = minim.save();

    if(poi) {
        std::cout << "final poi parameters" << std::endl;
        m_cSBModel->GetParametersOfInterest()->Print("v");
    }

    typedef std::numeric_limits< double > dbl;
    cout.precision(dbl::digits10);
    std::cout << "FINAL NLL = " << m_cNLL->getVal() << std::endl;
    m_fNLL = m_cNLL->getVal();

    if(poi) {
        cout << "PRINTING FIT RESULT " << poi->isConstant() << "\t" << poi->getVal() << endl;
    } else {
        cout << "PRINTING FIT RESULT " << endl;
    }
    r->Print();

    sw.Print();

    return r;
}

void WorkSpace::Fit(RooAbsData* cData) {
    if (!cData)
    {
        cData = m_cData;
    }
    // Information from workspace
    RooRealVar* cPOI = static_cast<RooRealVar*>(m_cPOIs->first());
    cPOI->setRange(m_cInfo->poi_range_low, m_cInfo->poi_range_high);
    auto sMinimizerType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
    Tools::println("POI [%] initial value is [% % %]", 
        cPOI->GetName(), cPOI->getVal(), cPOI->getErrorHi(), cPOI->getErrorLo());
    Tools::println("Set [%] as minimizer", sMinimizerType);
    // Remove constant paras
    RooArgSet cConstrainParas;
    cConstrainParas.add(*m_cNPs);
    RooStats::RemoveConstantParameters(&cConstrainParas);
    auto timeStart = steady_clock::now();

    RooFitResult* cRes = nullptr;
    // unique_ptr<RooFitResult> cRes = nullptr;

    // >>> core of fitting <<< START
    switch (m_cInfo->fit_func)
    {
    case FitFunction::ONELINE:
        cRes = OneLinerFit(cConstrainParas, *cData);
        break;
    case FitFunction::CUSTOM:
        cRes = CustomizedFit(cConstrainParas, *cData);
        break;
    case FitFunction::CST:
        cRes = CommonStatToolFit(cConstrainParas, *cData);
        break;
    case FitFunction::FCC:
        cRes = FCCFit(cConstrainParas, *cData);
        break;
    default:
        cRes = CommonStatToolFit(cConstrainParas, *cData);
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

    // if (cRes)
    // {
    //     delete cRes;
    // }
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
    SetConstantNP(sPara, fFixedVal);
    Fit();
    UpdateMapNPsFinal(m_cNPs);
    UpdateMapPOIsFitted(m_cPOIs);
}

void WorkSpace::FitWithFixedNP(const string& sPara, double fValue)
{
    SetConstantNP(sPara, fValue);
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
    SetConstantPOI(mu_0);
    Fit();
    UpdateMapNPsFinal(m_cNPs);
    UpdateMapPOIsFitted(m_cPOIs);
}

void WorkSpace::SetConstantPOI(double mu)
{
    static_cast<RooRealVar*>(m_cPOIs->first())->setVal(mu);
    static_cast<RooRealVar*>(m_cPOIs->first())->setConstant(true);
}

void WorkSpace::SetConstantNP(const string& sPara, double fValue)
{
    RooRealVar* cNP = (RooRealVar*)m_cNPs->find(sPara.c_str());
    cNP->setVal(fValue);
    cNP->setConstant(true);
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
    cPLLTestStat->SetStrategy(1);
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

    // ProofConfig cProofConfig(*m_cWs, nWorkers, "", false);
    // cToySampler->SetProofConfig(&cProofConfig); // enable proof

    RooArgSet* cAllParams = new RooArgSet();
    cAllParams->add(*m_cSBModel->GetParametersOfInterest());
    cAllParams->add(*m_cSBModel->GetNuisanceParameters());

    SamplingDistribution *cSamplingDist = cToySampler->GetSamplingDistribution(*cAllParams);
    SamplingDistPlot cSamplingPlot;
    cSamplingPlot.AddSamplingDistribution(cSamplingDist);
    
    TCanvas *cCanvas = new TCanvas("Sampling", "", 1200, 900);
    cCanvas->SetLogy();
    TH1F* cHistPlot = (TH1F*)cSamplingPlot.GetTH1F(cSamplingDist)->Clone("plotting");

    double xmin = std::max(cHistPlot->GetXaxis()->GetXmin(), 0.);
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
    
    cHistPlot->SetLineWidth(2);
    cHistPlot->SetMarkerStyle(kFullCircle);
    cHistPlot->Draw("E0 SAME");

    TF1 *cChi2 = new TF1("f", "2 * ROOT::Math::chisquared_pdf(2 * x, 1)", xmin, xmax);
    cChi2->SetLineColor(kRed+1);
    cChi2->SetLineWidth(1);
    // cChi2->Draw("SAME");

    TH1F* cHistHalfChi2 = (TH1F*)cHistPlot->Clone("half_chi2");
    for (int i = 0; i < cHistHalfChi2->GetNbinsX() + 2; ++i)
    {
        cHistHalfChi2->SetBinContent(i, 0);
        cHistHalfChi2->SetBinError(i, 0);
    }

    for (int i = 0; i < nToys; ++i)
    {
        int binary = rand() % 2;
        if (binary == 0)
        {
            cHistHalfChi2->Fill(0);
        }
        else
        {
            cHistHalfChi2->Fill(cChi2->GetRandom());
        }
    }
    double fSampleDist = cHistPlot->Integral();
    Tools::println("Integral of sampling distribution: %", fSampleDist);
    cHistHalfChi2->Scale(1. / (cHistHalfChi2->Integral() * cHistHalfChi2->GetBinWidth(1)));

    cHistHalfChi2->SetLineColor(kRed+1);
    cHistHalfChi2->SetMarkerColor(kRed+1);
    cHistHalfChi2->SetMarkerStyle(kOpenCircle);
    cHistHalfChi2->SetLineWidth(2);
    cHistHalfChi2->Draw("E0 SAME");

    // TLine* cLineObsData = new TLine(fTestStatData, cHistPlot->GetMinimum(), fTestStatData, cHistPlot->GetMaximum());
    // cLineObsData->SetLineColor(kBlue+1);
    // cLineObsData->SetLineStyle(2);
    // cLineObsData->SetLineWidth(2);
    // cLineObsData->Draw("SAME");

    double p = cHistHalfChi2->KolmogorovTest(cHistPlot);
    char p_s[64];
    sprintf(p_s, "p_{K-S test} = %.3f", p);
    Tools::println("K-S test result: %", p);

    TLegend* legend = new TLegend(0.40, 0.78, 0.90, 0.92);
    legend->SetNColumns(2);
    legend->SetTextFont(42);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.036);
    legend->SetTextAlign(12);
    legend->AddEntry(cHistPlot, "Pseudo experiment", "lep");
    legend->AddEntry(cHistHalfChi2, "Approximation", "lep");
    legend->AddEntry("", p_s, "");
    legend->Draw();

    cCanvas->SetLogy();
    cCanvas->Update();
    string sOutputPath = \
        "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/pll_test_stat_distribution." + m_cInfo->output_tag + ".pdf";
    cCanvas->SaveAs(sOutputPath.c_str());

    // delete cLineObsData;
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
