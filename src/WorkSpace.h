#ifndef WorkSpace_H
#define WorkSpace_H

#include <iostream>
#include <chrono>
#include <memory>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <thread>
#include <mutex>
#include <algorithm>
#include <atomic>
#include <fstream>
#include <sstream>
#include <filesystem>

#include "CommonInclude.h"
#include "CommonInclude_WS.h"
#include "Utils_WS.h"
#include "Utils.h"

using std::string;
using std::map;
using std::set;
using std::pair;
using std::tuple;
using std::cout;
using std::clog;
using std::cerr;
using std::endl;
using std::for_each;
using std::thread;
using std::mutex;
using std::lock_guard;
using std::ifstream;
using std::istringstream;

using std::chrono::steady_clock;
using std::chrono::milliseconds;
using std::chrono::duration_cast;

namespace fs = std::filesystem;

// Macros
#define MAP_UPDATE_PARA(map, v) \
    auto v1 = v->getErrorHi(); \
    auto v2 = v->getErrorLo(); \
    /*if (v1 < 0) { auto tmp = v1; v1 = v2; v2 = tmp; }*/ \
    map[v->GetName()] = std::make_tuple(v->getVal(), v1, v2);

#define PRINT_PARA_MAP(map) \
    Tools::println("%:", #map);\
    for (const auto& pp : map) \
    { \
        Tools::println("Name: %  =  %  [ %  % ]", \
        pp.first, std::get<0>(pp.second), std::get<1>(pp.second), std::get<2>(pp.second)); \
    }

#define ITX(f) _IterateAndDo(cParaSet, f);

struct Parameter
{
    Parameter() = default;
    string name;
    double val;
    double err_hi;
    double err_lo;
};

struct WorkspaceInfo
{
    WorkspaceInfo() = default;
    string path;
    string workspace_name;
    string config_name = "ModelConfig";
    string data_name = "obsData";
    double mu_asimov = 1.0;
    double tolerance = 1e-5;
    int8_t logLevel = -1;
    bool use_asimov = true;
    bool use_oneline_fit = true;
    bool use_relative_tolerance = false;
};

class WorkspaceLoad
{
public:
    static WorkspaceLoad& instance() 
    {
        static WorkspaceLoad singleton;
        return singleton;
    }

    WorkspaceLoad(WorkspaceLoad const&) = delete;
    WorkspaceLoad& operator=(WorkspaceLoad const&) = delete;

    void LoadRootFile(const WorkspaceInfo* cInfo)
    {
        if (!fs::exists(fs::path(cInfo->path)))
        {
            throw std::runtime_error("Input does not exist");
        }
        m_cFile = TFile::Open(cInfo->path.c_str(), "READ");
        m_cWs = static_cast<RooWorkspace*>(m_cFile->Get(cInfo->workspace_name.c_str()));
        if (!m_cWs)
        {
            throw std::runtime_error("Workspace does not exist");
        }
    }

    RooWorkspace* GetWorkspace() { return m_cWs; }

private:
    WorkspaceLoad() = default;
    ~WorkspaceLoad() 
    {
        m_cFile->Close();
        m_cFile = nullptr;
        m_cWs = nullptr;
    }

private:
    TFile* m_cFile;
    RooWorkspace* m_cWs;
};

class WorkSpace
{
public:
    WorkSpace(const WorkspaceInfo* info)
        : m_cInfo(info)
    {
        RooMsgService::instance().setGlobalKillBelow(ERROR);
        Load();
    };

    ~WorkSpace()
    {
        // this reduce half memory usage, but why?
        delete m_cWs;
    }

    WorkSpace(WorkSpace& old) = delete;

    WorkSpace& operator=(WorkSpace& old) = delete;

public:
    void Check() {
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

    void CheckParaMaps()
    {
        PRINT_PARA_MAP(m_mapNPsInit);
        PRINT_PARA_MAP(m_mapNPsFitted);
        PRINT_PARA_MAP(m_mapPOIsFitted);
    }

    string NameOfPOI()
    {
        return m_sPOIName;
    }

    enum class ePOI : uint8_t { VALUE, ERRORUP, ERRORDOWN };

private:
    void Load() {
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

protected:
    const WorkspaceInfo* m_cInfo;
    RooWorkspace* m_cWs;
    ModelConfig* m_cSBModel;
    const RooArgSet* m_cNPs;
    const RooArgSet* m_cPOIs;
    string m_sPOIName;
    RooAbsData* m_cData;
    map<string, tuple<double, double, double>> m_mapNPsInit; // Name, Val, Hi, Lo
    map<string, tuple<double, double, double>> m_mapNPsFitted; // Name, Val, Hi, Lo
    map<string, tuple<double, double, double>> m_mapPOIsFitted; // Name, Val, Hi, Lo
    set<string> m_setStrNPs;

// ============================================================================
// Useful Functions
// ============================================================================
private:
    template<typename Func>
    void _IterateAndDo(const RooArgSet*& cParaSet, Func pFunc)
    {
        TIterator* itParas = cParaSet->createIterator();
        RooRealVar* cPara = nullptr;
        while((cPara = (RooRealVar*)itParas->Next()))
            pFunc(cPara);
    }

    template<typename Func, typename ... Targs>
    void _IterateAndDo(const RooArgSet*& cParaSet, Func pFunc, Targs ... args)
    {
        TIterator* itParas = cParaSet->createIterator();
        RooRealVar* cPara = nullptr;
        while((cPara = (RooRealVar*)itParas->Next()))
            pFunc(cPara, args ...);
    }

    void _PrintParameters(const RooArgSet* cParaSet)
    {
        ITX([](const RooRealVar* v){
            Tools::println("Name: %  =  %  [ %  % ]",
                    v->GetName(), v->getVal(), v->getErrorHi(), v->getErrorLo());
        });
    }

    void UpdateMapNPsInit(const RooArgSet* cParaSet)
    {
        ITX([this](const RooRealVar* v){ MAP_UPDATE_PARA(m_mapNPsInit, v); });
    }

    void UpdateMapNPsFinal(const RooArgSet* cParaSet)
    {
        ITX([this](const RooRealVar* v){ MAP_UPDATE_PARA(m_mapNPsFitted, v); });
    }

    void UpdateMapPOIsFitted(const RooArgSet* cParaSet)
    {
        ITX([this](const RooRealVar* v){ MAP_UPDATE_PARA(m_mapPOIsFitted, v); });
    }

    void PrintNuisanceParameters()
    {
        _PrintParameters(m_cNPs);
    }

    void PrintParametersOfInterest()
    {
        _PrintParameters(m_cPOIs);
    }

    void PrintObservables(const RooArgSet* cParaSet)
    {
        ITX([this](const RooAbsArg* v){ v->Print(); });
    }

    void GetSetOfNPs(const RooArgSet* cParaSet)
    {
        ITX([this](const RooRealVar* v){ m_setStrNPs.insert(string(v->GetName())); });
    }

// ============================================================================
// Fitting
// ============================================================================
private:
    RooFitResult* OneLinerFit(RooArgSet& cConstrainParas, RooAbsData& cData)
    {
        std::cout << "Fitting with: \n"
        "o------------------------------------------------------------o\n"
        "| [One liner version]                                        |\n"
        "o------------------------------------------------------------o\n"
        "| This is quicker than the customized version (*)            |\n"
        "| I think it still makes sense to still use this for ranking |\n"
        "| (*) need to make more tests                                |\n"
        "o------------------------------------------------------------o\n";

        RooFitResult* cRes = m_cSBModel->GetPdf()->fitTo(
            cData, InitialHesse(false), Minos(false), Minimizer("Minuit", "Migrad"),
            Strategy(1), PrintLevel(m_cInfo->logLevel), Constrain(cConstrainParas), Save(true),
            Offset(RooStats::IsNLLOffset()));

        return cRes;
    }

    RooFitResult* CustomizedFit(RooArgSet& cConstrainParas, RooAbsData& cData)
    {
        std::cout << "Fitting with: \n"
        "o------------------------------------------------------------o\n"
        "| [Customized version]                                       |\n"
        "o------------------------------------------------------------o\n"
        "| The speed depends on the tolarence                         |\n"
        "| It is set to be 1e-5 x initial NLL value                   |\n"
        "| This might not be optimal                                  |\n"
        "o------------------------------------------------------------o\n";

        RooMsgService::instance().setGlobalKillBelow(ERROR);

        // Define NLL
        RooAbsReal *cNLL = m_cSBModel->GetPdf()->createNLL(cData, 
            Constrain(cConstrainParas), 
            GlobalObservables(*(m_cSBModel->GetGlobalObservables())), 
            Offset(RooStats::IsNLLOffset()), NumCPU(4));

        // Staring NLL value
        Tools::println("Starting NLL value = %", cNLL->getVal());

        // Define Minimizer
        RooMinimizer cMinimizer(*cNLL);
        // Print level
        cMinimizer.setPrintLevel(m_cInfo->logLevel);
        // Strategy
        const int nStrategy = ROOT::Math::MinimizerOptions::DefaultStrategy();
        cMinimizer.setStrategy(nStrategy);
        // double fTolerance = ROOT::Math::MinimizerOptions::DefaultTolerance();
        double fTolerance = m_cInfo->use_relative_tolerance ? m_cInfo->tolerance * cNLL->getVal() : m_cInfo->tolerance;
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
        Tools::println("Final NLL value = %", cNLL->getVal());

        RooFitResult* cRes = cMinimizer.save();

        delete cResPre;

        return cRes;
    }

    void Fit(RooAbsData* cData=nullptr) {
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

public:
    void FitAll() {
        // prefit paramters (not useful)
        UpdateMapNPsInit(m_cNPs);
        // Fitting
        Fit();
        // Store results in maps and take snapshot
        UpdateMapNPsFinal(m_cNPs);
        UpdateMapPOIsFitted(m_cPOIs);
    }

    /**
     * @param nMode
     * 1 corresponds to original up side, while -1 is the down side
     */
    void FitWithFixedNP(const string& sPara, 
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

    /**
     * @param nDirection
     * 1 corresponds to all going positive, while -1 is all going negative
     */
    void FitWithAllNPFixed(const map<string, tuple<double, double, double>>& mapNPsFromFitAll,
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

    /**
     * @param mu_0
     * mu_0 = 0 -> background only fit
     * mu_0 = x -> sig = x hypothesis fit
     */
    void FitWithFixedMu(double mu_0)
    {
        static_cast<RooRealVar*>(m_cPOIs->first())->setVal(mu_0);
        static_cast<RooRealVar*>(m_cPOIs->first())->setConstant(true);
        Fit();
        UpdateMapNPsFinal(m_cNPs);
        UpdateMapPOIsFitted(m_cPOIs);
    }

public:
    map<string, tuple<double, double, double>> GetFittedNPs()
    {
        return m_mapNPsFitted;
    }

    map<string, tuple<double, double, double>> GetFittedPOIs()
    {
        return m_mapPOIsFitted;
    }

    set<string> GetNPs()
    {
        return m_setStrNPs;
    }

    map<WorkSpace::ePOI, double> GetCache(const string& nm)
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

public:
    bool bFitted = false;

};

#endif