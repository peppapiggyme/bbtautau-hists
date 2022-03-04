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

enum class FitFunction
{
    ONELINE = 0x0,
    CUSTOM = 0x1,
    CST = 0x2,
    FCC = 0x4,
    NONE = 0xFF
};

struct WorkspaceInfo
{
    WorkspaceInfo() = default;
    string path;
    string workspace_name;
    string config_name = "ModelConfig";
    string data_name = "obsData";
    string output_tag = "TAG";
    string output_path = "";
    double mu_asimov = 1.0;
    double tolerance = 1e-5;
    FitFunction fit_func = FitFunction::CST;
    int8_t logLevel = -1;
    bool use_asimov = true;
    bool use_minos = true;
    bool use_relative_tolerance = false;
    double poi_range_low = -40.;
    double poi_range_high = 40.;
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

    inline RooWorkspace* GetWorkspace() { return m_cWs; }

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
    void Check();

    void CheckParaMaps();

    inline string NameOfPOI() { return m_sPOIName; }

    enum class ePOI : uint8_t { VALUE, ERRORUP, ERRORDOWN };

private:
    void Load();

protected:
    const WorkspaceInfo* m_cInfo;
    RooWorkspace* m_cWs;
    ModelConfig* m_cSBModel;
    const RooArgSet* m_cNPs;
    const RooArgSet* m_cPOIs;
    string m_sPOIName;
    RooAbsData* m_cData;
    RooAbsReal* m_cNLL;
    map<string, tuple<double, double, double>> m_mapNPsInit; // Name, Val, Hi, Lo
    map<string, tuple<double, double, double>> m_mapNPsFitted; // Name, Val, Hi, Lo
    map<string, tuple<double, double, double>> m_mapPOIsFitted; // Name, Val, Hi, Lo
    set<string> m_setStrNPs;
    double m_fNLL;

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
    /// title bait, can you fit bbtautau with one line of code 
    RooFitResult* OneLinerFit(RooArgSet& cConstrainParas, RooAbsData& cData);

    /// still need to be improved, not stable
    RooFitResult* CustomizedFit(RooArgSet& cConstrainParas, RooAbsData& cData);

    /// from Exotic CommomStatTool
    RooFitResult* CommonStatToolFit(RooArgSet& cConstrainParas, RooAbsData& cData);
    
    /// from WSMaker FitCrossCheckForLimits
    RooFitResult* FCCFit(RooArgSet& cConstrainParas, RooAbsData& cData);

    void Fit(RooAbsData* cData=nullptr);

public:
    void FitAll();

    /**
     * @param nMode
     * 1 corresponds to original up side, while -1 is the down side
     */
    void FitWithFixedNP(
        const string& sPara, 
        const map<string, tuple<double, double, double>>& mapNPsFromFitAll,
        double nMode);

    void FitWithFixedNP(
        const string& sPara, 
        double fValue);

    /**
     * @param nDirection
     * 1 corresponds to all going positive, while -1 is all going negative
     */
    void FitWithAllNPFixed(
        const map<string, tuple<double, double, double>>& mapNPsFromFitAll, 
        double nDirection);

    /**
     * @param mu_0
     * mu_0 = 0 -> background only fit
     * mu_0 = x -> sig = x hypothesis fit
     */
    void FitWithFixedMu(double mu_0);

    void SetConstantPOI(double mu);
    void SetConstantNP(const string& sPara, double fValue);
    void SetStatOnly();
    void SetFloatOnly(const set<string>& setNormFactorStr);

public:
    map<WorkSpace::ePOI, double> GetCache(const string& nm);

    /// Getters
    inline map<string, tuple<double, double, double>> GetFittedNPs() { return m_mapNPsFitted; }

    inline map<string, tuple<double, double, double>> GetFittedPOIs() { return m_mapPOIsFitted; }

    inline set<string> GetNPs() { return m_setStrNPs; }

    inline double GetNLL() const { return m_fNLL; }

    inline RooWorkspace* GetRooWorkspace() { return m_cWs; }  // cannot be const limited by RooFit

    inline const ModelConfig* GetModelConfig() const { return m_cSBModel; }

    inline RooAbsData* GetRooData() { return m_cData; } // cannot be const limited by RooFit

    inline const RooArgSet* GetRooNPs() const { return m_cNPs; }

    inline const RooArgSet* GetRooPOIs() const { return m_cPOIs; }

    inline const RooAbsReal* GetRooNLL() const { return m_cNLL; }

    inline int GetMinimizationStatus() const { return m_nStatus; }

public:
    bool bFitted = false;
    int m_nStatus = -1;


// High level usage that can go into stattools
public:
    void DrawProfiledLikelihoodTestStatDist(double fMu, int nToys, int nWorkers=4) const;

    void DrawProfiledLogLikelihood();
};

#endif