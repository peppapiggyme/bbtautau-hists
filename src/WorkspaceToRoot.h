#ifndef WORKSPACETOROOT_H
#define WORKSPACETOROOT_H

#include "StatisticTool.h"
#include "WorkSpace.h"

#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <memory>

using std::ofstream;
using std::vector;
using std::map;
using std::set;
using std::unique_ptr;

class W2RInfo
{
public:
    bool prefit = true;
    bool verbose = true; // show fit results
    /// @todo false -> ad-hoc postfit
    double mu = 1.; // signal scale
    bool postfit_from_file = true; 
    std::string postfit_result_file;
    std::string postfit_result_directory;
};

struct Yield
{
    double val = 0.;
    double err = 0.;
};

class WorkspaceToRoot : public StatisticTool
{
public:
    WorkspaceToRoot(
        const WorkspaceInfo* info,
        const W2RInfo* w2rInfo, 
        const string& name="WorkspaceToRoot");
    
    virtual ~WorkspaceToRoot() override;
    
    WorkspaceToRoot(WorkspaceToRoot& old) = delete;
    
    WorkspaceToRoot& operator=(WorkspaceToRoot& old) = delete;

public:
    virtual void Execute() override;
    void WriteToRootfile(const std::string& sOutName);
    
    /**
     * Function to recognise signal and control region (no signal)
     * @brief only sensitive to bbtautau signals
     */
    bool IsSignal(const string& sCompName);
    bool IsCR(const string& sCateName);

private:
    // m_mapContent is filled with components, RooProdcts and RooAdditions
    void FetchContents();
    void MakeYield();
    Yield GetIntegralAndError(const string& sCategory, RooAbsArg* cComponent) const;
    void GetFitResult();

private:
    const WorkspaceInfo* m_cInfo;
    const W2RInfo* m_cW2RInfo;
    unique_ptr<WorkSpace> m_cWs;
    unique_ptr<RooFitResult> m_cPreFitResult;
    unique_ptr<RooFitResult> m_cPostFitResult;
    /// { category_name : { component_name : RooAbsArg_object, ... }, ... }
    map<string, map<string, RooAbsArg*>> m_mapContent;
    map<string, RooAbsPdf*> m_mapPdf;
    map<string, RooAbsArg*> m_mapBinWidth;
    map<string, RooRealVar*> m_mapObs;
    map<string, map<string, Yield>> m_mapYields;
    map<string, TH1F*> m_mapHistPrefit;
};

#endif