#include "ExamplesInclude_WS.h"

#include "WorkSpace.h"

void test_ws_info(const std::string& filename) {

    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->fit_func = FitFunction::CST;
    info->use_minos = false;
    info->use_relative_tolerance = false;
    info->tolerance = 1e-3;

    WorkSpace* wst = new WorkSpace(info);
    wst->SetConstantPOI(0.);
    wst->Check();
    wst->FitAll();
    double bestfitVal = std::get<0>(wst->GetFittedNPs().at("alpha_SysFFVarrQCD"));
    double minNLL = wst->GetRooNLL()->getVal();

    Tools::println("Best fit value for the NP is [%]", bestfitVal);

    TFile* fout = TFile::Open("llh_scan.root", "RECREATE");
    TTree* t = new TTree("params", "params");
    t->SetDirectory(0);

    double alpha_SysFFVarrQCD;
    double PLL;
    double LL;
    int status;
    
    t->Branch("alpha_SysFFVarrQCD", &alpha_SysFFVarrQCD, "alpha_SysFFVarrQCD/D");
    t->Branch("PLL", &PLL, "PLL/D");
    t->Branch("LL", &LL, "LL/D");
    t->Branch("status", &status, "status/I");
    
    // vector<double> scan{
    //     -2., -1.75, -1.5, -1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.15, -0.1, -0.05, 
    //     0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.5, 1.75, 2.,
    // };

    vector<double> scan{
        -2., -1.5, -1., -0.5, -0.3, -0.2, -0.1, -0.05, 
        0, 0.05, 0.0814763, 0.1, 0.2, 0.3, 0.5, 1., 1.5, 2.,
    };

    for (auto val : scan)
    {
        wst->SetConstantNP("alpha_SysFFVarrQCD", val);
        LL = wst->GetRooNLL()->getVal() - minNLL;

        WorkSpace* wst_scan = new WorkSpace(info);
        wst_scan->SetConstantPOI(0.);
        wst_scan->FitWithFixedNP("alpha_SysFFVarrQCD", val);
        
        alpha_SysFFVarrQCD = val;
        PLL = wst_scan->GetRooNLL()->getVal() - minNLL;
        status = wst_scan->GetMinimizationStatus();
        t->Fill();
        
        delete wst_scan;
    }

    fout->cd();
    t->Write();

    // wst->Check();
    delete t;
    delete info;
    delete wst;
}
