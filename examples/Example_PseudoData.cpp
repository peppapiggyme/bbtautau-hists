#include "ExamplesInclude_WS.h"
#include "WorkSpace.h"

#include <vector>
#include <regex>

static void test_pseudodata_helper(const std::string& filename,
    double& obs_q0, double& obs_sig, double& mu_hat_uncon, int& status_con, int& status_uncon, int& sign);

void test_pseudodata(const std::vector<std::string>& filenames, const std::string& outfile)
{
    TFile* fout = TFile::Open(outfile.c_str(), "RECREATE");
    TTree* t = new TTree("params", "params");
    t->SetDirectory(0);

    double obs_q0 = -1;
    double obs_sig = 0;
    double obs_pvalue = 1;
    double mu_hat_uncon = 0;
    int status_con = -1;
    int status_uncon = -1;
    int sign = 0;
    int pseudo_data_idx = -1;

    t->Branch("obs_q0", &obs_q0, "obs_q0/D");
    t->Branch("obs_sig", &obs_sig, "obs_sig/D");
    t->Branch("obs_pvalue", &obs_pvalue, "obs_pvalue/D");
    t->Branch("status_con", &status_con, "status_con/I");
    t->Branch("status_uncon", &status_uncon, "status_uncon/I");
    t->Branch("mu_hat_uncon", &mu_hat_uncon, "mu_hat_uncon/D");
    t->Branch("sign", &sign, "sign/I");
    t->Branch("pseudo_data_idx", &pseudo_data_idx, "pseudo_data_idx/I");
    
    for (const auto& filename : filenames)
    {
        std::regex re("(.*)/([0-9]+).([0-9]+).root");
        std::smatch match;
        if (std::regex_search(filename, match, re)) 
        {
            pseudo_data_idx = std::stoi(match[3]);
        }
        else 
        {
            throw std::runtime_error("no match!");
        }

        test_pseudodata_helper(filename, obs_q0, obs_sig, mu_hat_uncon, status_con, status_uncon, sign);
        obs_pvalue = Utils_WS::significanceToPvalue(obs_sig);
        t->Fill();
    }

    fout->cd();
    t->Write();
    fout->ls();
    fout->Close();
    delete t;
    fout = nullptr;
}

void test_pseudodata_helper(const std::string& filename, 
    double& obs_q0, double& obs_sig, double& mu_hat_uncon, int& status_con, int& status_uncon, int& sign)
{
    const int doUncap = 1;

    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->fit_func = FitFunction::CST;
    info->tolerance = 1e-6;
    info->logLevel = -1;
    info->use_minos = false;
    info->poi_range_low = -40.;
    info->poi_range_high = 40.;
    
    WorkSpace* wst = new WorkSpace(info);

    // wst->GetRooWorkspace()->loadSnapshot("snapshot_paramsVals_initial");

    ((RooRealVar*)wst->GetRooPOIs()->first())->setVal(0);
    ((RooRealVar*)wst->GetRooPOIs()->first())->setConstant(1);
    wst->FitAll();
    status_con = wst->GetMinimizationStatus();

    double obs_nll_cond = ((RooRealVar*)wst->GetRooNLL())->getVal();

    // wst->GetRooWorkspace()->loadSnapshot("snapshot_paramsVals_initial");
    ((RooRealVar*)wst->GetRooPOIs()->first())->setConstant(0);
    wst->FitAll();
    status_uncon = wst->GetMinimizationStatus();
    mu_hat_uncon = ((RooRealVar*)wst->GetRooPOIs()->first())->getVal();

    double obs_nll_min = ((RooRealVar*)wst->GetRooNLL())->getVal();

    obs_q0 = 2 * (obs_nll_cond - obs_nll_min);
    if (doUncap && ((RooRealVar*)wst->GetRooPOIs()->first())->getVal() < 0) obs_q0 = -obs_q0;

    sign = int(obs_q0 == 0 ? 0 : obs_q0 / fabs(obs_q0));
    if (!doUncap && ((obs_q0 < 0 && obs_q0 > -0.1) || ((RooRealVar*)wst->GetRooPOIs()->first())->getVal() < 0.001))
        obs_sig = 0;
    else
        obs_sig = sign * sqrt(fabs(obs_q0));

    delete info;
    delete wst;
}
