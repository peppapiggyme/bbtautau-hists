#include "ExamplesInclude_WS.h"
#include "WorkSpace.h"

#include <vector>
#include <regex>

static void test_pseudodata_helper(const std::string& filename, const std::string& mass_point, 
    double& obs_q0, double& obs_sig, double& mu_hat_uncon, int& status_con, int& status_uncon, int& sign);

void test_pseudodata(const std::vector<std::string>& filenames, const std::string& outfile)
{
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

    std::string mass_point = "";

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
            mass_point = match[2];
            pseudo_data_idx = std::stoi(match[3]);
        }
        else 
        {
            throw std::runtime_error("no match!");
        }

        test_pseudodata_helper(filename, mass_point, 
            obs_q0, obs_sig, mu_hat_uncon, status_con, status_uncon, sign);
        obs_pvalue = Utils_WS::significanceToPvalue(obs_sig);
        t->Fill();
    }

    TFile* fout = TFile::Open(outfile.c_str(), "RECREATE");
    fout->cd();
    t->Write();
    fout->ls();
    fout->Close();
    delete t;
    fout = nullptr;
}

void test_pseudodata_helper(const std::string& filename, const std::string& mass_point, 
    double& obs_q0, double& obs_sig, double& mu_hat_uncon, int& status_con, int& status_uncon, int& sign)
{
    const int doUncap = 1;

    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->fit_func = FitFunction::FCC;
    info->tolerance = 1e-6;
    info->logLevel = -1;
    info->use_minos = false;
    info->poi_range_low = -40.;
    info->poi_range_high = 40.;
    
    WorkSpace* wst = new WorkSpace(info);

    std::vector<std::pair<std::string, std::string>> tn_and_prefix = {
        {"globs_hadhad", "nom_gamma_stat_Region_BMin0_incJet1_distPNN"+mass_point+"_J2_Y2015_DLLOS_T2_SpcTauHH_L0_bin_"},
        {"globs_ZCR", "nom_gamma_stat_Region_BMin0_incJet1_Y2015_DZllbbCR_T2_L2_distmLL_J2_bin_"},
    };
    wst->SetAlphaGlobalObservablesByFile(
        "/scratchfs/atlas/bowenzhang/ext/bbtt_toy_generation/ws_inputs/toy_globs_"+mass_point+".root", 
        "globs_alphas", 0);
    wst->SetGammaGlobalObservablesByFile(
        "/scratchfs/atlas/bowenzhang/ext/bbtt_toy_generation/ws_inputs/toy_globs_"+mass_point+".root", 
        tn_and_prefix, 0);

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
