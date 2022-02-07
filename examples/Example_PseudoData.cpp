#include "ExamplesInclude_WS.h"
#include "WorkSpace.h"

#include <vector>

static Float_t test_pseudodata_helper(const std::string& filename);

void test_pseudodata(const std::vector<std::string>& filenames, const std::string& outfile)
{
    TFile* fout = TFile::Open(outfile.c_str(), "RECREATE");
    TTree* t = new TTree("params", "params");

    Float_t muHat;
    t->Branch("muHat", &muHat, "muHat/F");
    for (const auto& filename : filenames)
    {
        muHat = test_pseudodata_helper(filename);
        cout << muHat << endl;
        t->Fill();
    }

    fout->cd();
    t->Write();
}

Float_t test_pseudodata_helper(const std::string& filename)
{
    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->fit_func = FitFunction::CST;

    WorkSpace* wst = new WorkSpace(info);
    // wst->Check();
    wst->FitAll();
    auto pois = wst->GetFittedPOIs();

    delete info;
    delete wst;

    return std::get<0>(pois.begin()->second);
}
