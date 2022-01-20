#include "WorkspaceToRoot.h"
#include "Utils.h"

using namespace std;

void test_w2r(const string& filename, const string& postfit_result_file)
{
    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->use_oneline_fit = false;
    info->use_relative_tolerance = false;
    info->tolerance = 1e-8;

    W2RInfo* w2rInfo = new W2RInfo();
    w2rInfo->verbose = true;
    // w2rInfo->prefit = true;
    w2rInfo->prefit = false; 
    w2rInfo->postfit_from_file = true;
    w2rInfo->mu = 1.0;
    w2rInfo->postfit_result_file = postfit_result_file;
    w2rInfo->postfit_result_directory = Utils::bkgOnlyFitResultString();

    WorkspaceToRoot* w2r = new WorkspaceToRoot(info, w2rInfo);
    w2r->Execute();
    w2r->WriteToRootfile("test.root");

    delete w2r;
    delete info;
    delete w2rInfo; 
}
