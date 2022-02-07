#include "ExamplesInclude_WS.h"

#include "WorkSpace.h"

using namespace std;

void test_samplingdist(const std::string& filename)
{
    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->fit_func = FitFunction::CST;
    info->tolerance = 1e-8;
    info->output_tag = "v13atMu0p4";

    WorkSpace* wst = new WorkSpace(info);
    wst->SetStatOnly();
    wst->Check();

    wst->DrawProfiledLikelihoodTestStatDist(0.4, 10000);

    delete info;
    delete wst;
}

