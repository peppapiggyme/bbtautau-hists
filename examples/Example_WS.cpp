#include "ExamplesInclude_WS.h"

#include "WorkSpace.h"

void test_ws_info(const std::string& filename) {

    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->use_oneline_fit = false;
    info->use_relative_tolerance = true;
    info->tolerance = 1e-6;

    WorkSpace* wst = new WorkSpace(info);
    wst->Check();
    // wst->FitWithFixedMu(0.);

    // wst->SetStatOnly();
    // wst->FitAll();

    wst->DrawProfiledLogLikelihood();

    wst->Check();

    delete info;
    delete wst;
}
