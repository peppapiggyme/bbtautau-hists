#include "ExamplesInclude_WS.h"

#include "WorkSpace.h"

void test_ws_info(const std::string& filename) {

    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->use_oneline_fit = false;

    WorkSpace* wst = new WorkSpace(info);
    wst->Check();
    wst->FitAll();
    // wst->CheckParaMaps();

    info->tolerance = 1e-8;

    WorkSpace* wst_up = new WorkSpace(info);
    wst_up->FitWithAllParaFixed(wst->GetFittedNPs(), 1.0f);

    WorkSpace* wst_down = new WorkSpace(info);
    wst_down->FitWithAllParaFixed(wst->GetFittedNPs(), -1.0f);

    delete info;
    delete wst;
    delete wst_up;
    delete wst_down;
}
