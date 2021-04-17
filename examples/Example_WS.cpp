#include "ExamplesInclude_WS.h"

#include "WorkSpace.h"

void test_ws_info(const std::string& filename) {

    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;

    WorkSpace* wst = new WorkSpace(info);
    wst->Check();
    wst->FitAll();
    wst->CheckParaMaps();

    delete info;
    delete wst;
}
