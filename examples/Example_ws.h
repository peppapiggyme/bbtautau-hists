#include "WorkspaceTool.h"

void test_ws_info(const std::string& filename) {
    WorkspaceTool *wst = new WorkspaceTool();
    wst->load(filename, "combined", "ModelConfig", "obsData");





    delete wst;
}
