#include "WorkspaceToRoot.cpp"

using namespace std;

void test_w2r(const string& filename)
{
    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->use_oneline_fit = false;
    info->use_relative_tolerance = false;
    info->tolerance = 1e-8;

    W2RInfo* w2rInfo = new W2RInfo();
    w2rInfo->make_yield = true;
    w2rInfo->prefit = true;

    WorkspaceToRoot* w2r = new WorkspaceToRoot(info, w2rInfo);
    w2r->Execute();
    w2r->WriteToRootfile("test.root");

    delete w2r;
    delete info;
    delete w2rInfo; 
}
