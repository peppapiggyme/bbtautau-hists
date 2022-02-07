#include "WorkspaceToRoot.h"
#include "Utils.h"

using namespace std;

void test_morphing(const string& filename)
{
    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->fit_func = FitFunction::CST;
    info->use_relative_tolerance = false;
    info->tolerance = 1e-8;

    W2RInfo* w2rInfo = new W2RInfo();
    w2rInfo->verbose = true;
    w2rInfo->prefit = true;
    w2rInfo->postfit_from_file = true;
    w2rInfo->mu = 1.0;

    WorkspaceToRoot* w2r = new WorkspaceToRoot(info, w2rInfo);
    w2r->Execute();
    w2r->DrawMorphing("alpha_SysFFVarrQCD", 21, -2., 2.);

    delete w2r;
    delete info;
    delete w2rInfo; 
}
