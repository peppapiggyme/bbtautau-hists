#include "ExamplesInclude_WS.h"
#include "WorkSpace.h"

#include <set>

using namespace std;

void test_samplingdist(const std::string& filename, const std::string& tag)
{
    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;
    info->fit_func = FitFunction::CST;
    info->tolerance = 1e-6;
    info->output_tag = tag;
    info->output_path = "/scratchfs/atlas/bowenzhang/bbtautau-hists/output/forThesis/SamplingDist/";

    WorkSpace* wst = new WorkSpace(info);
    wst->Check();
    wst->DrawProfiledLikelihoodTestStatDist(4., 100000, 4); // for boosted analysis

    delete info;
    delete wst;
}

