#include "ExamplesInclude_WS.h"

#include "PullsTool.h"
#include <iostream>

using namespace std;

int test_pulls(const std::string& filename, const std::string& outname)
{
    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    // info->workspace_name = "combWS";
    // info->data_name = "asimovData_SB_SM";
    info->use_asimov = false;
    info->use_oneline_fit = false;

    auto timeStart = steady_clock::now();

    PullsEngine* egn = new PullsEngine(info);
    egn->Execute();
    egn->WriteToTxt(outname);

    delete info;
    delete egn;

    auto timeEnd = steady_clock::now();
    Tools::println("Spent [%ms]", duration_cast<milliseconds>(timeEnd-timeStart).count());

    return EXIT_SUCCESS;
}

int test_pulls_plot(const std::string& in, const std::string& out)
{
    PullsPlotter* plt = new PullsPlotter();
    plt->LoadFromTxt(in.c_str());
    plt->Draw(out.c_str());

    delete plt;

    return EXIT_SUCCESS;
}
