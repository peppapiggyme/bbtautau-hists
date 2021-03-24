#include "Pulls.h"
#include <iostream>

using namespace std;

int test_pulls(const std::string& filename, const std::string& outname)
{
    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;

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

int test_pulls_plot(const char* in, const char* out)
{
    PullsPlotter* plt = new PullsPlotter();
    plt->LoadFromTxt(in);
    plt->Draw(out);

    delete plt;

    return EXIT_SUCCESS;
}
