#include "Ranking.h"
#include <iostream>

using namespace std;

int test_ranking(const std::string& filename, const std::string& outname)
{
    WorkspaceInfo* info = new WorkspaceInfo();
    info->path = filename;
    info->workspace_name = "combined";
    info->use_asimov = false;

    auto timeStart = steady_clock::now();

    RankingEngine* egn = new RankingEngine(info);
    egn->Execute();
    egn->WriteToTxt(outname);

    delete info;
    delete egn;

    auto timeEnd = steady_clock::now();
    Tools::println("Spent [%ms]", duration_cast<milliseconds>(timeEnd-timeStart).count());

    return EXIT_SUCCESS;
}

int test_ranking_plot(const char* in, const char* out)
{
    RankingPlotter* plt = new RankingPlotter();
    plt->LoadFromTxt(in);
    plt->Draw(out);

    delete plt;

    return EXIT_SUCCESS;
}
