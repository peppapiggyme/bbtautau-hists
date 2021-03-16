#include <iostream>
#include <fstream>
#include <iomanip>
#include "RankingTool.h"

using namespace std;

#define OUTPUT_TABLE_4(A, B, C, D) \
    std::left << setw(60) << A \
    << std::left << setw(15) << B \
    << std::left << setw(15) << C \
    << std::left << setw(15) << D << endl;

int test_ranking()
{
    ofstream output("/tmp/zhangbw/Ranking.txt");

    Engine* cEgn = new Engine();
    cEgn->RunRanking();

    auto val = cEgn->GetFittedPOI();
    auto res = cEgn->GetFittedPOIWithFixedNPs();

    output << OUTPUT_TABLE_4("FIXED NP", "POI_VALUE", "POI_ERROR_UP", "POI_ERROR_DOWN");

    Tools::println("% ->      VALUE %      ERROR_UP %      ERROR_DOWN %",
            "No fixed NP",
            val[FitResult::ePOI::VALUE],
            val[FitResult::ePOI::ERRORUP],
            val[FitResult::ePOI::ERRORDOWN]);
    output << OUTPUT_TABLE_4("None",
            val[FitResult::ePOI::VALUE],
            val[FitResult::ePOI::ERRORUP],
            val[FitResult::ePOI::ERRORDOWN])

    for (auto& p : res) {
        Tools::println("% Fixed ->      VALUE %      ERROR_UP %      ERROR_DOWN %",
                p.first,
                p.second[FitResult::ePOI::VALUE],
                p.second[FitResult::ePOI::ERRORUP],
                p.second[FitResult::ePOI::ERRORDOWN]);
        output << OUTPUT_TABLE_4(p.first,
                p.second[FitResult::ePOI::VALUE],
                p.second[FitResult::ePOI::ERRORUP],
                p.second[FitResult::ePOI::ERRORDOWN]);
    }

    delete cEgn;
    
    return EXIT_SUCCESS;
}
