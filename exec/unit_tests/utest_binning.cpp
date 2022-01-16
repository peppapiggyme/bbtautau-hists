#define CATCH_CONFIG_MAIN

#include "catch.h"
#include "Utils.h"

using namespace std;
using BU = BinningUtils;

vector<double> getBinningFromFileHadHad(const std::string& fn, bool BDT=false)
{ 
    return BU::intToDoubleBinEdgesForMVAInverse(BU::readBinningFromFile<int>(fn), 1000, BDT);
}

vector<double> getBinningFromFileLepHad(const std::string& fn) 
{ 
    return BU::intToDoubleBinEdgesForMVAInverse_LepHad1090Version(BU::readBinningFromFile<int>(fn));
}

bool equal(double a, double b, bool debug=false)
{
    if (debug)
        cout << "comparing " << (int)(a*10000) << " with " << (int)(b*10000) << endl;
    return (int)(a*10000) == (int)(b*10000);
}

TEST_CASE("hadhad pnn")
{
    auto binning = getBinningFromFileHadHad( // HadHad
        "/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_2HDM1000.txt"
    );

    REQUIRE(binning.size() == 4);
    REQUIRE(equal(binning[0], 0.));
    REQUIRE(equal(binning[1], 0.877));
    REQUIRE(equal(binning[2], 0.993));
    REQUIRE(equal(binning[3], 1.));
}

TEST_CASE("hadhad bdt")
{
    auto binning = getBinningFromFileHadHad( // HadHad
        "/scratchfs/atlas/bowenzhang/bbtautau-hists/data/Binning_Trafo14_SMBDT.txt", 
        true  // BDT = true
    );

    REQUIRE(binning.size() == 15);
    REQUIRE(equal(binning[0], -1.));
    REQUIRE(equal(binning[2], -0.954));
    REQUIRE(equal(binning[9], 0.75));
    REQUIRE(equal(binning[14], 1.));
}

TEST_CASE("lephad nns")
{
    auto binning = getBinningFromFileLepHad( // LepHad
        "/scratchfs/atlas/bowenzhang/bbtautau-hists/data/v9/Binning_LepHad_SLT_2HDM1100.v9.txt"
    );

    REQUIRE(binning.size() == 5);
    REQUIRE(equal(binning[0], 0.));
    REQUIRE(equal(binning[1], 0.001));
    REQUIRE(equal(binning[3], 0.9996));
    REQUIRE(equal(binning[4], 1.));
}

TEST_CASE("lephad special cases")
{
    auto binning = BU::intToDoubleBinEdgesForMVAInverse_LepHad1090Version( // LepHad
        {1091, 992, 991, 990, 1}
    );

    REQUIRE(binning.size() == 5);
    REQUIRE(equal(binning[0], 0.));
    REQUIRE(equal(binning[1], 0.989));
    REQUIRE(equal(binning[2], 0.990));
    REQUIRE(equal(binning[3], 0.9901));
    REQUIRE(equal(binning[4], 1.));
}