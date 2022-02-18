#include "ExamplesInclude_WS.h"
#include "TROOT.h"

#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    gROOT->SetBatch(1);

    if (argc < 2) 
    {
        cerr << "Usage: test <input> \n";
        return 0;
    }

    // test_ws_info(std::string(argv[1]));
    test_samplingdist(std::string(argv[1]));
    // test_w2r(std::string(argv[1]), std::string(argv[2]), std::string(argv[3]));
    // test_morphing(std::string(argv[1]));
}