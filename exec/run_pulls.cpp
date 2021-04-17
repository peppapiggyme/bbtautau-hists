#include "ExamplesInclude_WS.h"
#include "TROOT.h"

#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    gROOT->SetBatch(1);

    if (argc < 3) 
    {
        cerr << "Usage: run-pulls <input> <output>\n";
        return 0;
    }

    test_pulls(std::string(argv[1]), std::string(argv[2])+".txt");
    test_pulls_plot(std::string(argv[2])+".txt", std::string(argv[2])+".png");
}