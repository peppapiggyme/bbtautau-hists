#include "ExamplesInclude_WS.h"
#include "TROOT.h"

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        cerr << "Usage: run-pd <list.txt> <output>\n";
        return 0;
    }

    gROOT->SetBatch(1);

    vector<string> filenames;

    ifstream f(argv[1]);
    string line;
    
    while (getline(f, line))
    {
        istringstream iss(line);
        string x;

        while (iss >> x)
        {
            if (x[0] == '#') continue;
            filenames.push_back(x);
            cout << x << endl;
        }
    }

    test_pseudodata(filenames, argv[2]);

    return 0;
}