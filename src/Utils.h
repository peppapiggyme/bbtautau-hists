#ifndef UTILS_H
#define UTILS_H

#include "Processes.h"
#include "Regions.h"
#include "Variables.h"
#include "Systematics.h"

#include "TH1.h"
#include "TFile.h"
#include "TDirectory.h"

#include <iostream>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::pair;

namespace Utils {

    void histAssign(TH1* h, ProcessInfo* p, RegionInfo* r, VariableInfo* v);
    void histAssignSyst(TH1* h, ProcessInfo *p, const std::string& systname);

    /**
     * @note: so far only support CxAODReader like naming
     */
    string histString(const ProcessInfo* p, const RegionInfo* r, const VariableInfo* v);
    string systString(const SystematicInfo* s);
    string histStringSyst(const ProcessInfo* p, const RegionInfo* r, const VariableInfo* v, const SystematicInfo* s);

    static vector<pair<unsigned, unsigned>> paletteSysts = {         
        {kViolet, kAzure},
        {kViolet-1, kAzure-1},
        {kViolet-2, kAzure-2},
        {kViolet-3, kAzure-3},
        {kViolet-4, kAzure-4},
        {kViolet+1, kAzure+1},
        {kViolet+2, kAzure+2},
        {kViolet+3, kAzure+3},
        {kViolet+4, kAzure+4},
    };
}

#endif // UTILS_H
