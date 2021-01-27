#include "HistTool.h"

#include <exception>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <iomanip>

using namespace std;

# define FIVE_COLUMN_TABLE(A, B, C, D, E) \
       std::left << setw(25) << A \
    << std::left << setw(15) << B \
    << std::left << setw(15) << C \
    << std::left << setw(15) << D \
    << std::left << setw(15) << E << endl

bool HistTool::check(const Config* c) const 
{
    vector<ProcessInfo*>* ps = c->processes->content();
    clog << "removing nullptrs\n";
    ps->erase(remove_if(ps->begin(), ps->end(), [](const ProcessInfo* p) {
        return !p->histogram; }), ps->end());

    if (ps->size() < 1)
    {
        cerr << "FAIL: empty input\n";
        return false;
    }

    return true;
}

void HistTool::manipulate(const Config* c) 
{
    vector<ProcessInfo*>* ps_in_c = c->processes->content();
    clog << "merging\n";
    map<eProcess, vector<ProcessInfo*>> procs;
    for_each(ps_in_c->begin(), ps_in_c->end(), [&procs](ProcessInfo* p) {
        procs[p->process].push_back(p); });

    for_each(procs.begin(), procs.end(), [&ps_in_c](pair<const eProcess, vector<ProcessInfo*>>& pair) {
        // TODO: basically copy, copy assignment?
        ProcessInfo* front = pair.second.front();
        ProcessInfo* merged = new ProcessInfo(
            front->process_name, front->process_name, front->type, 
            front->process, front->process_name, front->color);
        
        merged->histogram = (TH1*)front->histogram->Clone();

        for_each(pair.second.begin() + 1, pair.second.end(), [&merged](const ProcessInfo* p) { 
            merged->histogram->Add(p->histogram); });

        merged->isMerged = true;
        merged->current_region = front->current_region;
        merged->current_variable = front->current_variable;
        
        ps_in_c->emplace_back(merged);
    });

    ps_in_c->erase(remove_if(ps_in_c->begin(), ps_in_c->end(), [](const ProcessInfo* p) {
        return !p->isMerged; }), ps_in_c->end());

    sort(ps_in_c->begin(), ps_in_c->end(), [](const ProcessInfo* p1, const ProcessInfo* p2) {
        return p1->type < p2->type; });

    // if do draw stack, the bkg and sig must be sorted!
}

void HistTool::makeYield(const Config* c, const string& fn) const
{
    ofstream fout(fn);
    vector<ProcessInfo*>* ps = c->processes->content();

    // total backgrouds
    int entriesBkg = 0;
    double sumBkg  = 0.0;
    double errBkg  = 0.0;

    fout << FIVE_COLUMN_TABLE("Process", "Entries", "Yield", "Error", "Rel.Err.");
    for (ProcessInfo* p : *ps)
    {
        double error;
        int from = 0;
        int to = p->histogram->GetNbinsX() + 1;
        int nentries = p->histogram->GetEntries();
        double integral = p->histogram->IntegralAndError(from, to, error, "");
        double eOverI = *(long*)(&integral) ? error / integral : 0.;
        fout << FIVE_COLUMN_TABLE(p->name, nentries, integral, error, eOverI);

        if (p->type == eProcessType::BKG)
        {
            entriesBkg += nentries;
            sumBkg += integral; 
            errBkg += error * error;
        }
    }
    errBkg = sqrt(errBkg);
    fout << FIVE_COLUMN_TABLE("Total Bkg", entriesBkg, sumBkg, errBkg, errBkg / sumBkg); 

    clog << "Yields saved in " << fn << '\n';
}
