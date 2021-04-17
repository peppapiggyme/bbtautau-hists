#include "HistTool.h"

#include <exception>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <iomanip>

using namespace std;

# define FOUR_COLUMN_TABLE(A, B, C, D) \
       std::left << setw(36) << A \
    << std::left << setw(15) << B \
    << std::left << setw(30) << C \
    << std::left << setw(15) << D << endl

# define FIVE_COLUMN_TABLE(A, B, C, D, E) \
       std::left << setw(40) << A \
    << std::left << setw(15) << B \
    << std::left << setw(15) << C \
    << std::left << setw(15) << D \
    << std::left << setw(15) << E << endl

bool HistTool::check(const Config* c) const 
{
    vector<ProcessInfo*>* ps = c->processes->content();
    clog << "INFO: removing nullptrs\n";
    ps->erase(remove_if(ps->begin(), ps->end(), [](const ProcessInfo* p) {
        return !p->histogram; }), ps->end());

    if (ps->size() < 1)
    {
        cerr << "FAIL: empty input\n";
        return false;
    }

    return true;
}

void HistTool::manipulate(Config* c) 
{
    c->setManipulated(true);
    vector<ProcessInfo*>* ps_in_c = c->processes->content();
    clog << "INFO: merging\n";
    map<eProcess, vector<ProcessInfo*>> procs;
    for_each(ps_in_c->begin(), ps_in_c->end(), [&procs](ProcessInfo* p) {
        procs[p->process].push_back(p); });

    for_each(procs.begin(), procs.end(), [&ps_in_c, &c](pair<const eProcess, vector<ProcessInfo*>>& pair) {
        // TODO: basically copy, copy assignment?
        ProcessInfo* front = pair.second.front();
        ProcessInfo* merged = new ProcessInfo(
            front->process_name, front->process_name, front->type, 
            front->process, front->process_name, front->color);
        
        merged->histogram = (TH1*)front->histogram->Clone();

        for_each(pair.second.begin() + 1, pair.second.end(), [&merged](const ProcessInfo* p) { 
            merged->histogram->Add(p->histogram); });

        if (!(front->systematic_histograms.empty()))
        {
            for (auto & pp : front->systematic_histograms)
            {
                std::cout << pp.first << std::endl;
                merged->systematic_histograms[pp.first] = (TH1*)front->systematic_histograms[pp.first]->Clone();
                for_each(pair.second.begin() + 1, pair.second.end(), [&merged, &pp](const ProcessInfo* p) {
                    std::cout << p->name << std::endl;
                    if (!p->systematic_histograms.empty()) {
                        merged->systematic_histograms[pp.first]->Add(p->systematic_histograms.at(pp.first));
                    }
                });
            }
        }

        merged->isMerged = true;
        /// @note: make sure the process are not mixed
        merged->norm_factor = front->norm_factor;
        merged->current_region = front->current_region;
        merged->current_variable = front->current_variable;
        
        c->current_region = merged->current_region;
        c->current_variable = merged->current_variable;

        ps_in_c->emplace_back(merged);
    });

    ps_in_c->erase(remove_if(ps_in_c->begin(), ps_in_c->end(), [](const ProcessInfo* p) {
        return !p->isMerged; }), ps_in_c->end());

    sort(ps_in_c->begin(), ps_in_c->end(), [](const ProcessInfo* p1, const ProcessInfo* p2) {
        return p1->type < p2->type; });

    // if do draw stack, the bkg and sig must be sorted!
}

void HistTool::rebin(const Config* c, eRebinOption opt) const
{
    switch (opt)
    {
    case eRebinOption::Self:
        rebin(c);
        break;
    case eRebinOption::N_Rebin:
        HistToolHelper::rebinByNRebin(c);
        break;
    case eRebinOption::Array:
        HistToolHelper::rebinByArray(c);
        break;
    default:
        /* never happen */
        break;
    }
}

void HistTool::makeYield(const Config* c, const std::string& tag) const
{   
    ostringstream oss;
    oss << output_path << "/" 
        << c->current_region->name << "_" 
        << tag << "_" 
        << c->current_variable->name << ".txt";

    ofstream fout(oss.str());
    vector<ProcessInfo*>* ps = c->processes->content();

    // total backgrouds
    bool hasBkg = false;
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
        double eOverI = integral > (double)0. ? error / integral : 0.;
        fout << FIVE_COLUMN_TABLE(p->name, nentries, integral * p->norm_factor, error * p->norm_factor, eOverI);

        for (auto& pp : p->systematic_histograms) {
            auto systEntries = pp.second->GetEntries();
            auto systInt = pp.second->Integral();
            fout << " |- " << FOUR_COLUMN_TABLE(pp.first, systEntries, systInt, systInt / integral - 1.f);
        }

        if (p->type == eProcessType::BKG)
        {
            hasBkg = true;
            entriesBkg += nentries;
            sumBkg += integral * p->norm_factor; 
            errBkg += error * error * p->norm_factor * p->norm_factor;
        }
    }
    if (hasBkg) {
        errBkg = sqrt(errBkg);
        fout << FIVE_COLUMN_TABLE("Total Bkg", entriesBkg, sumBkg, errBkg, errBkg / sumBkg); 
    }

    clog << "INFO: Yields saved in " << oss.str() << '\n';
}


bool HistToolHelper::check(const Config* c) 
{
    vector<ProcessInfo*>* ps = c->processes->content();
    clog << "INFO: removing nullptrs\n";
    ps->erase(remove_if(ps->begin(), ps->end(), [](const ProcessInfo* p) {
        return !p->histogram; }), ps->end());

    if (ps->size() < 1)
    {
        cerr << "FAIL: empty input\n";
        return false;
    }

    return true;
}

void HistToolHelper::rebinByNRebin(const Config* c)
{
    vector<ProcessInfo*>* ps = c->processes->content();
    for_each(ps->begin(), ps->end(), [&c](ProcessInfo* p) {
        p->histogram->Rebin(c->current_variable->n_rebin);
        for (auto& pp : p->systematic_histograms)
        {
            pp.second->Rebin(c->current_variable->n_bins);
        }
    });
}

void HistToolHelper::rebinByArray(const Config* c)
{
    /// @todo refactor this if statement
    vector<ProcessInfo*>* ps = c->processes->content();
    if (c->current_variable->binning)
    {
        for_each(ps->begin(), ps->end(), [&c](ProcessInfo* p) {
            TH1* rebinned = p->histogram->Rebin(c->current_variable->n_bins, p->histogram->GetName(), c->current_variable->binning);
            p->histogram = (TH1*)rebinned->Clone();
            for (auto& pp : p->systematic_histograms)
            {
                TH1* rebinned_pp = pp.second->Rebin(c->current_variable->n_bins, pp.second->GetName(), c->current_variable->binning);
                pp.second = (TH1*)rebinned_pp->Clone();
            }
        });
    }
    else
    {
        clog << "WARN: Binning not given, will rebin based on rebinByNRebin(c) instead\n";
        HistToolHelper::rebinByNRebin(c);
    }
}
