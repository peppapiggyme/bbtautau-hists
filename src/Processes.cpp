#include "Processes.h"

#include <algorithm>

using namespace std;

Processes::Processes() noexcept
{
    m_procs = make_unique<vector<ProcessInfo*>>();
}

Processes::~Processes()
{
    for_each(content()->begin(), content()->end(), [](ProcessInfo* p) { delete p; });
    content()->clear();
}

void Processes::add(const string& nm, const string& nmtex, eProcessType tp, 
                    eProcess proc, const string& nmproc, int col) const
{
    m_procs->emplace_back(new ProcessInfo(nm, nmtex, tp, proc, nmproc, col));
    if (assign_norm_factors) 
    {
        m_procs->back()->norm_factor = normFactors(m_procs->back());
    }
}

void Processes::setScaleFactorsMap(const map<eProcess, double>& scaleFactors)
{
    m_scaleFactorsMap = scaleFactors;
}

double Processes::normFactors(ProcessInfo* p) const
{
    if (!m_scaleFactorsMap.empty()) 
    {
        if (m_scaleFactorsMap.count(p->process) == 0) 
        {
            return 1.0;
        } 
        else 
        {
            return m_scaleFactorsMap.at(p->process);
        }
    }

    // cross section info updated 2021-03-14 with 
    // https://gitlab.cern.ch/CxAODFramework/CxAODOperations/-/tree/698a65d6d73188e68252ec3aaf9c0d1ef41b2aea     
    
    switch (p->process)
    {

    case eProcess::ZJETSMG:
    case eProcess::ZllLFMG:
    case eProcess::ZllHFMG:
        return 2091.16 / 2097.184;
    case eProcess::STOPT_AMC:
        return 70.426 / 64.326;
    case eProcess::STOPT_HW7:
        return 70.426 / 70.431;
    case eProcess::STOPWT_AMC:
        return 71.699 / 71.700;
    case eProcess::STOPWT_HW7:
        return 71.699 / 71.697;
    case eProcess::STOPWT_DS:
        return 71.699 / 71.701;
    default:
        return 1.00;
    }

    // for bbll ZllCR
    // Yanlin Mar 31
    // switch (p->process)
    // {
    // case eProcess::ZllHF:
    //     return 1.35;
    // case eProcess::TOP:
    //     return 0.97;
    // case eProcess::ZllLF:
    //     return 1.00;
    // default:
    //     return 1.00;
    // }

    return 1.0;
}
