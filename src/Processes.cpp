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

double Processes::normFactors(ProcessInfo* p) const
{
    // cross section info updated 2021-03-14 with 
    // https://gitlab.cern.ch/CxAODFramework/CxAODOperations/-/tree/698a65d6d73188e68252ec3aaf9c0d1ef41b2aea     
    
    if (p->name == "ggZHbb" || p->name == "qqZHbb")
    {
        return 0.9475481609231877;
    }
    else if (p->name == "ggZHtautau" || p->name == "qqZHtautau")
    {
        return 0.9483127795014143;
    }
    else if (p->name == "ggFHtautau")
    {
        return 1.2153682446311185;
    }

    switch (p->process)
    {
    case eProcess::ZllHF:
        return 1.2308196024245677;
    case eProcess::ZtautauHF:
        return 1.1808929113641928;
    case eProcess::ZllLF:
        return 0.8657669838456886;
    case eProcess::ZtautauLF:
        return 0.805916366803583;
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
    case eProcess::TTBAR:
        return 0.8662766660691217;
    case eProcess::TTBARFAKE:
        return 0.8367832792819666;
    case eProcess::TTBARTRUE:
    case eProcess::TTBARTT:
        return 1.00;
    case eProcess::HHKL3FROM1:
        return 3.00;
    case eProcess::HHKL0FROM1:
        return 0.50;
    case eProcess::FAKE:
        return 0.9508223744401554;
    default:
        return 1.00;
    }

    // for bbll ZllCR
    // Yanlin Mar 31
    // switch (p->process)
    // {
    // case eProcess::ZllHF:
    //     return 1.22;
    // case eProcess::TOP:
    //     return 0.98;
    // case eProcess::ZllLF:
    //     return 0.89;
    // default:
    //     return 1.00;
    // }
}
