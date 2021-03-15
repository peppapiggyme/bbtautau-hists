#ifndef PROCESSES_H
#define PROCESSES_H

#include "Regions.h"
#include "Variables.h"

#include "TColor.h"
#include "TH1.h"

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

using std::string;
using std::vector;
using std::shared_ptr;

// ENUMS
enum class eProcessType { DATA, BKG, SIG };

enum class eProcess {
    DATA, OTHERS, BKG, 
    TTBAR, TTBARTRUE, TTBARFAKE, TTBARTT, TTBARTF, TTBARFT, TTBARFF,
    STOP, STOPT, STOPS, STOPWT,
    STOPT_AMC, STOPT_HW7, 
    STOPWT_AMC, STOPWT_HW7, STOPWT_DS, 
    ZJETS, ZtautauJETS, ZtautauHF, ZtautauLF, ZllJETS, ZllLF, ZllHF,
    WJETS, WtauvJETS, WtauvHF, WtauvLF, WlvJETS, WlvLF, WlvHF, ttV,
    H, VH, WH, ZH, ttH, ggH, VBFH,
    DIBOSON, WW, WZ, ZZ, FAKE, QCD, MULTIJET,
    ZJETSMG, ZllLFMG, ZllHFMG,
    SIG, SMHH, HH, XtoHH, StoHH,
    HHKL1, HHKL3FROM1, HHKL0FROM1, HHKL10, HHKL10FROM1,
    HHKLXFROM1, HHKLXFROM10
};

// ENTRY
class ProcessInfo
{
public:
    ProcessInfo(const string& nm, const string& nmtex, eProcessType tp, 
                eProcess proc, const string& nmproc, int col) noexcept
        : name(nm), name_tex(nmtex), type(tp)
        , process(proc), process_name(nmproc), color(col) {}
public:
    string name; // same as histogram perfix
    string name_tex;
    eProcessType type;
    eProcess process;
    string process_name;
    int color; // EColor
    int rbg = 0xFFFFFF; // TODO: master of color platte
    double norm_factor = 1.0; // fitted norm
    TH1* histogram = nullptr; // depends on region and variable (will be set in Config)
    std::unordered_map<std::string, TH1*> systematic_histograms;
    bool isMerged = false;
    RegionInfo* current_region = nullptr;
    VariableInfo* current_variable = nullptr;
    
public:
    bool operator< (const ProcessInfo& p) const 
    { 
        return histogram->Integral() < p.histogram->Integral();
    }

    bool operator> (const ProcessInfo& p) const 
    { 
        return histogram->Integral() > p.histogram->Integral();
    }
};


// COLLECTION
class Processes
{
public:
    Processes() noexcept;
    ~Processes();
    Processes(Processes& ps) = delete;
    Processes& operator=(Processes& ps) = delete;

public:
    void add(const string& nm, const string& nmtex, eProcessType tp, 
             eProcess proc, const string& nmproc, int col) const;
    inline vector<ProcessInfo*>* content() const { return m_procs.get(); }

public:
    /**
     * @brief return current fitted normalisation factor
     */
    double normFactors(ProcessInfo* p) const;
    /**
     * @brief if true, will modify p->norm_factor, else it will keep as 1.0
     */
    bool assign_norm_factors = true;

private:
    unique_ptr<vector<ProcessInfo*>> m_procs;
};

#endif // PROCESSES_H
