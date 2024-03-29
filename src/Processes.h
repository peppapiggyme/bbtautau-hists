#ifndef PROCESSES_H
#define PROCESSES_H

#include "Regions.h"
#include "Variables.h"

#include "TColor.h"
#include "TH1.h"

#include <memory>
#include <string>
#include <vector>
#include <map>

using std::string;
using std::vector;
using std::map;
using std::shared_ptr;

// ENUMS
enum class eProcessType { DATA, BKG, SIG };

enum class eProcess {
    DATA, BKG, OTHERS, // 2 
    TOP, TTBAR, TTBARTRUE, TTBARFAKE, TTBARTT, TTBARTF, TTBARFT, TTBARFF, // 10
    TTBARPH7, TTBARAMCP8, TTBARPP8RADHI, // 13
    STOP, STOPT, STOPS, STOPWT, // 17
    STOPT_AMC, STOPT_HW7, // 19
    STOPWT_AMC, STOPWT_HW7, STOPWT_DS, // 22 
    ZJETS, ZtautauJETS, ZtautauHF, ZtautauLF, ZllJETS, ZllLF, ZllHF, // 29
    WJETS, WtauvJETS, WtauvHF, WtauvLF, WlvJETS, WlvLF, WlvHF, ttV, // 37
    H, VH, WH, ZH, ttH, ggH, VBFH, // 44
    DIBOSON, WW, WZ, ZZ, FAKE, QCD, MULTIJET, //51
    ZJETSMG, ZllLFMG, ZllHFMG, // 54
    SIG, SMHH, HH, XtoHH, StoHH, // 59
    HHKL1, HHKL3FROM1, HHKL0FROM1, HHKL10, HHKL10FROM1, // 64
    HHKLXFROM1, HHKLXFROM10, //66
    TTBARHH, TTBARLH, TTBARZCR, // 69
    X800, X900, X1000, X1100, //74
    TTBAR_NOM, TTBAR_AMC, TTBAR_HW7, TTBAR_HDAMP,

    // For sorting
    P1, P2, P3, P4, P5, P6, P7, P8, 
    P9, P10, P11, P12, P13, P14, P15, P16,
    P101, P102, P103,
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
    std::map<std::string, TH1*> systematic_histograms;
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
     * @brief user defined the process scale factor mapping rules
     * @note set before adding procs to be used
     */
    void setScaleFactorsMap(const map<eProcess, double>& scaleFactors);
    /**
     * @brief return current fitted normalisation factor
     * if the (user defined) map is not empty, then use the values in the map.
     * else use the default values
     */
    double normFactors(ProcessInfo* p) const;
    /**
     * @brief if true, will modify p->norm_factor, else it will keep as 1.0
     */
    bool assign_norm_factors = true;

private:
    unique_ptr<vector<ProcessInfo*>> m_procs;
    map<eProcess, double> m_scaleFactorsMap;
};

#endif // PROCESSES_H
