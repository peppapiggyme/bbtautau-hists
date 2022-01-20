#ifndef CONFIG_H
#define CONFIG_H

#include "Processes.h"
#include "Regions.h"
#include "Variables.h"
#include "Systematics.h"
#include "Utils.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"

#include <string>
#include <memory>

using std::string;
using std::unique_ptr;

class BasicInfo
{
public:
    BasicInfo(const string& ecm, const string& lumi) noexcept;

public:
    string ecm;
    string luminosity;
    Utils::NameConvention name_convention = Utils::NameConvention::CxAODReader;
    bool make_systs_shape_only = false; // used in comptool, but might also be useful in draw stack
};


class Config
{
public:
    Config(const BasicInfo* b, const Processes* ps, 
           const Regions* rs, const Variables* vs, const Systematics* ss=nullptr) noexcept;
    ~Config() noexcept;

    Config(Config& ht) = delete;
    Config& operator=(Config& ht) = delete;

public:
    const BasicInfo* basic;
    const Processes* processes;
    const Regions* regions;
    const Variables* variables;
    const Systematics* systematics;

public:
    void load(const string& fn, const string& dir="");
    void updateHistogramPtr(RegionInfo* r, VariableInfo* v);
    inline void setManipulated(bool m) { m_manipulated = m; }
    inline const unique_ptr<TFile>& getInputTFile() const { return m_fin; }
    inline const std::string& getInputDirStr() const { return m_dir; }
public:
    RegionInfo* current_region;
    VariableInfo* current_variable;

protected:
    unique_ptr<TFile> m_fin;
    std::string m_dir;
    bool m_loaded;
    bool m_manipulated;
};

#endif // CONFIG_H
