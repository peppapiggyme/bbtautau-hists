#ifndef HISTTOOL_H
#define HISTTOOL_H

#include "Config.h"
#include "TFile.h"

using std::string;

enum class eRebinOption
{
    Self,       /// By running the re-implementation of `rebin(c)` function, e.g. self is AutoBinningTool
    N_Rebin,    /// By n_rebin of current_variable
    Array       /// By a given binning array that describes the bin edges
};

// Use histograms
class HistTool
{
public:
    HistTool() noexcept {}
    virtual ~HistTool() noexcept {};

public:
    virtual bool check(const Config* c) const;
    /// @todo const
    virtual void manipulate(Config* c);
    virtual void rebin(const Config* c) const { (void)c; }
    virtual void rebin(const Config* c, eRebinOption opt) const;
    virtual void makeYield(const Config* c, const std::string& tag="TAG") const;

public:
    virtual void paint(const Config* c) const = 0;
    virtual void run(const Config* c) const = 0;

public:
    string output_path;
};

/**
 * Sometimes the static function is handy
 */
class HistToolHelper
{
public:
    static bool check(const Config* c);
    static void rebinByNRebin(const Config* c);
    static void rebinByArray(const Config* c);
};

#endif // HISTTOOL_H