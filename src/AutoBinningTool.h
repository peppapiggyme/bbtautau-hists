#ifndef AUTOBINNINGTOOL_H
#define AUTOBINNINGTOOL_H

#include "HistTool.h"

#include <vector>
#include <map>
#include <queue>

using std::vector;

// ----------------------------------------------------------------------------
// Info
// ----------------------------------------------------------------------------
class AutoBinningInfo
{
public:
    unsigned n_bins = 8;
    double min_bkg = 5; // all bkg
    double min_mcstats = 0.2; // all bkg
    double required_mcstats = 0.2; // individual component (leading half)
    double eff_factor = 0.5; // effective factor
    double significance_delta = 1e-5; // stop adding bins if delta significance < 1e-5
    const char* atlas_label = "Simulation Internal";
    std::string parameter = "TAG";
};

// ----------------------------------------------------------------------------
// Bin Content Data Structure
// ----------------------------------------------------------------------------
class BinContent
{
public:
    explicit BinContent(const Config* c) noexcept { fill(c); }
    ~BinContent() noexcept { free(); }
    BinContent(BinContent& bc) = delete;
    BinContent& operator=(BinContent& bc) = delete;

    std::size_t size() const noexcept;

private:
    void fill(const Config* c) noexcept;
    void free() noexcept;

public:
    /**
     * Helper functions
     */
    std::pair<double*, double*> getArrayFromContentAndError(ProcessInfo* p) const;
    void showData() const;
    bool checkRangeXY(std::size_t x, std::size_t y) const;
    double sumFromXtoY(double* arr, std::size_t x, std::size_t y) const;
    double errSumFromXtoY(double* arr, std::size_t x, std::size_t y) const;
    // asympotitic significance
    // Z = sqrt(2*(s+b)log(1+s/b)-2*s)
    double significance(double fSig, double fBkg) const;
    // by algebra, the asymptotic significance can be added in quadraure
    // pass by value so that the queue can be reused
    double significanceBinByBin(std::priority_queue<std::size_t, std::vector<std::size_t>, std::greater<std::size_t>> pq) const; 

public:
    /**
     * @brief 
     * Case 1 : 
     * from right to left, scan for bin edges that satisfy:
     * - # Bkg > MinBkg (typically _5_) 
     *   -> validity of asymptotic test statistic formula
     * - MC stat. error of total bkg > MinTotalErr (typically _20%_)
     *   -> generally reasonable mc stats
     * - MC stat. error of leading two components (typically _20%_)
     *   -> enough stats to show systematic variation of major bkg components
     * 
     * Q: any cons?
     */
    bool passCaseOneFromXtoY(double fBkg, double fBkgErr, double fReqFactor, std::size_t x, std::size_t y) const;
    /**
     * @brief
     * Case 2 : 
     * from right to left, scan for bin edges that satisfy:
     * - # Bkg > MinBkg (typically _5_) 
     *   -> validity of asymptotic test statistic formula
     * - MC stat. error of total bkg > min of MinTotalErr (typically _20%_) 
     *   and EffectiveFactor * (# sig_i / # sig_tot) (typically _1_)
     *   -> generally reasonable mc stats
     *   -> (sometimes sharply) monotonic behaviour of bkg (and sig) distribution
     * - MC stat. error of leading two components (typically _20%_)
     *   -> enough stats to show systematic variation of major bkg components
     */
    bool passCaseTwoFromXtoY(double fBkg, double fBkgErr, double fEffErr, double fReqFactor, std::size_t x, std::size_t y) const;

public:
    std::map<std::pair<eProcessType, eProcess>, std::pair<double*, double*>> data;
    // two pointers -> current interval
    std::size_t curX;
    std::size_t curY;

private:
    std::size_t m_size;
};

// ---------------------------------------------------------------------------
// Method Enum
// ---------------------------------------------------------------------------
enum class BinningCriteria 
{
    None,
    CaseOne,
    CaseTwo
};


// ----------------------------------------------------------------------------
// Auto Binning Engines
// ----------------------------------------------------------------------------
class AutoBinningTool : public HistTool
{
public:
    AutoBinningTool(const AutoBinningInfo* info);
    virtual ~AutoBinningTool() override { delete m_binEdges; }

public:
    inline vector<double> binning() const { return *m_binEdges; }
    /**
     * @brief re-bin based on `m_binEdges`
     */
    virtual void rebin(const Config* c) const override;
    /**
     * @brief this function just for completeness, will always use Self option not matter which is provided
     * i.e. still re-bin based on `m_binEdges`
     */
    virtual void rebin(const Config* c, eRebinOption opt) const override { opt = eRebinOption::Self; HistTool::rebin(c, opt); }
    /**
     * @brief binning tool does not do the art
     */
    virtual void paint(const Config* c) const override { (void)c; }

protected:
    const AutoBinningInfo* m_info;
    vector<double>* m_binEdges;
};

/**
 * @brief
 * v1: binning for MVA score, optimise sensitivity
 * testing with hadhad PNNs and BDT
 */
class AutoBinningTool_v1 : public AutoBinningTool
{
public:
    AutoBinningTool_v1(const AutoBinningInfo* info, BinningCriteria bc);

public: 
    virtual bool check(const Config* c) const override;
    virtual void run(const Config* c) const override;

private:
    BinningCriteria m_bc;
};

/**
 * @brief
 * v2: binning for MVA score or kinematics
 * flatten the distribution by splitting into equal quantiles
 * 
 * @note
 * The bin edges may not fit the original binning
 * in that case ROOT will complain, but it should automatically round to the nearest valid edge
 * So the final binning not perfectly flat  
 */
class AutoBinningTool_v2 : public AutoBinningTool
{
public:
    AutoBinningTool_v2(const AutoBinningInfo* info, eProcess proc);

public:
    virtual void run(const Config* c) const override;

private:
    eProcess m_proc;

};


#endif