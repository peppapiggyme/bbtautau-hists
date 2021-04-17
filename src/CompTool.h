#ifndef COMPTOOL_H
#define COMPTOOL_H

#include "HistTool.h"

#include <string>
#include <vector>
#include <iterator>

using std::string;
using std::vector;

class CompInfo
{
public:
    bool logx = false;
    bool logy = false;
    bool atlas = true;
    bool shape_only = false;
    bool auto_binning = false;
    double ratio_high = 1.24;
    double ratio_low = 0.76;
    const char* atlas_label = "Simulation Internal";
    std::string parameter;
    std::string ratio_tex = "Ratio";
    // add new configs here
};

/**
 * @brief compare with the first entry
 */
class CompTool : public HistTool
{
public:
    explicit CompTool(const CompInfo* info);
    virtual ~CompTool() override;

public:
    virtual bool check(const Config* c) const override;
    // colorful
    virtual void paint(const Config* c) const;
    // drawing
    virtual void run(const Config* c) const override;

protected:
    const CompInfo* m_info;
};

#endif // COMPTOOL_H
