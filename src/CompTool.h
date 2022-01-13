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
    bool save_ratio = false;
    bool do_envelope = false;
    double ratio_high = 1.24;
    double ratio_low = 0.76;
    double legend_scaling_horizontal = 1.0;
    double legend_scaling_vertical = 1.0;
    double xmin_for_test = DBL_MIN;
    double xmax_for_test = DBL_MAX;
    double ymax_ratio_nolog = 1.4;
    const char* atlas_label = "Simulation Internal";
    std::string parameter;
    std::string ratio_tex = "Ratio";
    std::string output_format = "png";
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
