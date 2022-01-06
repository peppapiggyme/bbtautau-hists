#ifndef DRAWSTACKTOOL_H
#define DRAWSTACKTOOL_H

#include "HistTool.h"

#include <string>
#include <vector>
#include <iterator>

using std::string;
using std::vector;

class DrawStackInfo
{
public:
    bool blind = true;
    bool logx = false;
    bool logy = false;
    bool atlas = true;
    bool show_scaling = true;
    bool draw_overflow = true;
    int legend_ncolumns = 2;
    int signal_linestyle = 1;
    int signal_linewidth = 2;
    double ratio_high = 1.24;
    double ratio_low = 0.76;
    double signal_scale = 100.;
    double legend_scaling_horizontal = 1.0;
    double legend_scaling_vertical = 1.0;
    const char* atlas_label = "Internal";
    std::string parameter = "TAG";
    std::string output_format = "png";
    // add new configs here
};

class DrawStackTool : public HistTool
{
public:
    explicit DrawStackTool(const DrawStackInfo* info);
    virtual ~DrawStackTool() override;

public:
    virtual bool check(const Config* c) const override;
    // merge and sort bkgs
    virtual void manipulate(Config* c) override;
    // colorful
    virtual void paint(const Config* c) const;
    // drawing
    virtual void run(const Config* c) const override;

protected:
    vector<ProcessInfo*>::iterator m_it_data;
    vector<ProcessInfo*>::iterator m_it_bkg;
    vector<ProcessInfo*>::iterator m_it_sig;
    vector<ProcessInfo*>::iterator m_it_end;
    const DrawStackInfo* m_info;
};

#endif // DRAWSTACKTOOL_H
