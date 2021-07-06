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
    double ratio_high = 1.24;
    double ratio_low = 0.76;
    const char* atlas_label = "Internal";
    double signal_scale = 100.;
    std::string parameter = "TAG";
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
