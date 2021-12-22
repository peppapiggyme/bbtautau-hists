#ifndef DISTRIBUTIONTOOL_H
#define DISTRIBUTIONTOOL_H

#include "CompTool.h"

#include <string>
#include <vector>
#include <iterator>

using std::string;
using std::vector;

class DistributionInfo : public CompInfo
{
public:
};

/**
 * @brief only draw/compare distributions, no ratio 
 */
class DistributionTool : public CompTool
{
public:
    explicit DistributionTool(const DistributionInfo* info);
    virtual ~DistributionTool() override;

public:
    // drawing
    virtual void run(const Config* c) const override;
    virtual void run_overfit(const Config* c_train, const Config* c_test) const;

protected:
    const DistributionInfo* m_info;
};

#endif // COMPTOOL_H
