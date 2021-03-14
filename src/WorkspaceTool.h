#ifndef WORKSPACETOOL_H
#define WORKSPACETOOL_H

#include "CommonInclude_WS.h"
#include "Utils_WS.h"

#include <memory>

struct NPInfo
{
  std::string name;
  std::string fromFit;
  double value;
  double error_high;
  double error_low;
};

class WorkspaceTool
{
public:
  WorkspaceTool() = default;

  void load(const std::string &fn, const std::string &ws_name, 
            const std::string &modelConfig_name, const std::string &obsData_name);

protected:
  std::shared_ptr<TFile> m_fin;
  RooWorkspace *m_ws;
  ModelConfig *m_mc;
  RooAbsData *m_data;
};

#endif