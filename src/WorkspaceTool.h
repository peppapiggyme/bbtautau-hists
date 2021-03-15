/**
 * Copy-paste a lot from FitCrossCheckforLimits.C
 */

#ifndef WORKSPACETOOL_H
#define WORKSPACETOOL_H

#include "CommonInclude_WS.h"
#include "Utils_WS.h"

#include <memory>
#include <vector>

enum class VerboseLevel { DEBUG, INFO };

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
  inline void setVerboseLevel(VerboseLevel lvl) { m_lvl = lvl; } 

public:
  bool isSimultaneousPdfOK();
  bool isChannelNameOK();
  void getNominalValueOfNP();

public:
  void printModelObservables();
  void printNPs();
  void printSubChannels();
  void printNumberOfEvents(RooAbsPdf *pdf);

public:
  static constexpr double LumiRelError = 0.017;

protected:
  VerboseLevel m_lvl = VerboseLevel::INFO;
  std::shared_ptr<TFile> m_fin;
  RooWorkspace *m_ws;
  ModelConfig *m_mc;
  RooAbsData *m_data;
  std::vector<NPInfo> m_allNPs_afterFit;
  std::unordered_map<std::string, double> m_nomNPs;
};

#endif