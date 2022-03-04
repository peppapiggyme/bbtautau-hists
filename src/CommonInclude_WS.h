// C++
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>

// Root
#include "Math/Minimizer.h"
#include "Math/MinimizerOptions.h"

// RooFit
#include "RooWorkspace.h"
#include "RooWorkspaceHandle.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooAbsData.h"
#include "RooHist.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooRealSumPdf.h"
#include "Roo1DTable.h"
#include "RooConstVar.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "RooRandom.h"
#include "TStopwatch.h"
#include "RooNLLVar.h"
#include "RooMsgService.h"
#include "TMatrixDSymEigen.h"

// RooStat
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileInspector.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ToyMCImportanceSampler.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/MinNLLTestStat.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
 
// Huirun's pdf in ext/
#include "RooTwoSidedCBShape.h"

#include "CommonInclude.h"

using namespace RooFit;
using namespace RooStats;