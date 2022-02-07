#ifndef __EXOSTATS_MINIMIZATION_H__
#define __EXOSTATS_MINIMIZATION_H__

/// \file
/// Tools to perform likelihood minimization

#include <TString.h>

class RooAbsData;
class RooAbsPdf;
class RooAbsReal;
class RooArgSet;
class RooFitResult;
class RooNLLVar;
class RooWorkspace;
namespace RooStats {
class ModelConfig;
}

namespace EXOSTATS {

/// Minimize a RooNLLVar
Int_t minimize(RooNLLVar *nll, Int_t maxRetries = 3, RooWorkspace *w = nullptr,
               TString mu0Snapshot = "conditionalNuis_0", TString nominalSnapshot = "nominalNuis", Int_t debugLevel = 2,
               Bool_t saveFitResult = kFALSE, RooFitResult **fitResult = nullptr, Bool_t doMinos = kFALSE);

/// Minimize a RooAbsReal
Int_t minimize(RooAbsReal *fcn, Int_t maxRetries = 3, RooWorkspace *w = nullptr,
               TString mu0Snapshot = "conditionalNuis_0", TString nominalSnapshot = "nominalNuis", Int_t debugLevel = 2,
               Bool_t saveFitResult = kFALSE, RooFitResult **fitResult = nullptr, Bool_t doMinos = kFALSE);

/// Create a NLL from a ModelConfig's p.d.f.
RooNLLVar *createNLL(RooStats::ModelConfig *modelConfig, RooAbsData *data, Int_t numCPU = 4);

/// Create a NLL from a generic p.d.f.
RooNLLVar *createNLL(RooAbsPdf *pdf, RooAbsData *data, const RooArgSet *nuis = nullptr, Int_t numCPU = 4);

/// Fit the default p.d.f of a ModelConfig
Int_t fitModelConfig(RooStats::ModelConfig *modelConfig, RooAbsData *data, Bool_t doMinos = kFALSE, Int_t numCPU = 4);

/// Fit a generic p.d.f. and retrieve the fit status
Int_t fitPdf(RooAbsPdf *pdf, RooAbsData *data, Bool_t doMinos = kFALSE, const RooArgSet *nuis = nullptr,
             Int_t numCPU = 4);

/// Fit a generic p.d.f. and retrieve the full fit results
RooFitResult *fitPdfRes(RooAbsPdf *pdf, RooAbsData *data, Bool_t doMinos = kFALSE, const RooArgSet *nuis = nullptr,
                        Int_t numCPU = 4);

} // namespace EXOSTATS

#endif
