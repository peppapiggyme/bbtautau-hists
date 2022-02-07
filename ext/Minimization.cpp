/// \file
/// Tools to perform likelihood minimization
#include <Math/MinimizerOptions.h>
#include <RooMinimizerFcn.h>
#include <RooMinimizer.h>
#include <RooWorkspace.h>
#include <RooNLLVar.h>
#include <RooAbsReal.h>
#include <RooStats/ModelConfig.h>

#include "Minimization.h"

using namespace std;

////////////////////////////////////////////////////
/// See minimize(RooAbsReal, ...) for details.
int EXOSTATS::minimize(RooNLLVar *nll, Int_t maxRetries, RooWorkspace *w, TString mu0Snapshot, TString nominalSnapshot,
                       Int_t debugLevel, Bool_t saveFitResult, RooFitResult **fitResult, Bool_t doMinos)
{
   RooAbsReal *fcn = (RooAbsReal *)nll;
   return EXOSTATS::minimize(fcn, maxRetries, w, mu0Snapshot, nominalSnapshot, debugLevel, saveFitResult, fitResult,
                             doMinos);
}

////////////////////////////////////////////////////
/// \param[in] nll pointer to the RooNLLVar to be minimized
/// \param[in] maxRetries maximum number of retries before giving up
/// \param[in] w pointer to the RooWorkSpace containing snapshots for fit retrials (see below)
/// \param[in] mu0Snapshot name of the snapshot of likelihood parameters after the background-only fit
/// \param[in] nominalSnapshot name of the snapshot of likelihood parameters at their nominal values
/// \param[in] debugLevel debug level (0 = verbose, 1 = debug, 2 = warning, 3 = error, 4 = fatal, 5 = silent)
/// \param[in] saveFitResult save fit results
/// \param[in] RooFitResult pointer to a pointer to a RooFitResult object, where fit results will be stored
/// \param[in] doMinos run Minos for error calculation
/// \param[out] status status of the minimization
///
/// The minimization is attempted in a iterative way, by playing with the Minuit version and strategy, until convergence
/// is reached or the maximum number of retrials has been tried.
///
/// If \c w is different from \c nullptr, the two specified parameter snapshots will also be used as a desperate measure
/// to improve fit stability.
///
/// Fit results can be optionally saved by activating the option \c saveFitResult and passing a pointer to a pointer to
/// RooFitResult.
int EXOSTATS::minimize(RooAbsReal *fcn, Int_t maxRetries, RooWorkspace *w, TString mu0Snapshot, TString nominalSnapshot,
                       Int_t debugLevel, Bool_t saveFitResult, RooFitResult **fitResult, Bool_t doMinos)
{
   static int nrItr = 0;
   if (debugLevel == 0) {
      cout << "Starting minimization" << endl;
   }

   int              printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
   RooFit::MsgLevel msglevel   = RooMsgService::instance().globalKillBelow();
   if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

   int          strat      = ROOT::Math::MinimizerOptions::DefaultStrategy();
   int          save_strat = strat;
   RooMinimizer minim(*fcn);
   minim.setStrategy(strat);
   minim.setPrintLevel(printLevel);
   minim.optimizeConst(2);

   int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                               ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());

   // up the strategy
   if (status != 0 && status != 1 && strat < 2) {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                              ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
   }

   if (status != 0 && status != 1 && strat < 2) {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                              ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
   }

   // cout << "status is " << status << endl;

   // //switch minuit version and try again
   if (status != 0 && status != 1) {
      string minType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
      string newMinType;
      if (minType == "Minuit2")
         newMinType = "Minuit";
      else
         newMinType = "Minuit2";

      cout << "Switching minuit type from " << minType << " to " << newMinType << endl;

      ROOT::Math::MinimizerOptions::SetDefaultMinimizer(newMinType.c_str());
      strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
      minim.setStrategy(strat);

      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                              ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());

      if (status != 0 && status != 1 && strat < 2) {
         strat++;
         cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
         minim.setStrategy(strat);
         status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                                 ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
      }

      if (status != 0 && status != 1 && strat < 2) {
         strat++;
         cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
         minim.setStrategy(strat);
         status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),
                                 ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
      }

      ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minType.c_str());
   }

   if (status != 0 && status != 1) {
      nrItr++;
      if (nrItr > maxRetries) {
         nrItr = 0;
         cout << "WARNING::Fit failure unresolved with status " << status << endl;
         return status;
      } else {
         if (nrItr == 0) { // retry with mu=0 snapshot
            if (w)
               w->loadSnapshot(mu0Snapshot);
            else
               cout << "WARNING: workspace not provided, unable to set mu=0 snapshot; will simply retry as is" << endl;
            return minimize(fcn);
         } else if (nrItr == 1) { // retry with nominal snapshot
            if (w)
               w->loadSnapshot(nominalSnapshot);
            else
               cout << "WARNING: workspace not provided, unable to set nominal NP snapshot; will simply retry as is"
                    << endl;
            return minimize(fcn);
         }
      }
   }

   if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(msglevel);
   ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);

   if (nrItr != 0) cout << "Successful fit" << endl;
   nrItr = 0;

   if (doMinos) {
      minim.minos();
   }

   // save
   if (saveFitResult) {
      *fitResult = minim.save();
   }

   return status;
}

/// \param[in] modelConfig pointer to the \c ModelConfig to take the p.d.f. from
/// \param[in] data dataset
/// \param[in] numCPU number of CPUs to use in the likelihood calculation
/// \param[out] nll negative log-likelihood
RooNLLVar *EXOSTATS::createNLL(RooStats::ModelConfig *modelConfig, RooAbsData *data, Int_t numCPU)
{
   return EXOSTATS::createNLL(modelConfig->GetPdf(), data, modelConfig->GetNuisanceParameters(), numCPU);
}

/// \param[in] pdf pointer to the p.d.f
/// \param[in] data dataset
/// \param[in] nuis nuisance parameters (if any)
/// \param[in] numCPU number of CPUs to use in the likelihood calculation
/// \param[out] nll negative log-likelihood
RooNLLVar *EXOSTATS::createNLL(RooAbsPdf *pdf, RooAbsData *data, const RooArgSet *nuis, Int_t /*numCPU*/)
{
   // VI: I commented out the number of CPU option due to this bug:
   // https://root-forum.cern.ch/t/numcpu-crash-with-very-simple-example/28532
   RooNLLVar *nll = nullptr;
   if (nuis != nullptr)
      nll = (RooNLLVar *)pdf->createNLL(*data, RooFit::Constrain(*nuis), // RooFit::NumCPU(numCPU, 3),
                                        RooFit::Optimize(2), RooFit::Offset(true));
   else
      nll =
         (RooNLLVar *)pdf->createNLL(*data, /*RooFit::NumCPU(numCPU, 3),*/ RooFit::Optimize(2), RooFit::Offset(true));
   return nll;
}

/// \param[in] modelConfig pointer to the \c ModelConfig to take the p.d.f. from
/// \param[in] data dataset
/// \param[in] doMinos if true, Minos errors are calculated
/// \param[in] numCPU number of CPUs to use in the likelihood calculation
/// \param[out] status Minuit status
Int_t EXOSTATS::fitModelConfig(RooStats::ModelConfig *modelConfig, RooAbsData *data, Bool_t doMinos, Int_t numCPU)
{
   return EXOSTATS::fitPdf(modelConfig->GetPdf(), data, doMinos, modelConfig->GetNuisanceParameters(), numCPU);
}

/// \param[in] pdf pointer to the p.d.f
/// \param[in] data dataset
/// \param[in] doMinos if true, Minos errors are calculated
/// \param[in] nuis nuisance parameters (if any)
/// \param[in] numCPU number of CPUs to use in the likelihood calculation
/// \param[out] status Minuit status
Int_t EXOSTATS::fitPdf(RooAbsPdf *pdf, RooAbsData *data, Bool_t doMinos, const RooArgSet *nuis, Int_t numCPU)
{
   auto nll = createNLL(pdf, data, nuis, numCPU);
   return EXOSTATS::minimize(nll, 3, nullptr, "", "", 2, kFALSE, nullptr, doMinos);
}

/// \param[in] pdf pointer to the p.d.f
/// \param[in] data dataset
/// \param[in] doMinos if true, Minos errors are calculated
/// \param[in] nuis nuisance parameters (if any)
/// \param[in] numCPU number of CPUs to use in the likelihood calculation
/// \param[out] result pointer to fit result
RooFitResult *EXOSTATS::fitPdfRes(RooAbsPdf *pdf, RooAbsData *data, Bool_t doMinos, const RooArgSet *nuis, Int_t numCPU)
{
   auto          nll = createNLL(pdf, data, nuis, numCPU);
   RooFitResult *result;
   EXOSTATS::minimize(nll, 3, nullptr, "", "", 2, kTRUE, &result, doMinos);
   return result;
}
