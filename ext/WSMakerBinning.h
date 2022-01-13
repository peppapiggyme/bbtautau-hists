#include <iostream>
#include <cmath>
#include <algorithm>
#include <forward_list>

#include <TMath.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphSmooth.h>
#include <TMathBase.h>
#include <memory>
#include <vector>

namespace WSMakerBinning
{

void smoothSyst(TH1* hnom, TH1* hsys);

std::vector<int> getLocalExtremaBinning(TH1* hnom, TH1* hsys, unsigned int nmax);

TH1* getRatioHist(TH1* hnom, TH1* hsys, const std::vector<int>& bins);

void getRatioHist(TH1* hnom, TH1* hsys, const std::vector<int>& bins, TH1* res);

std::vector<int> findExtrema(TH1* h) ;

void mergeBins(int lo, int hi, std::vector<int>& bins);

void smoothHistoRebin(TH1* hnom, TH1* hsys, const std::vector<int>& bins, bool smooth);

int findSmallerChi2(TH1* hnom, TH1* hsys, const std::vector<int>& extrema);

float computeChi2(TH1* hnom, TH1* hsys, int beg, int end);

float statError(TH1* hnom, int beg, int end);

double chi2(TH1 *const h_nom, TH1 *const h_smooth);

}