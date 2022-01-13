#include "WSMakerBinning.h"
#include "Utils.h"

namespace WSMakerBinning
{

void smoothSyst(TH1* hnom, TH1* hsys) {
  const auto& bins = getLocalExtremaBinning(hnom, hsys, 0);
  Tools::printVector(bins);
  // Call the rebinning function
  smoothHistoRebin(hnom, hsys, bins, false);
}

std::vector<int> getLocalExtremaBinning(TH1* hnom, TH1* hsys, unsigned int nmax) {
  // This implementation is iterative.
  // A faster (say analytic) implementation is possible if this one proves to be
  // too slow. This one is however easier to write and read
  std::vector<int> res;
  double err = 0;
  float sum = hnom->IntegralAndError(0, hnom->GetNbinsX()+1, err);
  // too large stat unc: no shape
  if(fabs(err / sum) > 0.05) {
    res.push_back(1);
    res.push_back(hnom->GetNbinsX()+1);
    return res;
  }

  // normal case. Then, beginning with no rebinning
  for(int i = 1 ; i < hnom->GetNbinsX()+2 ; i++) {
    res.push_back(i);
  }
  // Second pass first seems nicer
  // In practice, more noise in the fit

  // first pass start ?
  TH1* ratio = getRatioHist(hnom, hsys, res);
  std::vector<int> extrema = findExtrema(ratio);

  while(extrema.size() > nmax+2) {
    int pos = findSmallerChi2(hnom, hsys, extrema);
    mergeBins(extrema[pos], extrema[pos+1], res);
    getRatioHist(hnom, hsys, res, ratio);
    extrema = findExtrema(ratio);
  }
  // end of the first pass ?

  // second pass to avoid bins with too large stat uncertainty
  std::vector<int>::iterator fst = res.end();
  std::vector<int>::iterator lst = res.end();
  std::vector<int> to_remove;
  --lst;
  --fst;
  while(fst != res.begin()) {
    if(fst == lst) {
      --fst;
    }
    else {
      float statE = statError(hnom, *fst, *lst);
      if(statE > 0.05) {
        to_remove.push_back(fst - res.begin());
        --fst;
      }
      else {
        lst=fst;
      }
    }
  }
  for(unsigned int i = 0; i < to_remove.size(); i++) {
    res.erase(res.begin()+to_remove[i]);
  }
  // end of the second pass

  delete ratio;
  return res;
}


TH1* getRatioHist(TH1* hnom, TH1* hsys, const std::vector<int>& bins) {
  TH1* res = (TH1*) hsys->Clone();
  getRatioHist(hnom, hsys, bins, res);
  return res;
}

void getRatioHist(TH1* hnom, TH1* hsys, const std::vector<int>& bins, TH1* res) {
  for(unsigned int iRefBin = 0; iRefBin<bins.size()-1; iRefBin++) {
    float nomInt = hnom->Integral(bins.at(iRefBin), bins.at(iRefBin+1) - 1);
    float varInt = hsys->Integral(bins.at(iRefBin), bins.at(iRefBin+1) - 1);
    for(int b=bins.at(iRefBin); b<bins.at(iRefBin+1); b++) {
      if(nomInt != 0) {
        res->SetBinContent( b, varInt / nomInt );
      }
      else {
        res->SetBinContent( b, 0 );
      }
    }
  }
}

std::vector<int> findExtrema(TH1* h) {
  std::vector<int> res;
  res.push_back(1);
  int status = 0; // 1: potential max, -1: potential min
  int k = 1;
  for(int i=2; i<h->GetNbinsX()+1; i++) {
    // special rule for bins with 0 stat. Keep going on, until one finds another bin to compare to
    if(h->GetBinContent(i) < 1.e-6) {
      continue;
    }
    if(status == 1 && h->GetBinContent(i) < h->GetBinContent(k) - 1.e-6) {
      res.push_back(i-1);
      status = -1;
    }
    if(status == -1 && h->GetBinContent(i) > h->GetBinContent(k) + 1.e-6) {
      res.push_back(i-1);
      status = 1;
    }
    if(status == 0 && h->GetBinContent(i) < h->GetBinContent(k) - 1.e-6) {
      status = -1;
    }
    if(status == 0 && h->GetBinContent(i) > h->GetBinContent(k) + 1.e-6) {
      status = 1;
    }
    k = i;
  }
  res.push_back(h->GetNbinsX());

  return res;
}

// inclusive in lo and hi
void mergeBins(int lo, int hi, std::vector<int>& bins) {
  std::vector<int>::iterator beg = std::upper_bound(bins.begin(), bins.end(), lo);
  // +1 because inclusive merge
  std::vector<int>::iterator last = std::lower_bound(bins.begin(), bins.end(), hi+1);
  bins.erase(beg, last);
}

void smoothHistoRebin(TH1* hnom, TH1* hsys, const std::vector<int>& bins, bool smooth) {
  float norm_init = hsys->Integral();
  TH1* ratio = getRatioHist(hnom, hsys, bins);
  if(smooth && ratio->GetNbinsX() > 2) {
    std::vector<float> vals(ratio->GetNbinsX() - 2);
    for(int i=2 ; i<ratio->GetNbinsX(); i++) {
      vals[i-2] = (2. * ratio->GetBinContent(i) + ratio->GetBinContent(i-1) + ratio->GetBinContent(i+1)) / 4.;
    }
    for(int i=2 ; i<ratio->GetNbinsX(); i++) {
      ratio->SetBinContent(i, vals[i-2]);
    }
  }

  for(int i=1; i<hsys->GetNbinsX()+1; i++) {
    // float smoothed = (hnom->GetBinContent(i) != 0) ? ratio->GetBinContent(i) * hnom->GetBinContent(i) : 0.;
    if(hnom->GetBinContent(i) != 0) {
      hsys->SetBinContent(i, ratio->GetBinContent(i) * hnom->GetBinContent(i));
    }
    else {
      hsys->SetBinContent(i, 0);
    }
  }
  hsys->Scale(norm_init / hsys->Integral());
  delete ratio;
  // set bin errors to 0 for systematics. Easier later when doing chi2 tests
  for(int i=0 ; i<hsys->GetNbinsX()+2; i++) {
    hsys->SetBinError(i, 0);
  }
}

int findSmallerChi2(TH1* hnom, TH1* hsys, const std::vector<int>& extrema) {
  int pos = 0;
  float minval = 99999;
  for(unsigned int i=0; i<extrema.size()-1; i++) {
    float chi2 = computeChi2(hnom, hsys, extrema[i], extrema[i+1]);
    if(chi2 < minval) {
      pos = i;
      minval = chi2;
    }
  }
  return pos;
}

float computeChi2(TH1* hnom, TH1* hsys, int beg, int end) {
  float ratio = hsys->Integral(beg, end) / hnom->Integral(beg, end);
  float chi2 = 0;
  for(int i = beg; i < end+1; i++) {
    if(hnom->GetBinContent(i) != 0) {
    float iratio = hsys->GetBinContent(i) / hnom->GetBinContent(i);
    float err = hnom->GetBinError(i) / hnom->GetBinContent(i);
    chi2 += ((iratio-ratio)/err) * ((iratio-ratio)/err);
    }
  }
  return chi2;
}

float statError(TH1* hnom, int beg, int end) { // end is excluded
  double err = 0;
  float nomInt = hnom->IntegralAndError(beg, end-1, err);
  return fabs(err / nomInt);
}

double chi2(TH1 *const h_nom, TH1 *const h_smooth) {

  double chi2 = 0.0;
  for (int i=1; i < h_nom->GetNbinsX()+1; ++i) {
    chi2 += std::pow(h_smooth->GetBinContent(i) - h_nom->GetBinContent(i),2);
  }

  return chi2;

}

}