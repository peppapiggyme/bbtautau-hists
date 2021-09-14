#include "Utils_WS.h"

#include "TMath.h"

#include <Math/DistFunc.h>
#include <Math/QuantFuncMathCore.h>

namespace Utils_WS
{

double pvalueToSignificance(double p)
{
    // return TMath::Sqrt(2.0) * TMath::ErfInverse(1.0 - 2.0 * p);
    return TMath::NormQuantile(1.0 - p);
}

double significanceToPvalue(double z)
{
    // return 0.5 * (1.0 - TMath::Erf(z / TMath::Sqrt(2.0)));
    return 1.0 - ROOT::Math::gaussian_cdf(z);
}

double q0ToSignificance(double q0)
{
    return TMath::Sqrt(q0);
}

double nllToQ0(double nllMu0, double nllMuHat, double muHat)
{
    double q0 = (muHat < 0) ? 0 : (2 * (nllMu0 - nllMuHat));
    if (q0 < 0.) q0 = 0;
    return q0;
}


}