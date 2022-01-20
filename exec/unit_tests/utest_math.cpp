#define CATCH_CONFIG_MAIN

#include "catch.h"
#include "Utils.h"
#include "CommonInclude.h"

using namespace std;

double GaussianPdf(double p)
{
    return ROOT::Math::gaussian_cdf(ROOT::Math::gaussian_quantile(p, 1.), 1.);
}

TEST_CASE("Gaussain p.d.f.")
{
    REQUIRE(TMath::Abs(GaussianPdf(0.1) - 0.1) < 1e-6);
    REQUIRE(TMath::Abs(GaussianPdf(0.2) - 0.2) < 1e-6);
    REQUIRE(TMath::Abs(GaussianPdf(0.3) - 0.3) < 1e-6);
    REQUIRE(TMath::Abs(GaussianPdf(0.4) - 0.4) < 1e-6);
    REQUIRE(TMath::Abs(GaussianPdf(0.5) - 0.5) < 1e-6);
    REQUIRE(TMath::Abs(GaussianPdf(0.6) - 0.6) < 1e-6);
    REQUIRE(TMath::Abs(GaussianPdf(0.7) - 0.7) < 1e-6);
    REQUIRE(TMath::Abs(GaussianPdf(0.8) - 0.8) < 1e-6);
    REQUIRE(TMath::Abs(GaussianPdf(0.9) - 0.9) < 1e-6);
}
