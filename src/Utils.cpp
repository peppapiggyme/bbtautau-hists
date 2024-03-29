#include "Utils.h"

using namespace std;

namespace Utils
{

void histAssign(TH1* h, ProcessInfo *p, RegionInfo* r, VariableInfo* v)
{
    TH1* h_clone = (TH1*)h->Clone();
    h_clone->SetDirectory(0);
    p->histogram = std::move(h_clone);
    p->current_region = r;
    p->current_variable = v;
}

void histAssignSyst(TH1* h, ProcessInfo *p, const std::string& systname)
{
    TH1* h_clone = (TH1*)h->Clone();
    h_clone->SetDirectory(0);
    p->systematic_histograms[systname] = std::move(h_clone);
}

string histString(const ProcessInfo* p, const RegionInfo* r, const VariableInfo* v, const NameConvention nc)
{
    switch(nc)
    {
        case NameConvention::CxAODReader:
            return p->name + "_" + r->name + "_" + v->name;
        case NameConvention::WSMaker:
            return p->name;
        default:
            return p->name + "_" + r->name + "_" + v->name;
    }
}

string systString(const SystematicInfo* s, const NameConvention nc)
{
    (void)nc;
    return s->name_tex;
}

string histStringSyst(const ProcessInfo* p, const RegionInfo* r, const VariableInfo* v, const SystematicInfo* s, const NameConvention nc)
{
    return histString(p, r, v, nc) + "_Sys" + s->name;
}

string systStringShort(const string& sSyst, const NameConvention nc)
{
    (void)nc;
    const string gamma("gamma_");
    const string alpha("alpha_");
    const string alphaSys("alpha_Sys");
    const string mcStatHadHadSR("Y2015_DLLOS_T2_SpcTauHH");
    /// @todo not tested
    /// maybe this is SLT, LTT
    const string mcStatLepHadSLTSR("Y2015_DLLOS_T2_SpcTauLH_Y2015_LTT0");
    const string mcStatLepHadLTTSR("Y2015_DLLOS_T2_SpcTauLH_Y2015_LTT1");
    const string mcStatZCR("Y2015_DZllbbCR_T2_L2");

    if (sSyst.find(gamma) != std::string::npos)
    {
        string sShort = sSyst.substr(gamma.length());
        string sBin = sShort.substr(sShort.find("_bin_"));
        if (sShort.find(mcStatHadHadSR) != std::string::npos)
        {
            return "HadHad_SR_MVAScore" + sBin;
        }
        else if (sShort.find(mcStatLepHadSLTSR) != std::string::npos)
        {
            return "LepHad_SLT_SR_MVAScore" + sBin;
        }
        else if (sShort.find(mcStatLepHadLTTSR) != std::string::npos)
        {
            return "LepHad_LTT_SR_MVAScore" + sBin;
        }
        else if (sShort.find(mcStatZCR) != std::string::npos)
        {
            return "ZCR_MLL" + sBin;
        }
    }
    else if (sSyst.find(alpha) != std::string::npos)
    {
        if (sSyst.find(alphaSys) != std::string::npos)
        {
            return sSyst.substr(alphaSys.length());
        }

        return sSyst.substr(alpha.length());
    }

    return sSyst;
}

string signalTypeName(const string& sSigName)
{
    string name_SMHH("hhttbb");
    string name_2HDM("Hhhbbtautau");
    string sName(sSigName);

    if (sSigName == name_SMHH)
    {
        sName = "SM";
    }
    else if (sSigName.find(name_2HDM) != string::npos)
    {
        sName = "2HDM" + sSigName.substr(string(name_2HDM).length());
    }

    return sName;
}

void properties_copy(TH1* h1, TH1* h2)
{
    h1->SetLineColor(h2->GetLineColor());
    h1->SetLineWidth(h2->GetLineWidth());
    h1->SetLineStyle(h2->GetLineStyle());
    h1->SetFillColor(h2->GetFillColor());
    h1->SetFillStyle(h2->GetFillStyle());
    h1->SetFillColor(h2->GetFillColor());
    h1->SetMarkerColor(h2->GetMarkerColor());
    h1->SetMarkerStyle(h2->GetMarkerStyle());
    h1->SetMarkerSize(h2->GetMarkerSize());
}

}
