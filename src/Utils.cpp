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

string histString(const ProcessInfo* p, const RegionInfo* r, const VariableInfo* v)
{
    return p->name + "_" + r->name + "_" + v->name;
}

string systString(const SystematicInfo* s)
{
    return s->name;
}

string histStringSyst(const ProcessInfo* p, const RegionInfo* r, const VariableInfo* v, const SystematicInfo* s)
{
    return histString(p, r, v) + "_Sys" + s->name;
}

string systStringShort(const string& sSyst)
{
    const string gamma("gamma_");
    const string alpha("alpha_");
    const string alphaSys("alpha_Sys");
    const string mcStatHadHadSR("Y2015_DLLOS_T2_SpcTauHH");
    /// @todo not tested
    /// maybe this is SLT, LTT
    const string mcStatLepHadSR("Y2015_DLLOS_T2_SpcTauLH");
    const string mcStatZCR("Y2015_DZllbbCR_T2_L2");

    if (sSyst.find(gamma) != std::string::npos)
    {
        string sShort = sSyst.substr(gamma.length());
        string sBin = sShort.substr(sShort.find("_bin_"));
        if (sShort.find(mcStatHadHadSR) != std::string::npos)
        {
            return "HadHadSR_MVAScore" + sBin;
        }
        else if (sShort.find(mcStatLepHadSR) != std::string::npos)
        {
            return "LepHadSR_MVAScore" + sBin;
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

}
