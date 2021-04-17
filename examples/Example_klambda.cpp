#include "ExamplesInclude.h"

#include "Config.h"
#include "Utils.h"
#include "CompTool.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>

using std::cout; 
using std::endl;
using std::clog;

void test_hadhad_klambda(const std::string& filename)
{
    BasicInfo* b = new BasicInfo("#sqrt{s} = 13 TeV", "L = 139 fb^{-1}");

    Regions* rs = new Regions();
    rs->add("2tag2pjet_0ptv_LL_OS",     "2 b-tag, 2 loose #tau, OS",        eRegionType::SR);

    Variables* vs = new Variables();
    vs->add("mBB",                  "m_{BB} [GeV]",                                      10);
    vs->add("mMMC",                 "m_{#tau#tau} (MMC) [GeV]",                          12);
    vs->add("mHH",                  "m_{HH} [GeV]",                                      12);
    vs->add("mHHScaled",            "m_{HH} (Scaled) [GeV]",                             12);
    vs->add("dRBB",                 "#Delta R(B,B)",                                     4 );
    vs->add("pTBB",                 "p_{T}^{b,b}",                                       10);
    vs->add("dRTauTau",             "#Delta R(#tau_{had},#tau_{had})",                   4 );
    vs->add("pTTauTau",             "p_{T}^{#tau,#tau}",                                 10);
    vs->add("dPhiBBTauTau",         "#Delta #phi{bb,#tau#tau}",                          4 );
    vs->add("dRBBTauTau",           "#Delta R{bb,#tau#tau}",                             5 );
    vs->add("pTBalance",            "p_{T}^{b,b} / p_{T}^{#tau,#tau}",                   5 );
    vs->add("MET",                  "E_{T}^{miss} [GeV]",                                5 );
    vs->add("nJets",                "Number of Jets",                                    1 );
    vs->add("Yield",                "Event Yield",                                       1 );
    vs->add("metSig",               "MET significance",                                  2 );
    vs->add("metSigPU",             "MET significance (pile-up)",                        2 );
    vs->add("Jet0Pt",               "Leading b-jet p_{T} (before corr.) [GeV]",          10);
    vs->add("Jet1Pt",               "Sub-leading b-jet p_{T} (before corr.) [GeV]",      10);
    vs->add("pTB0",                 "Leading b-jet p_{T} [GeV]",                         10);
    vs->add("pTB1",                 "Sub-leading b-jet p_{T} [GeV]",                     10);
    vs->add("LeadJetPt",            "Leading jet p_{T} [GeV]",                           10);
    vs->add("SubleadJetPt",         "Sub-leading jet p_{T} [GeV]",                       10);
    vs->add("LeadJetEta",           "Leading jet #eta",                                  2 );
    vs->add("SubleadJetEta",        "Sub-leading jet #eta",                              2 );
    vs->add("Jet0Eta",              "Leading b-jet #eta (before corr.)",                 2 );
    vs->add("Jet1Eta",              "Sub-leading b-jet #eta (before corr.)",             2 );
    vs->add("Tau0Pt",               "Leading #tau_{had} p_{T} [GeV]",                    10);
    vs->add("Tau1Pt",               "Sub-leading #tau_{had} p_{T} [GeV]",                10);
    vs->add("Tau0Eta",              "Leading #tau_{had} #eta",                           2 );
    vs->add("Tau1Eta",              "Sub-leading #tau_{had} #eta",                       2 );
    vs->add("Tau0Ntrk",             "Leading #tau_{had} N_{tracks}",                     1 );
    vs->add("Tau1Ntrk",             "Sub-leading #tau_{had} N_{tracks}",                 1 );
    vs->add("dPhiTauTau",           "#Delta #phi (#tau_{had},#tau_{had})",               4 );

    CompInfo* info = new CompInfo();
    info->ratio_high = 1.38;
    info->ratio_low = 0.62;

    for (VariableInfo* v : *(vs->content()))
    {
        {
            Processes* ps = new Processes();
            ps->assign_norm_factors = false;
            ps->add("hhttbbKL10p0",         "HH (#kappa_{#lambda}=10) generated",   eProcessType::SIG,  eProcess::HHKL10,       "HH (#kappa_{#lambda}=10) generated",  kBlue+1);
            ps->add("hhttbbKL10p0from1p0",  "HH (#kappa_{#lambda}=10) reweighted",  eProcessType::SIG,  eProcess::HHKLXFROM1,   "HH (#kappa_{#lambda}=10) reweighted", kRed+1);

            Config* c = new Config(b, ps, rs, vs);
            c->load(filename, "Preselection");
            info->parameter = "from1p0";
            c->updateHistogramPtr(rs->content()->front(), v);
            CompTool* ct = new CompTool(info);
            ct->output_path = "/tmp/bowenzhang/hhttbbKL10p0from1p0";
            if (ct->check(c))
            {
                ct->rebin(c, eRebinOption::N_Rebin);
                ct->paint(c);
                ct->run(c);
            }
            else 
            {
                clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
            }

            delete ps;
            delete ct;
            delete c;
        }

        {
            Processes* ps = new Processes();
            ps->assign_norm_factors = false;
            ps->add("hhttbbKL1p0",          "HH (#kappa_{#lambda}=1) generated",   eProcessType::SIG,  eProcess::HHKL10,        "HH (#kappa_{#lambda}=10) generated",  kBlue+1);
            ps->add("hhttbbKL1p0from10p0",  "HH (#kappa_{#lambda}=1) reweighted",  eProcessType::SIG,  eProcess::HHKLXFROM10,   "HH (#kappa_{#lambda}=10) reweighted", kRed+1);

            Config* c = new Config(b, ps, rs, vs);
            c->load("/scratchfs/atlas/bowenzhang/CxAODReaderSemiBoosted/run/hist-klambda-v2.root", "Preselection");
            info->parameter = "from10p0";
            c->updateHistogramPtr(rs->content()->front(), v);
            CompTool* ct = new CompTool(info);
            ct->output_path = "/tmp/bowenzhang/hhttbbKL1p0from10p0";
            if (ct->check(c))
            {
                ct->rebin(c, eRebinOption::N_Rebin);
                ct->paint(c);
                ct->run(c);
            }
            else 
            {
                clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
            }

            delete ps;
            delete ct;
            delete c;
        }

        {
            for (const std::string& klambda : {"n20p0", "n10p0", "n5p0", "0p0", 
                                               "3p0", "5p0", "20p0"})
            {
                auto nice_kl = [&klambda]() {
                    std::string ret = klambda;
                    if (ret[0] == 'n') { ret[0] = '-'; }
                    ret[ret.length()-2] = '.'; 
                    return ret; };
                Processes* ps = new Processes();
                ps->assign_norm_factors = false;
                ps->add("hhttbbKL"+klambda+"from1p0",   
                        "HH (#kappa_{#lambda}="+nice_kl()+") from #kappa_{#lambda}=1",   
                        eProcessType::SIG,  eProcess::HHKLXFROM1,  
                        "HH (#kappa_{#lambda}="+nice_kl()+") from #kappa_{#lambda}=1",  kBlue+1);
                ps->add("hhttbbKL"+klambda+"from10p0",  
                        "HH (#kappa_{#lambda}="+nice_kl()+") from #kappa_{#lambda}=10",  
                        eProcessType::SIG,  eProcess::HHKLXFROM10, 
                        "HH (#kappa_{#lambda}="+nice_kl()+") from #kappa_{#lambda}=10", kRed+1);

                Config* c = new Config(b, ps, rs, vs);
                c->load("/scratchfs/atlas/bowenzhang/CxAODReaderSemiBoosted/run/hist-klambda-v2.root", "Preselection");
                info->parameter = "to"+klambda;
                c->updateHistogramPtr(rs->content()->front(), v);
                CompTool* ct = new CompTool(info);
                ct->output_path = "/tmp/bowenzhang/hhttbbKLxfrom1and10";
                if (ct->check(c))
                {
                    ct->rebin(c, eRebinOption::N_Rebin);
                    ct->paint(c);
                    ct->run(c);
                }
                else 
                {
                    clog << "Can not draw " << c->current_region->name << " " << c->current_variable->name << '\n';
                }

                delete ps;
                delete ct;
                delete c;
            }
        }
    }

    delete b;
    delete rs;
    delete vs;
    delete info;

}