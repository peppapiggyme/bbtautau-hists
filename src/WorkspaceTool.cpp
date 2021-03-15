#include "WorkspaceTool.h"

#include <iostream>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

void WorkspaceTool::load(const std::string &fn, const std::string &ws_name, 
    const std::string &modelConfig_name, const std::string &obsData_name)
{
    if (!fs::exists(fs::path(fn)))
    {
        throw std::runtime_error("Input does not exist");
    }
    m_fin.reset(TFile::Open(fn.c_str(), "read"));

    m_ws = (RooWorkspace*)m_fin->Get(ws_name.c_str());
    if (!m_ws)
    {
        throw std::runtime_error("Workspace does not exist");
    }

    m_mc = (ModelConfig*)m_ws->obj(modelConfig_name.c_str());
    m_data = m_ws->data(obsData_name.c_str());

    // save snapshot before any fit has been done
    RooSimultaneous* pdf = (RooSimultaneous*)m_ws->pdf("simPdf");
    if (!pdf) 
    {
        clog << "simPdf does not exist, use combPdf\n";
        pdf = (RooSimultaneous*)m_ws->pdf("combPdf");
    }

    RooArgSet *params = (RooArgSet*)pdf->getParameters(*m_data);
    if (!m_ws->loadSnapshot("snapshot_paramVals_initial"))
    {
        m_ws->saveSnapshot("snapshot_paramVals_initial", *params);
    }
    else
    {
        clog << "Snapshot \"snapshot_paramVals_initial\" already exists, will not override\n";
    }

    // m_ws->Print();

    delete params;

    if (!m_data || !m_mc) 
    {
        throw runtime_error("data or ModelConfig was not found");
    }

    if (!isSimultaneousPdfOK() || !isChannelNameOK()) 
    {
        throw runtime_error("Sanity check failed!");
    }

    getNominalValueOfNP();
    m_allNPs_afterFit.clear();

    if (m_lvl >= VerboseLevel::DEBUG) 
    {
        printModelObservables();
        printNPs();
        if (m_lvl == VerboseLevel::DEBUG) {
            printSubChannels();
        }
    }
}

bool WorkspaceTool::isSimultaneousPdfOK()
{
    const std::string pdfClass{m_mc->GetPdf()->ClassName()};
    if (pdfClass == "RooSimultaneous") 
    {
        return true;
    }

    cerr << "no Simultaneous PDF was found, class name is " << pdfClass << '\n';
    return false;
}

bool WorkspaceTool::isChannelNameOK()
{
    if (!isSimultaneousPdfOK()) return false;

    RooSimultaneous* simPdf = (RooSimultaneous*)(m_mc->GetPdf());
    RooCategory* channelCat = (RooCategory*) (&simPdf->indexCat());
    TIterator* iter = channelCat->typeIterator();
    RooCatType* tt = nullptr;
    while ((tt = (RooCatType*)iter->Next())) 
    {
        std::string channelName = tt->GetName();
        if (channelName.find("/") != std::string::npos){
        cerr << "\nOne of the channel name contain a caracter \"/\" : \n";
        cerr << "  - "  << channelName << '\n';
        cerr << "This is mis-intrepreted by roofit in the reading of the workspace.\n";
        cerr << "Please change the channel name in the xml file to run this code.\n\n";
        return false;
      }
    }

    return true;
}

void WorkspaceTool::getNominalValueOfNP()
{
    TIterator *it = m_mc->GetNuisanceParameters()->createIterator();
    RooRealVar *var = nullptr;
    if (!m_nomNPs.empty()) 
    {
        m_nomNPs.clear();
    }
    std::clog << "\nNuisance parameter names and values:\n";
    while ((var = (RooRealVar*)it->Next())){
      const double value = var->getVal();
      m_nomNPs[var->GetName()] = value;
    }
}

void WorkspaceTool::printModelObservables()
{
    RooArgSet* obses = (RooArgSet*)m_mc->GetObservables();
    std::clog << "\nList of model observables : ";
    std::clog << "\n----------------------------\n";
    for (auto& obs : *obses) 
    {
        obs->Print();
    }
}

void WorkspaceTool::printNPs()
{
    RooArgSet* nps = (RooArgSet*)m_mc->GetNuisanceParameters();
    std::clog << "\nList of nuisance parameters : ";
    std::clog << "\n----------------------------\n";
    for (auto& np : *nps) 
    {
        np->Print();
    }

}

void WorkspaceTool::printSubChannels()
{
    // ?
    RooMsgService::instance().setGlobalKillBelow(ERROR);
    if (!isSimultaneousPdfOK())
        return;
    
    RooSimultaneous* simPdf = (RooSimultaneous*)(m_mc->GetPdf());
    RooCategory* channelCat = (RooCategory*)(&simPdf->indexCat());
    TIterator* iter = channelCat->typeIterator();
    RooCatType* tt = nullptr;

    while((tt = (RooCatType*)iter->Next()) ){

      RooAbsPdf  *pdftmp  = simPdf->getPdf(tt->GetName());
      RooAbsData *datatmp = m_data->reduce(Form("%s==%s::%s", channelCat->GetName(), channelCat->GetName(), tt->GetName()));

      clog << "\nDetails on channel " << tt->GetName() << " : ";
      clog << "\n----------------------------------------------------------\n";
      datatmp->Print();
      pdftmp->Print();
      printNumberOfEvents(pdftmp);
    }
}

void WorkspaceTool::printNumberOfEvents(RooAbsPdf *pdf) 
{
    RooRealVar* firstPOI = (RooRealVar*)m_mc->GetParametersOfInterest()->first();
    constexpr double val_sym{1.};
    clog
      << Form(" %3s |","")
      << Form(" %-32s |","Nuisance Parameter")
      << Form(" %18s |","Signal events")
      << Form(" %18s |","% Change (+1sig)")
      << Form(" %18s |","% Change (-1sig)")
      << Form(" %18s |","Background events")
      << Form(" %18s |","% Change (+1sig)")
      << Form(" %18s |","% Change (-1sig)")
      << '\n';

    int inuis = -1;
    /// @note: two functions with same name different letter case..
    RooArgSet *obstmp = pdf->getObservables(*m_mc->GetObservables());
    RooRealVar *myobs = ((RooRealVar*)obstmp->first());

    RooArgSet nps = *m_mc->GetNuisanceParameters();
    TIterator* itr = nps.createIterator();
    RooRealVar* arg = nullptr;
    while ((arg = (RooRealVar*)itr->Next())) {
      if (!arg) continue;
      ++inuis;

      double val_hi = val_sym;
      double val_lo = -val_sym;
      double val_nom = arg->getVal();
      if (std::string(arg->GetName()) == "Lumi"){
        val_nom = m_ws->var("nominalLumi")->getVal();
        val_hi  = m_ws->var("nominalLumi")->getVal() * (1 + WorkspaceTool::LumiRelError);
        val_lo  = m_ws->var("nominalLumi")->getVal() * (1 - WorkspaceTool::LumiRelError);
      }
      //
      arg->setVal(val_hi);
      firstPOI->setVal(0);
      double b_hi = pdf->expectedEvents(*myobs);
      firstPOI->setVal(1);
      double s_hi = pdf->expectedEvents(*myobs) - b_hi;
      //
      arg->setVal(val_lo);
      firstPOI->setVal(0);
      double b_lo = pdf->expectedEvents(*myobs);
      firstPOI->setVal(1);
      double s_lo = pdf->expectedEvents(*myobs) - b_lo;
      //
      arg->setVal(val_nom);
      firstPOI->setVal(0);
      double b_nom = pdf->expectedEvents(*myobs);
      firstPOI->setVal(1);
      double s_nom = pdf->expectedEvents(*myobs) - b_nom;
      //
      double x_nom = s_nom ;
      double x_hi  = 0; if (s_nom) x_hi = (s_hi - s_nom) / s_nom;
      double x_lo  = 0; if (s_nom) x_lo = (s_lo - s_nom) / s_nom;
      double y_nom = b_nom ;
      double y_hi  = 0; if (b_nom) y_hi = (b_hi - b_nom) / b_nom;
      double y_lo  = 0; if (b_nom) y_lo = (b_lo - b_nom) / b_nom;

      std::clog
        << Form(" %3d |", inuis)
        << Form(" %-32s |", arg->GetName())
        << Form(" %18.2f |", x_nom)
        << Form(" %18.2f |", 100*x_hi)
        << Form(" %18.2f |", 100*x_lo)
        << Form(" %18.2f |", y_nom)
        << Form(" %18.2f |", 100*y_hi)
        << Form(" %18.2f |", 100*y_lo)
        << '\n';
    }
}

