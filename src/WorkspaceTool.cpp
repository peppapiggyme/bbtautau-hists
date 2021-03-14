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

    m_ws->Print();

    delete params;
}

