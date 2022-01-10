
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooBukinPdf.h"
#include "RooChi2Var.h"
#include "RooMinimizer.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
using namespace RooFit;

/**
 * Usage:
 * root -b -l -q mBBreso.C
 */

void mBBreso()
{
    // Declare observable
    RooRealVar x("mBB","m_{bb}",20000,220000) ;
    RooRealVar y("mBB_raw","m_{bb}",20000,220000) ;
    RooRealVar w("weight","weight",-10,10) ;
    x.setBins(50) ;
    y.setBins(50) ;

    TFile* f = TFile::Open("/scratchfs/atlas/bowenzhang/bbtautau-hists/data/MVATree_NonRes.root");
    TTree* tree = (TTree*)f->Get("mytree");
    // tree->Scan("mBB:mBB_raw");

    RooDataSet* data = new RooDataSet("data", "data", RooArgSet(x,y,w), Import(*tree));
    data->Print("V") ;

    RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,w.GetName()) ;
    wdata.Print("V") ;

    // for mBB corrected
    RooRealVar Xp_x("Xp_x","Xp_x", 120000.0, 20000.0, 220000.0) ;
    RooRealVar sigp_x("sigp_x","sigp_x",10000.0, 5000.0, 50000.0) ;
    RooRealVar xi_x("xi_x","xi_x",0.0, -1.0, 1.0) ;
    RooRealVar rho1_x("rho1_x","rho1_x",-0.1, -1.0, 0.0) ;
    RooRealVar rho2_x("rho2_x","rho2_x",0.1, 0.0, 1.0) ;
    RooBukinPdf bukin_x("bukin_x","bukin_x", x, Xp_x, sigp_x, xi_x, rho1_x, rho2_x) ;
    // for mBB raw
    RooRealVar Xp_y("Xp_y","Xp_y", 120000.0, 20000.0, 220000.0) ;
    RooRealVar sigp_y("sigp_y","sigp_y",10000.0, 5000.0, 50000.0) ;
    RooRealVar xi_y("xi_y","xi_y",0.0, -1.0, 1.0) ;
    RooRealVar rho1_y("rho1_y","rho1_y",-0.1, -1.0, 0.0) ;
    RooRealVar rho2_y("rho2_y","rho2_y",0.1, 0.0, 1.0) ;
    RooBukinPdf bukin_y("bukin_y","bukin_y", y, Xp_y, sigp_y, xi_y, rho1_y, rho2_y) ;


    RooFitResult* r_ml_wgt_corr_x = bukin_x.fitTo(wdata,Save(),SumW2Error(kTRUE)) ;
    RooFitResult* r_ml_wgt_corr_y = bukin_y.fitTo(wdata,Save(),SumW2Error(kTRUE)) ;

    RooPlot* xframe = x.frame(Title("")) ;
    RooPlot* yframe = y.frame(Title("")) ;

    wdata.plotOn(xframe, DataError(RooAbsData::SumW2), MarkerColor(TColor::GetColor(196, 110, 57)), LineColor(TColor::GetColor(196, 110, 57))) ;
    wdata.plotOn(yframe, DataError(RooAbsData::SumW2), MarkerColor(TColor::GetColor(73, 163, 166)), LineColor(TColor::GetColor(73, 163, 166))) ;

    bukin_x.plotOn(xframe, LineColor(TColor::GetColor(120, 46, 0)), MarkerColor(TColor::GetColor(120, 46, 0))) ;
    bukin_y.plotOn(yframe, LineColor(TColor::GetColor(19, 83, 84)), MarkerColor(TColor::GetColor(19, 83, 84))) ;

    TH1F* dummy_data_x = new TH1F();
    TH1F* dummy_data_y = new TH1F();
    TH1F* dummy_x = new TH1F();
    TH1F* dummy_y = new TH1F();

    dummy_data_x->SetLineColor(TColor::GetColor(196, 110, 57));
    dummy_data_x->SetMarkerColor(TColor::GetColor(196, 110, 57));
    dummy_data_y->SetLineColor(TColor::GetColor(73, 163, 166));
    dummy_data_y->SetMarkerColor(TColor::GetColor(73, 163, 166));

    dummy_x->SetLineColor(TColor::GetColor(120, 46, 0));
    dummy_x->SetMarkerColor(TColor::GetColor(120, 46, 0));
    dummy_y->SetLineColor(TColor::GetColor(19, 83, 84));
    dummy_y->SetMarkerColor(TColor::GetColor(19, 83, 84));

    dummy_x->SetLineWidth(2);
    dummy_y->SetLineWidth(2);
    dummy_data_x->SetLineWidth(2);
    dummy_data_y->SetLineWidth(2);

    r_ml_wgt_corr_x->Print() ;
    r_ml_wgt_corr_y->Print() ;

    gROOT->SetStyle("ATLAS");
    gStyle->SetErrorX(0.5);

    TCanvas* c = new TCanvas("c1","c1",1200,900) ;
    gPad->SetLeftMargin(0.15) ; 
    xframe->GetYaxis()->SetRangeUser(0, 0.70) ;
    xframe->GetYaxis()->SetTitleOffset(1.3) ;
    xframe->GetYaxis()->SetTitle("Event") ;
    xframe->GetXaxis()->SetTitle("m_{bb} [GeV]") ;
    xframe->SetTitle("") ;
    xframe->Draw("SAME") ;
    yframe->Draw("SAME") ;

    TLegend* legend = new TLegend(0.20, 0.82, 0.80, 0.92);
    legend->SetNColumns(2);
    legend->SetTextFont(42);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.036);
    legend->SetTextAlign(12);

    legend->AddEntry(dummy_data_y, "Original", "lep");
    legend->AddEntry(dummy_y, "Peak: 118 GeV, FWHM: 15.6 GeV", "l");
    legend->AddEntry(dummy_data_x, "Corrected ", "lep");
    legend->AddEntry(dummy_x, "Peak: 123 GeV, FWHM: 14.0 GeV", "l");

    legend->Draw();

    c->SaveAs("/scratchfs/atlas/bowenzhang/bbtautau-hists/output/mbbreso.pdf");

    delete c;
    delete data;
}