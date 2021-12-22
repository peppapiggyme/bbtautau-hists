#include "ExamplesInclude_WS.h"
#include "CommonInclude.h"

#include "TRandom.h"
#include "TRandom3.h"

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        cerr << "Usage:\n gen-pd-zcr <filename> <outname> <nToys>\n";
    }

    TFile* f = TFile::Open(argv[1]);

    vector<TH1F*> histograms;
    histograms.push_back((TH1F*)f->Get("Zcc")); histograms.back()->Scale(1.2308196024245677);
    histograms.push_back((TH1F*)f->Get("VBFHtautau"));
    histograms.push_back((TH1F*)f->Get("Zbl")); histograms.back()->Scale(0.8657669838456886);
    histograms.push_back((TH1F*)f->Get("Zcl")); histograms.back()->Scale(0.8657669838456886);
    histograms.push_back((TH1F*)f->Get("VHtautau"));
    //histograms.push_back((TH1F*)f->Get("data"));
    histograms.push_back((TH1F*)f->Get("DY"));
    histograms.push_back((TH1F*)f->Get("ttH"));
    histograms.push_back((TH1F*)f->Get("Zl")); histograms.back()->Scale(0.8657669838456886);
    histograms.push_back((TH1F*)f->Get("DYtt"));
    histograms.push_back((TH1F*)f->Get("ZZ"));
    // histograms.push_back((TH1F*)f->Get("hhWWbb"));
    histograms.push_back((TH1F*)f->Get("Zbb")); histograms.back()->Scale(1.2308196024245677);
    histograms.push_back((TH1F*)f->Get("stopWt"));
    histograms.push_back((TH1F*)f->Get("Zbc")); histograms.back()->Scale(1.2308196024245677);
    histograms.push_back((TH1F*)f->Get("W"));
    histograms.push_back((TH1F*)f->Get("HtoVV"));
    // histograms.push_back((TH1F*)f->Get("hhZZbb"));
    histograms.push_back((TH1F*)f->Get("ggHtautau")); histograms.back()->Scale(1.2153682446311185);
    histograms.push_back((TH1F*)f->Get("ttbar")); histograms.back()->Scale(0.8662766660691217);

    for (int i = 1; i < (int)histograms.size(); i++)
    {
        histograms[0]->Add(histograms[i]);
    }

    TRandom3 r;
    TFile* fout = TFile::Open(argv[2], "recreate");

    int nToys = stoi(argv[3]);
    for (int i = 0; i < nToys; i++)
    {
        if (i % 1000 == 0 ) cout << "Generated " << i << "toys\n";
        string histName("PseudoData");
        histName += to_string(i);
        TH1F *h = (TH1F*)histograms[0]->Clone();
        h->SetName(histName.c_str());
        h->SetTitle(histName.c_str());
        for (int ibin = 1; ibin <= h->GetNbinsX(); ibin++)
        {
            h->SetBinContent(ibin, r.Poisson(histograms[0]->GetBinContent(ibin)));
        }
        h->SetDirectory(fout);
        h->Write();
        delete h;
        h = nullptr;
    }
}