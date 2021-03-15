#include "ExamplesInclude.h"
#include "TROOT.h"

#include <iostream>

using namespace std;

int main(/*int argc, char *argv[]*/)
{
    gROOT->SetBatch(1);
    
    cout << "Hello!" << endl;

    // test_hadhad("/scratchfs/atlas/zhangbw/CxAODReaderSemiBoosted/run/hist-klambda-v2.root");

    // test_hadhad_klambda("/scratchfs/atlas/zhangbw/CxAODReaderSemiBoosted/run/hist-klambda-v2.root");

    // test_hadhad_yield("/scratchfs/atlas/zhangbw/CxAODReaderSemiBoosted/run/hist-klambda-v2.root");

    // test_hadhad_ZtautauMG("/scratchfs/atlas/zhangbw/CxAODReaderSemiBoosted/run/hist-Ztautau-v1.root");

    // test_hadhad_WtDS("/scratchfs/atlas/zhangbw/CxAODReaderForSystStudy/run/hist-top-v1.root");

    // test_hadhad_WtGen("/scratchfs/atlas/zhangbw/CxAODReaderForSystStudy/run/hist-top-v1.root");

    // test_hadhad_WtOTF("/scratchfs/atlas/zhangbw/CxAODReaderForSystStudy/run/hist-top-v1.root");
    
    // test_hadhad_NonRes("/scratchfs/atlas/zhangbw/CxAODReaderSemiBoosted/run/hist-NonRes-v0.root");

    // hadhadsr_v15("/scratchfs/atlas/zhangbw/CxAODReaderSemiBoosted/run/hist-HadHad-v1.root");

    // hadhadsr_v15_yield("/scratchfs/atlas/zhangbw/CxAODReaderSemiBoosted/run/hist-HadHadWSI-v1.all.root");

    cout << "Ciao!" << endl;

    return 0;
}
