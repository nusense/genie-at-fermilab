#include "Numerical/RandomGen.h"
#include "FluxDrivers/GSimpleNtpFlux.h"
#include "Utils/UnitUtils.h"

#include "TSystem.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <set>

using namespace std;
using namespace genie;
using namespace genie::flux;

//========================================================================
// main routine

// assumes a single input file (meta-data handling)


void dump_gsimple_file(string fname="gsimple_output.root",
                       long int nentries=10)
{
    string fluxfname(gSystem->ExpandPathName(fname.c_str()));
  flux::GSimpleNtpFlux* gsimplein = new GSimpleNtpFlux();
  gsimplein->LoadBeamSimData(fluxfname,"<<no-offset-index>>");
  gsimplein->SetEntryReuse(1); // don't reuse entries when dumping
  gsimplein->GenerateWeighted(true); // don't deweight, we want to see file entries
  gsimplein->SetUpstreamZ(-3e38);  // leave ray on flux window

  const string sepline = 
    "========================================================================";

  cout << sepline << endl << flush;
  gsimplein->PrintConfig();
  cout << sepline << endl << flush;

  GFluxI* fdriver = dynamic_cast<GFluxI*>(gsimplein);

  // if unspecified do all the entries in the input file
  if (nentries <= 0) {
    // nentries = 2147483647;
    nentries = gsimplein->GetFluxTChain()->GetEntries();
  }

  UInt_t last_metakey = 0;
  for ( long int ientry = 0; ientry < nentries; ++ientry ) {
    fdriver->GenerateNext();

    cout << *(gsimplein->GetCurrentEntry())
         << *(gsimplein->GetCurrentAux())
         << *(gsimplein->GetCurrentNuMI())
         << endl << flush;

    const genie::flux::GSimpleNtpMeta* fmeta_in = gsimplein->GetCurrentMeta();
    if ( last_metakey != fmeta_in->metakey ) {
      cout << *fmeta_in << endl << flush;
      last_metakey = fmeta_in->metakey;
    }
 
  }

  cout << "=========================== Complete " << endl;

  cout << "Generated/Dumped " << nentries << " entries"
       << endl
       << gsimplein->UsedPOTs() << " POTs " 
       << ", pulled NFluxNeutrinos " << gsimplein->NFluxNeutrinos()
       << endl;
  
}

