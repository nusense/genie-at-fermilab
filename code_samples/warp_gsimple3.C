//========================================================================
//
//  This is a script for creating a modified GSimpleNtpFlux file from
//  an original file, but modifying entries.  Assumes it will be passed
//  the name of a single file.
//
//  User must supply two functions "warp_entry()" and "warp_meta()"
//
//========================================================================
#include "Framework/Numerical/RandomGen.h"
#include "Tools/Flux/GSimpleNtpFlux.h"
#include "Framework/Utils/UnitUtils.h"

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
#include <algorithm>

using namespace std;
using namespace genie;
using namespace genie::flux;

// function prototype
double pick_energy(double, double);

bool gDoWarp = true;

//========================================================================
// this is what the USER must supply

void warp_entry(GSimpleNtpEntry* entry,
                GSimpleNtpAux*   aux,
                GSimpleNtpNuMI*  numi)
{
  if ( ! gDoWarp ) return;
  // this is a dumb weighting scheme
  //if ( entry->E > 1.5 ) {
    //// don't set absolute weight in case we started with weighted entries
    //// scale whatever is there from GSimpleNtpFlux driver
    //entry->wgt *= 0.5;
  //}
  //// change some flavors
  //if ( gRandom->Rndm() < 0.25 ) {
    //entry->pdg = ( entry->pdg > 0 ) ? 16 : -16;
  //}

  // record the 4 momentum before warp
  TLorentzVector fmb4 = TLorentzVector(entry->px, entry->py, entry->pz,
                                       entry->E);
  // sample an energy value from the 1/E distribution
  // min energy 0.01 GeV, where GENIE spline starts
  // max energy 100 GeV, by eye
  entry->E = pick_energy(0.01, 20);
  // scale momentum so that the new momentum is in the same direction as the
  // old one but obeys the energy-momentum relation
  double k = sqrt(1+(entry->E*entry->E-fmb4.E()*fmb4.E())/fmb4.Vect().Mag2());
  entry->px *= k;
  entry->py *= k;
  entry->pz *= k;

  // use all branches so there are no complaints from compiler
  static double trash = 0;
  if ( numi ) trash = numi->entryno;
  if ( aux  ) trash = aux->auxint.size();
  if ( trash == -1 ) cout << "trash was -1" << endl;
}

void warp_meta(GSimpleNtpMeta* meta, string fnamein)
{
  meta->infiles.push_back("THIS IS MDC WARPED FLUX FROM:");
  meta->infiles.push_back(fnamein);
  if ( ! gDoWarp ) {
    meta->infiles.push_back("NO ACTUAL WARPING APPLIED");
  } else {
    string msg ="USING WARP FUNCTION CODENAMED WHISKEY-TANGO-FOXTROT";
    meta->infiles.push_back(msg);
  }
}

double pick_energy(double emin, double emax)
{
  // pick a distribution of shape 1/E
  // between the values [emin, emax]
  // must have low energy cut to avoid infrared divergence
  double c = log(emin);
  double A = 1. / (log(emax) - c);
  double r = gRandom->Rndm(); // flat variate between (0,1]
  return exp(r/A + c);
}

void update_flavors(GSimpleNtpMeta* meta, std::set<int>& gFlavors);

//========================================================================
// main routine

// assumes a single input file (meta-data handling)


void warp_gsimple3(string fnameout="gsimple_output.root",
                   string fnamein="gsimple_input.root",
                   bool dowarp=true)
{

  gDoWarp = dowarp;

  string fnameinlong(gSystem->ExpandPathName(fnamein.c_str()));
  TFile* fin = TFile::Open(fnameinlong.c_str(),"READONLY");
  TTree* etreein = (TTree*)fin->Get("flux");
  TTree* mtreein = (TTree*)fin->Get("meta");
  genie::flux::GSimpleNtpEntry* entry_in = new genie::flux::GSimpleNtpEntry;
  genie::flux::GSimpleNtpNuMI*  numi_in  = new genie::flux::GSimpleNtpNuMI;
  genie::flux::GSimpleNtpAux*   aux_in   = new genie::flux::GSimpleNtpAux;
  genie::flux::GSimpleNtpMeta*  meta_in  = new genie::flux::GSimpleNtpMeta;

  long int nentries = etreein->GetEntries();

  int sba_status[4] = { -999, -999, -999, -999 };
  sba_status[0] = etreein->SetBranchAddress("entry",&entry_in);
  sba_status[1] = etreein->SetBranchAddress("numi",&numi_in);
  sba_status[2] = etreein->SetBranchAddress("aux",&aux_in);
  sba_status[3] = mtreein->SetBranchAddress("meta",&meta_in);
  cout << "SetBranchAddress results "
       << sba_status[0] << ","
       << sba_status[1] << ","
       << sba_status[2] << ","
       << sba_status[3]
       << endl;
  bool donumi = ( sba_status[1] == 0 );
  bool doaux  = ( sba_status[2] == 0 );

  int nindices = mtreein->BuildIndex("metakey"); // tie metadata to entry
  cout << "saw " << nindices << " metakey indices" << endl;

  cout << "Creating:    " << fnameout << endl;
  cout << "Input file:  " << fnamein << endl;
  cout << "Branches:    entry "
       << ( (doaux) ? "aux ":"" )
       << ( (donumi) ? "numi ":"" )
       << endl;
  if ( ! gDoWarp ) cout << "++++++++ NO ACTUAL WARP APPLIED" << endl;

  TFile* fout = TFile::Open(fnameout.c_str(),"RECREATE");
  TTree* fluxntp = new TTree("flux","a simple flux n-tuple");
  TTree* metantp = new TTree("meta","metadata for flux n-tuple");
  genie::flux::GSimpleNtpEntry* fentry = new genie::flux::GSimpleNtpEntry;
  genie::flux::GSimpleNtpAux*   faux   = new genie::flux::GSimpleNtpAux;
  genie::flux::GSimpleNtpNuMI*  fnumi  = new genie::flux::GSimpleNtpNuMI;
  genie::flux::GSimpleNtpMeta*  fmeta  = new genie::flux::GSimpleNtpMeta;

  fluxntp->Branch("entry",&fentry);
  if ( doaux  ) fluxntp->Branch("aux",&faux);
  if ( donumi ) fluxntp->Branch("numi",&fnumi);
  metantp->Branch("meta",&fmeta);

  cout << "=========================== Start " << endl;

  TStopwatch sw;
  sw.Start();

  const double large = 1.0e300;
  double minwgt = large, maxwgt = -large, maxenergy = -large;
  std::set<int> gFlavors;

  for ( long int ientry = 0; ientry < nentries; ++ientry ) {

    // reset what's been read in anticipation of a new entry
    entry_in->Reset();
    aux_in->Reset();
    numi_in->Reset();

    // read the next entry, get metadata if it's different
    //int nbytes =
    etreein->GetEntry(ientry);
    UInt_t metakey = entry_in->metakey;
    if ( fmeta->metakey != metakey ) {
      // UInt_t oldkey = meta_in->metakey;
      int nbmeta = mtreein->GetEntryWithIndex(metakey);
      cout << "on entry " << ientry << " fetched metadata w/ key "
           << metakey << " read " << nbmeta << " bytes" << endl;
    }

    // reset what we're writing
    fentry->Reset();
    fnumi->Reset();
    faux->Reset();

    // copy read in objects to output objects
    *fentry = *entry_in;
    *fnumi  = *numi_in;
    *faux   = *aux_in;

    // mess with the values to your hearts content
    warp_entry(fentry,faux,fnumi);
    // keep track of any weights for the meta data
    if ( fentry->wgt < minwgt ) minwgt = fentry->wgt;
    if ( fentry->wgt > maxwgt ) maxwgt = fentry->wgt;
    // user might have changed an energy larger than any in original file
    if ( fentry->E > maxenergy ) maxenergy = fentry->E;
    // user might have changed a flavor keep track of all we see
    gFlavors.insert(fentry->pdg);

    // process currently held metadata after transition to metadata key
    if ( fmeta->metakey != meta_in->metakey ) {
      // new meta data found
      cout << "new meta found " << fmeta->metakey << "/" << meta_in->metakey << endl;
      if ( fmeta->metakey != 0 ) {
        // have metadata that needs adjustment and writing
        // before processing new metadata
        warp_meta(fmeta,fnamein);
        fmeta->minWgt = minwgt;
        fmeta->maxWgt = maxwgt;
        fmeta->maxEnergy = maxenergy;
        update_flavors(fmeta,gFlavors);
        gFlavors.clear();
        cout << "metantp->Fill() " << *fmeta << endl;
        metantp->Fill();
        // next meta-data restarts weight range (and energy max)
        minwgt    =  large;
        maxwgt    = -large;
        maxenergy = -large;
      }
      // get a copy of metadata for accumulating adjustments
      *fmeta = *meta_in;
      //cout << " copy meta_in " << *meta_in << endl << *fmeta << endl;
    }

    fluxntp->Fill();
  }

  // write last set of meta-data after all is done
  if ( fmeta->metakey != 0 ) {
    cout << "final meta flush " << fmeta->metakey << endl;
    // have metadata that needs adjustment and writing
    warp_meta(fmeta,fnamein);
    fmeta->minWgt = minwgt;
    fmeta->maxWgt = maxwgt;
    fmeta->maxEnergy = maxenergy;
    update_flavors(fmeta,gFlavors);
    gFlavors.clear();
    cout << "metantp->Fill() " << *fmeta << endl;
    metantp->Fill();
  }

  cout << "=========================== Complete " << endl;

  sw.Stop();
  cout << "Generated " << nentries
       << endl
       << "Time to generate: " << endl;
  sw.Print();

  fout->cd();

  // write ntuples out
  fluxntp->Write();
  metantp->Write();
  fout->Close();

  cout << endl << endl;
}

void update_flavors(GSimpleNtpMeta* meta, std::set<int>& gFlavors)
{
  // add any new flavors here, if necessary
  std::set<int>::const_iterator flitr = gFlavors.begin();
  std::vector<Int_t>& flvlist = meta->pdglist;
  for ( ; flitr != gFlavors.end(); ++flitr ) {
    int  seen_pdg = *flitr;
    std::vector<Int_t>::iterator knwitr =
      std::find(flvlist.begin(),flvlist.end(),seen_pdg);
    bool known_pdg = ( knwitr != flvlist.end() );
    if ( ! known_pdg ) meta->pdglist.push_back(seen_pdg);
  }
}
