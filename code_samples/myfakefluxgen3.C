// a sample script for creating/filling a GSimpleNtpFlux file
// use:   $ genie -b -q myfakefluxgen.C+
//
// R. Hatcher <rhatcher@fnal.gov> 2014-05-07

// only for R-3_
#include "Framework/Conventions/GVersion.h"
#include "Framework/Numerical/RandomGen.h"
#include "Tools/Flux/GSimpleNtpFlux.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/AppInit.h"

#include "TSystem.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TFile.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <set>
#include <stdlib.h>  // for strtol, putenv, unsetenv

// my own laziness:
using namespace std;
using namespace genie;
using namespace genie::flux;

// forward declaration
void getInfo(TLorentzVector& p4u, TLorentzVector& x4u,
             int& pdg, double& wgt, double& decayDist);

/// main routine
void myfakefluxgen3(std::string fnameout="myfakeflux.gsimple.root",
                   long int nentries = 100000,
                   long int seed = 12345)
{
  /// Generate a GSimpleNtpFlux file
  /// ...out of nothing

  genie::utils::app_init::RandGen(seed);

  /// create the file, trees and branches
  TFile* file = TFile::Open(fnameout.c_str(),"RECREATE");
  TTree* fluxntp = new TTree("flux","a simple flux n-tuple");
  TTree* metantp = new TTree("meta","metadata for flux n-tuple");
  genie::flux::GSimpleNtpEntry* fentry = new genie::flux::GSimpleNtpEntry;
  genie::flux::GSimpleNtpMeta*  fmeta  = new genie::flux::GSimpleNtpMeta;
  fluxntp->Branch("entry",&fentry);
  metantp->Branch("meta",&fmeta);

  TLorentzVector p4u, x4u;
  double wgt, decayDist;
  int pdg;

  double POTs=0, minwgt=1.0e30, maxwgt=0, maxe=0;

  /// use metakey to associate entries w/ metadata
  UInt_t metakey = TString::Hash(fnameout.c_str(),strlen(fnameout.c_str()));

  /// generate & store entries
  long int ngen = 0;
  while ( ngen < nentries ) {
    fentry->Reset();
    ngen++;

    // make up an entry ...
    // p4u and x4u are the 4-momentum and 4-position in the user
    // (i.e. detector) frame
    getInfo(p4u,x4u,pdg,wgt,decayDist);

    // in our fake setup each entry represents
    //  3.1415 protons-on-target on average
    POTs += 3.1415;

    // fill the entry in the tree
    fentry->metakey = metakey;
    fentry->pdg     = pdg;
    fentry->wgt     = wgt;

    // flux positions are expected to be in SI units (i.e. meters)
    fentry->vtxx    = x4u.X();
    fentry->vtxy    = x4u.Y();
    fentry->vtxz    = x4u.Z();
    fentry->dist    = decayDist;

    fentry->px      = p4u.Px();
    fentry->py      = p4u.Py();
    fentry->pz      = p4u.Pz();
    fentry->E       = p4u.E();

    fluxntp->Fill();

    // accumulate meta data
    fmeta->AddFlavor(fentry->pdg);
    minwgt = TMath::Min(minwgt,fentry->wgt);
    maxwgt = TMath::Max(maxwgt,fentry->wgt);
    maxe   = TMath::Max(maxe,fentry->E);

  }

  fmeta->maxEnergy = maxe;
  fmeta->minWgt    = minwgt;
  fmeta->maxWgt    = maxwgt;
  fmeta->protons   = POTs;

  fmeta->seed    = seed;
  fmeta->metakey = metakey;
  metantp->Fill();

  // flush/write ntuples out
  fluxntp->Write();
  metantp->Write();
  file->Close();

}

// forward declaration
void getInfo(TLorentzVector& p4u, TLorentzVector& x4u,
             int& pdg, double& wgt, double& decayDist)
{

  // GENIE's pseudo-random number central station
  RandomGen* rg  = RandomGen::Instance();

  // here's the pseudo-random # engine we'll use
  // make sure to use a _reference_ here
  // otherwise you're creating a copy and it will not reflect changing
  // internal state (and thus you'll get the same values each call to getInfo)
  // Ah, CINT barfs at this use of a reference ... don't know why
  // ACLiC compiles it fine ... whatever!  just use this expression everywhere
  // TRandom3& rndm = rg->RndFlux();

  // pick a flavor
  double r = rg->RndFlux().Uniform(1.0);
  pdg = 14;  // nu_mu
  if ( r > 0.75 ) pdg = -14;  // nu_mu_bar
  if ( r > 0.95 ) pdg =  12;  // nu_e
  if ( r > 0.99 ) pdg = -12;  // nu_e_bar

  // energy spectrum ... how about Landau
  // make sure it sensible, though ... no zero or negative energies
  // nothing crazy high either (cut it at 100GeV)
  double e = 0;
  while ( e <= 0 || e > 100 ) {
    e = rg->RndFlux().Landau(2.0,1.0);
  }
  // just have it come face on +z direction
  p4u = TLorentzVector(0,0,e,e);

  // flux coming off a plane at z=0
  // guassian intensity in x-y, centered at x=0.5 w/ sigmax=4, sigmay=3
  // but cut so nothing |x|>15 or |y|>15
  bool okay_xy = false;
  double a, b, x, y, z=0, t=0;
  while ( ! okay_xy ) {
    rg->RndFlux().Rannor(a,b);
    x = 0.5 + a*4.0;
    y = b*3.0;
    if ( fabs(x) < 15 || fabs(y) < 15 ) okay_xy = true;
  }
  x4u = TLorentzVector(x,y,z,t);

  // decayDist is distance from decay of particle giving the neutrino
  // to the x4u position ... for this test just put 735e3
  decayDist = 735.e3;

  // unweighted rays are easiest
  wgt = 1.0;

}
