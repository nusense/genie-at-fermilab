#include "Numerical/RandomGen.h"
#include "FluxDrivers/GNuMIFlux.h"
#include "FluxDrivers/GSimpleNtpFlux.h"
#include "Utils/UnitUtils.h"
#include "Utils/AppInit.h"

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

using namespace std;
using namespace genie;
using namespace genie::flux;

// forward declaration
void encode_transform(flux::GNuMIFlux* gnumi, string& rot, string& pos);

// main routine
void gnumi2simple_basic(string fnameout="fluxntgenie.root",
                        string fluxpatt="",
                        string cfg="MINOS-Near", 
                        long int nentries=0, 
                        double pots=1.0e30,
                        double enumin=0.0,
                        bool doaux=true,
                        string srcpatt="",
                        bool weighted=false)
{
  
  cout << "Creating:    " << fnameout << endl;
  cout << "Input files: " << fluxpatt << endl;
  cout << "Config:      " << cfg << endl;
  cout << "NEntries:    " << nentries << endl;
  cout << "POTs:        " << pots << endl;
  if ( enumin > 0.0 ) 
    cout << "EnuMin:      " << enumin << endl;
  cout << "AuxTree:     " << doaux << endl;
  if ( srcpatt != "" )
    cout << "SrcPatt:     " << srcpatt << endl;
  if ( weighted ) cout << "Generate weighted entries" << endl;
  else            cout << "Generate unweighted entries" << endl;

  // post R-2_8_0 one can NOT use GSEED
  const char* gseedstr = std::getenv("GSEED");
  if ( ! gseedstr ) {
    cout << "no GSEED set ... exit" << endl;
    return;
  }
  long int seed = strtol(gseedstr,NULL,0);
  ::unsetenv("GSEED");  // must be unset
  genie::utils::app_init::RandGen(seed);
  //old pre R-2_8_0 // RandomGen::SetSeed(long int seed)

  string fluxfname(gSystem->ExpandPathName(fluxpatt.c_str()));
  flux::GNuMIFlux* gnumi = new GNuMIFlux();
  gnumi->LoadBeamSimData(fluxfname,cfg);
  gnumi->SetEntryReuse(1); // don't reuse entries when reformatting

  if ( weighted ) gnumi->GenerateWeighted(weighted);

  //gnumi->PrintConfig();
  //cout << " change to cm:" << endl;
  //double cm_unit = genie::utils::units::UnitFromString("cm");
  //gnumi->SetLengthUnits(cm_unit);
  //gnumi->PrintConfig();
  //cout << " change to m:" << endl;
  ///////RWH: the following 2 lines should work ...
  //double m_unit = genie::utils::units::UnitFromString("m");
  //gnumi->SetLengthUnits(m_unit);
  gnumi->PrintConfig();
  //double scaleusr = 100.;

  gnumi->SetEntryReuse(1);     // reuse entries
  gnumi->SetUpstreamZ(-3e38);  // leave ray on flux window
  
  GFluxI* fdriver = dynamic_cast<GFluxI*>(gnumi);

  if (nentries == 0) nentries = 2147483647;

  // so as not to include scan time generate 1 nu
  fdriver->GenerateNext();

  TFile* file = TFile::Open(fnameout.c_str(),"RECREATE");
  TTree* fluxntp = new TTree("flux","a simple flux n-tuple");
  TTree* metantp = new TTree("meta","metadata for flux n-tuple");
  genie::flux::GSimpleNtpEntry* fentry = new genie::flux::GSimpleNtpEntry;
  genie::flux::GSimpleNtpNuMI*  fnumi  = new genie::flux::GSimpleNtpNuMI;
  genie::flux::GSimpleNtpAux*   faux   = new genie::flux::GSimpleNtpAux;
  genie::flux::GSimpleNtpMeta*  fmeta  = new genie::flux::GSimpleNtpMeta;
  fluxntp->Branch("entry",&fentry);
  fluxntp->Branch("numi",&fnumi);
  if (doaux) fluxntp->Branch("aux",&faux);
  metantp->Branch("meta",&fmeta);

  TLorentzVector p4u, x4u;
  long int ngen = 0;

  set<int> pdglist;
  double maxe = 0;
  double minwgt = +1.0e10;
  double maxwgt = -1.0e10;

  UInt_t metakey = TString::Hash(fnameout.c_str(),strlen(fnameout.c_str()));
  cout << "metakey " << metakey << endl;
  // ensure that we don't get smashed by UInt_t vs Int_t
  metakey &= 0x7FFFFFFF;
  cout << "metakey " << metakey << " after 0x7FFFFFFF" << endl;
  cout << "=========================== Start " << endl;

  TStopwatch sw;
  sw.Start();
  int nok = 0, nwrite = 0;
  while ( nok < nentries && gnumi->UsedPOTs() < pots ) {
    fdriver->GenerateNext();
    //if ( nok < 10 ) cout << gnumi->PassThroughInfo() << endl;

    ngen++;
    fentry->Reset();
    fnumi->Reset();
    faux->Reset();

    fentry->metakey = metakey;
    fentry->pdg     = fdriver->PdgCode();
    fentry->wgt     = gnumi->Weight();
    x4u     = fdriver->Position();
    fentry->vtxx    = x4u.X();
    fentry->vtxy    = x4u.Y();
    fentry->vtxz    = x4u.Z();
    fentry->dist    = gnumi->GetDecayDist();
    p4u     = fdriver->Momentum();
    fentry->px      = p4u.Px();
    fentry->py      = p4u.Py();
    fentry->pz      = p4u.Pz();
    fentry->E       = p4u.E();

    fnumi->run      = gnumi->PassThroughInfo().run;
    fnumi->evtno    = gnumi->PassThroughInfo().evtno;
    fnumi->entryno  = gnumi->GetEntryNumber();

    fnumi->tpx      = gnumi->PassThroughInfo().tpx;
    fnumi->tpy      = gnumi->PassThroughInfo().tpy;
    fnumi->tpz      = gnumi->PassThroughInfo().tpz;
    fnumi->vx       = gnumi->PassThroughInfo().vx;
    fnumi->vy       = gnumi->PassThroughInfo().vy;
    fnumi->vz       = gnumi->PassThroughInfo().vz;

    fnumi->pdpx     = gnumi->PassThroughInfo().pdpx;
    fnumi->pdpy     = gnumi->PassThroughInfo().pdpy;
    fnumi->pdpz     = gnumi->PassThroughInfo().pdpz;

    double apppz = gnumi->PassThroughInfo().pppz;
    fnumi->pppx     = gnumi->PassThroughInfo().ppdxdz * apppz;
    fnumi->pppy     = gnumi->PassThroughInfo().ppdydz * apppz;
    fnumi->pppz     = apppz;

    fnumi->ndecay   = gnumi->PassThroughInfo().ndecay;
    fnumi->ptype    = gnumi->PassThroughInfo().ptype;
    fnumi->ppmedium = gnumi->PassThroughInfo().ppmedium;
    fnumi->tptype   = gnumi->PassThroughInfo().tptype;


    if ( doaux ) {
      faux->auxint.push_back(gnumi->PassThroughInfo().tgen);
      faux->auxdbl.push_back(gnumi->PassThroughInfo().fgXYWgt);
      faux->auxdbl.push_back(gnumi->PassThroughInfo().nimpwt);
    }

    if ( fentry->E > enumin ) ++nok;
    fluxntp->Fill();
    ++nwrite;

    // accumulate meta data
    pdglist.insert(fentry->pdg);
    minwgt = TMath::Min(minwgt,fentry->wgt);
    maxwgt = TMath::Max(maxwgt,fentry->wgt);
    maxe   = TMath::Max(maxe,fentry->E);

  }
  cout << "=========================== Complete " << endl;

  fmeta->pdglist.clear();
  set<int>::const_iterator setitr = pdglist.begin();
  for ( ; setitr != pdglist.end(); ++setitr)  fmeta->pdglist.push_back(*setitr);
 
  fmeta->maxEnergy = maxe;
  fmeta->minWgt    = minwgt;
  fmeta->maxWgt    = maxwgt;
  fmeta->protons   = gnumi->UsedPOTs();
  TVector3 p0, p1, p2;
  gnumi->GetFluxWindow(p0,p1,p2);
  TVector3 d1 = p1 - p0;
  TVector3 d2 = p2 - p0;
  fmeta->windowBase[0] = p0.X();
  fmeta->windowBase[1] = p0.Y();
  fmeta->windowBase[2] = p0.Z();
  fmeta->windowDir1[0] = d1.X();
  fmeta->windowDir1[1] = d1.Y();
  fmeta->windowDir1[2] = d1.Z();
  fmeta->windowDir2[0] = d2.X();
  fmeta->windowDir2[1] = d2.Y();
  fmeta->windowDir2[2] = d2.Z();
  if ( doaux ) {
    fmeta->auxintname.push_back("tgen");
    fmeta->auxdblname.push_back("fgXYWgt");
    fmeta->auxdblname.push_back("nimpwt");
  }

  std::string rotinfo, posinfo;
  encode_transform(gnumi,rotinfo,posinfo);
  fmeta->infiles.push_back(rotinfo);
  fmeta->infiles.push_back(posinfo);

  if ( srcpatt != "" ) {
    string srcstring = "srcpatt=" + srcpatt;
    fmeta->infiles.push_back(srcstring);
  }
  string fname = "flux_pattern=" + fluxfname;
  fmeta->infiles.push_back(fname); // should get this expanded from gnumi
  vector<string> flist = gnumi->GetFileList();
  for (size_t i = 0; i < flist.size(); ++i) {
    fname = "flux_infile=" + flist[i];
    fmeta->infiles.push_back(fname);
  }

  fmeta->seed    = RandomGen::Instance()->GetSeed();
  fmeta->metakey = metakey;

  metantp->Fill();

  sw.Stop();
  cout << "Generated " << nwrite << " (" << nok << " w/ E >" << enumin << " ) "
       << " ( request " << nentries << " ) "
       << endl
       << gnumi->UsedPOTs() << " POTs " << " ( request " << pots << " )"
       << endl
       << " pulled NFluxNeutrinos " << gnumi->NFluxNeutrinos()
       << endl
       << "Time to generate: " << endl;
  sw.Print();
  
  cout << "===================================================" << endl;

  cout << "Last GSimpleNtpEntry: " << endl
       << *fentry
       << endl
       << "Last GSimpleNtpMeta: " << endl
       << *fmeta
       << endl;

  file->cd();

  // write meta as simple object at top of file
  //fmeta->Write();

  // write ntuples out
  fluxntp->Write();
  metantp->Write();
  file->Close();

  cout << endl << endl;
  gnumi->PrintConfig();
}

void encode_transform(flux::GNuMIFlux* gnumi, string& rot, string& pos)
{
  TRotation rot3x3 = gnumi->GetBeamRotation();
  TVector3  posv3  = gnumi->GetBeamCenter();

  std::ostringstream rotstr;
  std::ostringstream posstr;

  const int p=15, w=20;
  
  rotstr << std::setprecision(p);
  rotstr << "<beamdir type=\"newxyz\"> \n"
         << "[ "
         << std::setw(w) << rot3x3.XX() << " "
         << std::setw(w) << rot3x3.XY() << " "
         << std::setw(w) << rot3x3.XZ() << " ]\n"
         << "[ "
         << std::setw(w) << rot3x3.YX() << " "
         << std::setw(w) << rot3x3.YY() << " "
         << std::setw(w) << rot3x3.YZ() << " ]\n"
         << "[ "
         << std::setw(w) << rot3x3.ZX() << " "
         << std::setw(w) << rot3x3.ZY() << " "
         << std::setw(w) << rot3x3.ZZ() << " ]\n"
         << "    </beamdir>" << std::endl;
  posstr << std::setprecision(p);
  posstr << "<beampos> ( "
         << std::setw(w) << posv3.X() << " "
         << std::setw(w) << posv3.Y() << " "
         << std::setw(w) << posv3.Z() << " ) "
         << "</beampos>" << std::endl;

  rot = rotstr.str();
  pos = posstr.str();
}

