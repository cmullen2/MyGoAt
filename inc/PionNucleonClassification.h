#ifndef __PionNucleonClassification_h__
#define __PionNucleonClassification_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "GTreeManager.h"
#include "chrisPPhysics.h"
#include "GTreeTrack.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVector3.h"
#include "TGraph.h"
#include "THSParticle.h"
#include "THSEventInfo.h"
#include "TROOT.h"


using namespace std;

class	PionNucleonClassification  : public chrisPPhysics
{
 private:

  TTree *treePi0;
  Int_t	   usePeriodMacro;
  Int_t    period;
  Int_t    nEventsWritten;
  //Branches in the trees
  vector<THSParticle> Particles;
  vector<THSParticle> Generated;
  THSParticle part;
  THSParticle Gen;
  Double_t TruthPLab=0;
  Double_t TruthELab=0;
  Double_t TruthTheta=0;
  Double_t DetectedClusterE=0;
  Double_t DetectedParticleE=0;
  Double_t DeltaE=0;
  Double_t DetectedTheta=0;
  THSEventInfo fEventInfo;
  Double_t X =0;
  Double_t X2=0;
  Double_t Y =0;
  Double_t Y2=0;
  Double_t Z=0;
  Double_t Numerator=-1;
  TLorentzVector fcbPhoton;
  TLorentzVector fglasgowTaggerPhoton;
  TLorentzVector beam;
  TLorentzVector frootino;
  Double_t DE = 1.2;
//Mott Measurements runs closest 14579,14673,14777,14893,15029,15122,15426,15584,15949
  Double_t fbeamHelicity;
  Double_t ftaggedTime;
  Double_t feventNo;
  Double_t fenergyBeam;
  Double_t fenergySum;
  Double_t generatedPDGs[3]={2112, 211, 2112}; //For Neutron Pi0 to be tested (2212=proton according to root website , 2112 = Neutron)
  Int_t fmultiplicity;
  Int_t fedgePlane;
  Int_t ftaggChannel;
  Int_t PidHitIndex;
  Int_t multiplicity;
  Int_t EventNumber;
  Int_t NPidhits;
  Int_t particleindex=-1;
  Int_t rootindex=-1;
  Int_t photindex=-1;
  Int_t photindex2=-1;
  TVector3 targetPosition;
  Double_t flinPol;
  Double_t ePol; 
  Double_t Pcirc; 
  std::string fileNo;
  
 protected:
  virtual Bool_t  Start();
  virtual void    ProcessEvent();
  virtual void	ProcessScalerRead();
  virtual Bool_t    Write();
  
 public:
  PionNucleonClassification();
  virtual ~PionNucleonClassification();
  virtual Bool_t  Init();
  
};
#endif
