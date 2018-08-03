#ifndef __chrisForTHSProject_h__
#define __chrisForTHSProject_h__

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

class	chrisForTHSProject  : public chrisPPhysics
{
 private:

  TTree *treePi0;
 
  Int_t	   usePeriodMacro;
  Int_t    period;
  Int_t    nEventsWritten;
  
  //Branches in the trees
  vector<THSParticle> Particles;
  vector<THSParticle> Generated;
  THSEventInfo fEventInfo;


  TLorentzVector fcbPhoton;
  TLorentzVector fglasgowTaggerPhoton;
  TLorentzVector beam;
  TLorentzVector frootino;



  //Mott Measurements runs closest 14579,14673,14777,14893,15029,15122,15426,15584,15949
  Double_t MottMeas[9]={0.7515, 0.7629915, 0.762434, 0.7662845, 0.769435, 0.773674, 0.7759165,0.7736475, 0.7282685    };
  Double_t MCMottNeg = 1;
  Double_t MCMottPos = -1;
  Double_t MCMottFla = 0;
  
  Double_t fbeamHelicity;
  Double_t ftaggedTime;
  Double_t feventNo;
  Double_t fenergyBeam;
  Double_t fenergySum;

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

  Double_t generatedPDGs[6]={2212, -211, 2212, 22, 22, -22}; //For Proton Pi0PiM
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
  chrisForTHSProject();
  virtual ~chrisForTHSProject();
  virtual Bool_t  Init();
  
};
#endif
