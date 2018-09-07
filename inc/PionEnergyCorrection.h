#ifndef __PionEnergyCorrection_h__
#define __PionEnergyCorrection_h__

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

class	PionEnergyCorrection  : public chrisPPhysics
{
 private:

  TTree *treePi0;
  Int_t	   usePeriodMacro;
  Int_t    period;
  Int_t    nEventsWritten;
  //Branches in the trees
  vector<THSParticle> Particles;
  vector<THSParticle> Generated;
  Double_t TruthPLab;
  Double_t TruthELab;
  Double_t TruthTheta;
  Double_t DetectedClusterE;
  Double_t DetectedParticleE;
  Double_t DeltaE;
  Double_t DetectedTheta;
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
  PionEnergyCorrection();
  virtual ~PionEnergyCorrection();
  virtual Bool_t  Init();
  
};
#endif
