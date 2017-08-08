#ifndef __chrisForTHSProject_h__
#define __chrisForTHSProject_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<THSParticle*>+;//want to make tree branch
#pragma link C++ class vector<THSParticle>+;//want to make tree branch Can get rid all pragma links here since they are in THSParticle and have no effect here.
#pragma link C++ class vector<TVector3*>+;//want to make tree branch
#endif
#include "GTreeManager.h"
#include "chrisPPhysics.h"
#include "GTreeTrack.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVector3.h"
#include "TGraph.h"
#include "THSParticle.h"
#include "TROOT.h"
#pragma link C++ class vector<TVector3*>+;//want to make tree branch
#pragma link C++ class vector<TVector3>+;//want to make tree branch

using namespace std;



class	chrisForTHSProject  : public chrisPPhysics
{
private:



  TTree *treePi0;
 
  Int_t	   usePeriodMacro;
  Int_t    period;
  Int_t    nEventsWritten;


  //PID plan create an array of length 24 to store the elements of the PID phi. DATA version
  Double_t PIDElemPhi[24] = {8.574,22.974,37.379,51.784,66.188,80.593,94.997,109.402,123.806,138.211,152.615,167.02,-178.93,-163.16,-147.39,-131.62,-115.85,-100.08,-84.31,-68.54,-52.77,-37.01,-21.24,-5.47};
  //Simulation Version
  //Double_t PIDElemPhi[24] = { 172.5, 157.5, 142.5, 127.5, 112.5, 97.5, 82.5, 67.5, 52.5, 37.5, 22.5, 7.5,  -7.5, -22.5, -37.5, -52.5, -67.5, -82.5, -97.5, -112.5, -127.5, -142.5, -157.5, -172.5};

//Mott Measurements runs closest 14579,14673,14777,14893,15029,15122,15426,15584,15949
Double_t MottMeas[9]={0.7515, 0.7629915, 0.762434, 0.7662845, 0.769435, 0.773674, 0.7759165,0.7736475, 0.7282685    };

Double_t MCMottNeg = 1;
Double_t MCMottPos = -1;


  //Branches in the trees
  vector<THSParticle*> Particles;

  
  Double_t fbeamHelicity;
  Double_t ftaggedTime;
  Double_t feventNo;
  Double_t fenergyBeam;
  Double_t fenergySum;
  Int_t fmultiplicity;
  Double_t fpidPhi;
  Int_t fpidIndex;
  TVector3 fchamber1Vec;
  TVector3 fchamber2Vec;
  Double_t fpidRootinoPhi;
  Int_t fedgePlane;
  Double_t frootinoPhi;

  //Other parameters used in .cc

  vector<THSParticle*> * fReadParticles=nullptr;

  THSParticle particle1;
  //THSParticle *DummyProton;
  //THSParticle *fballPhotons;
  //THSParticle *ftaggPhotons;
  THSParticle ftaggPhotons;
  THSParticle fballPhotons;
  THSParticle fcbRootino;
  

  TLorentzVector fcbPhoton;
  TLorentzVector fglasgowTaggerPhoton;
  TLorentzVector missingp4;
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector particle;
  TLorentzVector frootino;

  
  Int_t PidHitIndex;
  Int_t multiplicity;
  Int_t EventNumber;
  Int_t NPidhits;
  Int_t fNin;
  Double_t energySum;
  Double_t fchamber1VecPhi;
  Double_t fphidiff;
  TVector3 targetPosition;
  Double_t flinPol;
  Double_t ePol; 
  Double_t Pcirc; 
  std::string fileNo;




//Dwelete me
Double_t  rootinoClustE;
Double_t  rootinoTheta;
Double_t  rootinoPhi;
Double_t  rootinoTime;
Double_t  rootinoClustS;
Double_t  rootinoClustC;
Double_t  rootinoVetoC;
Double_t  rootinoDet;
Double_t  rootinoVetoE; 
Double_t  rootinoCham1E;
Double_t  rootinoCham2E;
Int_t testcounter=0;





  
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
