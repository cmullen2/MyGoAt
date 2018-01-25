#ifndef __chrisChargedPions_h__
#define __chrisChargedPions_h__

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



class	chrisChargedPions  : public chrisPPhysics
{
private:


  //Output Trees
  TTree *treePi;

  //Initialisation 
  Int_t	   usePeriodMacro;
  Int_t    period;
  Int_t    nEventsWritten;


  //PID plan create an array of length 24 to store the elements of the PID phi. DATA version
  Double_t PIDElemPhi[24] = {8.574,22.974,37.379,51.784,66.188,80.593,94.997,109.402,123.806,138.211,152.615,167.02,-178.93,-163.16,-147.39,-131.62,-115.85,-100.08,-84.31,-68.54,-52.77,-37.01,-21.24,-5.47};
  //Simulation Version
  //Double_t PIDElemPhi[24] = { 172.5, 157.5, 142.5, 127.5, 112.5, 97.5, 82.5, 67.5, 52.5, 37.5, 22.5, 7.5,  -7.5, -22.5, -37.5, -52.5, -67.5, -82.5, -97.5, -112.5, -127.5, -142.5, -157.5, -172.5};

  //Mott Measurements runs closest 14579,14673,14777,14893,15029,15122,15426,15584,15949
  Double_t MottMeas[9]={0.7515, 0.7629915, 0.762434, 0.7662845, 0.769435, 0.773674, 0.7759165,0.7736475, 0.7282685};
  Double_t MCMottNeg = 1;
  Double_t MCMottPos = -1;
  Double_t MCMottFla = 0;


  //Branches in the trees
  vector<THSParticle*> Particles;
  vector<THSParticle*> Generated;

  //vectors for number of thsparticles per event 
  vector<THSParticle*> * fReadParticles=nullptr;
  vector<THSParticle*> * fGenParticles=nullptr;

  TVector3 fchamber2Vec;
  TVector3 fchamber1Vec;
  Double_t fbeamHelicity;
  Double_t ftaggedTime;
  Double_t feventNo;
  Double_t fenergyBeam;
  Double_t fenergySum;
  Double_t fpidPhi;
  Double_t fpidRootinoPhi;
  Double_t frootinoPhi;
  Int_t fmultiplicity;
  Int_t fpidIndex;
  Int_t fedgePlane;
  Int_t ftaggChannel;

  //Other parameters used in .cc

  TLorentzVector fcbPhoton;
  TLorentzVector fglasgowTaggerPhoton;
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector frootino;

  TVector3 targetPosition;
  Int_t PidHitIndex;
  Int_t multiplicity;
  Int_t EventNumber;
  Int_t NPidhits;
  Int_t fNin;
  Int_t particleindex=-1; //Track index for rootinos 
  Int_t fphotindex=-1;	//Track index for photons
  //Double_t generatedPDGs[5]={2112, 2212, 22 ,22 ,-22}; //p pi0
  //Double_t generatedPDGs[5]={2212, 2112, 22 ,22 ,-22};//n pi0
  Double_t generatedPDGs[4]={2112, 2112, -211 ,-22};//p pi-
  //Double_t generatedPDGs[4]={2212, 2212, 211 ,-22};//n pi+
  Double_t energySum;
  Double_t fchamber1VecPhi;
  Double_t fphidiff;
  Double_t flinPol;
  Double_t ePol; 
  Double_t Pcirc; 
  Double_t taggUpRange;
  Double_t taggLowRange;
  std::string fileNo;


  //For testing channel at goat stage
  Int_t ftopology=0;
  TLorentzVector fcrystalphot;
  Int_t fcrystalphotindex;
  Double_t fcrystalphotmass;
  Double_t fcrystalphotphi;
  Double_t fcrystalphottheta;
  Double_t fcrystalphotpidE;
  Double_t fcrystalphotmwpc0E;
  Double_t fcrystalphotmwpc1E;

  TLorentzVector fcrystalroot;
  Int_t fcrystalrootindex;
  Double_t fcrystalrootmass;
  Double_t fcrystalrootphi;
  Double_t fcrystalroottheta;
  Double_t fcrystalrootpidE;
  Double_t fcrystalrootmwpc0E;
  Double_t fcrystalrootmwpc1E;

  Double_t fcrystalcoplan;


  TLorentzVector fcrystalroot1;
  Int_t fcrystalrootindex1;
  Double_t fcrystalrootmass1;
  Double_t fcrystalrootphi1;
  Double_t fcrystalroottheta1;
  Double_t fcrystalrootpidE1;
  Double_t fcrystalrootmwpc0E1;
  Double_t fcrystalrootmwpc1E1;

  TLorentzVector fcrystalroot2;
  Int_t fcrystalrootindex2;
  Double_t fcrystalrootmass2;
  Double_t fcrystalrootphi2;
  Double_t fcrystalroottheta2;
  Double_t fcrystalrootpidE2;
  Double_t fcrystalrootmwpc0E2;
  Double_t fcrystalrootmwpc1E2;

  Double_t fcrystalcoplan2;





 
  
 protected:
  virtual Bool_t  Start();
  virtual void    ProcessEvent();
  virtual void	ProcessScalerRead();
  virtual Bool_t    Write();
  
 public:
  chrisChargedPions();
  virtual ~chrisChargedPions();
  virtual Bool_t  Init();
  
};
#endif
