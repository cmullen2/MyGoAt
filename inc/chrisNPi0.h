#ifndef __chrisNPi0_h__
#define __chrisnPi0_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<THSParticle*>+;//want to make tree branch
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

using namespace std;



class	chrisNPi0  : public chrisPPhysics
{
private:



  TTree *treeselected;

  Int_t	   usePeriodMacro;
  Int_t    period;
  Int_t    nEventsWritten;


//PID plan create an array of length 24 to store the elements of the PID phi. DATA version
//  Double_t PIDElemPhi[24] = {8.574,22.974,37.379,51.784,66.188,80.593,94.997,109.402,123.806,138.211,152.615,167.02,-178.93,-163.16,-147.39,-131.62,-115.85,-100.08,-84.31,-68.54,-52.77,-37.01,-21.24,-5.47};
//Simulation Version
Double_t PIDElemPhi[24] = { 172.5, 157.5, 142.5, 127.5, 112.5, 97.5, 82.5, 67.5, 52.5, 37.5, 22.5, 7.5,  -7.5, -22.5, -37.5, -52.5, -67.5, -82.5, -97.5, -112.5, -127.5, -142.5, -157.5, -172.5};



//Selection Criteria 
Double_t MMPi0Criteria;
Double_t IMPi0Criteria;
Double_t MMuppercut;
Double_t MMlowercut;
Double_t invlowercut;
Double_t invuppercut;



//Testing THSParticle
vector<THSParticle*> fparticleList;  //A vector of type THSParticle, will push back to this.

THSParticle particle1;
THSParticle *DummyProton;
THSParticle *DummyProton2;



//Splots branch parameters
Double_t Coplanarity;
Double_t OpeningAngle;
Double_t EventNo;


  TLorentzVector missingp4;
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector particle;
  TLorentzVector PionCan;
  TLorentzVector NeutronCan;



  Double_t BeamHelicity;
  Double_t time_beam;
 
  TLorentzVector Mp41;
  TLorentzVector Mp42;
  TLorentzVector Mp43;
  Double_t MM1;
  Double_t MM2;
  Double_t MM3;
  Double_t MMDiff1;
  Double_t MMDiff2;
  Double_t MMDiff3;

  TLorentzVector cluster1_pi;
  TLorentzVector cluster2_pi;
  Double_t cluster_pi_ang_diff;


  Int_t PidHitIndex;
  Int_t multiplicity;
  Double_t energySum;
  Int_t EventNumber;
  

  Double_t inv_M_value;
  Double_t energy_beam;
  Double_t MissingM; 
  Double_t NeutronMass;

  Double_t ProtonEnergy;

  Int_t NPidhits;
  TVector3 targetPosition;
  

  Int_t    l_g1;
  Int_t    l_g2;

  TLorentzVector g1_Vec1;
  TLorentzVector g2_Vec1;
  TLorentzVector g1_Vec2;
  TLorentzVector g2_Vec2;
  TLorentzVector g1_Vec3;
  TLorentzVector g2_Vec3;
  TLorentzVector C1;
  TLorentzVector C2;
  TLorentzVector C3;
  TLorentzVector P1;
  TLorentzVector P2;
  TLorentzVector P3;


  Double_t IMDiff1;
  Double_t IMDiff2;
  Double_t IMDiff3;
  Double_t CopDiff1;
  Double_t CopDiff2;
  Double_t CopDiff3;
  Double_t SumDiff1;
  Double_t SumDiff2;
  Double_t SumDiff3;



  
 protected:
  virtual Bool_t  Start();
  virtual void    ProcessEvent();
  virtual void	ProcessScalerRead();
  virtual Bool_t    Write();
  
 public:
  chrisNPi0();
  virtual ~chrisNPi0();
  virtual Bool_t  Init();
  
};
#endif