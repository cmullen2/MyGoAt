#ifndef __chrisPPi0Example_h__
#define __chrisPPi0Example_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>

#include "GTreeManager.h"
#include "chrisPPhysics.h"
#include "GTreeTrack.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVector3.h"
#include "TGraph.h"

class	chrisPPi0Example  : public chrisPPhysics
{
private:

//    GH1*	nChamberHitsin1; //17/5/16 Chris Mullen 
//    GH1*	Chamber1X; // 17/5/16 Chris Mullen 
//    GH1*	Chamber1Y; // 17/5/16 Chris Mullen 
//    GH1*	Chamber1Z; // 17/5/16 Chris Mullen 
//    GH1*	nChamberHitsin2; // 17/5/16 Chris Mullen 
//    GH1*	Chamber2X; // 17/5/16 Chris Mullen 
//    GH1*	Chamber2Y; // 17/5/16 Chris Mullen 
//    GH1*	Chamber2Z; // 17/5/16 Chris Mullen 


  TTree *treeselected;

  Int_t	   usePeriodMacro;
  Int_t    period;
  Int_t    nEventsWritten;


//PID plan create an array of length 24 to store the elements of the PID phi.
  Double_t PIDElemPhi[24] = {8.574,22.974,37.379,51.784,66.188,80.593,94.997,109.402,123.806,138.211,152.615,167.02,-178.93,-163.16,-147.39,-131.62,-115.85,-100.08,-84.31,-68.54,-52.77,-37.01,-21.24,-5.47};






//First Branch Parameters 1 rootino branch

  TLorentzVector cutmissingp4;
  TLorentzVector missingp4;
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector particle;
  TLorentzVector PionCan;
  TLorentzVector ProtonCan;
  TLorentzVector ReconProton;
  TLorentzVector PionBoost;
  TLorentzVector ProtonBoost;
  TLorentzVector ReconProtonBoost;

  TLorentzVector NewEnergy;
  Int_t detectornum;
  Int_t PromptRegion;
  Int_t RandomRegion;

 Int_t test;
 
  Double_t Chamber1_VecPhi;
  Double_t Phidiff;
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



  TVector3 TestChamber;

  Int_t PidHitIndex;
  Int_t multiplicity;
  Double_t energySum;
  Int_t EventNumber;
  

  Double_t inv_M_value;
  Double_t energy_beam;
  Double_t MissingM;
  Double_t phiPi;
  Double_t phiProton;
  Double_t coplanarityMeasure;
  Double_t thetaPi;
  Double_t thetaProton;
  Double_t thetadiff;
  Double_t* MWPC_Cham1_X;

  Double_t cutphiProton;  
  Double_t cutthetaProton;  
  Double_t cutphiPion;  
  Double_t cutthetaPion;  
  Double_t cutinv_M_value;  
  Double_t cutcoplanarity;  
  Double_t cutenergy_beam;  
  Double_t cutMissingM;  
    
  Double_t cutthetadiff;  

  Double_t phiReconProton; 
  Double_t thetaReconProton; 
  Double_t coplanarityRecon; 
  Double_t Reconthetdiff;

  Double_t ProtonEnergy;


  Double_t ChamberTest;
  Double_t ChamberTest2;
  Double_t ChamberTest3;


  Double_t PIDPhi;
  Int_t NPidhits;


  TVector3 Chamber1_Vec;
  TVector3 Chamber2_Vec;
  TVector3 targetPosition;
  TVector3 PVector1;
  TVector3 PVector2;

//Second Branch Parameters 3 Gamma branch

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



  // Center of mass boost vector //
  TLorentzVector CMVector(TLorentzVector vec , TLorentzVector vec_t , TLorentzVector vec_b)
  {
    TLorentzVector vec_cm=vec;
    TVector3 beta=((vec_t+vec_b).BoostVector())*(-1.0);
    vec_cm.Boost(beta);
    return vec_cm;
  }
  
 protected:
  virtual Bool_t  Start();
  virtual void    ProcessEvent();
  virtual void	ProcessScalerRead();
  virtual Bool_t    Write();
  
 public:
  chrisPPi0Example();
  virtual ~chrisPPi0Example();
  virtual Bool_t  Init();
  
};
#endif
