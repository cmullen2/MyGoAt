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


  

    GH1*	nChamberHitsin1; //17/5/16 Chris Mullen 
//    GH1*	Chamber1X; // 17/5/16 Chris Mullen 
//    GH1*	Chamber1Y; // 17/5/16 Chris Mullen 
//    GH1*	Chamber1Z; // 17/5/16 Chris Mullen 
    GH1*	nChamberHitsin2; // 17/5/16 Chris Mullen 
//    GH1*	Chamber2X; // 17/5/16 Chris Mullen 
//    GH1*	Chamber2Y; // 17/5/16 Chris Mullen 
//    GH1*	Chamber2Z; // 17/5/16 Chris Mullen 


  TTree *tree1r;
  TTree *tree3g;
  TTree *treetest;
  TTree *treetrigger;
  TTree *treescatter;
  TTree *treeselected;
  TTree *treepid;



  Int_t	   usePeriodMacro;
  Int_t    period;
  Int_t    nEventsWritten;

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


 Int_t test;
 
  Double_t Chamber1_VecPhi;
  Double_t Phidiff;
  Double_t BeamHelicity;


  TLorentzVector Mp41;
  TLorentzVector Mp42;
  TLorentzVector Mp43;
  Double_t MM1;
  Double_t MM2;
  Double_t MM3;
  Double_t MMDiff1;
  Double_t MMDiff2;
  Double_t MMDiff3;





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

//  Double_t ReconphiProton;  
//  Double_t Reconcoplanarity;  
//  Double_t ReconthetaProton;  
//  Double_t Reconthetadiff;  

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



//Currently unused or commented out Parameters

//  TLorentzVector Rootino4Vec;

//  Double_t theta1;
//  Double_t theta2;
//  Double_t phi1;
//  Double_t phi2;
//  Double_t timev;
//  Double_t timeT;
//  Double_t timep;
//  Double_t phi3;
//  Double_t phi_para_pos;
//  Double_t phi_para_neg;
//  Double_t phi_perp_pos;
//  Double_t phi_perp_neg;
//  Double_t pt;
//  Double_t pb;
//  Double_t pb2;
//  Double_t cosine_CMS;
//  Double_t Energy_CMS;

//  Int_t edgeplane;
//  Double_t edge_setting;
//  Double_t edge_position;


//  Double_t time;
//  Double_t T;
//  Double_t T1;
//  Double_t T2;
//  Double_t T3;
//  Double_t C1_M_value;
//  Double_t C2_M_value;
//  Double_t C3_M_value;
//  Double_t inv_M_value;
//  Double_t energy_beam;
//  Double_t phi1;
//  Double_t phi2;
//  Double_t phiP;
//  Double_t theta1;
//  Double_t theta2;
//  Double_t thetaP;
//  Double_t thetadiff;
//  Double_t thetadiff2;
//  Double_t coplanarity;
//  Double_t coplanarity2;
//  Double_t pt;
//  Double_t pb;
//  Double_t cosine_CMS;
//  Double_t Energy_CMS;
//  Double_t MissingM_C1;
//  Double_t MissingM_C2;
//  Double_t MissingM_C3;
//  Double_t MissingM;


//  Double_t coplanarity23;
//  Double_t phiMissingp4;
//  Double_t CoplanMissingp4;
//  Double_t Zerocoplanarity;
//  Double_t ZeromissingCoplan;


  

//  Double_t AngDiff1;
//  Double_t AngDiff2;
//  Double_t AngDiff3;

//  Int_t    edgeplane;
//  Double_t edge_setting;
//  Double_t edge_position;
  


  

//  TLorentzVector boost_4vec;



//  TLorentzVector beam;
//  TLorentzVector target;
//  TLorentzVector particle;
//  TLorentzVector C1missingp4;
//  TLorentzVector C2missingp4;
//  TLorentzVector C3missingp4;
//  TLorentzVector missingp4;
//  TLorentzVector boost_4vec;
//  TLorentzVector boost_MMvec;
//  TLorentzVector boost_G3vec;

 
//  TLorentzVector meson_Vec;
//  TLorentzVector meson;
//  TLorentzVector P;


//  TLorentzVector ZeroPion;
//  TLorentzVector ZeroProton;
//  TLorentzVector PionTLorentz;
//  TLorentzVector Zeromissingp4;
  


  
//  TVector3 Pi0_vec;
//  TVector3 P_vec;

/*
      // Center of mass boost vector //
      TLorentzVector CMVector(TLorentzVector vec,TLorentzVector vec_t,TLorentzVector vec_b)
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
*/

/* // Target polarisation loop //
  Double_t targetpol(TFile *f){
    
    TString* filename = new TString(f->GetPath());
    TString* path1 = new TString(filename->Tokenize("_")->At(filename->Tokenize("_")->GetEntries()-1)->GetName());
    path1->Resize(path1->Length()-7);
    stringstream ss(path1->Data());
    Int_t runnumber1;
    ss >> runnumber1;
    string plane;
    string poledge;
    
    ifstream in2;
    in2.open("/home/roddym/Documents/a2GoAT/allinfos_recentbeamtimes.txt");
    Double_t runnumber;
    Double_t pol;
    Double_t pol1;
    Int_t N = 0;
    
    if (in2.is_open()) {
      while (1) {
	
	in2 >> runnumber >> plane >> poledge >> pol;
	//       cout << runnumber << endl;
	if (!in2.good()){pol1=500.;break;}
	if(runnumber==runnumber1){pol1=pol;break;}
	
      }
      N++;
    }
    else{cout << "File does not open !" << endl;}
    return pol1/100;
  }
 */ 
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
