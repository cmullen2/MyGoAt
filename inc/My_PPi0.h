#ifndef __My_PPi0_h__
#define __My_PPi0_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"


class	My_PPi0  : public PPhysics
{
 private:
  
  TTree *tree;
  
  Double_t timev;
  Double_t timeT;
  Double_t timep;
  Double_t inv_M_value;
  Double_t energy_beam;
  Double_t phi1;
  Double_t phi2;
  Double_t phi3;
  Double_t phi_para_pos;
  Double_t phi_para_neg;
  Double_t phi_perp_pos;
  Double_t phi_perp_neg;
  Double_t theta1;
  Double_t theta2;
  Double_t thetadiff;
  Double_t coplanarity;
  Double_t pt;
  Double_t pb;
  Double_t pb2;
  Double_t cosine_CMS;
  Double_t Energy_CMS;
  Double_t MissingM;

  Int_t edgeplane;
  Double_t edge_setting;
  Double_t edge_position;
  
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector particle;
  TLorentzVector missingp4;
  TLorentzVector boost_4vec;


  // Target polarisation loop //
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

      // Center of mass boost vector //
      TLorentzVector CMVector(TLorentzVector vec,TLorentzVector 
			      vec_t,TLorentzVector vec_b)
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
    My_PPi0();
    virtual ~My_PPi0();
    virtual Bool_t  Init();

};
#endif
