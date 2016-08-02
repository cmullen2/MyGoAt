#ifndef __chrisPPhysics_h__
#define __chrisPPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GH1.h"
#include "GConfigFile.h"
#include "TGraph.h"



class	chrisPPhysics : virtual public GTreeManager
{
private:

	Double_t targetmass;
	Double_t Prompt_low;
	Double_t Prompt_high;
	Double_t Random_low1;
	Double_t Random_high1;
	Double_t Random_low2;
	Double_t Random_high2;
	
	Double_t PvR_ratio;
		
	TLorentzVector beam;
	TLorentzVector target;
	TLorentzVector particle;
	TLorentzVector missingp4;
	TLorentzVector PRecon; //CAM 18/5/16
	TLorentzVector PNormal; //CAM 18/5/16
	
	Double_t time;
	Bool_t 	Prompt;
	Bool_t 	Random;

	Int_t TC_cut_min;
	Int_t TC_cut_max;

	Int_t TC_scaler_min;
	Int_t TC_scaler_max;
	
	
protected:


public:
    chrisPPhysics();
    virtual ~chrisPPhysics();

    virtual Bool_t	Init();
	virtual void	Analyse() {;}
	virtual void	Reconstruct();
    virtual Bool_t	Write();


	//Chris Mullen CAM 29/6/16 to plot the cylindrical nature of the wire chamber hits using TGraph
	void	CylindricalChamber(const GTreeParticle& treePi ,const GTreeMWPCHit& tree, TGraph* gGraph,Int_t ChamberNo, Int_t PlotNum, Bool_t TaggerBinning = kFALSE);
	void CylindricalChamber(const GTreeMWPCHit& tree,  Int_t particle_index, Int_t tagger_index, TGraph* gGraph, Int_t ChamberNo,  Int_t PlotNum, Bool_t TaggerBinning = kFALSE);

	//Chris Mullen CAM 6/4/16 added the Fill P* and FillE
	void	FillE(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void	FillE(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillE(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);

	void	FillP(const GTreeParticle& tree, GH1* gHist, Int_t VecElem ,Bool_t TaggerBinning = kFALSE);
	void	FillP(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Int_t VecElem,  Bool_t TaggerBinning = kFALSE);
	void 	FillP(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Int_t VecElem, Bool_t TaggerBinning = kFALSE);


	//Chris Mullen CAM 7/4/16 trying to add the wire chamber hit positions

	void FillnChamber(const GTreeMWPCHit& tree, GH1* gHist, Int_t ChamberNo);
	void FillChamber(const GTreeParticle& treePi ,const GTreeMWPCHit& tree, GH1* gHist, Int_t ChamberNo, Int_t ChamberAxis, Bool_t TaggerBinning = kFALSE);
	void FillChamber(const GTreeMWPCHit& tree,  Int_t particle_index, Int_t tagger_index, GH1* gHist,  Int_t ChamberNo, Int_t ChamberAxis,Bool_t TaggerBinning = kFALSE);


	//Chris Mullen CAM 18/5/16 trying to add the angle reconstruction

	void	FillProtonMM(const GTreeParticle& treePi,const GTreeParticle& treePro, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillProtonMM(const GTreeParticle& treePi, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);


	void	FillProton(const GTreeParticle& treePi,const GTreeParticle& treePro, GH1* gHist, Int_t ProtonElem, Bool_t TaggerBinning = kFALSE);
	void 	FillProton(const GTreeParticle& treePi, Int_t particle_index, Int_t tagger_index, GH1* gHist, Int_t ProtonElem, Bool_t TaggerBinning = kFALSE);

	void	FillProtonPhi(const GTreeParticle& treePi,const GTreeParticle& treePro, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillProtonPhi(const GTreeParticle& treePi, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);

	void	FillProtonTheta(const GTreeParticle& treePi,const GTreeParticle& treePro, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillProtonTheta(const GTreeParticle& treePi, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);

	void	FillProtonAngleDiffPhi(const GTreeParticle& treePi,const GTreeParticle& treePro, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillProtonAngleDiffPhi(const GTreeParticle& treePi, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);






	void	FillProtonNormal(const GTreeParticle& treePro,const GTreeParticle& treePi, GH1* gHist, Int_t ProtonElem, Bool_t TaggerBinning = kFALSE);
	void 	FillProtonNormal(const GTreeParticle& treePro, Int_t particle_index, Int_t tagger_index, GH1* gHist, Int_t ProtonElem, Bool_t TaggerBinning = kFALSE);


	void	FillProtonNormalMass(const GTreeParticle& treePro,const GTreeParticle& treePi, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillProtonNormalMass(const GTreeParticle& treePro, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);


	void	FillProtonNormalPhi(const GTreeParticle& treePro,const GTreeParticle& treePi, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillProtonNormalPhi(const GTreeParticle& treePro, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);

	void	FillProtonNormalTheta(const GTreeParticle& treePro,const GTreeParticle& treePi, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillProtonNormalTheta(const GTreeParticle& treePro, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);

	void	FillProtonNormalAngleDiffTheta(const GTreeParticle& treePro,const GTreeParticle& treePi, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillProtonNormalAngleDiffTheta(const GTreeParticle& treePro, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);





	void	FillMissingMass(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void	FillMissingMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);

    Double_t CalcMissingMass(const GTreeParticle &tree, Int_t particle_index, Int_t tagger_index);
    Double_t CalcMissingEnergy(const GTreeParticle &tree, Int_t particle_index, Int_t tagger_index);
    TLorentzVector CalcMissingP4(const GTreeParticle &tree, Int_t particle_index, Int_t tagger_index);

	void 	FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning = kFALSE, Double_t MM_min = -100000, Double_t MM_max = 100000);
	void 	FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE, Double_t MM_min = -100000, Double_t MM_max = 100000);

	void 	FillTime(const GTreeParticle& tree, GH1* gHist);
	void 	FillTime(const GTreeParticle& tree, Int_t particle_index, GH1* gHist);
	void 	FillTimeCut(const GTreeParticle& tree, GH1* gHist);
	void 	FillTimeCut(const GTreeParticle& tree, Int_t particle_index, GH1* gHist);

	void 	FillMass(const GTreeParticle& tree, GH1* gHist);
	void 	FillMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist);
				
	void	SetTarget(Double_t mass) {target = TLorentzVector(0.,0.,0.,mass);}
	TLorentzVector GetTarget() {return target;}

	void 	SetTC_cut(Int_t cut_min, Int_t cut_max) { TC_cut_min = cut_min; TC_cut_max = cut_max; }
	Int_t 	GetTC_cut_min() { return TC_cut_min;}
	Int_t 	GetTC_cut_max() { return TC_cut_max;}

	void 	SetTC_scalers(Int_t sc_min, Int_t sc_max) { TC_scaler_min = sc_min; TC_scaler_max = sc_max; }
	Int_t 	GetTC_scaler_min() { return TC_scaler_min;}
	Int_t 	GetTC_scaler_max() { return TC_scaler_max;}	

	// TH1 routines
//Chris Mullen CAM 6/4/16 added the Fill P* and FillE Note I have changed to simply P for all three now for GH1's but TH1's are not needed yet could simply delete these


// CAM 17/5/16 to add wire chamber stuff
 
//Chris Mullen CAM 7/4/16 trying to add the wire chamber hit positions

// CAM 17/5/16 to add angle reconstruction stuff



	void FillMissingMass(const GTreeParticle& tree, TH1* Hprompt, TH1* Hrandom);
	void FillMissingMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom);
	void FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom);



	void FillTime(const GTreeParticle& tree, TH1* Hist);
	void FillTime(const GTreeParticle& tree, Int_t particle_index, TH1* Hist);
	void FillTimeCut(const GTreeParticle& tree, TH1* Hist);
	void FillTimeCut(const GTreeParticle& tree, Int_t particle_index, TH1* Hist);

	void FillMass(const GTreeParticle& tree, TH1* Hist);
	void FillMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hist);

	void FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max);
	void FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max);

	Double_t CalcCoplanarity(const GTreeParticle& tree1, Int_t particle_index1, const GTreeParticle& tree2, Int_t particle_index2);

	void FillScalers(Int_t low_scaler_number, Int_t high_scaler_number, TH1* hist);

	Bool_t InitBackgroundCuts();
	Bool_t InitTargetMass();
	Bool_t InitTaggerChannelCuts();
	Bool_t InitTaggerScalers();
};
#endif
