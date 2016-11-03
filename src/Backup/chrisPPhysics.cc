#ifndef __CINT__

#include "chrisPPhysics.h"

chrisPPhysics::chrisPPhysics() 
{
	TC_cut_min = 0;
	TC_cut_max = 352;
}

chrisPPhysics::~chrisPPhysics()
{
}

Bool_t	chrisPPhysics::Init()
{
	return kTRUE;
}

void	chrisPPhysics::Reconstruct()
{
}

// ----------------------------------------------------------------------------------------
// TH1 routines
// ----------------------------------------------------------------------------------------
void chrisPPhysics::FillScalers(Int_t low_scaler_number, Int_t high_scaler_number, TH1* hist)
{
	Int_t nFillScalers = high_scaler_number - low_scaler_number + 1;

	if( nFillScalers < hist->GetNbinsX())
	{
	    cout << "Error: FillScalers - histogram has insufficient bins for range" << endl;
	    return;
	}
	
	// To properly accumulate, create a histogram for this scaler read
	// cloning input histogram means the axis will be equivalent
	TH1* hist_current_SR = (TH1D*) hist->Clone();
	hist_current_SR->Reset();

    // Loop over scaler range, don't pull anything higher than the real # GetScalers()
	if (low_scaler_number  < 0)
	{
		cout << "FillScalers given scaler number outside range: " << low_scaler_number << endl;
		cout << "Setting lower limit to zero and continuing" << endl;
		low_scaler_number = 0;
	}
    if (high_scaler_number > GetScalers()->GetNScalers())
	{
		cout << "FillScalers given scaler number outside range: " << high_scaler_number << endl;
		cout << "Setting upper limit to "<< high_scaler_number << " and continuing" << endl;	
        high_scaler_number = GetScalers()->GetNScalers();
	}

    for (Int_t i = low_scaler_number; i <= high_scaler_number; i++)
	{
		Int_t bin = i - low_scaler_number;
        hist_current_SR->SetBinContent(bin,GetScalers()->GetScaler(i));
	}

	// Add to accumulated
	hist->Add(hist_current_SR);
}	


//EDIT Chris Mullen CAM 6/4/16 added the function to fill P4 px,py,pz,E 12 functions as follows the format below
// CAM adding functions for wire chamber stuff


//the following are functions for the TH1's which can be left out atm so I can compile they will be filled in later
//*****************************************************************************************************************************************


//END OF CHAMBER FUNCTIONS ****************************************************************************************************************


//ANGLE RECON





//******************************************************************************************************************************************





void chrisPPhysics::FillMissingMass(const GTreeParticle& tree, TH1* Hprompt, TH1* Hrandom)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        	FillMissingMass(tree, i, j, Hprompt, Hrandom);
	}
    }
}

void chrisPPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        	FillMissingMass(tree, particle_index, i, Hprompt, Hrandom);
	}
}

void chrisPPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom)
{
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

	if (GHistBGSub::IsPrompt(time)) Hprompt->Fill(missingp4.M());
	if (GHistBGSub::IsRandom(time)) Hrandom->Fill(missingp4.M());						
}

void chrisPPhysics::FillTime(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			Hist->Fill(time);
		}
	}
}

void chrisPPhysics::FillTime(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
    // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;
	
        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
	Hist->Fill(time);
	}
}

void chrisPPhysics::FillTimeCut(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                    time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) Hist->Fill(time);
		}
	}
}

void chrisPPhysics::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) Hist->Fill(time);
	}
}

void chrisPPhysics::FillMass(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        Hist->Fill(tree.GetMass(i));
	}
}

void chrisPPhysics::FillMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    Hist->Fill(tree.GetMass(particle_index));
}

void chrisPPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        FillBeamAsymmetry(tree, particle_index, i, Hprompt, Hrandom, MM_min, MM_max);
    }

}

void chrisPPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max)
{
    // Is tagger channel rejected by user?
//    cout << tagger->GetTaggedChannel(tagger_index) << " " << TC_cut_min << " " << TC_cut_max << endl;
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
//    cout << "time " << time << endl; 
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
 //   cout << "MM " << missingp4.M() << endl;     
    if((missingp4.M() < MM_min) || (missingp4.M() > MM_max)) return;
   
   if (GHistBGSub::IsPrompt(time)) Hprompt->Fill(tree.GetPhi(particle_index)); //cout << "prompt" << endl;}
   if (GHistBGSub::IsRandom(time)) Hrandom->Fill(tree.GetPhi(particle_index));	//cout << "random" << endl;}
}

Double_t chrisPPhysics::CalcCoplanarity(const GTreeParticle& tree1, Int_t particle_index1, const GTreeParticle& tree2, Int_t particle_index2)
{
   Double_t phi1 = tree1.GetPhi(particle_index1);
   Double_t phi2 = tree2.GetPhi(particle_index2);
   Double_t phidiff = TMath::Abs(phi1 - phi2);

   return phidiff;
}

//-----------------------------------------------------------------------------------------------
// Wire Chamber Cylinder (x vs y etc.)
//------------------------------------------------------------------------------------------------

void chrisPPhysics::CylindricalChamber(const GTreeParticle& treePi ,const GTreeMWPCHit& tree, TGraph* gGraph, Int_t ChamberNo, Int_t PlotNum, Bool_t TaggerBinning)
{

	//gGraph = new TGraph();
	for (Int_t i = 0; i < treePi.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			CylindricalChamber(tree, i, j, gGraph, ChamberNo, PlotNum,  TaggerBinning);
		}
	}




}

void chrisPPhysics::CylindricalChamber(const GTreeMWPCHit& tree,  Int_t particle_index, Int_t tagger_index, TGraph* gGraph, Int_t ChamberNo, Int_t PlotNum, Bool_t TaggerBinning)
{

	if( ChamberNo == 1 )
	{	if(PlotNum == 1)
		{	Double_t ChamberPosIX = GetMWPCHitsChris()->GetMWPCChamber1X(particle_index);
			Double_t ChamberPosIY = GetMWPCHitsChris()->GetMWPCChamber1Y(particle_index);
			gGraph->SetPoint(35000, ChamberPosIX, ChamberPosIY);
			//gHist->Fill(ChamberPosI);
		}
		else if(PlotNum == 2)
		{	//Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber1Y(particle_index);
			//gHist->Fill(ChamberPosI);
		}
		else if(PlotNum == 3)
		{	//Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber1Z(particle_index);
			//gHist->Fill(ChamberPosI);
		}
		else
		{
			std::cout << " Incorrect Chamber Axis Provided for chamber 1 " <<std::endl;
		}

	}

	else if( ChamberNo == 2 )
	{	if(PlotNum == 1)
		{	//Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber2X(particle_index);
			//gHist->Fill(ChamberPosI);
		}
		else if(PlotNum == 2)
		{	//Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber2Y(particle_index);
			//gHist->Fill(ChamberPosI);
		}
		else if(PlotNum == 3)
		{	//Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber2Z(particle_index);
			//gHist->Fill(ChamberPosI);
		}

		else
		{
			std::cout << " Incorrect Chamber Axis Provided for chamber 2 " <<std::endl;
		}

	}
	else 
	{
		std::cout << " Incorrect Chamber Number Provided " <<std::endl;
	}


	//gGraph->Write("MyGraph"); 

}

// ----------------------------------------------------------------------------------------
// GH1 routines
// ----------------------------------------------------------------------------------------


//EDIT Chris Mullen CAM 6/4/16 added the function to fill P4 px,py,pz,E 12 functions as follows the format below
// CAM 17/5/16 added functions for wire chamber
// CAM 18/5/16 added functions for reconstructing angles


void chrisPPhysics::FillProtonMM(const GTreeParticle& treePi,const GTreeParticle& treePro, GH1* gHist, Bool_t TaggerBinning)
{
	

	for (Int_t i = 0; i < treePi.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillProtonMM(treePi, i, j, gHist, TaggerBinning);
		}
	}

//	missingP4pion = CalcMissingP4(treePi, particle_index,tagger_index); 
//	gHist->Fill();

}

/*
void chrisPPhysics::FillAngle(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillPx(tree, particle_index, i, gHist, TaggerBinning);
	}
}

*/

void chrisPPhysics::FillProtonMM(const GTreeParticle& treePi, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - treePi.GetTime(particle_index);
    
    // calc missing p4
    PRecon = CalcMissingP4(treePi, particle_index,tagger_index);

   // Fill GH1
   if(TaggerBinning)   gHist->Fill(PRecon.M(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(PRecon.M(),time);

}



void chrisPPhysics::FillProton(const GTreeParticle& treePi,const GTreeParticle& treePro, GH1* gHist, Int_t ProtonElem,Bool_t TaggerBinning)
{
	
	for (Int_t i = 0; i < treePi.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillProton(treePi, i, j, gHist, ProtonElem, TaggerBinning);
		}
	}
}



void chrisPPhysics::FillProton(const GTreeParticle& treePi, Int_t particle_index, Int_t tagger_index, GH1* gHist,Int_t ProtonElem,Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - treePi.GetTime(particle_index);
    
    // calc missing p4
    PRecon = CalcMissingP4(treePi, particle_index,tagger_index);

   // Fill GH1

   if(ProtonElem ==1)
	{
  	 	if(TaggerBinning)   gHist->Fill(PRecon.X(),time, GetTagger()->GetTaggedChannel(tagger_index));
   		else gHist->Fill(PRecon.X(),time);

	}

   else if(ProtonElem==2)
	{

   		if(TaggerBinning)   gHist->Fill(PRecon.Y(),time, GetTagger()->GetTaggedChannel(tagger_index));
   		else gHist->Fill(PRecon.Y(),time);
		
	}

   else if(ProtonElem==3)
	{

   		if(TaggerBinning)   gHist->Fill(PRecon.Z(),time, GetTagger()->GetTaggedChannel(tagger_index));
   		else gHist->Fill(PRecon.Z(),time);


	}
   else
	{
		std::cout<<"The proton vector element provided does not exist. Please select between 1,2 or 3 for x,y or z respectively."<<std::endl;
	}

}





void chrisPPhysics::FillProtonPhi(const GTreeParticle& treePi,const GTreeParticle& treePro, GH1* gHist, Bool_t TaggerBinning)
{
	
	for (Int_t i = 0; i < treePi.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillProtonPhi(treePi, i, j, gHist, TaggerBinning);
		}
	}
}



void chrisPPhysics::FillProtonPhi(const GTreeParticle& treePi, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - treePi.GetTime(particle_index);
    
    // calc missing p4
    PRecon = CalcMissingP4(treePi, particle_index,tagger_index);
    Double_t ph = PRecon.Phi();
    Double_t phi = ph  * TMath::RadToDeg();
   
    

   // Fill GH1
   if(TaggerBinning)   gHist->Fill(PRecon.Phi(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(phi,time);

}







void chrisPPhysics::FillProtonTheta(const GTreeParticle& treePi,const GTreeParticle& treePro, GH1* gHist, Bool_t TaggerBinning)
{
	
	for (Int_t i = 0; i < treePi.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillProtonTheta(treePi, i, j, gHist, TaggerBinning);
		}
	}
}



void chrisPPhysics::FillProtonTheta(const GTreeParticle& treePi, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - treePi.GetTime(particle_index);
    
    // calc missing p4
    PRecon = CalcMissingP4(treePi, particle_index,tagger_index);
    Double_t th = PRecon.Theta();
    Double_t Theta = th  * TMath::RadToDeg();



   // Fill GH1
   if(TaggerBinning)   gHist->Fill(PRecon.Theta(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(Theta,time);

}





void chrisPPhysics::FillProtonAngleDiffPhi(const GTreeParticle& treePi,const GTreeParticle& treePro, GH1* gHist, Bool_t TaggerBinning)
{
	
	for (Int_t i = 0; i < treePi.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillProtonAngleDiffPhi(treePi, i, j, gHist, TaggerBinning);
		}
	}
}



void chrisPPhysics::FillProtonAngleDiffPhi(const GTreeParticle& treePi, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - treePi.GetTime(particle_index);
    
    // calc missing p4
    PRecon = CalcMissingP4(treePi, particle_index,tagger_index);

   // Fill GH1
   //if(TaggerBinning)   gHist->Fill(PRecon.Angle(),time, GetTagger()->GetTaggedChannel(tagger_index));
   //else gHist->Fill(PRecon.Angle(),time);

}



//******************************************************************************************************************************************



void chrisPPhysics::FillProtonNormal(const GTreeParticle& treePro,const GTreeParticle& treePi, GH1* gHist, Int_t ProtonElem, Bool_t TaggerBinning)
{

	for (Int_t i = 0; i < treePro.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillProtonNormal(treePro, i, j, gHist, ProtonElem, TaggerBinning);
		}
	}



}


void chrisPPhysics::FillProtonNormal(const GTreeParticle& treePro, Int_t particle_index, Int_t tagger_index, GH1* gHist, Int_t ProtonElem, Bool_t TaggerBinning)
{
 
				    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - treePro.GetTime(particle_index);
    
    // calc missing p4
    PNormal = treePro.Particle(particle_index);//CalcMissingP4(treePi, particle_index,tagger_index);

   // Fill GH1

  if(ProtonElem ==1)
	{
  	 	if(TaggerBinning)   gHist->Fill(PNormal.X(),time, GetTagger()->GetTaggedChannel(tagger_index));
   		else gHist->Fill(PNormal.X(),time);

	}

   else if(ProtonElem==2)
	{

   		if(TaggerBinning)   gHist->Fill(PNormal.Y(),time, GetTagger()->GetTaggedChannel(tagger_index));
   		else gHist->Fill(PNormal.Y(),time);
		
	}

   else if(ProtonElem==3)
	{

   		if(TaggerBinning)   gHist->Fill(PNormal.Z(),time, GetTagger()->GetTaggedChannel(tagger_index));
   		else gHist->Fill(PNormal.Z(),time);


	}
   else
	{
		std::cout<<"The proton vector element provided does not exist. Please select between 1,2 or 3 for x,y or z respectively."<<std::endl;
	}


}






void chrisPPhysics::FillProtonNormalMass(const GTreeParticle& treePro,const GTreeParticle& treePi, GH1* gHist, Bool_t TaggerBinning)
{

	for (Int_t i = 0; i < treePro.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillProtonNormalMass(treePro, i, j, gHist, TaggerBinning);
		}
	}



}


void chrisPPhysics::FillProtonNormalMass(const GTreeParticle& treePro, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
 
		    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - treePro.GetTime(particle_index);
    
    // calc missing p4
    PNormal = treePro.Particle(particle_index);//CalcMissingP4(treePi, particle_index,tagger_index);

   // Fill GH1
   if(TaggerBinning)   gHist->Fill(PNormal.M(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(PNormal.M(),time);


}





void chrisPPhysics::FillProtonNormalPhi(const GTreeParticle& treePro,const GTreeParticle& treePi, GH1* gHist, Bool_t TaggerBinning)
{

	for (Int_t i = 0; i < treePro.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillProtonNormalPhi(treePro, i, j, gHist, TaggerBinning);
		}
	}



}


void chrisPPhysics::FillProtonNormalPhi(const GTreeParticle& treePro, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
 
				    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - treePro.GetTime(particle_index);
    
    // calc missing p4
    PNormal = treePro.Particle(particle_index);//CalcMissingP4(treePi, particle_index,tagger_index);
    Double_t ph = PNormal.Phi();
    Double_t phi = ph  * TMath::RadToDeg();
   // Fill GH1
   if(TaggerBinning)   gHist->Fill(PNormal.Phi(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(phi,time);


}




void chrisPPhysics::FillProtonNormalTheta(const GTreeParticle& treePro,const GTreeParticle& treePi, GH1* gHist, Bool_t TaggerBinning)
{

	for (Int_t i = 0; i < treePro.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillProtonNormalTheta(treePro, i, j, gHist, TaggerBinning);
		}
	}



}


void chrisPPhysics::FillProtonNormalTheta(const GTreeParticle& treePro, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
 
				    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - treePro.GetTime(particle_index);
    
    // calc missing p4
    PNormal = treePro.Particle(particle_index);//CalcMissingP4(treePi, particle_index,tagger_index);
    Double_t th = PNormal.Theta();
    Double_t Theta = th  * TMath::RadToDeg();

   // Fill GH1
   if(TaggerBinning)   gHist->Fill(PNormal.Theta(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(Theta,time);


}




void chrisPPhysics::FillProtonNormalAngleDiffTheta(const GTreeParticle& treePro,const GTreeParticle& treePi, GH1* gHist, Bool_t TaggerBinning)
{

	for (Int_t i = 0; i < treePro.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillProtonNormalAngleDiffTheta(treePro, i, j, gHist, TaggerBinning);
		}
	}



}


void chrisPPhysics::FillProtonNormalAngleDiffTheta(const GTreeParticle& treePro, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
 
				    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - treePro.GetTime(particle_index);
    
    // calc missing p4
    PNormal = treePro.Particle(particle_index);//CalcMissingP4(treePi, particle_index,tagger_index);

   // Fill GH1
   //if(TaggerBinning)   gHist->Fill(PNormal.Angle(),time, GetTagger()->GetTaggedChannel(tagger_index));
   //else gHist->Fill(PNormal.Angle(),time);


}




//END OF ANgLE FUNCTIONS *****************************************************************************************





void chrisPPhysics::FillnChamber(const GTreeMWPCHit& tree, GH1* gHist, Int_t ChamberNo)
{

   if (ChamberNo ==1){

   Int_t ChamberHits = GetMWPCHitsChris()->GetNMWPCHitsChrisChamber1();
   gHist->Fill(ChamberHits);

}

   else if (ChamberNo == 2){

   Int_t ChamberHits = GetMWPCHitsChris()->GetNMWPCHitsChrisChamber2();
   gHist->Fill(ChamberHits);

}
 
}


void chrisPPhysics::FillChamber(const GTreeParticle& treePi,const GTreeMWPCHit& tree, GH1* gHist, Int_t ChamberNo, Int_t ChamberAxis,Bool_t TaggerBinning)
{

	for (Int_t i = 0; i < treePi.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillChamber(tree, i, j, gHist, ChamberNo, ChamberAxis, TaggerBinning);
		}
	}

}


void chrisPPhysics::FillChamber(const GTreeMWPCHit& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Int_t ChamberNo, Int_t ChamberAxis, Bool_t TaggerBinning)
{
 
	if( ChamberNo == 1 )
	{	if(ChamberAxis == 1)
		{	Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber1X(particle_index);
			gHist->Fill(ChamberPosI);
		}
		else if(ChamberAxis == 2)
		{	Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber1Y(particle_index);
			gHist->Fill(ChamberPosI);
		}
		else if(ChamberAxis == 3)
		{	Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber1Z(particle_index);
			gHist->Fill(ChamberPosI);
		}
		else
		{
			std::cout << " Incorrect Chamber Axis Provided for chamber 1 " <<std::endl;
		}

	}

	else if( ChamberNo == 2 )
	{	if(ChamberAxis == 1)
		{	Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber2X(particle_index);
			gHist->Fill(ChamberPosI);
		}
		else if(ChamberAxis == 2)
		{	Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber2Y(particle_index);
			gHist->Fill(ChamberPosI);
		}
		else if(ChamberAxis == 3)
		{	Double_t ChamberPosI = GetMWPCHitsChris()->GetMWPCChamber2Z(particle_index);
			gHist->Fill(ChamberPosI);
		}

		else
		{
			std::cout << " Incorrect Chamber Axis Provided for chamber 2 " <<std::endl;
		}

	}
	else 
	{
		std::cout << " Incorrect Chamber Number Provided " <<std::endl;
	}



}







//END OF CHAMBER FUNCTIONS *****************************************************************************************




void chrisPPhysics::FillP(const GTreeParticle& tree, GH1* gHist, Int_t VecElem,  Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillP(tree, i, j, gHist, VecElem, TaggerBinning);
		}
	}
}


void chrisPPhysics::FillP(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Int_t VecElem,  Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillP(tree, particle_index, i, gHist, VecElem, TaggerBinning);
	}
}



void chrisPPhysics::FillP(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Int_t VecElem, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

   // Fill GH1

   if ( VecElem==1 ) //case for Px
	{
	   if(TaggerBinning)   gHist->Fill(missingp4.Px(),time, GetTagger()->GetTaggedChannel(tagger_index));
	   else gHist->Fill(missingp4.Px(),time);
	}

   else if (VecElem==2) //case for Py
	{
	   if(TaggerBinning)   gHist->Fill(missingp4.Py(),time, GetTagger()->GetTaggedChannel(tagger_index));
	   else gHist->Fill(missingp4.Py(),time);
	
	}

   else if (VecElem==3) //case for Pz
	{
	   if(TaggerBinning)   gHist->Fill(missingp4.Pz(),time, GetTagger()->GetTaggedChannel(tagger_index));
	   else gHist->Fill(missingp4.Pz(),time);

	}

   else 
	{
		std::cout << " The vector element provided does not exist. Please select between 1,2 or 3 for x,y or z respectively. "<< std::endl;
	}
}


//END OF CHAMBER FUNCTIONS *****************************************************************************************


void chrisPPhysics::FillE(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillE(tree, i, j, gHist, TaggerBinning);
		}
	}
}


void chrisPPhysics::FillE(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillE(tree, particle_index, i, gHist, TaggerBinning);
	}
}



void chrisPPhysics::FillE(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

   // Fill GH1
   if(TaggerBinning)   gHist->Fill(missingp4.E(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(missingp4.E(),time);

}




void chrisPPhysics::FillMissingMass(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillMissingMass(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void chrisPPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillMissingMass(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void chrisPPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

   // Fill GH1
   if(TaggerBinning)   gHist->Fill(missingp4.M(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(missingp4.M(),time);

}

Double_t chrisPPhysics::CalcMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.M();
}

Double_t chrisPPhysics::CalcMissingEnergy(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree,particle_index, tagger_index);

	return missingp4.T();
}

TLorentzVector chrisPPhysics::CalcMissingP4(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    particle	= tree.Particle(particle_index);
    beam 		= TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));
	missingp4 	= beam + target - particle;						

	return missingp4;
}

void chrisPPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning, Double_t MM_min, Double_t MM_max)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        FillBeamAsymmetry(tree, particle_index, i, gHist, TaggerBinning);
    }

}

void chrisPPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning, Double_t MM_min, Double_t MM_max)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
    if((missingp4.M() < MM_min) || (missingp4.M() > MM_max)) return;
   
   if(TaggerBinning)   gHist->Fill(tree.GetPhi(particle_index),time,GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(tree.GetPhi(particle_index),time);

}

void chrisPPhysics::FillTime(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			gHist->Fill(time);
		}
	}
}

void chrisPPhysics::FillTime(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;
	
        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		gHist->Fill(time);
	}
}

void chrisPPhysics::FillTimeCut(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                    time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
		}
	}
}

void chrisPPhysics::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
	}
}

void chrisPPhysics::FillMass(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        gHist->Fill(tree.GetMass(i));
	}
}

void chrisPPhysics::FillMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    gHist->Fill(tree.GetMass(particle_index));
}

Bool_t 	chrisPPhysics::Write()
{
	return kTRUE;
}

// Some common initialisation stuff
Bool_t 	chrisPPhysics::InitBackgroundCuts()
{
	// Set background cuts
	Double_t p1, p2, r1, r2;
	string config = ReadConfig("Set-Prompt-Cut");
	if(strcmp(config.c_str(), "nokey") == 0) 
		cout << "No BG subtraction - At least 1 prompt and random cut required" << endl;
	else if(sscanf( config.c_str(), "%lf %lf\n", &p1, &p2) == 2)
	{
	   config = ReadConfig("Add-Random-Cut",0);
	   if(strcmp(config.c_str(), "nokey") == 0) 
	   	cout << "No BG subtraction - At least 1 random cut required" << endl;
	   else if(sscanf( config.c_str(), "%lf %lf\n", &r1, &r2) == 2)
	   {
		cout << "Init BG cuts:" << endl;
		cout << "prompt(" << p1 << "," << p2 << ") ";
		cout << "random(" << r1 << "," << r2 << ") " << endl;

		GHistBGSub::InitCuts(p1,p2,r1,r2);

		// Look for additional random windows
		Int_t instance = 1;
		do
		{
			config = ReadConfig("Add-Random-Cut",instance);
			if(sscanf( config.c_str(), "%lf %lf\n", &r1, &r2) == 2)
			{
				cout << "Adding random cuts: ";
				cout << "random(" << r1 << "," << r2 << ") " << endl;

				GHistBGSub::AddRandCut(r1,r2);
			}		
			instance++;
		} while (strcmp(config.c_str(), "nokey") != 0);
	   }
	   else {cout << "Random window not set correctly" << endl; return kFALSE;}
	}
	else {cout << "Prompt window not set correctly" << endl; return kFALSE;}

	cout << endl;
	return kTRUE;

}

Bool_t 	chrisPPhysics::InitTargetMass()
{
	Double_t mass;
	string config = ReadConfig("Target-Mass");
	if(strcmp(config.c_str(), "nokey") == 0)
	{
		cout << "Target mass unknown!" << endl;
	}
	else if(sscanf( config.c_str(), "%lf\n", &mass) == 1)
	{
		cout << "Setting Target mass: " << mass << " MeV" << endl;
		SetTarget(mass);		
	}
	else 
	{
		cout << "Target Mass not set correctly" << endl; 
		return kFALSE;
	}

	cout << endl;
	return kTRUE;

}

Bool_t 	chrisPPhysics::InitTaggerChannelCuts()
{
	Double_t tc1, tc2;
	string config = ReadConfig("Tagger-Channel-Cut");
	if(sscanf( config.c_str(), "%lf %lf\n", &tc1, &tc2) == 2)
	{
		if ((tc1 < 0) || (tc1 > 352))
		{
           cout << "Invalid tagger channel cut: " << tc1 << endl;
		   return kFALSE;
		}
		else if ((tc2 < 0) || (tc2 > 352))
		{
           cout << "Invalid tagger channel cut: " << tc2 << endl;
		   return kFALSE;
		}
		
        cout << "Setting cut on tagger channels: " << tc1 << " to " << tc2 << endl;
		SetTC_cut(tc1,tc2);
	}
	else if(strcmp(config.c_str(), "nokey") != 0)
	{
		cout << "Tagger Channel cut not set correctly" << endl; 
		return kFALSE;
	}

	cout << endl;
	return kTRUE;

}

Bool_t 	chrisPPhysics::InitTaggerScalers()
{
	Int_t sc1, sc2;
	string config = ReadConfig("Tagger-Scalers");
	if(sscanf( config.c_str(), "%d %d\n", &sc1, &sc2) == 2)
	{
		cout << "Setting Tagger scaler channels: " << sc1 << " to " << sc2 << endl;
        SetTC_scalers(sc1,sc2);
	}
	else if(strcmp(config.c_str(), "nokey") != 0)
	{
        cout << "Tagger Channel GetScalers() not set correctly" << endl;
		return kFALSE;
	}

	cout << endl;
	return kTRUE;

}
#endif
