#include "chrisPPi0Example.h"

chrisPPi0Example::chrisPPi0Example()
{ 
    nChamberHitsin1 = new GH1("nChamberHitsin1", 	"nChamberHitsin1",	10,	0, 9);
    nChamberHitsin2 = new GH1("nChamberHitsin2", 	"Chamber2",	18,	0, 17);

  treeselected =  new TTree("Selected","Event selection tree"); // selected tree   

//  treeselected->Branch("Pion.",&PionCan);
//  treeselected->Branch("Chamber1.",&Chamber1_Vec);
//  treeselected->Branch("Chamber2.",&Chamber2_Vec);
//  treeselected->Branch("Phidiff.",&Phidiff);       //Wire chamber phi minus PID phi
//  treeselected->Branch("Ebeam.",&energy_beam);
//  treeselected->Branch("ESum.",&energySum);
//  treeselected->Branch("Multiplicity.",&multiplicity);
  treeselected->Branch("Inv_Mass_Pion.",&inv_M_value);
//  treeselected->Branch("PidIndex.",&PidHitIndex);
//  treeselected->Branch("Eventno.",&EventNumber);
//  treeselected->Branch("PIDPhi.",&PIDPhi);
//  treeselected->Branch("Chamber1Phi.",&Chamber1_VecPhi);
//  treeselected->Branch("ProtonCandidate.",&ProtonCan);
//  treeselected->Branch("MissingMass.",&MissingM);
//  treeselected->Branch("BeamHelicity.",&BeamHelicity);
//  treeselected->Branch("TaggedTime.",&time_beam);
//  treeselected->Branch("TestChamber.",&TestChamber);
//  treeselected->Branch("PromptRegion.",&PromptRegion);
//  treeselected->Branch("RandomRegion.",&RandomRegion);

}

chrisPPi0Example::~chrisPPi0Example()
{
}

Bool_t	chrisPPi0Example::Init()
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

	target.SetXYZM(0.0,0.0,0,1875.613);	//NEEDS CHANGING currently deuteron
	if(!InitBackgroundCuts()) return kFALSE;
	if(!InitTargetMass()) return kFALSE;
	if(!InitTaggerChannelCuts()) return kFALSE;
	if(!InitTaggerScalers()) return kFALSE;
	cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	chrisPPi0Example::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();

    return kTRUE;
}

void	chrisPPi0Example::ProcessEvent()
{

  
  if(usePeriodMacro == 1)
    {
      if(GetEventNumber() % period == 0)
	cout << "Events: " << GetEventNumber() << "  Events Accepted: " << nEventsWritten << endl;
    }

    FillnChamber(*GetMWPCHitsChris(),nChamberHitsin1,1);
    FillnChamber(*GetMWPCHitsChris(),nChamberHitsin2,2);
  // Loop over tagged events 
  for(Int_t i=0; i<=GetTagger()->GetNTagged() ;i++)
    {
     

      if((GetRootinos()->GetNParticles()==1 && GetPhotons()->GetNParticles()==2)  ) 
 //     if((GetPhotons()->GetNParticles() != 3 )  ) 
      
	{
								


} //closing rootino = 1 if.

       else if ( (GetPhotons()->GetNParticles()==3))
	{
	energy_beam = GetTagger()->GetTaggedEnergy(i);
	beam.SetXYZM(0.,0.,energy_beam,0.);  
	time_beam = GetTagger()->GetTaggedTime(i);

     // Loop over the different combination of photons and set the different vectors for each
            for(Int_t l=0;l<3;l++){
		if(l==0) { l_g1 = 1; l_g2 = 2;
	  	g1_Vec1 = GetPhotons()->Particle(l_g1);
	  	g2_Vec1 = GetPhotons()->Particle(l_g2);
	  	C1 = g1_Vec1 + g2_Vec1;  //pion candidate
	  	IMDiff1 = TMath::Abs(C1.M() - 134.97)/134.97;
	  	P1 = GetPhotons()->Particle(0); // proton candidate
	  	//CopDiff1 = TMath::Abs(C1.TLorentzVector::DeltaPhi(-P1)*TMath::RadToDeg() )/180;
		Mp41 = beam + target - C1;
		MM1 = Mp41.M();
		MMDiff1 = TMath::Abs(MM1 - 1912)/1912; //What should the missing mass value be?
	  	SumDiff1 = (IMDiff1) + (MMDiff1);

		}
		if(l==1) { l_g1 = 2; l_g2 = 0;
	  	g1_Vec2 = GetPhotons()->Particle(l_g1);
	  	g2_Vec2 = GetPhotons()->Particle(l_g2);
	  	C2 = g1_Vec2 + g2_Vec2;
	  	IMDiff2 = TMath::Abs(C2.M() - 134.97)/134.97;
	  	P2 = GetPhotons()->Particle(1);
	  	//CopDiff2 = TMath::Abs(C2.TLorentzVector::DeltaPhi(-P2)*TMath::RadToDeg() )/180;
		Mp42 = beam + target - C2;
		MM2 = Mp42.M();
		MMDiff2 = TMath::Abs(MM2 - 1912)/1912; //What should the missing mass value be? should i divide the difference by the peak position to get them both scaled percentage wise.
	  	SumDiff2 = (IMDiff2) + (MMDiff2) ;

		}
		if(l==2) { l_g1 = 0; l_g2 = 1;
	  	g1_Vec3 = GetPhotons()->Particle(l_g1);
	  	g2_Vec3 = GetPhotons()->Particle(l_g2);
	  	C3 = g1_Vec3 + g2_Vec3;
	  	IMDiff3 = TMath::Abs(C3.M() - 134.97)/134.97;
	  	P3 = GetPhotons()->Particle(2);
	  	//CopDiff3 = TMath::Abs(C3.TLorentzVector::DeltaPhi(-P3)*TMath::RadToDeg() )/180;
		Mp43 = beam + target - C3;
		MM3 = Mp43.M();
		MMDiff3 = TMath::Abs(MM3 - 1912)/1912; //What should the missing mass value be?
	  	SumDiff3 = (IMDiff3) + (MMDiff3) ;
		}
	
      	   } //closing for candidate loop
      
      		if( SumDiff1 < SumDiff2 && SumDiff1 < SumDiff3){
		PionCan.SetXYZM(C1.X(),C1.Y(),C1.Z(),C1.M());
		ProtonCan = P1;
      	  	}
      		else if(SumDiff2 < SumDiff1 && SumDiff2 < SumDiff3){
		PionCan.SetXYZM(C2.X(),C2.Y(),C2.Z(),C2.M());
		ProtonCan = P2;
	      }
	      else if(SumDiff3 < SumDiff1 && SumDiff3 < SumDiff2){
		PionCan.SetXYZM(C3.X(),C3.Y(),C3.Z(),C3.M());
		ProtonCan = P3;
	      }  
     
      	// Invariant Mass
      	inv_M_value = PionCan.M();
	treeselected->Fill();
		//}
		//}

	} //closing elif
	

} //closing for loop
} //closing function


void	chrisPPi0Example::ProcessScalerRead()
{
	// Fill Tagger Scalers
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	chrisPPi0Example::Write()
{

//    tree1r->Write();
//    tree1r->Reset();
	//treeselected->Reset();
//    tree3g->Write();
//    tree3g->Reset();
//    return 0;

    // Write all GH1's and TObjects defined in this class
    return GTreeManager::Write();
}


