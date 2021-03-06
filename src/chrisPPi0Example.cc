#include "chrisPPi0Example.h"

chrisPPi0Example::chrisPPi0Example()
{ 
 
//CAM 17/5/16 Adding wire chamber hits
    nChamberHitsin1 = new GH1("nChamberHitsin1", 	"nChamberHitsin1",	10,	0, 9);
    nChamberHitsin2 = new GH1("nChamberHitsin2", 	"Chamber2",	18,	0, 17);

   // Create Tree
  tree3g =  new TTree("Tree3gamma","Event selection tree"); //3 gamma tree 
  treetrigger =  new TTree("Trigger","Event selection tree"); //trigger tree
  treetest =  new TTree("TreeTest","Event selection tree"); // test tree 
  treescatter =  new TTree("Scatter","Event selection tree"); // scatter tree 

  treeselected =  new TTree("Selected","Event selection tree"); // selected tree   

  treepid = new TTree("Pid","Event selection tree"); //pid tree



  // Define branches of TTree 3 gamma
  tree3g->Branch("Beam_Energy",&energy_beam);       // Beam Energy distribution
  tree3g->Branch("inv_M",&inv_M_value);             // Invariant mass 
  tree3g->Branch("MissingM",&MissingM);             // Missing mass distribution
  tree3g->Branch("MissingP4",&missingp4);           // Missing mass
  tree3g->Branch("Pion_phi",&phiPi);                // Pi0 phi
  tree3g->Branch("Proton_phi",&phiProton);          // Missing Mass phi
  tree3g->Branch("CoplanarityMeasure",&coplanarityMeasure);       // Coplanarity between Proton and Pi0 
  tree3g->Branch("Pion_Theta",&thetaPi);            // Pi0 theta
  tree3g->Branch("Proton_Theta",&thetaProton);      // Missing Mass theta
  tree3g->Branch("Theta_Difference",&thetadiff);    // Difference between Pi0 and Missing Mass theta
  tree3g->Branch("Proton_Candidate",&ProtonCan);    // The Proton Candidate 4 vector
  tree3g->Branch("Pion_Candidate",&PionCan);	    // The Pion Candidate 4 vector
  tree3g-> Branch("Num_PID_Hits",&NPidhits);	    // Number of Pid Hits
  tree3g->Branch("PIDPhi",&PIDPhi);		    // The PID phi 
  tree3g->Branch("MWPCChamber1",&Chamber1_Vec);	    // The inner wire chamber 
  tree3g->Branch("MWPCChamber2",&Chamber2_Vec);     // The outer wire chamber
  tree3g->Branch("ReconProton",&ReconProton);	    // The Reconstructed Proton from the missingp4
  tree3g->Branch("Phi_ReconProton",&phiReconProton); // The Reconstructed Proton Phi
  tree3g->Branch("Theta_ReconProton",&thetaReconProton); // The Reconstructed Proton theta
  tree3g->Branch("Coplanarity_Recon",&coplanarityRecon); // Coplanarity between Pi0 and Recon Proton(missing p4)
  tree3g->Branch("ProtonEnergy", &ProtonEnergy);     // Proton Energy 
  tree3g->Branch("energySum",&energySum);		// Energy Sum from trigger
  tree3g->Branch("multiplicity",&multiplicity);		//multiplicity from trigger
  tree3g->Branch("Boosted_Proton",&ProtonBoost);	// Proton in COM frame
  tree3g->Branch("Boosted_Pion",&PionBoost);		// Pion in COM frame
  tree3g->Branch("Boosted_Recostructed_Proton",&ReconProtonBoost); // Recon Proton in COM frame
  tree3g->Branch("Pvector1",&PVector1); 
  tree3g->Branch("Pvector2",&PVector2); 
  tree3g->Branch("NewEnergy",&NewEnergy);
  tree3g->Branch("DetectorNumber",&detectornum);

  // Define branches of TTree trigger
  treetrigger->Branch("energySum",&energySum);
  treetrigger->Branch("multiplicity",&multiplicity);

  // Define branches of TTree test
  treetest->Branch("PIDPhi",&PIDPhi);
  treetest->Branch("NTracks", &test);


  //Define branches of TTree scatter
  treescatter->Branch("Pion",&PionCan);
  treescatter->Branch("Chamber1",&Chamber1_Vec);
  treescatter->Branch("Chamber2",&Chamber1_Vec);


  //Define branches of TTree selected
  treeselected->Branch("Pion.",&PionCan);
  treeselected->Branch("Chamber1.",&Chamber1_Vec);
  treeselected->Branch("Chamber2.",&Chamber2_Vec);
  treeselected->Branch("Phidiff.",&Phidiff);       //Wire chamber phi minus PID phi
  treeselected->Branch("Ebeam.",&energy_beam);
  treeselected->Branch("ESum.",&energySum);
  treeselected->Branch("Multiplicity.",&multiplicity);
  treeselected->Branch("Inv_Mass_Pion.",&inv_M_value);
  treeselected->Branch("PidIndex.",&PidHitIndex);
  treeselected->Branch("Eventno.",&EventNumber);
  treeselected->Branch("PIDPhi.",&PIDPhi);
  treeselected->Branch("Chamber1Phi.",&Chamber1_VecPhi);
  treeselected->Branch("ProtonCandidate.",&ProtonCan);
  treeselected->Branch("MissingMass.",&MissingM);
  treeselected->Branch("BeamHelicity.",&BeamHelicity);
  treeselected->Branch("TaggedTime.",&time_beam);

//  treeselected->Branch("PromptRegion.",&PromptRegion);
//  treeselected->Branch("RandomRegion.",&RandomRegion);

//tree pidhits
  treepid->Branch("PidHits",&PidHitIndex);


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
	//CAM 17/5/16 Fill wire chamber info.

//	For FillChamber the extra parameter is an integer for the axis of the wc 1=x 2=y 3=z. The first param is the chamber number.
//      For FillnChamber the extra parameter is an integer for the Chamber number 

   

    FillnChamber(*GetMWPCHitsChris(),nChamberHitsin1,1);
    FillnChamber(*GetMWPCHitsChris(),nChamberHitsin2,2);
  
	// Does this take account of multiple hits in the wire chamber in the same event  
    Int_t iii = 0; // for the wire chamber info


	//Constructing a wire chamber vector for the hit position(These won't fill fully as the fill command is in the if statement for 3gamma events)
								//(Need to create own tree or move to three gamma tree).


if (GetMWPCHitsChris()->GetNMWPCHitsChrisChamber1() != 0){

    Chamber1_Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber1X(iii),GetMWPCHitsChris()->GetMWPCChamber1Y(iii),GetMWPCHitsChris()->GetMWPCChamber1Z(iii));

}
else
{
	
    Chamber1_Vec.SetXYZ(-1000,-1000,-1000);

}


if (GetMWPCHitsChris()->GetNMWPCHitsChrisChamber2() != 0){


    Chamber2_Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber2X(iii),GetMWPCHitsChris()->GetMWPCChamber2Y(iii),GetMWPCHitsChris()->GetMWPCChamber2Z(iii));
}
else
{
	
    Chamber2_Vec.SetXYZ(-1000,-1000,-1000);

}


    targetPosition.SetXYZ(0.,0.,0);   //NEEDS CHANGING and the one at the top

if  (Chamber1_Vec.X()==-1000){ //Can I get rid of all events if theres no hit in wc1 cos I wont be able to form the vectors properly.


    PVector1.SetXYZ(-1000,-1000,-1000);
    PVector2.SetXYZ(-1000,-1000,-1000);


}
else{

    PVector1 = Chamber1_Vec - targetPosition;
    PVector2 = Chamber2_Vec - PVector1;
}


	//PID info
    NPidhits = GetDetectorHits()->GetNPIDHits();

   EventNumber =GetEventNumber() ;

//std::cout << GetDetectorHits()->GetMWPCHits(0)<< std::endl;

 

if (NPidhits>0){
   PidHitIndex = GetDetectorHits()->GetPIDHits(0);    //The parameter in the pidhits is the number of pid hits in the event.
   treepid->Fill();


//PID plan create an array of length 24 to store the elements of the PID phi.  Put this in the header later andd tidy header.


        Double_t PIDElemPhi[24] = {8.574,22.974,37.379,51.784,66.188,80.593,94.997,109.402,123.806,138.211,152.615,167.02,-178.93,-163.16,-147.39,-131.62,-115.85,-100.08,-84.31,-68.54,-52.77,-37.01,-21.24,-5.47}; 


        PIDPhi = PIDElemPhi[PidHitIndex]; //here

if  (Chamber1_Vec.X()==-1000){


	Chamber1_VecPhi =-10000 ;
}

else{

	Chamber1_VecPhi =Chamber1_Vec.Phi()*TMath::RadToDeg() ;

}

	PIDPhi = PIDPhi + 180;
	Chamber1_VecPhi =Chamber1_VecPhi +180;

	Phidiff = TMath::Abs( PIDPhi - Chamber1_VecPhi) ;   
        if (Phidiff >180){

	Phidiff = 360 -Phidiff;


} //closing phidiff
} //closingpihits

	else{

	Phidiff = -10;
        PIDPhi = -10;
        Chamber1_VecPhi = -10;


}


 	



	//Int_t PIDIndices = GetGTreeA2Geant()->GetPIDHitIndices();

	test = GetTracks()->GetNTracks();

//	std::cout <<test <<std::endl;

	treetest->Fill();



//*************************************************TRIGGER*************************************

	energySum =GetTrigger()->GetEnergySum();
	multiplicity = GetTrigger()->GetMultiplicity();




	treetrigger->Fill();
 
//  ************************************************************************************************************
  // Loop over tagged events 
  for(Int_t i=0; i<=GetTagger()->GetNTagged() ;i++)
    {
      


      // Select number and type of events. Event selection loop.

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
	  
     	
	    
	//Invariant mass cut here! that can be commented out to get both plots.
	if (inv_M_value > 100.0 && inv_M_value <170.0){    


//	energy_beam = GetTagger()->GetTaggedEnergy(i);
//	beam.SetXYZM(0.,0.,energy_beam,0.);
		    
	missingp4     = beam + target - PionCan;
	MissingM = missingp4.M();
	
	ReconProton = missingp4;

	if (MissingM > 1700.0 && MissingM <2100.0){    


	//Angular Distribution of particles		    

	phiPi= (PionCan.Phi()* TMath::RadToDeg()); 
	phiProton = (ProtonCan.Phi()* TMath::RadToDeg());
	phiReconProton = (ReconProton.Phi()* TMath::RadToDeg());
	coplanarityMeasure  = TMath::Abs((ProtonCan.TLorentzVector::DeltaPhi(PionCan))*TMath::RadToDeg()); //Changed to use abs, CHange to pid v pion
	coplanarityRecon  = TMath::Abs((ReconProton.TLorentzVector::DeltaPhi(PionCan))*TMath::RadToDeg()); //Changed to use abs
		  	
	//Method1
	thetaPi = (PionCan.Theta()* TMath::RadToDeg()); 
	thetaProton = (ProtonCan.Theta()* TMath::RadToDeg());
	thetaReconProton = (ReconProton.Theta()* TMath::RadToDeg());
	thetadiff  = ProtonCan.TLorentzVector::Angle(ReconProton.Vect())*TMath::RadToDeg();


  
	ProtonEnergy = ProtonCan.E();
//	Mygraph->SetPoint(GetEventNumber()coplanarityMeasure,ProtonEnergy,);


	//Boost Frame
	ProtonBoost = CMVector(ProtonCan,target,beam);
	PionBoost = CMVector(PionCan,target,beam);
	ReconProtonBoost = CMVector(ReconProton,target,beam);

	if (ProtonCan.X() ==0 && ProtonCan.Y() ==0 && ProtonCan.Z() ==0){

	NewEnergy = ProtonCan;
	detectornum = GetTracks()->GetDetectors(0);
//	std::cout << detectornum <<std::endl;
}

	else {
	NewEnergy.SetPxPyPzE(-10,-10,-10,-10);
}


	tree3g->Fill();
	treescatter->Fill();
	treeselected->Fill();
		//}
		//}
}//closing missing mass cut
} //Closing invariant mass cut

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

//    tree3g->Write();
//    tree3g->Reset();
//    return 0;

    // Write all GH1's and TObjects defined in this class
    return GTreeManager::Write();
}


