#include "chrisPPi0Example.h"

chrisPPi0Example::chrisPPi0Example()
{ 
 
//CAM 17/5/16 Adding wire chamber hits
    nChamberHitsin1 = new GH1("nChamberHitsin1", 	"nChamberHitsin1",	10,	0, 9);
    nChamberHitsin2 = new GH1("nChamberHitsin2", 	"Chamber2",	18,	0, 17);

   // Create Tree
  treeselected =  new TTree("Selected","Event selection tree"); // selected tree   
  treepid = new TTree("Pid","Event selection tree"); //pid tree
 
  //Define branches of TTree selected
  treeselected->Branch("Pion.",&PionCan); //Bad						****
  treeselected->Branch("Chamber1.",&Chamber1_Vec); //Bad				****
  treeselected->Branch("Chamber2.",&Chamber2_Vec); //Bad				****
//  treeselected->Branch("Phidiff.",&Phidiff); //Checked      //Wire chamber phi minus PID phi
//  treeselected->Branch("Ebeam.",&energy_beam); //Checked
//  treeselected->Branch("ESum.",&energySum);  //Checked
//  treeselected->Branch("Multiplicity.",&multiplicity); //Checked
//  treeselected->Branch("Inv_Mass_Pion.",&inv_M_value); //Checked
//  treeselected->Branch("PidIndex.",&PidHitIndex); //Checked
//  treeselected->Branch("Eventno.",&EventNumber); //Checked
//  treeselected->Branch("PIDPhi.",&PIDPhi); //Checked
//  treeselected->Branch("Chamber1Phi.",&Chamber1_VecPhi);
  treeselected->Branch("ProtonCandidate.",&ProtonCan); //Bad				****
  treeselected->Branch("MissingMass.",&MissingM); //Checked
//  treeselected->Branch("BeamHelicity.",&BeamHelicity); //Checked
//  treeselected->Branch("TaggedTime.",&time_beam); //Checked
  treeselected->Branch("TestChamber.",&TestChamber); //Bad 				****
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
//std::cout << GetOutputFile(1) <<std::endl;
    //Fill Beam Helicity
    Bool_t Helicity = GetTrigger()->GetHelicity();

    if(Helicity){

    //std::cout << Helicity <<  std::endl;
    BeamHelicity = 1;
}


if (GetMWPCHitsChris()->GetNMWPCHitsChrisChamber1() != 0){
    Chamber1_Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber1X(iii),GetMWPCHitsChris()->GetMWPCChamber1Y(iii),GetMWPCHitsChris()->GetMWPCChamber1Z(iii));

	TestChamber.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber1X(1),GetMWPCHitsChris()->GetMWPCChamber1Y(1),GetMWPCHitsChris()->GetMWPCChamber1Z(1));

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
treeselected->Write();
treeselected->Reset();
//    tree3g->Write();
//    tree3g->Reset();
    return 0;

    // Write all GH1's and TObjects defined in this class
//    return GTreeManager::Write();
}


