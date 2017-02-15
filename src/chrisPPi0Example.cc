#include "chrisPPi0Example.h"

chrisPPi0Example::chrisPPi0Example()
{ 

   // Create Tree,
  treeselected =  new TTree("Selected","Event selection tree"); // selected tree   
  
  //Define branches of TTree selected
  treeselected->Branch("Pion.",&PionCan);
  treeselected->Branch("Chamber1.",&Chamber1_Vec);
  treeselected->Branch("Chamber2.",&Chamber2_Vec);
  treeselected->Branch("Phidiff.",&Phidiff);       //Wire chamber phi minus PID phi
  treeselected->Branch("Ebeam.",&energy_beam);
  treeselected->Branch("ESum.",&energySum); // To get use 	energySum =GetTrigger()->GetEnergySum();
  treeselected->Branch("Multiplicity.",&multiplicity); // To get use 	multiplicity = GetTrigger()->GetMultiplicity();
  treeselected->Branch("Inv_Mass_Pion.",&inv_M_value);
  treeselected->Branch("PidIndex.",&PidHitIndex);
 // removed 17/10/16 treeselected->Branch("Eventno.",&EventNumber);
  treeselected->Branch("PIDPhi.",&PIDPhi);
  treeselected->Branch("Chamber1Phi.",&Chamber1_VecPhi);
  treeselected->Branch("ProtonCandidate.",&ProtonCan);
  treeselected->Branch("MissingMass.",&MissingM);
  treeselected->Branch("BeamHelicity.",&BeamHelicity);
  treeselected->Branch("TaggedTime.",&time_beam);
// removed 17/10/16  treeselected->Branch("TestChamber.",&TestChamber); //Not needed 
  treeselected->Branch("Cluster1Pi.",&cluster1_pi);
  treeselected->Branch("Cluster2Pi.",&cluster2_pi);
  treeselected->Branch("ClusterPiAngDiff.",&cluster_pi_ang_diff);
}

chrisPPi0Example::~chrisPPi0Example()
{
}

Bool_t	chrisPPi0Example::Init()
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;



	std::string config = ReadConfig("Period-Macro");
	if( sscanf(config.c_str(),"%d\n", &period) == 1 ) usePeriodMacro = 1;

	target.SetXYZM(0.0,0.0,0.0,1875.613);	//NEEDS CHANGING currently deuteron
//	target.SetXYZM(0.0,0.0,0.0,938.272);	//NEEDS CHANGING only TEMP for a SIM

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

  //Event Selection Parameters! Compare the reconstruction of pions to these to determine if a pion0
  IMPi0Criteria = 134.977; //pion mass
  MMPi0Criteria =target.M(); //dependent on the target used!


  targetPosition.SetXYZ(target.X(),target.Y(),target.Z());   



  Int_t iii = 0; // for the wire chamber info

  //Fill Beam Helicity
  Bool_t Helicity = GetTrigger()->GetHelicity();

  if(Helicity == 1){
    BeamHelicity = 1;
  }
  else{
    BeamHelicity = 0;
  }

  if (GetMWPCHitsChris()->GetNMWPCHitsChrisChamber1() != 0){
    Chamber1_Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber1X(iii),GetMWPCHitsChris()->GetMWPCChamber1Y(iii),GetMWPCHitsChris()->GetMWPCChamber1Z(iii));
  }
  else{
    Chamber1_Vec.SetXYZ(-1000,-1000,-1000);
  }

  if (GetMWPCHitsChris()->GetNMWPCHitsChrisChamber2() != 0){
    Chamber2_Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber2X(iii),GetMWPCHitsChris()->GetMWPCChamber2Y(iii),GetMWPCHitsChris()->GetMWPCChamber2Z(iii));
  }
  else{
    Chamber2_Vec.SetXYZ(-1000,-1000,-1000);
  }


  if  (Chamber1_Vec.X()==-1000){ //Can I get rid of all events if theres no hit in wc1 cos I wont be able to form the vectors properly.
    PVector1.SetXYZ(-1000,-1000,-1000);
    PVector2.SetXYZ(-1000,-1000,-1000);
  }
  else{
    PVector1 = Chamber1_Vec - targetPosition;
    PVector2 = Chamber2_Vec - PVector1;
  }

  //PID info(This section needs a serious revision over christmas)
  NPidhits = GetDetectorHits()->GetNPIDHits();

  if (NPidhits>0){
    PidHitIndex = GetDetectorHits()->GetPIDHits(0);    //The parameter in the Npidhits is the number of pid hits in the event while getPIDHits is the element number of a hit
    PIDPhi = PIDElemPhi[PidHitIndex]; //here

    if (Chamber1_Vec.X()==-1000){
      Chamber1_VecPhi =-10000 ;
    }

    else{
      Chamber1_VecPhi =Chamber1_Vec.Phi()*TMath::RadToDeg() ;
    }
//Need to justify this (why not simple difference!!)
    PIDPhi = PIDPhi + 180;
    Chamber1_VecPhi =Chamber1_VecPhi +180;
    Phidiff = TMath::Abs( PIDPhi - Chamber1_VecPhi) ;   

    if (Phidiff >180){
      Phidiff = 360 -Phidiff;
    } //closing phidiff

  } //closingpihits

  else{  //Changing from -10 to -1000
    Phidiff = -10;
    PIDPhi = -10;
    Chamber1_VecPhi = -10;
  }


  if (Phidiff< -100){
     Phidiff = -100;
  }


  // Loop over tagged events     Needs change from Roddy, Not = to
  for(Int_t i=0; i<=GetTagger()->GetNTagged() ;i++){
    
    if ( (GetPhotons()->GetNParticles()==3)){

      multiplicity = GetTrigger()->GetMultiplicity();
      energySum =GetTrigger()->GetEnergySum();
      energy_beam = GetTagger()->GetTaggedEnergy(i);
      beam.SetXYZM(0.,0.,energy_beam,0.);  
      time_beam = GetTagger()->GetTaggedTime(i);
      // Change the selection into a function
      // Loop over the different combination of photons and set the different vectors for each
      for(Int_t l=0;l<3;l++){
	if(l==0) {
	  l_g1 = 1;
	  l_g2 = 2;
	  g1_Vec1 = GetPhotons()->Particle(l_g1);
	  g2_Vec1 = GetPhotons()->Particle(l_g2);
	  C1 = g1_Vec1 + g2_Vec1;  //pion candidate
	  IMDiff1 = TMath::Abs(C1.M() - IMPi0Criteria)/IMPi0Criteria;
	  P1 = GetPhotons()->Particle(0); // proton candidate
	  Mp41 = beam + target - C1;
	  MM1 = Mp41.M();
	  MMDiff1 = TMath::Abs(MM1 - MMPi0Criteria)/MMPi0Criteria; //What should the missing mass value be?
	  SumDiff1 = (IMDiff1) + (MMDiff1);
	}

	if(l==1) {
	  l_g1 = 2;
	  l_g2 = 0;
	  g1_Vec2 = GetPhotons()->Particle(l_g1);
	  g2_Vec2 = GetPhotons()->Particle(l_g2);
	  C2 = g1_Vec2 + g2_Vec2;
	  IMDiff2 = TMath::Abs(C2.M() - IMPi0Criteria)/IMPi0Criteria;
	  P2 = GetPhotons()->Particle(1);
	  Mp42 = beam + target - C2;
	  MM2 = Mp42.M();
	  MMDiff2 = TMath::Abs(MM2 - MMPi0Criteria)/MMPi0Criteria; //What should the missing mass value be? should i divide the difference by the peak position to get them both scaled percentage wise.
	  SumDiff2 = (IMDiff2) + (MMDiff2) ;
	}

	if(l==2) { 
	  l_g1 = 0;
	  l_g2 = 1;
	  g1_Vec3 = GetPhotons()->Particle(l_g1);
	  g2_Vec3 = GetPhotons()->Particle(l_g2);
	  C3 = g1_Vec3 + g2_Vec3;
	  IMDiff3 = TMath::Abs(C3.M() - IMPi0Criteria)/IMPi0Criteria;
	  P3 = GetPhotons()->Particle(2);
	  Mp43 = beam + target - C3;
	  MM3 = Mp43.M();
	  MMDiff3 = TMath::Abs(MM3 - MMPi0Criteria)/MMPi0Criteria; //What should the missing mass value be?
	  SumDiff3 = (IMDiff3) + (MMDiff3) ;
	}
	

      } //closing for candidate loop
      
      if( SumDiff1 < SumDiff2 && SumDiff1 < SumDiff3){
	PionCan.SetXYZM(C1.X(),C1.Y(),C1.Z(),C1.M());
	ProtonCan = P1;
	cluster1_pi=GetPhotons()->Particle(1);
	cluster2_pi=GetPhotons()->Particle(2); 
	cluster_pi_ang_diff = cluster1_pi.TLorentzVector::Angle(cluster2_pi.Vect())*TMath::RadToDeg()/180 ;
     }
      	
      else if(SumDiff2 < SumDiff1 && SumDiff2 < SumDiff3){
	PionCan.SetXYZM(C2.X(),C2.Y(),C2.Z(),C2.M());
	ProtonCan = P2;
	cluster1_pi=GetPhotons()->Particle(2);
        cluster2_pi=GetPhotons()->Particle(0);   
	cluster_pi_ang_diff = cluster1_pi.TLorentzVector::Angle(cluster2_pi.Vect())*TMath::RadToDeg()/180 ;
      }
	  
      else if(SumDiff3 < SumDiff1 && SumDiff3 < SumDiff2){
	PionCan.SetXYZM(C3.X(),C3.Y(),C3.Z(),C3.M());
	ProtonCan = P3;
        cluster1_pi=GetPhotons()->Particle(0);
        cluster2_pi=GetPhotons()->Particle(1);   
	cluster_pi_ang_diff = cluster1_pi.TLorentzVector::Angle(cluster2_pi.Vect())*TMath::RadToDeg()/180 ;

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

    } //closing if 3 Particles
	

  } //closing for loop

  nEventsWritten++;


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
treeselected->Write(); //Added CAM 28/09/16
treeselected->Reset(); //Added CAM 28/09/16
//    tree3g->Write();
//    tree3g->Reset();
    return 0; //Added CAM 28/09/16

    // Write all GH1's and TObjects defined in this class
//    return GTreeManager::Write(); //Removed CAM 28/09/16
}


