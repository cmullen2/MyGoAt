#include "chrisNPi0.h"
#include "TROOT.h"

chrisNPi0::chrisNPi0()
{ 
  gROOT->ProcessLine("#include <vector>");
  // gROOT->ProcessLine(".x /home/chris/GoAT/TestVector.C++");
  //gROOT->LoadMacro("$HSANA/THSParticle.C+");
  // Create Tree,
  treeselected =  new TTree("Selected","Event selection tree"); // selected tree   
  
  //Define branches of TTree selected

  //Splots branches
  treeselected->Branch("Inv_Mass_Pion.",&inv_M_value);
  treeselected->Branch("MissingMass.",&MissingM);
  treeselected->Branch("TaggedTime.",&time_beam);
  treeselected->Branch("Coplanarity.",&Coplanarity);
 // treeselected->Branch("OpeningAngle.",&OpeningAngle);
  treeselected->Branch("EventNo.",&EventNo);

  //Non-Splots branches
  treeselected->Branch("Pion.",&PionCan);
  treeselected->Branch("Ebeam.",&energy_beam);
  treeselected->Branch("ESum.",&energySum); // To get use 	energySum =GetTrigger()->GetEnergySum();
  treeselected->Branch("Multiplicity.",&multiplicity); // To get use 	multiplicity = GetTrigger()->GetMultiplicity();
  treeselected->Branch("NeutronCandidate.",&NeutronCan);
  treeselected->Branch("NeutronMass.",&NeutronMass);
  treeselected->Branch("BeamHelicity.",&BeamHelicity);
  treeselected->Branch("Cluster1Pi.",&cluster1_pi);
  treeselected->Branch("Cluster2Pi.",&cluster2_pi);
  treeselected->Branch("ClusterPiAngDiff.",&cluster_pi_ang_diff);
  treeselected->Branch("fparticlelist","vector<THSParticle*>",&fparticleList);
  treeselected->Branch("particle1",&particle1);
  DummyProton=0;
  DummyProton=new THSParticle();
}

chrisNPi0::~chrisNPi0()
{
  delete DummyProton;
}

Bool_t	chrisNPi0::Init()
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;



	std::string config = ReadConfig("Period-Macro");
	if( sscanf(config.c_str(),"%d\n", &period) == 1 ) usePeriodMacro = 1;

	target.SetXYZM(0.0,0.0,-65.0,1875.613);	//NEEDS CHANGING currently deuteron
//	target.SetXYZM(0.0,0.0,-65.0,938.272);	//NEEDS CHANGING only TEMP for a SIM

	if(!InitBackgroundCuts()) return kFALSE;
	if(!InitTargetMass()) return kFALSE;
	if(!InitTaggerChannelCuts()) return kFALSE;
	if(!InitTaggerScalers()) return kFALSE;
	cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	chrisNPi0::Start()
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

void	chrisNPi0::ProcessEvent()
{

  if(usePeriodMacro == 1)
    {
      if(GetEventNumber() % period == 0)
	cout << "Events: " << GetEventNumber() << "  Events Accepted: " << nEventsWritten << endl;
    }


  std::string outFile = outputFile->GetName();

  std::size_t pos = outFile.find("_1");
  std::string filert = outFile.substr(pos);
  std::string fileNo =  filert.substr(2,4);
  Int_t tempEventno = GetEventNumber();
  std::string tevent = std::to_string(tempEventno);
  std::string eventName = fileNo + tevent;
  EventNo = std::stod(eventName);
  
 

  //Event Selection Parameters! Compare the reconstruction of pions to these to determine if a pion0
  IMPi0Criteria = 134.977; //pion mass
  MMPi0Criteria =target.M(); //dependent on the target used!
  MMlowercut = 1700.0;
  MMuppercut = 2100.0;
  //  MMlowercut = 800.0;
  //  MMuppercut = 1200.0;
  invlowercut = 100.0;
  invuppercut = 170.0;

  targetPosition.SetXYZ(target.X(),target.Y(),target.Z());   


  //Fill Beam Helicity
  Bool_t Helicity = GetTrigger()->GetHelicity();

  if(Helicity == 1){
    BeamHelicity = 1;
  }
  else{
    BeamHelicity = 0;
  }

  //Attempts to use THSParticle class to store information:
  // Single THSParticle storage works
  // Moving on to storing them in vectors using pointers

  TVector3 Vertex;
  Vertex.SetXYZ(1,2,3);
  TLorentzVector Party1;
  Party1.SetXYZM(4,5,6,7);
  particle1.SetP4(Party1);
  particle1.SetVertex(Vertex);
  particle1.SetTime(8);
  particle1.SetMeasMass(9);
  particle1.SetPDGcode(11);
      
  //std::vector<THSParticle> particleList;
  fparticleList.clear();
  //  TVector3 Vertex;
  //  Vertex.SetXYZ(1,2,3);
  //  TLorentzVector Party1;
  //  Party1.SetXYZM(4,5,6,7);
  //THSParticle DummyProton; // if doing pointer then these should be members of class too surely?
  DummyProton->SetP4(Party1);
  DummyProton->SetVertex(Vertex);
  DummyProton->SetTime(8);
  DummyProton->SetMeasMass(9);
  DummyProton->SetPDGcode(11);
  fparticleList.push_back(DummyProton);
  //  fparticleList.push_back(DummyProton); 
  //    fparticleList.push_back(particle1); 
   



  if ( (GetPhotons()->GetNParticles()==3)){


    //PID info(This section needs a serious revision over christmas)
    NPidhits = GetDetectorHits()->GetNPIDHits();
    if (NPidhits==0){

      //if (GetMWPCHitsChris()->GetNMWPCHitsChrisChamber1() == 0){

      // Loop over tagged events     Needs change from Roddy, Not = to
      for(Int_t i=0; i<GetTagger()->GetNTagged() ;i++){
    
	time_beam = GetTagger()->GetTaggedTime(i);

	if ( time_beam > 40 && time_beam <80 ) {

	  multiplicity = GetTrigger()->GetMultiplicity();
	  energySum =GetTrigger()->GetEnergySum();
	  energy_beam = GetTagger()->GetTaggedEnergy(i);
	  beam.SetXYZM(0.,0.,energy_beam,0.);  
      
	  Functionx3
	    if (SumDiff1 < SumDiff2 && SumDiff1 < SumDiff3 ){
	      AssignFunction(SumDiff1)

		}//Closing SumDiff1 

	  if (SumDiff2 < SumDiff1 && SumDiff2 < SumDiff3 ){
	    AssignFunction(SumDiff2)

	      }//Closing SumDiff1 


	  if (SumDiff3 < SumDiff2 && SumDiff3 < SumDiff1 ){
	    AssignFunction(SumDiff3)

	      }//Closing SumDiff3 
      
	  /*	
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
	  NeutronCan = P1;
	  cluster1_pi=GetPhotons()->Particle(1);
	  cluster2_pi=GetPhotons()->Particle(2); 
	  cluster_pi_ang_diff = cluster1_pi.TLorentzVector::Angle(cluster2_pi.Vect())*TMath::RadToDeg()/180 ;
	  }
      	
	  else if(SumDiff2 < SumDiff1 && SumDiff2 < SumDiff3){
	  PionCan.SetXYZM(C2.X(),C2.Y(),C2.Z(),C2.M());
	  NeutronCan = P2;
	  cluster1_pi=GetPhotons()->Particle(2);
	  cluster2_pi=GetPhotons()->Particle(0);   
	  cluster_pi_ang_diff = cluster1_pi.TLorentzVector::Angle(cluster2_pi.Vect())*TMath::RadToDeg()/180 ;
	  }
	  
	  else if(SumDiff3 < SumDiff1 && SumDiff3 < SumDiff2){
	  PionCan.SetXYZM(C3.X(),C3.Y(),C3.Z(),C3.M());
	  NeutronCan = P3;
	  cluster1_pi=GetPhotons()->Particle(0);
	  cluster2_pi=GetPhotons()->Particle(1);   
	  cluster_pi_ang_diff = cluster1_pi.TLorentzVector::Angle(cluster2_pi.Vect())*TMath::RadToDeg()/180 ;

	  }  
	  */    
	  // Invariant Mass
	  inv_M_value = PionCan.M();
	  NeutronMass = NeutronCan.M();
  
	  //Invariant mass cut here! that can be commented out to get both plots.
	  //CAM 16/03/17     if (inv_M_value > invlowercut && inv_M_value <invuppercut){   
		    
	  missingp4     = beam + target - PionCan;
	  MissingM = missingp4.M();

	  //CAM 16/03/17	if (MissingM > MMlowercut && MissingM <MMuppercut){    

	  //Angular Distribution of particles		    
	
	  Coplanarity = (PionCan).TLorentzVector::DeltaPhi(NeutronCan)*TMath::RadToDeg();
	  if (Coplanarity <0){
	    Coplanarity =Coplanarity +360;

	  }

	  treeselected->Fill();
	  //}
	  //}
	  //CAM 16/03/17	}//closing missing mass cut

	  //CAM 16/03/17      } //Closing invariant mass cut

	} //closing if TaggedTime
	

      } //closing for loop

      //} //closing NMWPCHits if

    } //closing NPIDHits if


    //Do Proton channel here 

    //***********************************************************************************************


    if (NPidhits > 0){
	
      for(Int_t i=0; i<GetTagger()->GetNTagged() ;i++){

	time_beam = GetTagger()->GetTaggedTime(i);

	if ( time_beam > 40 && time_beam <80 ) {

	  multiplicity = GetTrigger()->GetMultiplicity();
	  energySum =GetTrigger()->GetEnergySum();
	  energy_beam = GetTagger()->GetTaggedEnergy(i);
	  beam.SetXYZM(0.,0.,energy_beam,0.);  
      
	  Functionx3
	    if (SumDiff1 < SumDiff2 && SumDiff1 < SumDiff3 ){
	      AssignFunction(SumDiff1)

		}//Closing SumDiff1 

	  if (SumDiff2 < SumDiff1 && SumDiff2 < SumDiff3 ){
	    AssignFunction(SumDiff2)

	      }//Closing SumDiff1 


	  if (SumDiff3 < SumDiff2 && SumDiff3 < SumDiff1 ){
	    AssignFunction(SumDiff3)

	      }//Closing SumDiff3 


	  // Invariant Mass
	  inv_M_value = PionCan.M();
	  NeutronMass = NeutronCan.M();
  
		    
	  missingp4     = beam + target - PionCan;
	  MissingM = missingp4.M();
	  //Angular Distribution of particles		    
	
	  Coplanarity = (PionCan).TLorentzVector::DeltaPhi(NeutronCan)*TMath::RadToDeg();
	  if (Coplanarity <0){
	    Coplanarity =Coplanarity +360;

	  }

	  treeselected2->Fill();




	}//closing p pi0 Taggedtime

      } //clossing p pi0 for NTagged

    } //Closing p pi0 NPIDHITS



  } //closing if 3 particles

  nEventsWritten++;


} //closing function


void	chrisNPi0::ProcessScalerRead()
{
	// Fill Tagger Scalers
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	chrisNPi0::Write()
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


