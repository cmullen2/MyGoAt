#include "chrisForTHSProject.h"
#include "TROOT.h"

chrisForTHSProject::chrisForTHSProject()
{ 
  gROOT->ProcessLine("#include <vector>");

  treePi0 =  new TTree("Particles","Event selection tree"); // Neutron pi0 final state tree   
  //treeProtonPi0 =  new TTree("ProtonPi0","Event selection tree"); // Proton pi0 final state tree   


  //Neutron Tree branches
  //  treeNeutronPi0->Branch("fparticleList","vector<THSParticle*>",&fparticleList);
  //  treeNeutronPi0->Branch("ftaggedPhotons","vector<THSParticle*>",&ftaggedPhotons);
  treeNeutronPi0->Branch("fparticleList","vector<THSParticle>",&fparticleList);
//  treeNeutronPi0->Branch("ftaggedPhotons","vector<THSParticle>",&ftaggedPhotons);
  treeNeutronPi0->Branch("fbeamHelicity",&fbeamHelicity);
  treeNeutronPi0->Branch("ftaggedTime",&ftaggedTime);
  treeNeutronPi0->Branch("feventNo",&feventNo);
  treeNeutronPi0->Branch("fenergyBeam",&fenergyBeam);
  treeNeutronPi0->Branch("fenergySum",&fenergySum);
  treeNeutronPi0->Branch("fmultiplicity",&fmultiplicity);
  treeNeutronPi0->Branch("fpidPhi",&fpidPhi);
  treeNeutronPi0->Branch("fpidIndex",&fpidIndex);
  treeNeutronPi0->Branch("fchamber1Vec",&fchamber1Vec);
  treeNeutronPi0->Branch("fchamber2Vec",&fchamber2Vec);


  //Proton Tree branches
//  treeProtonPi0->Branch("fparticlelist","vector<THSParticle>",&fparticleList);
//  treeProtonPi0->Branch("fparticlelist","vector<THSParticle*>",&fparticleList);
//  treeProtonPi0->Branch("ftaggedPhotons","vector<THSParticle*>",&ftaggedPhotons);
  treeProtonPi0->Branch("fbeamHelicity",&fbeamHelicity);
  treeProtonPi0->Branch("ftaggedTime",&ftaggedTime);
  treeProtonPi0->Branch("feventNo",&feventNo);
  treeProtonPi0->Branch("fenergyBeam",&fenergyBeam);
  treeProtonPi0->Branch("fenergySum",&fenergySum);
  treeProtonPi0->Branch("fmultiplicity",&fmultiplicity);
  treeProtonPi0->Branch("fpidPhi",&fpidPhi);
  treeProtonPi0->Branch("fpidIndex",&fpidIndex);
  treeProtonPi0->Branch("fchamber1Vec",&fchamber1Vec);
  treeProtonPi0->Branch("fchamber2Vec",&fchamber2Vec);


  
  // gROOT->ProcessLine(".x /home/chris/GoAT/TestVector.C++");
  //gROOT->LoadMacro("$HSANA/THSParticle.C+");
  // Create Tree,
 //  treeselected =  new TTree("Selected","Event selection tree"); // selected tree   
  
//   //Define branches of TTree selected

//   //Splots branches
// //  treeselected->Branch("Inv_Mass_Pion.",&inv_M_value);
// //  treeselected->Branch("MissingMass.",&MissingM);
//   treeselected->Branch("TaggedTime.",&time_beam);
//   treeselected->Branch("Coplanarity.",&Coplanarity);
//  // treeselected->Branch("OpeningAngle.",&OpeningAngle);
//   treeselected->Branch("EventNo.",&EventNo);

//   //Non-Splots branches
//   treeselected->Branch("Pion.",&PionCan);
//   treeselected->Branch("Ebeam.",&energy_beam);
//   treeselected->Branch("ESum.",&energySum); // To get use 	energySum =GetTrigger()->GetEnergySum();
//   treeselected->Branch("Multiplicity.",&multiplicity); // To get use 	multiplicity = GetTrigger()->GetMultiplicity();
//   treeselected->Branch("NeutronCandidate.",&NeutronCan);
//   treeselected->Branch("NeutronMass.",&NeutronMass);
//   treeselected->Branch("BeamHelicity.",&BeamHelicity);
//   treeselected->Branch("Cluster1Pi.",&cluster1_pi);
//   treeselected->Branch("Cluster2Pi.",&cluster2_pi);
//   treeselected->Branch("ClusterPiAngDiff.",&cluster_pi_ang_diff);
//   treeselected->Branch("fparticlelist","vector<THSParticle*>",&fparticleList);
//   treeselected->Branch("particle1",&particle1);

// Setting up the THSParticles

  // fballPhotons=0;
  // fballPhotons=new THSParticle();

  // ftaggPhotons=0;
  // ftaggPhotons=new THSParticle();

  // DummyProton=0;
  // DummyProton=new THSParticle();
}

chrisForTHSProject::~chrisForTHSProject()
{
//   delete DummyProton;
//   delete fballPhotons;
//   delete ftaggPhotons;
 }

Bool_t	chrisForTHSProject::Init()
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

Bool_t	chrisForTHSProject::Start()
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

void	chrisForTHSProject::ProcessEvent()
{

  if(usePeriodMacro == 1)
    {
      if(GetEventNumber() % period == 0)
	cout << "Events: " << GetEventNumber() << "  Events Accepted: " << nEventsWritten << endl;
    }

	// Unique Event ID
  std::string outFile = outputFile->GetName();
  std::size_t pos = outFile.find("_1");
  std::string filert = outFile.substr(pos);
  std::string fileNo =  filert.substr(2,4);
  Int_t tempEventno = GetEventNumber();
  std::string tevent = std::to_string(tempEventno);
  std::string eventName = fileNo + tevent;
  feventNo = std::stod(eventName);

  //Clearing the vectors.
  fparticleList.clear();
  ftaggedPhotons.clear();

  

  //MWPC information
  targetPosition.SetXYZ(target.X(),target.Y(),target.Z());

  Int_t iii = 0;

  if (GetMWPCHitsChris()->GetNMWPCHitsChrisChamber1() != 0){
    fchamber1Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber1X(iii),GetMWPCHitsChris()->GetMWPCChamber1Y(iii),GetMWPCHitsChris()->GetMWPCChamber1Z(iii));
  }
  else{
    fchamber1Vec.SetXYZ(-1000,-1000,-1000);
  }

  if (GetMWPCHitsChris()->GetNMWPCHitsChrisChamber2() != 0){
    fchamber2Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber2X(iii),GetMWPCHitsChris()->GetMWPCChamber2Y(iii),GetMWPCHitsChris()->GetMWPCChamber2Z(iii));
  }
  else{
    fchamber2Vec.SetXYZ(-1000,-1000,-1000);
  }


  //PID information
  NPidhits = GetDetectorHits()->GetNPIDHits();
  
  if (NPidhits>0){
    fpidIndex = GetDetectorHits()->GetPIDHits(0);    //The parameter in the Npidhits is the number of pid hits in the event while getPIDHits is the element number of a hit
    fpidPhi = PIDElemPhi[fpidIndex]; //here

  } //closingpihits

  else{  
    
    fpidPhi = -10;
  }


  
  //Fill Beam Helicity
  Bool_t Helicity = GetTrigger()->GetHelicity();
  if(Helicity == 1){
    fbeamHelicity = 1;
  }
  else{
    fbeamHelicity = 0;
  }


  
//************************************************************************************************************************

  if ( (GetPhotons()->GetNParticles()==3)){

    //PID info
    if (NPidhits==0){
	
      for(Int_t j=0; j<GetPhotons()->GetNParticles(); j++){ 

		
	//pushback each of the three photons here to THSParticle.
	//Need a single THSParticle defined as dummyproton above.
	//This should take the p4,vertex,time,measmass,pdgcode(see THSParticle.h) of each j photon
	//then pushback to the particle list

	fcbPhoton = GetPhotons()->Particle(j); //TLorentzVector
	// fballPhotons->SetP4(fcbPhoton); //THSParticle
	// fballPhotons->SetPDGcode(22);
	fballPhotons.SetP4(fcbPhoton); //THSParticle
	fballPhotons.SetPDGcode(22);
	fparticleList.push_back(fballPhotons);

      } //Closing For NParticles 


      for(Int_t i=0; i<GetTagger()->GetNTagged() ;i++){
 	
	ftaggedTime = GetTagger()->GetTaggedTime(i);

	if ( ftaggedTime > 40 && ftaggedTime <80 ) {
	  
	  fmultiplicity = GetTrigger()->GetMultiplicity();
	  fenergySum =GetTrigger()->GetEnergySum();
	  fenergyBeam = GetTagger()->GetTaggedEnergy(i);
	  beam.SetXYZM(0.,0.,fenergyBeam,0.);  
	  
	  //pushback the tagged photons to THSParticle here using the pdg code 0
	  fglasgowTaggerPhoton.SetPxPyPzE(0,0,fenergyBeam, fenergyBeam);  //TLorentzVector
	  // ftaggPhotons->SetP4(fglasgowTaggerPhoton); //THSParticle
	  // ftaggPhotons->SetPDGcode(0);
	  // ftaggPhotons->SetTime(ftaggedTime);
	  ftaggPhotons.SetP4(fglasgowTaggerPhoton); //THSParticle
	  ftaggPhotons.SetPDGcode(0);
	  ftaggPhotons.SetTime(ftaggedTime);
	  // ftaggedPhotons.push_back(ftaggPhotons);
	  fparticleList.push_back(ftaggPhotons);
	  
	  
	} //closing if TaggedTime
	
      } //Closing for NTagged.

      
      treeNeutronPi0->Fill();
    } //closing NPIDHits if
    
    
    
    //Do Proton channel here 

    //********************************************************************************************************************************
    
    //PID info
    NPidhits = GetDetectorHits()->GetNPIDHits();
    if (NPidhits>0){
      
      //PID determination of proton.
      if (fchamber1Vec.X()==-1000){
	fchamber1VecPhi =-10000 ;
      }
      else{
	fchamber1VecPhi =fchamber1Vec.Phi()*TMath::RadToDeg() ;
      }
      
      fpidPhi = fpidPhi + 180;
      fchamber1VecPhi =fchamber1VecPhi +180;
      fphidiff = TMath::Abs( fpidPhi - fchamber1VecPhi) ;

      if (fphidiff >180){
	fphidiff = 360 -fphidiff;
      } //closing phidiff
    
    else{  //Changing from -10 to -1000
      fphidiff = -10;
      fpidPhi = -10;
      fchamber1VecPhi = -10;
    }

    
    if (fphidiff< -100){
      fphidiff = -100;
    }
    
    
    

    
    


	
    for(Int_t k=0; k<GetPhotons()->GetNParticles(); k++){ 
      
      
      //pushback each of the three photons here to THSParticle.
      if (fphidiff <15 && fphidiff >-1){
	
    	fcbPhoton = GetPhotons()->Particle(k); //TLorentzVector
  	fballPhotons.SetP4(fcbPhoton); //THSParticle
    	fballPhotons.SetPDGcode(2212);
    	fparticleList.push_back(fballPhotons);
	
      }
      else{
    	fcbPhoton = GetPhotons()->Particle(k); //TLorentzVector
    	// fballPhotons->SetP4(fcbPhoton); //THSParticle
    	// fballPhotons->SetPDGcode(22);
	fballPhotons.SetP4(fcbPhoton); //THSParticle
    	fballPhotons.SetPDGcode(22);
    	fparticleList.push_back(fballPhotons);
      }
      
    } //Closing For NParticles 
    
    
    
    for(Int_t l=0; l<GetTagger()->GetNTagged() ;l++){
      
      ftaggedTime = GetTagger()->GetTaggedTime(l);
      
    	if ( ftaggedTime > 40 && ftaggedTime <80 ) {
	  
    	  fmultiplicity = GetTrigger()->GetMultiplicity();
    	  fenergySum =GetTrigger()->GetEnergySum();
    	  fenergyBeam = GetTagger()->GetTaggedEnergy(l);
    	  beam.SetXYZM(0.,0.,fenergyBeam,0.);  
	  
    	  //pushback the tagged photons to THSParticle here using the pdg code 0
    	  fglasgowTaggerPhoton.SetPxPyPzE(0,0,fenergyBeam, fenergyBeam);  //TLorentzVector
    	  // ftaggPhotons->SetP4(fglasgowTaggerPhoton); //THSParticle
    	  // ftaggPhotons->SetPDGcode(0);
    	  // ftaggPhotons->SetTime(ftaggedTime);
	  ftaggPhotons.SetP4(fglasgowTaggerPhoton); //THSParticle
    	  ftaggPhotons.SetPDGcode(0);
    	  ftaggPhotons.SetTime(ftaggedTime);
    	  //ftaggedPhotons.push_back(ftaggPhotons);
    	  fparticleList.push_back(ftaggPhotons);
	  
	  
	  
	  
	  
    	} //closing if TaggedTime
	
    } //Closing for NTagged.
    
    
    treeProtonPi0->Fill();
  } //closing NPIDHits if
  
  
  //*************************************************************************************************************************************

  

} //closing if 3 particles

nEventsWritten++;

  
} //closing function


void	chrisForTHSProject::ProcessScalerRead()
{
	// Fill Tagger Scalers
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	chrisForTHSProject::Write()
{

  treeNeutronPi0->Write(); 
  treeProtonPi0->Write();
  treeNeutronPi0->Reset();
  treeProtonPi0->Reset();

    return 0; //Added CAM 28/09/16

    // Write all GH1's and TObjects defined in this class
//    return GTreeManager::Write(); //Removed CAM 28/09/16
}


