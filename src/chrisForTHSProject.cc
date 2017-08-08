#include "chrisForTHSProject.h"
#include "TROOT.h"


Double_t CircPol(Double_t , Double_t );


chrisForTHSProject::chrisForTHSProject()
{ 
  gROOT->ProcessLine("#include <vector>");

  treePi0 =  new TTree("HSParticles","Event selection tree"); // Neutron pi0 final state tree   
  //treeProtonPi0 =  new TTree("ProtonPi0","Event selection tree"); // Proton pi0 final state tree   


  //Neutron Tree branches
  // treePi0->Branch("Particles.","vector<THSParticle>",&Particles);
  treePi0->Branch("Particles",&Particles);
  treePi0->Branch("fbeamHelicity",&fbeamHelicity);
  treePi0->Branch("ftaggedTime",&ftaggedTime);
  treePi0->Branch("feventNo",&feventNo);
  treePi0->Branch("fenergyBeam",&fenergyBeam);
  treePi0->Branch("fenergySum",&fenergySum);
  treePi0->Branch("fmultiplicity",&fmultiplicity);
  treePi0->Branch("fpidPhi",&fpidPhi);
  treePi0->Branch("fpidIndex",&fpidIndex);
  //  treePi0->Branch("fpidRootinoPhi",&fpidRootinoPhi);
  treePi0->Branch("fchamber1Vec",&fchamber1Vec);
  treePi0->Branch("fchamber2Vec",&fchamber2Vec);
  treePi0->Branch("fedgeplane",&fedgePlane);
  treePi0->Branch("flinPol",&flinPol);
  treePi0->Branch("frootinoPhi",&frootinoPhi);
treePi0->Branch("rootinoClustE",&rootinoClustE);
treePi0->Branch("rootinoTheta",&rootinoTheta);
treePi0->Branch("rootinoPhi",&rootinoPhi);
treePi0->Branch("rootinoTime",&rootinoTime);
treePi0->Branch("rootinoClustS",&rootinoClustS);
treePi0->Branch("rootinoClustC",&rootinoClustC);
treePi0->Branch("rootinoVetoC",&rootinoVetoC);
treePi0->Branch("rootinoDet",&rootinoDet);
treePi0->Branch("rootinoVetoE",&rootinoVetoE);
treePi0->Branch("rootinoCham1E",&rootinoCham1E);
treePi0->Branch("rootinoCham2E",&rootinoCham2E);
 


}

chrisForTHSProject::~chrisForTHSProject()
{

}

Bool_t	chrisForTHSProject::Init()
{
  cout << "Initialising physics analysis..." << endl;
  cout << "--------------------------------------------------" << endl << endl;



  std::string config = ReadConfig("Period-Macro");
  if( sscanf(config.c_str(),"%d\n", &period) == 1 ) usePeriodMacro = 1;

  target.SetXYZM(0.0,0.0,-65.0,1875.613);	//NEEDS CHANGING currently deuteron
  //	target.SetXYZM(0.0,0.0,-65.0,938.272);	//NEEDS CHANGING only TEMP for a SIM

  //Following line may break due to GetNTagged not being defined yet, if so simply replace with a number and the while statement in main will correct it.
  fNin = (3+ (GetTagger()->GetNTagged()) ) ; //Estimate of Number of input particles(tagged+ball+rootinos etc.) Can and will calc this!(3+NTagged)
  fReadParticles=new vector<THSParticle*>;
  for(Int_t m=0; m<fNin; m++ ){
    fReadParticles->push_back(new THSParticle());
  }

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
	cout << "Events: " << GetEventNumber() << "  Events Accepted: " << nEventsWritten <<" testcounter= " <<testcounter <<  endl;
    }

Int_t mc =0;


if(!mc){

  // Unique Event ID
  std::string outFile = outputFile->GetName();
  std::size_t pos = outFile.find("_1");
  std::string filert = outFile.substr(pos);
  fileNo =  filert.substr(2,4);
  Int_t tempEventno = GetEventNumber();
  std::string tevent = std::to_string(tempEventno);
  std::string tevent2 = fileNo + "0000";
  Double_t eventName = std::stod(tevent2) + std::stod(tevent);
  feventNo = eventName;
//  feventNo = std::stod(eventName);
}
else{
//MC needs to deal with event no, ePol and linpol and where para or perp
cout << " Using MC data " << endl;

std::string outFile = outputFile->GetName();
std::string filert = outFile.substr( outFile.length() - 11  );
//cout <<" Substring is filert " << filert << endl;

std::string planesetting = filert.substr(0,3);
//cout << " plane setting " << planesetting <<endl;

if(planesetting=="Neg" ) fileNo = "11110000";
if(planesetting=="Pos" ) fileNo = "22220000";

Int_t tempEventno = GetEventNumber();
feventNo = std::stod(fileNo) + tempEventno;

}




  //Clearing the vectors.
  Particles.clear();

 //Mott measurements for electron beam polarisation
  Double_t runNo = std::stod(fileNo);
  if (runNo<=4626) ePol = MottMeas[0]; 
  if (runNo>4626 && runNo<=4725 ) ePol = MottMeas[1]; 
  if (runNo>4725 && runNo<=4835 ) ePol = MottMeas[2]; 
  if (runNo>4835 && runNo<=4961 ) ePol = MottMeas[3]; 
  if (runNo>4961 && runNo<=5076 ) ePol = MottMeas[4]; 
  if (runNo>5076 && runNo<=5274 ) ePol = MottMeas[5]; 
  if (runNo>5274 && runNo<=5505 ) ePol = MottMeas[6]; 
  if (runNo>5505 && runNo<=5767 ) ePol = MottMeas[7]; 
  if (runNo>5767 ) ePol = MottMeas[8]; 
  if (runNo==11110000 ) ePol = MCMottNeg; 
  if (runNo==22220000 ) ePol = MCMottPos; 

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
  
  // if (NPidhits>0){
  //   fpidIndex = GetDetectorHits()->GetPIDHits(0);    //The parameter in the Npidhits is the number of pid hits in the event while getPIDHits is the element number of a hit
  //   fpidPhi = PIDElemPhi[fpidIndex]; //here

  // } //closingpihits

  // else{  
    
  //   fpidPhi = -10;
  // }


  
  //Fill Beam Helicity (Two different determined states when helicity is 0 or 1, others are error codes) Chosen 1 and -1 for helicity from now on. 
  //Check for error codes here?? Need to decide if it is necessary? Should only be a small effect so come back to it once everything is working
  Bool_t Helicity = GetTrigger()->GetHelicity();
  if(Helicity == 1){
    fbeamHelicity = 1;
  }
  if(Helicity == 0){
    fbeamHelicity = -1;
  }

  //Linear Polarisation Information
  if(GetLinpol()->GetPolarizationPlane()==0) fedgePlane = -1; // Para
  if(GetLinpol()->GetPolarizationPlane()==1){ fedgePlane = 1;} // Perp
  else{
    fedgePlane=0; //Moeller or other
  }



  fNin= 3+ (GetTagger()->GetNTagged());
  Int_t countb4 = 0;
  Int_t countaft = 0;
  Int_t counter = 0;

  
  
  //************************************************************************************************************************
  //Methodology: Now using a rootino for proton so how  does this effect neutron will it be rootino or photon? photon right?
  if ( (GetPhotons()->GetNParticles()==3)){

    //PID info
    if (NPidhits==0){


      while(fNin>fReadParticles->size()){
	fReadParticles->push_back(new THSParticle());
      }

      
      for(Int_t j=0; j<GetPhotons()->GetNParticles(); j++){ 

	Particles.push_back(fReadParticles->at(j));
	fcbPhoton = GetPhotons()->Particle(j); //TLorentzVector
	Particles[j]->SetP4(fcbPhoton);
	Particles[j]->SetPDGcode(22);

      } //Closing For NParticles 

      for(Int_t i=0; i<GetTagger()->GetNTagged() ;i++){
 	
	ftaggedTime = GetTagger()->GetTaggedTime(i);
	countb4 = countb4+1;
	
	if ( ftaggedTime > 40 && ftaggedTime <80 ) {

	  countaft=countaft+1;
	  counter = countb4 - countaft;
	  
	  Particles.push_back(fReadParticles->at(i+3-counter));
	  fmultiplicity = GetTrigger()->GetMultiplicity();
	  fenergySum =GetTrigger()->GetEnergySum();
	  fenergyBeam = GetTagger()->GetTaggedEnergy(i);
	  beam.SetXYZM(0.,0.,fenergyBeam,0.);  

	  //Linear polarisation
          flinPol = GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(i) );
	  //Circular polarisation
	  Pcirc = CircPol(fenergyBeam, ePol );	
	  
	  fglasgowTaggerPhoton.SetPxPyPzE(0,0,fenergyBeam, fenergyBeam);  //TLorentzVector
	  Particles[i+3-counter]->SetP4(fglasgowTaggerPhoton);
	  Particles[i+3-counter]->SetPDGcode(-22);
	  Particles[i+3-counter]->SetTime(ftaggedTime);
	  if (flinPol>0)  Particles[i+3-counter]->SetVertex(flinPol,0,Pcirc*fbeamHelicity);
	  if (flinPol<0)  Particles[i+3-counter]->SetVertex(0,flinPol,Pcirc*fbeamHelicity);
	  
	} //closing if TaggedTime
	
      } //Closing for NTagged.
      
    } //closing NPIDHits if
    
  } //closing Number Protons ==3 
    
    //Do Proton channel here 

    //********************************************************************************************************************************
    
 
  if (GetPhotons()->GetNParticles()==2 && GetRootinos()->GetNParticles()==1){

    //PID info
    if (NPidhits>0){

      if (NPidhits>0){
        fpidIndex = GetDetectorHits()->GetPIDHits(0);    //The parameter in the Npidhits is the number of pid hits in the event while getPIDHits is the element number of a hit
        fpidPhi = PIDElemPhi[fpidIndex]; //here
      } //closingpihits

      else{  
    
        fpidPhi = -10;
      }
      

      while(fNin>fReadParticles->size()){
	fReadParticles->push_back(new THSParticle());
      }
         
      Particles.push_back(fReadParticles->at(0));
      frootino = GetRootinos()->Particle(0);
     // std::cout << "X= " <<frootino.X() << "   Y= " <<frootino.Y() <<"    Z= " <<frootino.Z() << std::endl;
      if (frootino.X()==0 && frootino.Y()==0)testcounter= testcounter +1 ;
      Particles[0]->SetP4(frootino);
      Particles[0]->SetPDGcode(2212);
      frootinoPhi= frootino.Phi();

	rootinoClustE =frootino.E();// GetTracks()->GetClusterEnergy(0);
	rootinoTheta =frootino.Theta();//GetTracks()->GetTheta(0);
	rootinoPhi =frootino.Phi();//GetTracks()->GetPhi(0);
	rootinoTime =0;//GetTracks()->GetTime(0);
	rootinoClustS =0;//GetTracks()->GetClusterSize(0);
	rootinoClustC =0;//GetTracks()->GetCentralCrystal(0);
	rootinoVetoC =0;//GetTracks()->GetCentralVeto(0);
	rootinoDet =0;//GetTracks()->GetDetectors(0);
	rootinoVetoE =0;//GetTracks()->GetVetoEnergy(0);
	rootinoCham1E =0;//GetTracks()->GetMWPC0Energy(0);
	rootinoCham2E =0;//GetTracks()->GetMWPC1Energy(0);



      for(Int_t k=0; k<GetPhotons()->GetNParticles(); k++){

	Particles.push_back(fReadParticles->at(k+1));
      	fcbPhoton = GetPhotons()->Particle(k); //TLorentzVector
	Particles[k+1]->SetP4(fcbPhoton);
	Particles[k+1]->SetPDGcode(22);
	      
      } //Closing For NParticles 
        
      for(Int_t l=0; l<GetTagger()->GetNTagged() ;l++){
      
	ftaggedTime = GetTagger()->GetTaggedTime(l);

	countb4 = countb4 +1;


       	if ( ftaggedTime > 40 && ftaggedTime <80 ) { 

	  countaft = countaft +1;
	  counter = countb4 - countaft;
	  Particles.push_back(fReadParticles->at(l+3-counter));
       	  fmultiplicity = GetTrigger()->GetMultiplicity();
    	  fenergySum =GetTrigger()->GetEnergySum();
    	  fenergyBeam = GetTagger()->GetTaggedEnergy(l);
    	  beam.SetXYZM(0.,0.,fenergyBeam,0.);


      	  fglasgowTaggerPhoton.SetPxPyPzE(0,0,fenergyBeam, fenergyBeam);  //TLorentzVector

	//Linear Polarisation
          flinPol = GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(l) );
	//Circular Polarisation
	  Pcirc = CircPol(fenergyBeam, ePol );	

	  
	  Particles[l+3-counter]->SetP4(fglasgowTaggerPhoton);
	  Particles[l+3-counter]->SetPDGcode(-22);
	  Particles[l+3-counter]->SetTime(ftaggedTime);
	  if (flinPol>0)  Particles[l+3-counter]->SetVertex(flinPol,0,Pcirc*fbeamHelicity);
	  if (flinPol<0)  Particles[l+3-counter]->SetVertex(0,flinPol,Pcirc*fbeamHelicity);
	  
	  
	} //closing if TaggedTime
	
      } //Closing for NTagged.
        
    } //closing NPIDHits if
  
  } //closing 2 photons and 1 rootino.  



  //*************************************************************************************************************************************
else{

	rootinoClustE =-111;  
        rootinoTheta =-111;
        rootinoPhi =-111;
        rootinoTime =-111;
        rootinoClustS =-111;
        rootinoClustC =-111;
        rootinoVetoC =-111;
        rootinoDet =-111;
        rootinoVetoE =-111;
        rootinoCham1E =-111;
        rootinoCham2E =-111;


}
  

  //} //closing if 3 particles


  treePi0->Fill();
  
  nEventsWritten++;

  
} //closing function


void	chrisForTHSProject::ProcessScalerRead()
{
  // Fill Tagger Scalers
  //FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	chrisForTHSProject::Write()
{

  treePi0->Write(); 
  treePi0->Reset();

  return 0; //Added CAM 28/09/16

  // Write all GH1's and TObjects defined in this class
  //    return GTreeManager::Write(); //Removed CAM 28/09/16
}



Double_t CircPol( Double_t Eg , Double_t ePol   ){
//function to calculate the circ photon polarisation produced from elec beam given elec beam pol, photon energy and elec beam energy 1557Mev

Double_t E0 = 1557;

Double_t Pg = ( (ePol*Eg)/E0 )* (  ( 1+ ((1/3)*(1-(Eg/E0))) )/(1 - ((2/3)*(1- (Eg/E0))) +  pow( (1 - (Eg/E0)), 2 )  ) );


return Pg;
}



