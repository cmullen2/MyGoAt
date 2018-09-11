#include "chrisForTHSProject.h"
#include "TROOT.h"


Double_t CircPol(Double_t , Double_t );


chrisForTHSProject::chrisForTHSProject()
{ 
  gROOT->ProcessLine("#include <vector>");

  treePi0 =  new TTree("HSParticles","Event selection tree"); //Proton/Neutron pi0 final state tree   


  //pi0 Tree branches
  treePi0->Branch("Particles",&Particles);
  treePi0->Branch("Generated",&Generated);
  treePi0->Branch("EventInfo",&fEventInfo);


}

chrisForTHSProject::~chrisForTHSProject()
{

}

Bool_t	chrisForTHSProject::Init()
{
  cout << "Initialising physics analysis..." << endl;
  cout << "--------------------------------------------------" << endl << endl;
  Int_t mcc=1; // 1 for simulation, 0 for  production data
  if(mcc){
    cout <<"MC File detected. Processing as MC file with Neutron spectator and Proton Participant" <<endl;
    cout << "Please check the histograms for Truth branch for Elab,Plab and dircos for each particle" <<endl;
    cout << "Truth->Draw('dircos[0][1]') Truth->Draw('truthElab[0]')  etc ..." <<endl;
  }

  std::string config = ReadConfig("Period-Macro");
  if( sscanf(config.c_str(),"%d\n", &period) == 1 ) usePeriodMacro = 1;
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

  Int_t mc =1; //0 for production data, non-zero for simulation

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

    //Linear polarisation plane setting
    if(GetLinpol()->GetPolarizationPlane()==0) fedgePlane = -1; // Para
    if(GetLinpol()->GetPolarizationPlane()==1) fedgePlane = 1; // Perp
    if(GetLinpol()->GetPolarizationPlane()==2) fedgePlane = 0; //Moeller or other


  }
  else{
    //MC needs to deal with event no, ePol and linpol and where para or perp. Also needs to deal with truth values by adding to Generated
    //IMPORTANT! A strict naming convention is applied to the files that requires Neg and Pos or Fla to be in the name for the polarisation plane. blah_Fla.root _Pos.root etc

    Generated.clear();
    std::string outFile = outputFile->GetName();
    std::string filert = outFile.substr(outFile.size() -8  );
    std::string planesetting = filert.substr(0,3);

    if(planesetting=="Neg" ) fileNo = "11110000";
    if(planesetting=="Pos" ) fileNo = "22220000";
    if(planesetting=="Fla" ) fileNo = "33330000";

    if(planesetting!="Fla" && planesetting!="Neg" && planesetting!="Pos" ){
      throw std::invalid_argument("Received incorrected FileName Ending!  You have ignored the strict filename requirements which describe the polarisation state. Please make sure the filenames end in Fla, Neg, or Pos. ");
    }

    Int_t tempEventno = GetEventNumber();
    feventNo = std::stod(fileNo) + tempEventno;


    //Linear Polarsation plane setting
    if(planesetting=="Neg" ) fedgePlane = -1; //Equating Neg file with Para polarisation state
    if(planesetting=="Pos" ) fedgePlane = 1;
    if(planesetting=="Fla" ) fedgePlane = 0;  //Flat simulation is equated with Amo runs

    //Linear Polarisation value
    //Use if !mc for the linPol parameter and assign the mc value here based on name 
    if(planesetting=="Neg" ) flinPol = 1; //Equating Neg file with Para polarisation with max value CAM CHANGE DUE TO linpol no longer being neg
    if(planesetting=="Pos" ) flinPol = 1;
    if(planesetting=="Fla" ) flinPol = 0; 


    for(Int_t j=1; j<(GetTruth()->GetfNMC()+1); j++){ //push back the 4 particles but not beam here?

       THSParticle Gen;
       if(j<GetTruth()->GetfNMC()){
	Gen.SetXYZT(1000 * (GetTruth()->GettruthPlab(j)) * (GetTruth()->Getdircos(0+(j*3)) ), 1000 *  (GetTruth()->GettruthPlab(j)) * (GetTruth()->Getdircos(1+(j*3))), 1000 * (GetTruth()->GettruthPlab(j)) * (GetTruth()->Getdircos(2+(j*3))), 1000*GetTruth()->GettruthElab(j));
      }
      else{
	Gen.SetXYZT(1000 * (GetTruth()->GettruthBeam(3)) * (GetTruth()->GettruthBeam(0)), 1000 *  (GetTruth()->GettruthBeam(3)) * (GetTruth()->GettruthBeam(1)), 1000 * (GetTruth()->GettruthBeam(3)) * (GetTruth()->GettruthBeam(2)), 1000*GetTruth()->GettruthBeam(4));
      }
      Gen.SetPDGcode(generatedPDGs[j-1]);
      Gen.SetVertex(GetTruth()->GettruthVertex(0),GetTruth()->GettruthVertex(1),GetTruth()->GettruthVertex(2));
      Generated.push_back(Gen);
    }

  }


  //Clearing the vectors.
  Particles.clear();

  //Mott measurements for electron beam polarisation  //CAM NEED more runs(sept)
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
  if (runNo==33330000 ) ePol = MCMottFla;


  //PID information
  NPidhits = GetDetectorHits()->GetNPIDHits();
  
  //Fill Beam Helicity (Two different determined states when helicity is 0 or 1, others are error codes) Chosen 1 and -1 for helicity from now on. 
  //Check for error codes here?? Need to decide if it is necessary? Should only be a small effect so come back to it once everything is working
  Bool_t Helicity = GetTrigger()->GetHelicity();
  if(Helicity == 1){
    fbeamHelicity = 1;
  }
  if(Helicity == 0){
    fbeamHelicity = -1;
  }
  
  if(GetTagger()->GetNTagged()<200){  
    if(GetTagger()->GetNTagged()>0){  
      //************************************************************************************************************************

      if ( (GetPhotons()->GetNParticles()==3)){

	//PID info
	if (NPidhits==0){
	  if (GetTagger()->GetNTagged()>100)cout<< GetTagger()->GetNTagged()<<endl; 
	  for(Int_t j=0; j<GetPhotons()->GetNParticles(); j++){ 

	    THSParticle part;
	    //Particles.push_back(fReadParticles->at(j));
	    fcbPhoton = GetPhotons()->Particle(j); //TLorentzVector
	    part.SetP4(fcbPhoton);
	    part.SetPDGcode(22);
	    part.SetTime(GetTracks()->GetTime(j));//Should this be using gettrack index then get time on the index?Yes but since the particles are all photons trackNum=photonNum,see below for alternative
	    part.SetDetector(GetTracks()->GetDetectors(j)); //DETECTOR_NONE = 0,DETECTOR_NaI = 1, DETECTOR_PID = 2, DETECTOR_MWPC = 4, DETECTOR_BaF2 = 8, DETECTOR_PbWO4 = 16, DETECTOR_Veto = 32,(Additive)
	    Particles.push_back(part);

	  } //Closing For NParticles 

	  fEventInfo.SetBeamHel(fbeamHelicity );
	  fEventInfo.SetTarPolDir(fedgePlane);

	  for(Int_t i=0; i<GetTagger()->GetNTagged() ;i++){
	    THSParticle part;

	    ftaggedTime = GetTagger()->GetTaggedTime(i);
	    fmultiplicity = GetTrigger()->GetMultiplicity();
	    fenergySum =GetTrigger()->GetEnergySum();
	    fenergyBeam = GetTagger()->GetTaggedEnergy(i);
	    if(fenergyBeam==0)cout << fenergyBeam << " BEam Energy " << GetTagger()->GetTaggedChannel(i) <<" Tagger Channel"<<endl;
	    beam.SetXYZM(0.,0.,fenergyBeam,0.);  
	    //Linear Polarisation
	    if(!mc)flinPol =( GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(i) ) );
	    fEventInfo.SetTarPol(flinPol);
	    //Circular polarisation
	    fEventInfo.SetBeamPol( CircPol(fenergyBeam, ePol ));
	    ftaggChannel = GetTagger()->GetTaggedChannel(i);
	    part.SetDetector(ftaggChannel);
	    fglasgowTaggerPhoton.SetPxPyPzE(0,0,fenergyBeam, fenergyBeam);  //TLorentzVector(HSLorentzVector)
	    part.SetP4(fglasgowTaggerPhoton);
	    part.SetPDGcode(-22);
	    part.SetTime(ftaggedTime);
	    Particles.push_back(part);
	  } //Closing for NTagged.
      
	} //closing NPIDHits if
    
      } //closing Number Protons ==3 
    
      //Do Proton channel here 

      //********************************************************************************************************************************
 
      if (GetPhotons()->GetNParticles()==2 && GetRootinos()->GetNParticles()==1){

	//PID info
	if (NPidhits>0){


	  THSParticle part;
	  //Particles.push_back(fReadParticles->at(0));
	  frootino = GetRootinos()->Particle(0);
	  particleindex=GetRootinos()->GetTrackIndex(0) ;//Set as -1 in header so will break if has any unexpected behaviour
	  part.SetP4(frootino);
	  part.SetPDGcode(2212);
	  part.SetTime(GetTracks()->GetTime(particleindex));
	  part.SetDetector(GetTracks()->GetDetectors(particleindex));
	  Particles.push_back(part);
	  part.Clear();//clear or reset??


	  for(Int_t k=0; k<GetPhotons()->GetNParticles(); k++){
	    THSParticle part2;
	    particleindex=GetPhotons()->GetTrackIndex(k);
	    fcbPhoton = GetPhotons()->Particle(k); //TLorentzVector
	    part2.SetP4(fcbPhoton);
	    part2.SetPDGcode(22);
	    part2.SetTime(GetTracks()->GetTime(particleindex));
	    part2.SetDetector(GetTracks()->GetDetectors(particleindex));
	    Particles.push_back(part2);  
	  } //Closing For NParticles 

	  fEventInfo.SetBeamHel(fbeamHelicity);
	  fEventInfo.SetTarPolDir(fedgePlane);

	  for(Int_t l=0; l<GetTagger()->GetNTagged() ;l++){
      
	    THSParticle part3;
	    ftaggedTime = GetTagger()->GetTaggedTime(l);
	    fmultiplicity = GetTrigger()->GetMultiplicity();
	    fenergySum =GetTrigger()->GetEnergySum();
	    fenergyBeam = GetTagger()->GetTaggedEnergy(l);
	    if(fenergyBeam==0)cout << fenergyBeam << " BEam Energy " << GetTagger()->GetTaggedChannel(l) <<" Tagger Channel"<<endl;
	    beam.SetXYZM(0.,0.,fenergyBeam,0.);
	    fglasgowTaggerPhoton.SetPxPyPzE(0,0,fenergyBeam, fenergyBeam);  //TLorentzVector
	    //Linear Polarisation
	    if(!mc) flinPol = ( GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(l) ) );
	    fEventInfo.SetTarPol(flinPol);
	    //Circular Polarisation
	    fEventInfo.SetBeamPol( CircPol(fenergyBeam, ePol ));	
	    //Tagger Channel
	    ftaggChannel  = GetTagger()->GetTaggedChannel(l) ;
	    part3.SetDetector(ftaggChannel);
	    part3.SetP4(fglasgowTaggerPhoton);
	    part3.SetPDGcode(-22);
	    part3.SetTime(ftaggedTime);
	    Particles.push_back(part3);

	  } //Closing for NTagged.

	  //Commented out for SIMS	  }//Closing if pseudoVertexZ
        
	} //closing NPIDHits if
  
      } //closing 2 photons and 1 rootino.  


      //*************************************************************************************************************************************
  
    }//closing if tagged exists(>0)
    //} //closing if 3 particles
  }//closing if <200 tagged particles

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

  return 0; 

 
}



Double_t CircPol( Double_t Eg , Double_t ePol   ){
  //function to calculate the circ photon polarisation produced from elec beam given elec beam pol, photon energy and elec beam energy 1557Mev

  Double_t E0 = 1557;

  Double_t Pg = ( (ePol*Eg)/E0 )* (  ( 1+ ((1/3)*(1-(Eg/E0))) )/(1 - ((2/3)*(1- (Eg/E0))) +  pow( (1 - (Eg/E0)), 2 )  ) );


  return Pg;
}



