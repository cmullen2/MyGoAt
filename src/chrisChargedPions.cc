#include "chrisChargedPions.h"
#include "TROOT.h"


Double_t CircPol(Double_t , Double_t );
TCutG* OpenCut(TString,TString);

chrisChargedPions::chrisChargedPions()
{ 
  gROOT->ProcessLine("#include <vector>");
  treePi =  new TTree("HSParticles","Event selection tree"); //Proton/Neutron pi-+ final state tree   
  //Charged Pion channel Tree branches
  treePi->Branch("Particles",&Particles);
  treePi->Branch("Generated",&Generated);
  treePi->Branch("EventInfo",&fEventInfo);  

}

chrisChargedPions::~chrisChargedPions()
{

}

Bool_t	chrisChargedPions::Init()
{
  cout << "Initialising physics analysis..." << endl;
  cout << "--------------------------------------------------" << endl << endl;
  Int_t mcc=0; // 1 for simulation, 0 for  production data
  if(mcc){

    cout <<"MC File detected. Processing as MC file with Neutron spectator and Proton Participant" <<endl;
    cout << "Please check the histograms for Truth branch for Elab,Plab and dircos for each particle" <<endl;
    cout << "Truth->Draw('dircos[0][1]') Truth->Draw('truthElab[0]')  etc ..." <<endl;
    cout << "Is the correct PID element array being used for the simulation????? Check the header " << endl;

  }

  std::string config = ReadConfig("Period-Macro");
  if( sscanf(config.c_str(),"%d\n", &period) == 1 ) usePeriodMacro = 1;

  target.SetXYZM(0.0,0.0,-65.0,1875.613);	//NEEDS CHANGING currently deuteron
  //	target.SetXYZM(0.0,0.0,-65.0,938.272);	//NEEDS CHANGING only TEMP for a SIM
  fNin = (3 + (GetTagger()->GetNTagged()) ) ; //Estimate of Number of input particles(tagged+ball+rootinos etc.) Can and will calc this!(3+NTagged)

  cout << "Target 4-Vector is (" << target.X() <<"," << target.Y() << ","<< target.Z() <<","<<target.M()<<")"  << endl;
  if(!InitBackgroundCuts()) return kFALSE;
  if(!InitTargetMass()) return kFALSE;
  if(!InitTaggerChannelCuts()) return kFALSE;
  if(!InitTaggerScalers()) return kFALSE;
  cout << "--------------------------------------------------" << endl;
  /*  cout << "Setting up Graphical cuts " << endl;
  //Graphical cuts
  ProtonCut = OpenCut("PProtonCut.root","PProtonCut");
  PipCut = OpenCut("NPipCut.root","NPipCut");
  PimCut = OpenCut("PPimCut.root","PPimCut");
  */
  return kTRUE;

}

Bool_t	chrisChargedPions::Start()
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

void	chrisChargedPions::ProcessEvent()
{

  if(GetEventNumber() == 10)str1 = clock();
  if(GetEventNumber() == 300000)
    {

      end1=clock();
      cout << "Time required for execution: "
	   << (Double_t)(end1-str1)/CLOCKS_PER_SEC
	   << " seconds." << "\n\n";
    }
  auto start1 = std::chrono::high_resolution_clock::now();

  if(usePeriodMacro == 1)
    {
      if(GetEventNumber() % period == 0)
	cout << "Events: " << GetEventNumber() << "  Events Accepted: " << nEventsWritten << endl;
    }

  Int_t mc =0; //0 for production data, non-zero for simulation

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


    if(usePeriodMacro == 1)
      {
	if(GetEventNumber() % period == 0){
	  cout << "Using Production Data var mc = " << mc << endl;
	  cout << "Edge plane setting " << fedgePlane << endl;
	}
      }


  }
  else{
    
    Generated.clear();
    std::string outFile = outputFile->GetName();
    std::string filert = outFile.substr(outFile.size() -8  ); //Pos.root etc is last 8 chars
    std::string planesetting = filert.substr(0,3);  // Isolate pos,fla etc
   
    if(planesetting=="Neg" ) fileNo = "11110000";
    if(planesetting=="Pos" ) fileNo = "22220000";
    if(planesetting=="Fla" ) fileNo = "33330000";

    if(planesetting!="Fla" && planesetting!="Neg" && planesetting!="Pos" )
      {
	throw std::invalid_argument("Received incorrect FileName Ending!  You have ignored the strict filename requirements which describe the polarisation state. Please make sure the filenames end in Fla, Neg, or Pos. ");

      }

    Int_t tempEventno = GetEventNumber();
    feventNo = std::stod(fileNo) + tempEventno;

    //Linear Polarsation plane setting
    if(planesetting=="Neg" ) fedgePlane = -1; //Equating Neg file with Para polarisation state
    if(planesetting=="Pos" ) fedgePlane = 1;
    if(planesetting=="Fla" ) fedgePlane = 0;  //Flat simulation is equated with Amo runs

    if(planesetting=="Neg" ) flinPol = 1; //Equating Neg file with Para polarisation with max value
    if(planesetting=="Pos" ) flinPol = 1;
    if(planesetting=="Fla" ) flinPol = 0; 


    for(Int_t j=0; j<(GetTruth()->GetfNMC()+1); j++){ //push back the 4 particles but not beam here? For the charged channels will be 3 particles + beam
      THSParticle Gen;
      if(j<GetTruth()->GetfNMC()){
	Gen.SetXYZT(1000 * (GetTruth()->GettruthPlab(j)) * (GetTruth()->Getdircos(0+(j*3)) ), 1000 *  (GetTruth()->GettruthPlab(j)) * (GetTruth()->Getdircos(1+(j*3))), 1000 * (GetTruth()->GettruthPlab(j)) * (GetTruth()->Getdircos(2+(j*3))), 1000*GetTruth()->GettruthElab(j));
      }

      else{ //push back the beam here
	Gen.SetXYZT(1000*(GetTruth()->GettruthBeam(3))*(GetTruth()->GettruthBeam(0)), 1000*(GetTruth()->GettruthBeam(3))*(GetTruth()->GettruthBeam(1)), 1000*(GetTruth()->GettruthBeam(3))*(GetTruth()->GettruthBeam(2)), 1000*GetTruth()->GettruthBeam(4));
      }

      Gen.SetPDGcode(generatedPDGs[j]);
      Gen.SetVertex(GetTruth()->GettruthVertex(0),GetTruth()->GettruthVertex(1),GetTruth()->GettruthVertex(2));//What is this line about???
      Generated.push_back(Gen);

    }


    if(usePeriodMacro == 1)
      {
	if(GetEventNumber() % period == 0){
	  cout <<"Using MC Data variable mc = " << mc << endl;
	  cout << "Edge plane setting " << fedgePlane << " Lin pol setting " << flinPol << endl;
	}
      }


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
  if (runNo==33330000 ) ePol = MCMottFla;


  //MWPC information
  targetPosition.SetXYZ(target.X(),target.Y(),target.Z());

  //PID information
  NPidhits = GetDetectorHits()->GetNPIDHits();
  //Helicity information
  Bool_t Helicity = GetTrigger()->GetHelicity();
  if(Helicity == 1){
    fbeamHelicity = 1;
  }
  if(Helicity == 0){
    fbeamHelicity = -1;
  }
  
  fNin= 3+ (GetTagger()->GetNTagged());

  if(GetTagger()->GetNTagged()<200){  
    if(GetTagger()->GetNTagged()>0){  
      //************************************************************************************************************************

      if ( (GetPhotons()->GetNParticles()==1 && GetRootinos()->GetNParticles()==1 )){
	ftopology=2;
	//PID info
	if (NPidhits>0){
	  for(Int_t j=0; j<GetPhotons()->GetNParticles(); j++){ 
	    THSParticle part;
	    fcbPhoton = GetPhotons()->Particle(j); //TLorentzVector
	    fphotindex =GetPhotons()->GetTrackIndex(j);
	    part.SetP4(fcbPhoton);
	    part.SetPDGcode(2112);
	    part.SetTime(GetTracks()->GetTime(fphotindex));
	    part.SetDetector(GetTracks()->GetDetectors(fphotindex));
	    part.SetEPid(GetTracks()->GetVetoEnergy(fphotindex));
	    part.SetEMWPC0(GetTracks()->GetMWPC0Energy(fphotindex));
	    part.SetEMWPC1(GetTracks()->GetMWPC1Energy(fphotindex));
	    Particles.push_back(part);

	  } //Closing For NParticles 
	
	  for(Int_t k=0; k<GetRootinos()->GetNParticles(); k++){
	    THSParticle part;
	    frootino = GetRootinos()->Particle(k); //TLorentzVector
	    particleindex=GetRootinos()->GetTrackIndex(k);
	    part.SetP4(frootino);
	    part.SetPDGcode(211);
	    part.SetTime(GetTracks()->GetTime(particleindex));
	    part.SetDetector(GetTracks()->GetDetectors(particleindex));
	    part.SetEPid(GetTracks()->GetVetoEnergy(particleindex));
	    part.SetEMWPC0(GetTracks()->GetMWPC0Energy(particleindex));  
	    part.SetEMWPC1(GetTracks()->GetMWPC1Energy(particleindex));
	    Particles.push_back(part);
	
	  } //Closing For NParticles 
	
	  fEventInfo.SetBeamHel(fbeamHelicity);  
	  fEventInfo.SetTarPolDir(fedgePlane);
	
	  for(Int_t i=0; i<GetTagger()->GetNTagged() ;i++){
	    THSParticle part;
	    ftaggedTime = GetTagger()->GetTaggedTime(i);
	    fmultiplicity = GetTrigger()->GetMultiplicity();
	    fenergySum =GetTrigger()->GetEnergySum();
	    fenergyBeam = GetTagger()->GetTaggedEnergy(i);
	    beam.SetXYZM(0.,0.,fenergyBeam,0.);  
	    //Linear Polarisation
	    if(!mc)flinPol =( GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(i) ) );
	    fEventInfo.SetTarPol(flinPol);
	    //Circular polarisation
	    fEventInfo.SetBeamPol( CircPol(fenergyBeam, ePol ));	
	    //Tagger Channel
	    ftaggChannel = GetTagger()->GetTaggedChannel(i);
	    // part.SetEdgePlane(fedgePlane);
	    part.SetDetector(ftaggChannel);
	    fglasgowTaggerPhoton.SetPxPyPzE(0,0,fenergyBeam, fenergyBeam);  //TLorentzVector
	    part.SetP4(fglasgowTaggerPhoton);
	    part.SetPDGcode(-22);
	    part.SetTime(ftaggedTime);
	    part.SetVertex(flinPol,0,Pcirc*fbeamHelicity);
	    Particles.push_back(part); 

	  } //Closing for NTagged.
	
	} //closing NPIDHits if
    
      } //closing Number rootinos and photons
  
      //Do Proton channel here 

      //********************************************************************************************************************************
      if (GetPhotons()->GetNParticles()==0 && GetRootinos()->GetNParticles()==2){
	ftopology=-2;
	//PID info
	if (NPidhits>0){



	  for(Int_t k=0; k<GetRootinos()->GetNParticles(); k++){
	    THSParticle part2;
	    frootino2 = GetRootinos()->Particle(k); //TLorentzVector
	    particleindex=GetRootinos()->GetTrackIndex(k);
	    part2.SetP4(frootino2);
	    part2.SetPDGcode(1E4);
	    part2.SetTime(GetTracks()->GetTime(particleindex));
	    part2.SetDetector(GetTracks()->GetDetectors(particleindex));
	    part2.SetEPid(GetTracks()->GetVetoEnergy(particleindex));
	    part2.SetEMWPC0(GetTracks()->GetMWPC0Energy(particleindex));  
	    part2.SetEMWPC1(GetTracks()->GetMWPC1Energy(particleindex));
	    Particles.push_back(part2);
	
	  } //Closing For NParticles 

	  fEventInfo.SetBeamHel(fbeamHelicity);
	  fEventInfo.SetTarPolDir(fedgePlane);

	  for(Int_t l=0; l<GetTagger()->GetNTagged() ;l++){
	    ftaggedTime = GetTagger()->GetTaggedTime(l);
	    THSParticle part;
	    fmultiplicity = GetTrigger()->GetMultiplicity();
	    fenergySum =GetTrigger()->GetEnergySum();
	    fenergyBeam = GetTagger()->GetTaggedEnergy(l);
	    beam.SetXYZM(0.,0.,fenergyBeam,0.);
	    fglasgowTaggerPhoton.SetPxPyPzE(0,0,fenergyBeam, fenergyBeam);  //TLorentzVector
	    //Linear Polarisation
	    if(!mc) flinPol =( GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(l) ) );
	    fEventInfo.SetTarPol(flinPol);
	    //Circular Polarisation
	    fEventInfo.SetBeamPol( CircPol(fenergyBeam, ePol ));
	    //Tagger Channel
	    ftaggChannel  = GetTagger()->GetTaggedChannel(l) ;
	    // part.SetEdgePlane(fedgePlane);
	    part.SetDetector(ftaggChannel);
	    part.SetP4(fglasgowTaggerPhoton);
	    part.SetPDGcode(-22);
	    part.SetTime(ftaggedTime);
	    part.SetVertex(flinPol,0,Pcirc*fbeamHelicity); //Perp
	    Particles.push_back(part);
	
	  } //Closing for NTagged.
 
	} //closing NPIDHits if
  
      } //closing 2 photons and 1 rootino.  

      //*************************************************************************************************************************************

    } // If tagged exists ()>0

  }// If tagged photons <200


  treePi->Fill();
  
  nEventsWritten++;

  auto finish = std::chrono::high_resolution_clock::now();  

  std::chrono::duration<double> elapsed = finish -start1;
  //cout << " Elapsed for function: " << elapsed.count() << endl;

} //closing function


void	chrisChargedPions::ProcessScalerRead()
{
  // Fill Tagger Scalers
  //FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	chrisChargedPions::Write()
{

  treePi->Write(); 
  treePi->Reset();

  return 0; 

 
}



Double_t CircPol( Double_t Eg , Double_t ePol   ){
  //function to calculate the circ photon polarisation produced from elec beam given elec beam pol, photon energy and elec beam energy 1557Mev

  Double_t E0 = 1557;

  Double_t Pg = ( (ePol*Eg)/E0 )* (  ( 1+ ((1/3)*(1-(Eg/E0))) )/(1 - ((2/3)*(1- (Eg/E0))) +  pow( (1 - (Eg/E0)), 2 )  ) );


  return Pg;
}



TCutG* OpenCut(TString file, TString cutname){

  TCutG *cut;

  TFile cutFile(file,"READ");
  cutFile.GetObject(cutname,cut);
  cutFile.Close();
 
  return cut;


}





