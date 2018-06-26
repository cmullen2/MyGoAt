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
  Int_t mcc=1; // 1 for simulation, 0 for  production data
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
  cout << "Setting up Graphical cuts " << endl;
  //Graphical cuts
  ProtonCut = OpenCut("PProtonCut.root","PProtonCut");
  PipCut = OpenCut("NPipCut.root","NPipCut");
  PimCut = OpenCut("PPimCut.root","PPimCut");

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


	  //Photons
	  fcrystalphot = GetPhotons()->Particle(0);
	  fcrystalphotindex =GetPhotons()->GetTrackIndex(0);
	  fcrystalphotmass= fcrystalphot.M();
	  fcrystalphotphi= fcrystalphot.Phi();
	  fcrystalphottheta= fcrystalphot.Theta();
	  fcrystalphotpidE = GetTracks()->GetVetoEnergy(fcrystalphotindex);
	  fcrystalphotmwpc0E = GetTracks()->GetMWPC0Energy(fcrystalphotindex);
	  fcrystalphotmwpc1E = GetTracks()->GetMWPC1Energy(fcrystalphotindex);
	  //Rootinos
	  fcrystalroot = GetRootinos()->Particle(0);
	  fcrystalrootindex =GetRootinos()->GetTrackIndex(0);
	  fcrystalrootmass= fcrystalroot.M();
	  fcrystalrootphi= fcrystalroot.Phi();
	  fcrystalroottheta= fcrystalroot.Theta();
	  fcrystalrootpidE = GetTracks()->GetVetoEnergy(fcrystalrootindex);
	  fcrystalrootmwpc0E = GetTracks()->GetMWPC0Energy(fcrystalrootindex);
	  fcrystalrootmwpc1E = GetTracks()->GetMWPC1Energy(fcrystalrootindex);

	  fcrystalcoplan=fcrystalroot.TLorentzVector::DeltaPhi(-fcrystalphot);


	  if( PipCut->IsInside(fcrystalroot.E(),fcrystalrootpidE) ){


	    fcrystalrootIn = GetRootinos()->Particle(0);
	    fcrystalrootindexIn = GetRootinos()->GetTrackIndex(0);
	    fcrystalrootmassIn = fcrystalroot.M();
	    fcrystalrootphiIn = fcrystalroot.Phi();
	    fcrystalrootthetaIn = fcrystalroot.Theta();
	    fcrystalrootpidEIn = GetTracks()->GetVetoEnergy(fcrystalrootindex);
	    fcrystalrootmwpc0EIn = GetTracks()->GetMWPC0Energy(fcrystalrootindex);
	    fcrystalrootmwpc1EIn = GetTracks()->GetMWPC1Energy(fcrystalrootindex);

	    // for(Int_t iii=0; iii<GetMWPCHitsChris()->GetNMWPCHitsChrisChamber1(); iii++){
	
	    //   fchamber1Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber1X(iii),GetMWPCHitsChris()->GetMWPCChamber1Y(iii),GetMWPCHitsChris()->GetMWPCChamber1Z(iii));
	    //   MWPC1Hits.push_back(fchamber1Vec);
	    // }


	    // for(Int_t jjj=0; jjj<GetMWPCHitsChris()->GetNMWPCHitsChrisChamber2(); jjj++){

	    //   fchamber2Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber2X(jjj),GetMWPCHitsChris()->GetMWPCChamber2Y(jjj),GetMWPCHitsChris()->GetMWPCChamber2Z(jjj));
	    //   MWPC2Hits.push_back(fchamber2Vec);
	    // }

      
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
	      //      part.SetMWPC0Hits(MWPC1Hits); //a vector of TVector3's.
	      //      part.SetMWPC1Hits(MWPC2Hits);
	      //part.Set(GetTracks()->Get(fphotindex));
	      Particles.push_back(part);
	    } //Closing For NParticles 



	    for(Int_t k=0; k<GetRootinos()->GetNParticles(); k++){

	      THSParticle part;
	      frootino = GetRootinos()->Particle(k); //TLorentzVector
	      particleindex=GetRootinos()->GetTrackIndex(k);
	      part.SetP4(fcbPhoton);
	      part.SetPDGcode(211);
	      part.SetTime(GetTracks()->GetTime(particleindex));
	      part.SetDetector(GetTracks()->GetDetectors(particleindex));
	      part.SetEPid(GetTracks()->GetVetoEnergy(particleindex));
	      part.SetEMWPC0(GetTracks()->GetMWPC0Energy(particleindex));
	      part.SetEMWPC1(GetTracks()->GetMWPC1Energy(particleindex));
	      //  part.SetMWPC0Hits(MWPC1Hits); //a vector of TVector3's.
	      //  part.SetMWPC1Hits(MWPC2Hits);
	      Particles.push_back(part);

	    } //Closing For NParticles 

	    MWPC1Hits.clear();
	    MWPC2Hits.clear();

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
	
	  }//Closing if IsInside graphical cut      

	} //closing NPIDHits if
    
      } //closing Number rootinos and photons
  
      else{
	//Photons
	fcrystalphot.SetXYZM(-10000,-10000,-10000,-10000);
	fcrystalphotindex=-1; 
	fcrystalphotmass=-10000; 
	fcrystalphotphi=-10000;
	fcrystalphottheta=-10000;
	fcrystalphotpidE =-10000;
	fcrystalphotmwpc0E =-10000; 
	fcrystalphotmwpc1E =-10000; 
	//Rootinos
	fcrystalroot.SetXYZM(-10000,-10000,-10000,-10000);
	fcrystalrootindex =-1;
	fcrystalrootmass=-10000;
	fcrystalrootphi= -10000;
	fcrystalroottheta= -10000;
	fcrystalrootpidE = -10000;
	fcrystalrootmwpc0E = -10000;
	fcrystalrootmwpc1E = -10000;

	fcrystalcoplan=-10000;

	
	fcrystalrootIn.SetXYZM(-10000,-10000,-10000,-10000);
	fcrystalrootindexIn =-1;
	fcrystalrootmassIn=-10000;
	fcrystalrootphiIn= -10000;
	fcrystalrootthetaIn= -10000;
	fcrystalrootpidEIn = -10000;
	fcrystalrootmwpc0EIn = -10000;
	fcrystalrootmwpc1EIn = -10000;

      }  
      //Do Proton channel here 

      //********************************************************************************************************************************
      //NEED TO ADD PARTICLES VECTORS TO THSFINALSTATE FOR NEUTRONS 2122 and for charged but dunno +or- currently using vec0.
      // How do I identify which pid hit and mwpc hit came from which particle, should be a track associated with each right so use position of phi hit in pid and phi in mwpc to check. What about scattered? Also Coplanarity should also be a good discriminator since back to back in pid and mwpc and crystal.
 
      if (GetPhotons()->GetNParticles()==0 && GetRootinos()->GetNParticles()==2){
	ftopology=-2;
	//PID info
	if (NPidhits>0){

	  if (NPidhits>0){
	    fpidIndex = GetDetectorHits()->GetPIDHits(0);    //The parameter in the Npidhits is the number of pid hits in the event while getPIDHits is the element number of a hit
	    fpidPhi = PIDElemPhi[fpidIndex]; //here
	  } //closingpihits

	  else{  
    
	    fpidPhi = -10;
	  }

	  //Rootinos, channel first look
	  fcrystalroot1 = GetRootinos()->Particle(0);
	  fcrystalrootindex1 =GetRootinos()->GetTrackIndex(0);
	  fcrystalrootmass1= fcrystalroot1.M();
	  fcrystalrootphi1= fcrystalroot1.Phi();
	  fcrystalroottheta1= fcrystalroot1.Theta();
	  fcrystalrootpidE1 = GetTracks()->GetVetoEnergy(fcrystalrootindex1);
	  fcrystalrootmwpc0E1 = GetTracks()->GetMWPC0Energy(fcrystalrootindex1);
	  fcrystalrootmwpc1E1 = GetTracks()->GetMWPC1Energy(fcrystalrootindex1);

	  fcrystalroot2 = GetRootinos()->Particle(1);
	  fcrystalrootindex2 =GetRootinos()->GetTrackIndex(1);
	  fcrystalrootmass2= fcrystalroot2.M();
	  fcrystalrootphi2= fcrystalroot2.Phi();
	  fcrystalroottheta2= fcrystalroot2.Theta();
	  fcrystalrootpidE2 = GetTracks()->GetVetoEnergy(fcrystalrootindex2);
	  fcrystalrootmwpc0E2 = GetTracks()->GetMWPC0Energy(fcrystalrootindex2);
	  fcrystalrootmwpc1E2 = GetTracks()->GetMWPC1Energy(fcrystalrootindex2);

	  fcrystalcoplan2=fcrystalroot2.TLorentzVector::DeltaPhi(-fcrystalroot1);

	  //Form two new variables one for each particle with the energy correction for proton applied. combos check 1p&2m(1ecorr) 2p&1m(2ecorr) 
          fcrystalroot1ECorr = ProtonELossCorrection(fcrystalroot1.Theta(), fcrystalroot1.E());
	  fcrystalroot2ECorr = ProtonELossCorrection(fcrystalroot2.Theta(), fcrystalroot2.E());

	  if(  ( (ProtonCut->IsInside(fcrystalroot1ECorr,fcrystalrootpidE1)) && (PimCut->IsInside(fcrystalroot2.E(),fcrystalrootpidE2)) )  || (    (PimCut->IsInside(fcrystalroot1.E(),fcrystalrootpidE1)) && (ProtonCut->IsInside(fcrystalroot2ECorr,fcrystalrootpidE2)) ) ){
	    //  if(  ( (ProtonCut->IsInside(fcrystalroot1.E(),fcrystalrootpidE1)) && (PimCut->IsInside(fcrystalroot2.E(),fcrystalrootpidE2)) )  || (    (PimCut->IsInside(fcrystalroot1.E(),fcrystalrootpidE1)) && (ProtonCut->IsInside(fcrystalroot2.E(),fcrystalrootpidE2)) ) ){

	    fcrystalroot1In = GetRootinos()->Particle(0);
	    fcrystalrootindex1In = GetRootinos()->GetTrackIndex(0);
	    fcrystalrootmass1In = fcrystalroot1.M();
	    fcrystalrootphi1In = fcrystalroot1.Phi();
	    fcrystalroottheta1In = fcrystalroot1.Theta();
	    fcrystalrootpidE1In = GetTracks()->GetVetoEnergy(fcrystalrootindex1);
	    fcrystalrootmwpc0E1In = GetTracks()->GetMWPC0Energy(fcrystalrootindex1);
	    fcrystalrootmwpc1E1In = GetTracks()->GetMWPC1Energy(fcrystalrootindex1);

	    fcrystalroot2In = GetRootinos()->Particle(1);
	    fcrystalrootindex2In = GetRootinos()->GetTrackIndex(1);
	    fcrystalrootmass2In = fcrystalroot2.M();
	    fcrystalrootphi2In = fcrystalroot2.Phi();
	    fcrystalroottheta2In = fcrystalroot2.Theta();
	    fcrystalrootpidE2In = GetTracks()->GetVetoEnergy(fcrystalrootindex2);
	    fcrystalrootmwpc0E2In = GetTracks()->GetMWPC0Energy(fcrystalrootindex2);
	    fcrystalrootmwpc1E2In = GetTracks()->GetMWPC1Energy(fcrystalrootindex2);


	    // for(Int_t iii=0; iii<GetMWPCHitsChris()->GetNMWPCHitsChrisChamber1(); iii++){
	
	    //   fchamber1Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber1X(iii),GetMWPCHitsChris()->GetMWPCChamber1Y(iii),GetMWPCHitsChris()->GetMWPCChamber1Z(iii));
	    //   MWPC1Hits.push_back(fchamber1Vec);
	    // }


	    // for(Int_t jjj=0; jjj<GetMWPCHitsChris()->GetNMWPCHitsChrisChamber2(); jjj++){

	    //   fchamber2Vec.SetXYZ(GetMWPCHitsChris()->GetMWPCChamber2X(jjj),GetMWPCHitsChris()->GetMWPCChamber2Y(jjj),GetMWPCHitsChris()->GetMWPCChamber2Z(jjj));
	    //   MWPC2Hits.push_back(fchamber2Vec);
	    // }



	    if( (ProtonCut->IsInside(fcrystalroot1ECorr,fcrystalrootpidE1)) && (PimCut->IsInside(fcrystalroot2.E(),fcrystalrootpidE2)) ){
	      //first rootino is proton, second is pi-

	      THSParticle part;
	      frootino = GetRootinos()->Particle(0);
	      particleindex=GetRootinos()->GetTrackIndex(0) ;
	      part.SetP4(frootino);
	      part.SetPDGcode(2212);
	      part.SetTime(GetTracks()->GetTime(particleindex));
	      part.SetDetector(GetTracks()->GetDetectors(particleindex));
	      part.SetEPid(GetTracks()->GetVetoEnergy(particleindex));
	      part.SetEMWPC0(GetTracks()->GetMWPC0Energy(particleindex));
	      part.SetEMWPC1(GetTracks()->GetMWPC1Energy(particleindex));
	      //	      part.SetMWPC0Hits(MWPC1Hits); //a vector of TVector3's.
	      //	      part.SetMWPC1Hits(MWPC2Hits);
	      Particles.push_back(part);

	      THSParticle part2;
	      frootino = GetRootinos()->Particle(1);
	      particleindex=GetRootinos()->GetTrackIndex(1) ;
	      part2.SetP4(frootino);
	      part2.SetPDGcode(-211);
	      part2.SetTime(GetTracks()->GetTime(particleindex));
	      part2.SetDetector(GetTracks()->GetDetectors(particleindex));
	      part2.SetEPid(GetTracks()->GetVetoEnergy(particleindex));
	      part2.SetEMWPC0(GetTracks()->GetMWPC0Energy(particleindex));
	      part2.SetEMWPC1(GetTracks()->GetMWPC1Energy(particleindex));
	      //	      part2.SetMWPC0Hits(MWPC1Hits); //a vector of TVector3's.
	      //	      part2.SetMWPC1Hits(MWPC2Hits);
	      Particles.push_back(part2);
	    }

	    if(       (PimCut->IsInside(fcrystalroot1.E(),fcrystalrootpidE1)) && (ProtonCut->IsInside(fcrystalroot2ECorr,fcrystalrootpidE2))        ){
	      //first rootino is pi-, second is proton

	      THSParticle part;
	      frootino = GetRootinos()->Particle(0);
	      particleindex=GetRootinos()->GetTrackIndex(0) ;//Set as -1 in header so will break if has any unexpected behaviour
	      part.SetP4(frootino);
	      part.SetPDGcode(-211);
	      part.SetTime(GetTracks()->GetTime(particleindex));
	      part.SetDetector(GetTracks()->GetDetectors(particleindex));
	      part.SetEPid(GetTracks()->GetVetoEnergy(particleindex));
	      part.SetEMWPC0(GetTracks()->GetMWPC0Energy(particleindex));
	      part.SetEMWPC1(GetTracks()->GetMWPC1Energy(particleindex));
	      //	      part.SetMWPC0Hits(MWPC1Hits); //a vector of TVector3's.
	      //	      part.SetMWPC1Hits(MWPC2Hits);
	      Particles.push_back(part);

	      THSParticle part2;
	      frootino = GetRootinos()->Particle(1);
	      particleindex=GetRootinos()->GetTrackIndex(1) ;//Set as -1 in header so will break if has any unexpected behaviour
	      part2.SetP4(frootino);
	      part2.SetPDGcode(2212);
	      part2.SetTime(GetTracks()->GetTime(particleindex));
	      part2.SetDetector(GetTracks()->GetDetectors(particleindex));
	      part2.SetEPid(GetTracks()->GetVetoEnergy(particleindex));
	      part2.SetEMWPC0(GetTracks()->GetMWPC0Energy(particleindex));
	      part2.SetEMWPC1(GetTracks()->GetMWPC1Energy(particleindex));
	      //	      part2.SetMWPC0Hits(MWPC1Hits); //a vector of TVector3's.
	      //	      part2.SetMWPC1Hits(MWPC2Hits);
	      Particles.push_back(part2);

	    }

	    MWPC1Hits.clear();
	    MWPC2Hits.clear();

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
        
	  } //Closing if graphical cut options

	} //closing NPIDHits if
  
      } //closing 2 photons and 1 rootino.  

      else{

	//Rootinos
	fcrystalroot1.SetXYZM(-10000,-10000,-10000,-10000);
	fcrystalrootindex1 =-1;
	fcrystalrootmass1=-10000;
	fcrystalrootphi1= -10000;
	fcrystalroottheta1= -10000;
	fcrystalrootpidE1 = -10000;
	fcrystalrootmwpc0E1 = -10000;
	fcrystalrootmwpc1E1 = -10000;

	fcrystalroot2.SetXYZM(-10000,-10000,-10000,-10000);
	fcrystalrootindex2 =-1;
	fcrystalrootmass2=-10000;
	fcrystalrootphi2= -10000;
	fcrystalroottheta2= -10000;
	fcrystalrootpidE2 = -10000;
	fcrystalrootmwpc0E2 = -10000;
	fcrystalrootmwpc1E2 = -10000;

	fcrystalcoplan2=-10000;



	fcrystalroot1In.SetXYZM(-10000,-10000,-10000,-10000);
	fcrystalrootindex1In =-1;
	fcrystalrootmass1In=-10000;
	fcrystalrootphi1In= -10000;
	fcrystalroottheta1In= -10000;
	fcrystalrootpidE1In = -10000;
	fcrystalrootmwpc0E1In = -10000;
	fcrystalrootmwpc1E1In = -10000;

	fcrystalroot2In.SetXYZM(-10000,-10000,-10000,-10000);
	fcrystalrootindex2In =-1;
	fcrystalrootmass2In=-10000;
	fcrystalrootphi2In= -10000;
	fcrystalroottheta2In= -10000;
	fcrystalrootpidE2In = -10000;
	fcrystalrootmwpc0E2In = -10000;
	fcrystalrootmwpc1E2In = -10000;




      }  

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





