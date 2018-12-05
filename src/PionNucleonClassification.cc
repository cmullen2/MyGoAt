#include "PionNucleonClassification.h"
#include "TROOT.h"

PionNucleonClassification::PionNucleonClassification()
{ 
  gROOT->ProcessLine("#include <vector>");
  treePi0 =  new TTree("HSParticles","Event selection tree"); //Proton/Neutron pi0 final state tree   
  treePi0->Branch("Particles",&Particles);
  treePi0->Branch("Generated",&Generated);
  //  treePi0->Branch("EventInfo",&fEventInfo);

}

PionNucleonClassification::~PionNucleonClassification()
{

}

Bool_t	PionNucleonClassification::Init()
{
  cout << "Initialising physics analysis..." << endl;
  cout << "--------------------------------------------------" << endl << endl;

  std::string config = ReadConfig("Period-Macro");
  if( sscanf(config.c_str(),"%d\n", &period) == 1 ) usePeriodMacro = 1;
  if(!InitBackgroundCuts()) return kFALSE;
  if(!InitTargetMass()) return kFALSE;
  if(!InitTaggerChannelCuts()) return kFALSE;
  if(!InitTaggerScalers()) return kFALSE;
  cout << "--------------------------------------------------" << endl;
  return kTRUE;
}

Bool_t	PionNucleonClassification::Start()
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

void	PionNucleonClassification::ProcessEvent()
{

  if(usePeriodMacro == 1)
    {
      if(GetEventNumber() % period == 0)
	cout << "Events: " << GetEventNumber() << "  Events Accepted: " << nEventsWritten << endl;
    }
  //********************************************************************************************************************************
  Generated.clear();
  Particles.clear();

  for(Int_t j=0; j<(GetTruth()->GetfNMC()+1); j++){ //push back the 4 particles but not beam here?

    THSParticle Gen;
    if(j<GetTruth()->GetfNMC()){
      Gen.SetXYZT(1000 * (GetTruth()->GettruthPlab(j)) * (GetTruth()->Getdircos(0+(j*3)) ), 1000 *  (GetTruth()->GettruthPlab(j)) * (GetTruth()->Getdircos(1+(j*3))), 1000 * (GetTruth()->GettruthPlab(j)) * (GetTruth()->Getdircos(2+(j*3))), 1000*GetTruth()->GettruthElab(j));
    }
    else{
      Gen.SetXYZT(1000 * (GetTruth()->GettruthBeam(3)) * (GetTruth()->GettruthBeam(0)), 1000 *  (GetTruth()->GettruthBeam(3)) * (GetTruth()->GettruthBeam(1)), 1000 * (GetTruth()->GettruthBeam(3)) * (GetTruth()->GettruthBeam(2)), 1000*GetTruth()->GettruthBeam(4));
    }
    Gen.SetPDGcode(generatedPDGs[j]);
    Gen.SetVertex(GetTruth()->GettruthVertex(0),GetTruth()->GettruthVertex(1),GetTruth()->GettruthVertex(2));
    Generated.push_back(Gen);
    Gen.Clear();
  }


  //Add in a generated bit similar to bit from PPi0 analysis. Add a particles bit for two rootinos and add new info like cluster time rms central cluster ncrystals in cluster etc. Separate by dcorrect in haspect
  //


  if (GetRootinos()->GetNParticles()==2){ //Need to know each particle!!!!
    for(Int_t i=0;i<GetRootinos()->GetNParticles();i++){
      THSParticle part;
      frootino = GetRootinos()->Particle(i);
      particleindex=GetRootinos()->GetTrackIndex(i);
     cout << "Particle " <<i << "  ClusterEnergy=" << GetRootinos()->GetClusterEnergy(i)<< " particleindex="<<particleindex<<endl;
      part.SetP4(frootino);
      part.SetPDGcode(1E4);
      part.SetEdep(GetRootinos()->GetClusterEnergy(i));
      part.SetDoca(GetRootinos()->GetClusterSize(i));
      part.SetDeltaE(GetRootinos()->GetCentralCrystal(i));
      part.SetEMWPC0(GetRootinos()->GetMWPC0Energy(i));
      part.SetEMWPC1(GetRootinos()->GetMWPC1Energy(i));
      part.SetPreE(GetRootinos()->GetCentralVeto(i));
      part.SetEPid(GetRootinos()->GetVetoEnergy(i));
      part.SetDetector(GetRootinos()->GetDetectors(i));
      //    GetRootinos()->GetPIDEnergy(i);
      Particles.push_back(part);
      part.Clear();
    }
  } //closing 2  rootinos.  
  
  //*************************************************************************************************************************************
 
  treePi0->Fill();
  nEventsWritten++;
} //closing function

void	PionNucleonClassification::ProcessScalerRead()
{
  // Fill Tagger Scalers
  //FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	PionNucleonClassification::Write()
{
  treePi0->Write(); 
  treePi0->Reset();
  return 0; 
}
