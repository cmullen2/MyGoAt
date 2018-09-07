#include "PionEnergyCorrection.h"
#include "TROOT.h"

PionEnergyCorrection::PionEnergyCorrection()
{ 
  gROOT->ProcessLine("#include <vector>");
  treePi0 =  new TTree("HSParticles","Event selection tree"); //Proton/Neutron pi0 final state tree   
  treePi0->Branch("Particles",&Particles);
  treePi0->Branch("Generated",&Generated);
  treePi0->Branch("EventInfo",&fEventInfo);

  /*  //pi0 Tree branches
      treePi0->Branch("TruthPLab",&TruthPLab);
      treePi0->Branch("TruthELab",&TruthELab);
      treePi0->Branch("TruthTheta",&TruthTheta);
      treePi0->Branch("DetectedClusterE",&DetectedClusterE);
      treePi0->Branch("DetectedParticleE",&DetectedParticleE);
      treePi0->Branch("DetectedTheta",&DetectedTheta);
      treePi0->Branch("DeltaE",&DeltaE);
  */
}

PionEnergyCorrection::~PionEnergyCorrection()
{

}

Bool_t	PionEnergyCorrection::Init()
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

Bool_t	PionEnergyCorrection::Start()
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

void	PionEnergyCorrection::ProcessEvent()
{

  if(usePeriodMacro == 1)
    {
      if(GetEventNumber() % period == 0)
	cout << "Events: " << GetEventNumber() << "  Events Accepted: " << nEventsWritten << endl;
    }

  //********************************************************************************************************************************
  Generated.clear();
  Particles.clear();

  if (GetRootinos()->GetNParticles()==1){

    //Truth Info
    TruthPLab = 1000*(GetTruth()->GettruthPlab(1));
    TruthELab = 1000*(GetTruth()->GettruthElab(1));
    //   TruthTheta = atan2((( (GetTruth()->GettruthPlab(1))*1000*(GetTruth()->Getdircos(3)) )^2 + ( (GetTruth()->GettruthPlab(1)) *1000* (GetTruth()->Getdircos(4)) )^2 )^0.5, GetTruth()->GettruthPlab(1)*1000*GetTruth()->Getdircos(5) ); //tan-1 ( (x^2+y^2)/z  )
    Double_t X = ( (GetTruth()->GettruthPlab(1))*1000*(GetTruth()->Getdircos(3)) ) ;
    Double_t X2= X*X;
    Double_t Y =( (GetTruth()->GettruthPlab(1)) *1000* (GetTruth()->Getdircos(4)) ) ;
    Double_t Y2= Y*Y;
    Double_t Z = (GetTruth()->GettruthPlab(1)*1000*GetTruth()->Getdircos(5)); //tan-1 ( (x^2+y^2)/z  )

    Double_t Numerator = sqrt(X2 + Y2)  ;
    TruthTheta = atan2( Numerator,Z ) ;

    frootino = GetRootinos()->Particle(0);
    particleindex=GetRootinos()->GetTrackIndex(0);//Set as -1 in header so will break if has any unexpected behaviour

    DetectedParticleE = frootino.E();
    DetectedTheta = frootino.Theta();
    DetectedClusterE = GetRootinos()->GetClusterEnergy(0);

    DeltaE = TruthELab - DetectedParticleE;




    THSParticle Gen;
    Gen.SetXYZT(X,Y,Z,TruthELab);
    //	Gen.SetVertex();  //Is This needed?
    Gen.SetPDGcode(-211);
    Generated.push_back(Gen);//Generated defined in header, need to add the libs to PionE in Cmake

    THSParticle part;
    part.SetP4(frootino);
    part.SetPDGcode(-211);
    part.SetDetector(GetTracks()->GetDetectors(particleindex ));
    Particles.push_back(part);
    
    //part.Clear();
    //Gen.Clear();

    //treePi0->Fill();
  } //closing 2 photons and 1 rootino.  
  //*************************************************************************************************************************************
   
  treePi0->Fill();
  nEventsWritten++;
} //closing function

void	PionEnergyCorrection::ProcessScalerRead()
{
  // Fill Tagger Scalers
  //FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	PionEnergyCorrection::Write()
{
  treePi0->Write(); 
  treePi0->Reset();
  return 0; 
}
