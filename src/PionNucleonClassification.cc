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
  //cout<<" tree ? "<<treePi0<<endl;
  //treePi0->SetDirectory(gDirectory);
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

  part.Clear();
  Gen.Clear();
  if (GetRootinos()->GetNParticles()==1){ //Need to know each particle!!!!

    //Truth Info
    TruthPLab = 1000*(GetTruth()->GettruthPlab(1));
    TruthELab = 1000*(GetTruth()->GettruthElab(1));
    //   TruthTheta = atan2((( (GetTruth()->GettruthPlab(1))*1000*(GetTruth()->Getdircos(3)) )^2 + ( (GetTruth()->GettruthPlab(1)) *1000* (GetTruth()->Getdircos(4)) )^2 )^0.5, GetTruth()->GettruthPlab(1)*1000*GetTruth()->Getdircos(5) ); //tan-1 ( (x^2+y^2)/z  )
    X = ( (GetTruth()->GettruthPlab(1))*1000*(GetTruth()->Getdircos(3)) ) ;
    X2= X*X;
    Y =( (GetTruth()->GettruthPlab(1)) *1000* (GetTruth()->Getdircos(4)) ) ;
    Y2= Y*Y;
    Z = (GetTruth()->GettruthPlab(1)*1000*GetTruth()->Getdircos(5)); //tan-1 ( (x^2+y^2)/z  )

    Numerator = sqrt(X2 + Y2)  ;
    TruthTheta = atan2( Numerator,Z ) ;

    frootino = GetRootinos()->Particle(0);
    particleindex=GetRootinos()->GetTrackIndex(0);//Set as -1 in header so will break if has any unexpected behaviour

    DetectedParticleE = frootino.E();
    DetectedTheta = frootino.Theta();
    DetectedClusterE = GetRootinos()->GetClusterEnergy(0);
    DeltaE = TruthELab - DetectedParticleE;

    Gen.SetXYZT(X,Y,Z,TruthELab);
    //	Gen.SetVertex();  //Is This needed?
    Gen.SetPDGcode(-211);
    Generated.push_back(Gen);
 
    part.SetP4(frootino);
    part.SetPDGcode(-211);
    part.SetDetector(GetTracks()->GetDetectors(particleindex ));
    Particles.push_back(part);
  } //closing 2 photons and 1 rootino.  
  
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
//cout << treePi0->GetDirectory()->GetName() << endl;
//treePi0->GetDirectory()->Print(); 
 treePi0->Write(); 
  //treePi0->FlushBaskets();
  treePi0->Reset();
 /* delete treePi0;
 treePi0=nullptr;
  treePi0 =  new TTree("HSParticles","Event selection tree"); //Proton/Neutron pi0 final state tree   
  treePi0->Branch("Particles",&Particles);
  treePi0->Branch("Generated",&Generated);
  treePi0->Branch("EventInfo",&fEventInfo);*/
  return 0; 
}
