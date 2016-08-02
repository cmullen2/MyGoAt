#include "My_PPi0.h"

My_PPi0::My_PPi0()
{ 
   // Create Tree
  tree = new TTree("Tree","Event selection tree");
  
  // Define branches of TTree
  tree->Branch("time",&timev);                   // creates time leaf under the tree directory
  tree->Branch("timeTagg",&timeT);                   // creates time leaf under the tree directory
  tree->Branch("timeP",&timep);                   // creates time leaf under the tree directory
  tree->Branch("Beam_Energy",&energy_beam);
  tree->Branch("CMS_Energy",&Energy_CMS);
  tree->Branch("inv_M",&inv_M_value);            // creates invariant mass leaf under the tree directory
  tree->Branch("MissingM",&MissingM);       
  tree->Branch("MM",&missingp4);                 // creates missing mass leaf under the tree directory
  tree->Branch("gamma_phi",&phi1);               // creates 2-photon phi leaf under the tree directory
  tree->Branch("Rootino_phi",&phi2);             // creates proton phi leaf under the tree directory
  tree->Branch("Coplanarity",&coplanarity);
  tree->Branch("cosTheta",&cosine_CMS);
  tree->Branch("gamma_Theta",&theta1);
  tree->Branch("Rootino_Theta",&theta2);
  tree->Branch("Theta_Difference",&thetadiff);
  tree->Branch("Lin_pol_tagged",&pb);
  tree->Branch("Lin_pol",&pb2);
  tree->Branch("Target_pol",&pt);
  tree->Branch("Edge_plane",&edgeplane);
  tree->Branch("Edge_setting",&edge_setting); // General edge setting
  tree->Branch("Edge_position",&edge_position); // Precise edge position
  
  //  tree->Print();
  
}

My_PPi0::~My_PPi0()
{
}

Bool_t	My_PPi0::Init()
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

	target.SetXYZM(0.0,0.0,0.0,938.0);   //Is this where we decide on the target? Is it not contained in any of the trees before this???
	if(!InitBackgroundCuts()) return kFALSE;
	if(!InitTargetMass()) return kFALSE;
	if(!InitTaggerChannelCuts()) return kFALSE;
	if(!InitTaggerScalers()) return kFALSE;
	cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	My_PPi0::Start()
{
  pt = My_PPi0::targetpol(inputFile);
  if(GetLinpol()->GetPolarizationPlane()==0) edgeplane = 0; // Para
  if(GetLinpol()->GetPolarizationPlane()==1) edgeplane = 1; // Perp
  
  cout << "File polarisation plane = " << GetLinpol()->GetPolarizationPlane() << endl;
  cout << "File edge position = " << edge_setting << endl;
  
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();

    return kTRUE;
}

void	My_PPi0::ProcessEvent()
{
  
  
  // Loop over tagged events
  for(Int_t i=0;i<=GetTagger()->GetNTagged();i++)
    {
      
      
      pb = GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(i));
      //      if(pb>0) cout << "Lin pol = " << pb << endl;
      if(pb == -1) continue;
      
      edge_position = GetLinpol()->GetEdge();
      edge_setting = GetLinpol()->GetEdgeSetting();
      
      
      // Select number and type of events. Event selection loop.
      if((GetRootinos()->GetNParticles()==1 && GetPhotons()->GetNParticles()==2))               // REMEMBER TO TURN ON ROOTINO CALLS WITHIN LOOP   -> 2 photons & 1 or fewer rootinos
	{
	  //if(GetPhotons()->GetNParticles()==2)                                                  // REMEMBER TO TURN OFF ROOTINO CALLS WITHIN LOOP  -> 2 photons, no rootinos
	  //{
	  pb2 = GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(i)); 
	  
	  // Timing calculation
	  timev = GetTagger()->GetTaggedTime(i) - 0.5*(GetPhotons()->GetTime(0) + GetPhotons()->GetTime(1));
	  timeT = GetTagger()->GetTaggedTime(i);
	  timep = 0.5*(GetPhotons()->GetTime(0) + GetPhotons()->GetTime(1)); // Check on photon time distribution
	  //	  if (i%10000==0) cout << timev << endl; // checks to see if timev is actually getting values.
	  
	  
	  // Invariant Mass
	  inv_M_value = (GetPhotons()->Particle(0) + GetPhotons()->Particle(1)).M();
	  energy_beam = GetTagger()->GetTaggedEnergy(i);
	  
	  
	  // Missing Mass calculation (in MeV)
	  // TLorentzVector target(0.0,0.0,0.0,938.0);
	  //beam.SetXYZM(0.,0.,GetTagger()->GetTaggedEnergy(i),GetTagger()->GetTaggedEnergy(i));
	  beam.SetXYZM(0.,0.,GetTagger()->GetTaggedEnergy(i),0.);
	  missingp4     = beam + target - (GetPhotons()->Particle(0) + GetPhotons()->Particle(1));
	  MissingM  = missingp4.M();
	  
	  // Angular distribution of particles
	  phi1 = ( (GetPhotons()->Particle(0) + GetPhotons()->Particle(1)).Phi()* TMath::RadToDeg() );
	  phi2 = (GetRootinos()->Particle(0) ).Phi()* TMath::RadToDeg();
	  
	  //theta1 = (GetPhotons()->Particle(0) + GetPhotons()->Particle(1)).Theta()* TMath::RadToDeg();
	  theta1 = (missingp4).Theta()* TMath::RadToDeg();
	  theta2 = (GetRootinos()->Particle(0) ).Theta()* TMath::RadToDeg();
	  thetadiff = (theta1 - theta2);
	 	  
	  
	  // Coplanarity of rootino and 2 photons
	  coplanarity = (phi1 - phi2);
	  if(coplanarity<0) coplanarity = coplanarity + 360; // Corrects values to one peak around 180 degrees
	  
	  
	  // Boost system //
	  boost_4vec = CMVector((GetPhotons()->Particle(0) + GetPhotons()->Particle(1)), target, beam);
	  cosine_CMS = boost_4vec.CosTheta();
	  Energy_CMS = boost_4vec.Energy();
	  
	  tree->Fill();
	  
	}
    }
}

void	My_PPi0::ProcessScalerRead()
{
  //  Fill Tagger Scalers
  //  FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}
 
Bool_t	My_PPi0::Write()
{
  
  
  tree->Write();
  tree->Reset();
  return 0;
  //  return GTreeManager::Write();
  
}
