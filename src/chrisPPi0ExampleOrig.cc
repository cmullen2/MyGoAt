#include "chrisPPi0Example.h"

chrisPPi0Example::chrisPPi0Example()
{ 
    time 	= new GH1("time", 	"time", 	1400, -700, 700);
    time_cut 	= new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

    time_2g 	= new GH1("time_2g",	"time_2g", 	1400, -700, 700);
    time_2g_cut = new GH1("time_2g_cut","time_2g_cut", 	1400, -700, 700);

    IM 		= new GH1("IM", 	"IM", 		400,   0, 400);
    IM_2g 	= new GH1("IM_2g", 	"IM_2g", 	400,   0, 400);
  
    MM		= new GH1("MM", 	"MM", 	 	400,   800, 1200);     
    MM_2g	= new GH1("MM_2g", 	"MM_2g", 	400,   800, 1200);

    TaggerAccScal = new TH1D("TaggerAccScal","TaggerAccScal",352,0,352);
// Edit 5/4/16 Chris Mullen need to add a process event and corresponding instruction in PPhysics.cc see PPhysics::CalcMissingMass (there already is a 4mom calced for all these just propagate it through print
    Px = new GH1("Px", 	"Px",	1400,	-700, 700);
    Py = new GH1("Py", 	"Py",	1400,	-700, 700);
    Pz = new GH1("Pz", 	"Pz",	1400,	-200, 1200);
    E = new GH1("E", 	"E",	800,	600, 1400);



//CAM 17/5/16 Adding wire chamber hits
    nChamberHitsin1 = new GH1("nChamberHitsin1", 	"nChamberHitsin1",	10,	0, 9);
    Chamber1X = new GH1("Chamber1X", 	"Chamber1X",	180,	-90, 90);
    Chamber1Y = new GH1("Chamber1Y", 	"Chamber1Y",	180,	-90, 90);
    Chamber1Z = new GH1("Chamber1Z", 	"Chamber1Z",	1400,	-700, 700);
    nChamberHitsin2 = new GH1("nChamberHitsin2", 	"Chamber2",	18,	0, 17);
    Chamber2X = new GH1("Chamber2X", 	"Chamber2X",	220,	-110, 110);
    Chamber2Y = new GH1("Chamber2Y", 	"Chamber2Y",	220,	-110, 110);
    Chamber2Z = new GH1("Chamber2Z", 	"Chamber2Z",	680,	-340, 340);


//CAM 18/5/16 Adding reconstructing angles 


    ProtonMM = new GH1("ReconMMusingPion0", 	"Recon MM using Pion0",	400,	800, 1200);
    ProtonX =  new GH1("ReconXusingPion0", 	"Recon X using Pion0",	1200,	-600, 610);
    ProtonY = new GH1("ReconYusingPion0", 	"Recon Y using Pion0",	1200,	-600, 610);
    ProtonZ = new GH1("ReconZusingPion0", 	"Recon Z using Pion0",	1200,	-600, 610);
    ProtonPhi = new GH1("ReconPhiusingPion0", 	"Recon Phi using Pion0",	400,	-200, 200);
    ProtonTheta = new GH1("ReconThetausingPion0", 	"Recon Theta using Pion0",	400,	-200, 200);
    ProtonAngleDiffPhi = new GH1("ProtonAngleDiffPhi", 	"Proton Angle vs Recon Proton Angle (Phi)",	12,	-6, 6);
     


    ProtonNormalX = new GH1("ProtonNormal4vecX", 	"Proton Normal X",	2000,	0, 000);
    ProtonNormalY = new GH1("ProtonNormal4vecY", 	"Proton Normal Y",	2000,	0, 2000);
    ProtonNormalZ = new GH1("ProtonNormal4vecZ", 	"Proton Normal Z",	2000,	0, 2000);
    ProtonNormalMass = new GH1("ProtonNormalMass", 	"Proton Normal Mass",	2000,	0, 2000);
    ProtonNormalPhi = new GH1("ProtonNormalPhi", 	"Proton Normal Phi",	400,	-200, 200);
    ProtonNormalTheta = new GH1("ProtonNormalTheta", 	"Proton Normal Theta",	400,	-200, 200);
    ProtonNormalAngleDiffTheta = new GH1("ProtonAngleDiffTheta", 	"Proton Angle vs Recon Proton Angle (Theta) ",	12,	-6, 6);

    CylinderChamberxvsy = new TGraph();

    //Py = new GH1("Py", 	"Py",	1400,	-700, 700);
    //Pz = new GH1("Pz", 	"Pz",	1400,	-200, 1200);
    //E = new GH1("E", 	"E",	800,	600, 1400);


   // Create Tree(Different method than the function based one for filling above)
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

chrisPPi0Example::~chrisPPi0Example()
{
}

Bool_t	chrisPPi0Example::Init()
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

	if(!InitBackgroundCuts()) return kFALSE;
	if(!InitTargetMass()) return kFALSE;
	if(!InitTaggerChannelCuts()) return kFALSE;
	if(!InitTaggerScalers()) return kFALSE;
	cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	chrisPPi0Example::Start()
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

void	chrisPPi0Example::ProcessEvent()
{
	// fill time diff (tagger - pi0), all pi0
    FillTime(*GetNeutralPions(),time);
    FillTimeCut(*GetNeutralPions(),time_cut);
	
	// fill missing mass, all pi0
    FillMissingMass(*GetNeutralPions(),MM);
	
	// fill invariant mass, all pi0
    FillMass(*GetNeutralPions(),IM);
	//Chris Mullen CAM 6/4/16 Fill Px,Py,Pz and E of all pi0
    FillE(*GetNeutralPions(),E);		//could put this into FillP and use option 4 for it? 
    FillP(*GetNeutralPions(),Px, 1);
    FillP(*GetNeutralPions(),Py, 2);
    FillP(*GetNeutralPions(),Pz, 3);


	//CAM 17/5/16 Fill wire chamber info.

//	For FillChamber the extra parameter is an integer for the axis of the wc 1=x 2=y 3=z. The first param is the chamber number.
//      For FillnChamber the extra parameter is an integer for the Chamber number 
  //Should have just used a for loop or 2 for this
    FillChamber(*GetNeutralPions(),*GetMWPCHitsChris(),Chamber1X,1,1);
    FillChamber(*GetNeutralPions(),*GetMWPCHitsChris(),Chamber1Y,1,2);
    FillChamber(*GetNeutralPions(),*GetMWPCHitsChris(),Chamber1Z,1,3);
    FillnChamber(*GetMWPCHitsChris(),nChamberHitsin1,1);
    FillChamber(*GetNeutralPions(),*GetMWPCHitsChris(),Chamber2X,2,1);
    FillChamber(*GetNeutralPions(),*GetMWPCHitsChris(),Chamber2Y,2,2);
    FillChamber(*GetNeutralPions(),*GetMWPCHitsChris(),Chamber2Z,2,3);
    FillnChamber(*GetMWPCHitsChris(),nChamberHitsin2,2);

//For cylinder chamber first number is chamber 1 or 2 and second number is plot number i.e. 1=x vs y 2 = x vs z and 3 = y vs z
    CylindricalChamber(*GetNeutralPions(), *GetMWPCHitsChris(),CylinderChamberxvsy,1,1);
    //CylindricalChamber(*GetNeutralPions(), *GetMWPCHitsChris(),CylinderChamberxvsz,1,2);
    //CylindricalChamber(*GetNeutralPions(), *GetMWPCHitsChris(),CylinderChamberyvsz,1,3);
    //CylindricalChamber(*GetNeutralPions(), *GetMWPCHitsChris(),CylinderChamberxvsz,2,1);
    //CylindricalChamber(*GetNeutralPions(), *GetMWPCHitsChris(),CylinderChamberxvsz,2,2);
    //CylindricalChamber(*GetNeutralPions(), *GetMWPCHitsChris(),CylinderChamberxvsz,2,3);
    //CylinderChamberxvsy->Write("MyGraph");
    

	//CAM 18/5/16 Angle reconstruction
	//The Proton from missing P4 i.e. the reconstructed one (1,2,3 are the elements of the 4vector X,y,z
    FillProton(*GetNeutralPions(),*GetProtons(),ProtonX, 1); 
    FillProton(*GetNeutralPions(),*GetProtons(),ProtonY, 2);
    FillProton(*GetNeutralPions(),*GetProtons(),ProtonZ, 3);

    FillProtonMM(*GetNeutralPions(),*GetProtons(),ProtonMM);
    FillProtonPhi(*GetNeutralPions(),*GetProtons(),ProtonPhi);
    FillProtonTheta(*GetNeutralPions(),*GetProtons(),ProtonTheta);
    FillProtonAngleDiffPhi(*GetNeutralPions(),*GetProtons(),ProtonAngleDiffPhi);

	//The Proton from identified protons by goat (as it stands but soon will be my own id procedure).
    FillProtonNormal(*GetProtons(),*GetNeutralPions(),ProtonNormalX, 1);
    FillProtonNormal(*GetProtons(),*GetNeutralPions(),ProtonNormalY, 2);
    FillProtonNormal(*GetProtons(),*GetNeutralPions(),ProtonNormalZ, 3);   

    FillProtonNormalMass(*GetProtons(),*GetNeutralPions(),ProtonNormalMass); 
    FillProtonNormalPhi(*GetProtons(),*GetNeutralPions(),ProtonNormalPhi);
    FillProtonNormalTheta(*GetProtons(),*GetNeutralPions(),ProtonNormalTheta);
    FillProtonNormalAngleDiffTheta(*GetProtons(),*GetNeutralPions(),ProtonNormalAngleDiffTheta);
   

		
    // Some neutral decays
    for (Int_t i = 0; i < GetNeutralPions()->GetNParticles(); i++)
    {
        // Fill MM for 2 photon decay
        if ((GetNeutralPions()->GetNSubParticles(i) == 2) & (GetNeutralPions()->GetNSubPhotons(i) == 2))
        {
		// fill time diff (tagger - pi0), this pi0
        FillTime(*GetNeutralPions(),i,time_2g);
        FillTimeCut(*GetNeutralPions(),i,time_2g_cut);
			
		// fill missing mass, this pi0
                FillMissingMass(*GetNeutralPions(),i,MM_2g);
            
		// fill invariant mass, this pi0
            FillMass(*GetNeutralPions(),i,IM_2g);
        }

    }



//  ************************************************************************************************************
  // Loop over tagged events 
  for(Int_t i=0; i<=GetTagger()->GetNTagged() ;i++)
    {
      
/*      
      pb = GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(i));
      //      if(pb>0) cout << "Lin pol = " << pb << endl;
      if(pb == -1) continue;
      
      edge_position = GetLinpol()->GetEdge();
      edge_setting = GetLinpol()->GetEdgeSetting();
*/      

      // Select number and type of events. Event selection loop.
      //if((GetRootinos()->GetNParticles()==1 && GetPhotons()->GetNParticles()==2))  // REMEMBER TO TURN ON ROOTINO CALLS WITHIN LOOP   -> 2 photons & 1 or fewer rootinos


      if((GetRootinos()->GetNParticles()==1 && GetPhotons()->GetNParticles()==2)  ) 
      //if(GetPhotons()->GetNParticles()==3)  
	{
	  //if(GetPhotons()->GetNParticles()==2)                                                  // REMEMBER TO TURN OFF ROOTINO CALLS WITHIN LOOP  -> 2 photons, no rootinos
	  //{
//	  pb2 = GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(i)); 

	  // Timing calculation
	  timev = GetTagger()->GetTaggedTime(i) - 0.5*(GetPhotons()->GetTime(0) + GetPhotons()->GetTime(1));
	  timeT = GetTagger()->GetTaggedTime(i);
	  timep = 0.5*(GetPhotons()->GetTime(0) + GetPhotons()->GetTime(1)); // Check on photon time distribution
	  //	  if (i%10000==0) cout << timev << endl; // checks to see if timev is actually getting values.
	  
	  
	  // Invariant Mass
	  inv_M_value = (GetPhotons()->Particle(0) + GetPhotons()->Particle(1)).M();
	  energy_beam = GetTagger()->GetTaggedEnergy(i);
	  
	  // Cut PID angle, cut inv mass, plot mm then loose cut after looking at mm on mm. Phi plots of pi0 and identify scattered vecs using recons.
	  // Phi asyms change pol and run a diff sim name for it. Read simons paper properly.  
	  // Missing Mass calculation (in MeV)
	  TLorentzVector target(0.0,0.0,0.0,938.0); //this line was commented out in Roddy's one
	  beam.SetXYZM(0.,0.,GetTagger()->GetTaggedEnergy(i),0.);
	  missingp4     = beam + target - (GetPhotons()->Particle(0) + GetPhotons()->Particle(1));
	  MissingM  = missingp4.M();
	  

	  // Angular distribution of particles
	  phi1 = ( (GetPhotons()->Particle(0) + GetPhotons()->Particle(1)).Phi()* TMath::RadToDeg() ); //what would this give me? The pion phi which should be back to back with the proton(rootino phi)
	  phi2 = ((GetRootinos()->Particle(0) ).Phi()* TMath::RadToDeg());
	  

	  if (phi2 == 0)
	  {  
	  	Double_t X2 = ((GetRootinos()->Particle(0) ).X());
	  	Double_t Y2 = ((GetRootinos()->Particle(0) ).Y());
	  	Double_t Z2 = ((GetRootinos()->Particle(0) ).Z());
	  	Double_t E2 = ((GetRootinos()->Particle(0) ).E());
	  	std::cout <<"rootino X " << X2 <<std::endl;
	  	std::cout <<"rootino Y " << Y2 <<std::endl;
	  	std::cout <<"rootino Z " << Z2 <<std::endl;
	  	std::cout <<"rootino E " << E2 <<std::endl;
	  }
	  //theta1 = (GetPhotons()->Particle(0) + GetPhotons()->Particle(1)).Theta()* TMath::RadToDeg();
	  theta1 = (missingp4).Theta()* TMath::RadToDeg();
	  theta2 = (GetRootinos()->Particle(0) ).Theta()* TMath::RadToDeg();
	  thetadiff = (theta1 - theta2);
	 	  
	  
	  // Coplanarity of rootino and 2 photons
	  coplanarity = (phi1 - phi2);
	  if(coplanarity<0) coplanarity = coplanarity + 360; // Corrects values to one peak around 180 degrees
	  
	  
	  // Boost system //
//	  boost_4vec = CMVector((GetPhotons()->Particle(0) + GetPhotons()->Particle(1)), target, beam);
//	  cosine_CMS = boost_4vec.CosTheta();
//	  Energy_CMS = boost_4vec.Energy();
	  
	  tree->Fill();



}

       else if ( (GetPhotons()->GetNParticles()==3))
	{
		
	  // Timing calculation
	  timev = GetTagger()->GetTaggedTime(i) - 0.5*(GetPhotons()->GetTime(0) + GetPhotons()->GetTime(1));
	  timeT = GetTagger()->GetTaggedTime(i);
	  timep = 0.5*(GetPhotons()->GetTime(0) + GetPhotons()->GetTime(1)); // Check on photon time distribution
	  //	  if (i%10000==0) cout << timev << endl; // checks to see if timev is actually getting values.
	  
	  
	  // Invariant Mass
	  inv_M_value = (GetPhotons()->Particle(0) + GetPhotons()->Particle(1)).M();
	  energy_beam = GetTagger()->GetTaggedEnergy(i);
	  
	  // Cut PID angle, cut inv mass, plot mm then loose cut after looking at mm on mm. Phi plots of pi0 and identify scattered vecs using recons.
	  // Phi asyms change pol and run a diff sim name for it. Read simons paper properly.  
	  // Missing Mass calculation (in MeV)
	  TLorentzVector target(0.0,0.0,0.0,938.0); //this line was commented out in Roddy's one
	  beam.SetXYZM(0.,0.,GetTagger()->GetTaggedEnergy(i),0.);
	  missingp4     = beam + target - (GetPhotons()->Particle(0) + GetPhotons()->Particle(1));
	  MissingM  = missingp4.M();
	  

	  // Angular distribution of particles
	  phi1 = ( (GetPhotons()->Particle(0) + GetPhotons()->Particle(1)).Phi()* TMath::RadToDeg() );
//	  phi2 = (GetRootinos()->Particle(0) ).Phi()* TMath::RadToDeg();
	  
	  //theta1 = (GetPhotons()->Particle(0) + GetPhotons()->Particle(1)).Theta()* TMath::RadToDeg();
	  theta1 = (missingp4).Theta()* TMath::RadToDeg();
//	  theta2 = (GetRootinos()->Particle(0) ).Theta()* TMath::RadToDeg();
	  thetadiff = (theta1 - theta2);
	 	  
	  
	  // Coplanarity of rootino and 2 photons
//	  coplanarity = (phi1 - phi2);
//	  if(coplanarity<0) coplanarity = coplanarity + 360; // Corrects values to one peak around 180 degrees
	  
	  
	  // Boost system //
//	  boost_4vec = CMVector((GetPhotons()->Particle(0) + GetPhotons()->Particle(1)), target, beam);
//	  cosine_CMS = boost_4vec.CosTheta();
//	  Energy_CMS = boost_4vec.Energy();
	  
	  tree->Fill();



}


}
}


void	chrisPPi0Example::ProcessScalerRead()
{
	// Fill Tagger Scalers
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	chrisPPi0Example::Write()
{

    CylinderChamberxvsy->Write("MyGraph");
    tree->Write();
    tree->Reset();
//    return 0;

    // Write all GH1's and TObjects defined in this class
    return GTreeManager::Write();
}


