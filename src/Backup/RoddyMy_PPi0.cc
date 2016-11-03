#include "My_PPi0.h"

My_PPi0::My_PPi0()
{ 
   // Create Tree
  tree = new TTree("Tree","Event selection tree");
  
  // Define branches of TTree
  tree->Branch("time",&time);                     // Time distribution for Pi0
  tree->Branch("Beam_Energy",&energy_beam);       // Beam Energy distribution
  tree->Branch("CMS_Energy",&Energy_CMS);         // Boosted Energy distribution
  //tree->Branch("C1_M",&C1_M_value);               // Invariant mass for photon combination C1
  //tree->Branch("C2_M",&C2_M_value);               // Invariant mass for photon combination C2
  //tree->Branch("C3_M",&C3_M_value);               // Invariant mass for photon combination C3
  tree->Branch("inv_M",&inv_M_value);             // Invariant mass 
  //tree->Branch("MissingM_C1",&MissingM_C1);       // Missing mass distribution for photon combination C1
  //tree->Branch("MissingM_C2",&MissingM_C2);       // Missing mass distribution for photon combination C2
  //tree->Branch("MissingM_C3",&MissingM_C3);       // Missing mass distribution for photon combination C3
  tree->Branch("MissingM",&MissingM);             // Missing mass distribution
  //tree->Branch("MM1",&C1missingp4);               // Missing mass for photon combination C1
  //tree->Branch("MM2",&C2missingp4);               // Missing mass for photon combination C2
  //tree->Branch("MM3",&C3missingp4);               // Missing mass for photon combination C3
  tree->Branch("MM",&missingp4);                  // Missing mass
  tree->Branch("gamma_phi",&phi1);                // Pi0 phi
  tree->Branch("MissingM_phi",&phi2);             // Missing Mass phi
  tree->Branch("Coplanarity",&coplanarity);       // Coplanarity between Pi0 and Missing Mass phi
  tree->Branch("CoPl_MMg3",&coplanarity2);       // Coplanarity between Pi0 and Missing Mass phi
  tree->Branch("cosTheta",&cosine_CMS);           // Boosted CosTheta distribution
  tree->Branch("gamma_Theta",&theta1);            // Pi0 theta
  tree->Branch("MissingM_Theta",&theta2);         // Missing Mass theta
  tree->Branch("Theta_Difference",&thetadiff);    // Difference between Pi0 and Missing Mass theta
  tree->Branch("Theta_Diff_MMg3",&thetadiff2);    // Difference between Pi0 and Missing Mass theta
  tree->Branch("Lin_pol",&pb);                    // Linear polarisation distribution
  tree->Branch("Target_pol",&pt);                 // Target Polarisation
  tree->Branch("Edge_plane",&edgeplane);          // Edge plane: para or perp
  tree->Branch("Edge_setting",&edge_setting);     // General edge setting
  tree->Branch("Edge_position",&edge_position);   // Precise edge position
  
  //  tree->Print();
  
}

My_PPi0::~My_PPi0()
{
}

Bool_t	My_PPi0::Init()
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

	std::string config = ReadConfig("Period-Macro");
	if( sscanf(config.c_str(),"%d\n", &period) == 1 ) usePeriodMacro = 1;

	target.SetXYZM(0.0,0.0,0.0,938.0);
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
  
  //   cout << "File polarisation plane = " << GetLinpol()->GetPolarizationPlane() << endl;
  //   cout << "File edge position = " << edge_setting << endl;
  
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
  
  if(usePeriodMacro == 1)
    {
      if(GetEventNumber() % period == 0)
	cout << "Events: " << GetEventNumber() << "  Events Accepted: " << nEventsWritten << endl;
    }
  if(GetPhotons()->GetNParticles()==3) // Three neutral particle analysis (use of PID information not in the latest config files
    {
      
      // Loop over the different combination of photons and set the different vectors for each
      for(Int_t l=0;l<3;l++){
	if(l==0) { l_g1 = 1; l_g2 = 2;
	  g1_Vec1 = GetPhotons()->Particle(l_g1);
	  g2_Vec1 = GetPhotons()->Particle(l_g2);
	  C1 = g1_Vec1 + g2_Vec1;
	  //IMDiff1 = TMath::Abs(C1.M() - 134.97)/13.497;
	  IMDiff1 = TMath::Abs(C1.M() - 134.97);
	  P1 = GetPhotons()->Particle(0);
	  CopDiff1 = TMath::Abs(C1.TLorentzVector::DeltaPhi(-P1)*TMath::RadToDeg() )/180;
	  AngDiff1 = TMath::Abs(C1.TLorentzVector::Angle(P1.Vect())*TMath::RadToDeg() )/180;
	  SumDiff1 = ((IMDiff1)*(IMDiff1)) + ((CopDiff1)*(CopDiff1)) + ((AngDiff1)*(AngDiff1));
	  T1 = 0.5*(GetPhotons()->GetTime(l_g1) + GetPhotons()->GetTime(l_g2));
	}
	if(l==1) { l_g1 = 2; l_g2 = 0;
	  g1_Vec2 = GetPhotons()->Particle(l_g1);
	  g2_Vec2 = GetPhotons()->Particle(l_g2);
	  C2 = g1_Vec2 + g2_Vec2;
	  IMDiff2 = TMath::Abs(C2.M() - 134.97);
	  P2 = GetPhotons()->Particle(1);
	  CopDiff2 = TMath::Abs(C2.TLorentzVector::DeltaPhi(-P2)*TMath::RadToDeg() )/180;
	  AngDiff2 = TMath::Abs(C2.TLorentzVector::Angle(P2.Vect())*TMath::RadToDeg() )/180;
	  SumDiff2 = ((IMDiff2)*(IMDiff2)) + ((CopDiff2)*(CopDiff2)) + ((AngDiff2)*(AngDiff2));
	  T2 = 0.5*(GetPhotons()->GetTime(l_g1) + GetPhotons()->GetTime(l_g2));
	}
	if(l==2) { l_g1 = 0; l_g2 = 1;
	  g1_Vec3 = GetPhotons()->Particle(l_g1);
	  g2_Vec3 = GetPhotons()->Particle(l_g2);
	  C3 = g1_Vec3 + g2_Vec3;
	  IMDiff3 = TMath::Abs(C3.M() - 134.97);
	  P3 = GetPhotons()->Particle(2);
	  CopDiff3 = TMath::Abs(C3.TLorentzVector::DeltaPhi(-P3)*TMath::RadToDeg() )/180;
	  AngDiff3 = TMath::Abs(C3.TLorentzVector::Angle(P3.Vect())*TMath::RadToDeg() )/180;
	  SumDiff3 = ((IMDiff3)*(IMDiff3)) + ((CopDiff3)*(CopDiff3)) + ((AngDiff3)*(AngDiff3));
	  T3 = 0.5*(GetPhotons()->GetTime(l_g1) + GetPhotons()->GetTime(l_g2));
	}
	
      }
      
      
      // Sort through different photon combinations and decide on the best pi0 candidates, also setting the time
 //      if( (TMath::Abs(C1.M() - 134.97) < TMath::Abs(C2.M() - 134.97)) && (TMath::Abs(C1.M() - 134.97) < TMath::Abs(C3.M() - 134.97))){
// 	meson_Vec.SetXYZM(C1.X(),C1.Y(),C1.Z(),C1.M());
// 	T = T1;
// 	P = P1;
//       }
//       else if(TMath::Abs(C2.M() - 134.97) < TMath::Abs(C3.M() - 134.97) && (TMath::Abs(C2.M() - 134.97) < TMath::Abs(C1.M() - 134.97))){
// 	meson_Vec.SetXYZM(C2.X(),C2.Y(),C2.Z(),C2.M());
// 	T = T2;
// 	P = P2;
//       }
//       else if(TMath::Abs(C3.M() - 134.97) < TMath::Abs(C2.M() - 134.97) && (TMath::Abs(C3.M() - 134.97) < TMath::Abs(C1.M() - 134.97))){
// 	meson_Vec.SetXYZM(C3.X(),C3.Y(),C3.Z(),C3.M());
// 	T = T3;
// 	P = P3;
//       }
      if( SumDiff1 < SumDiff2 && SumDiff1 < SumDiff3){
	meson_Vec.SetXYZM(C1.X(),C1.Y(),C1.Z(),C1.M());
	T = T1;
	P = P1;
      }
      else if(SumDiff2 < SumDiff1 && SumDiff2 < SumDiff3){
	meson_Vec.SetXYZM(C2.X(),C2.Y(),C2.Z(),C2.M());
	T = T2;
	P = P2;
      }
      else if(SumDiff3 < SumDiff1 && SumDiff3 < SumDiff2){
	meson_Vec.SetXYZM(C3.X(),C3.Y(),C3.Z(),C3.M());
	T = T3;
	P = P3;
      }
  
      meson.SetXYZM(meson_Vec.X(),meson_Vec.Y(),meson_Vec.Z(),134.97);
      
      // Invariant Mass
      C1_M_value = C1.M();
      C2_M_value = C2.M();
      C3_M_value = C3.M();
      inv_M_value = meson_Vec.M();
      

      // Loop over tagged events
      for(Int_t i=0;i<=GetTagger()->GetNTagged();i++)
	{
	  
	  time = GetTagger()->GetTaggedTime(i) - T;
	  
	  // Restrict tagger range to reduce file size (but keep a meaningful background range)
	  if(time > -50 && time < 50){

	    edge_setting  = GetLinpol()->GetEdgeSetting();
	    energy_beam   = GetTagger()->GetTaggedEnergy(i);
	    edge_position = GetLinpol()->GetEdge();
	    
	    pb = GetLinpol()->GetPolarizationDegree(GetTagger()->GetTaggedChannel(i));
	    if(pb == -1) continue;
	    
	    
	    
	    energy_beam = GetTagger()->GetTaggedEnergy(i);
	    beam.SetXYZM(0.,0.,GetTagger()->GetTaggedEnergy(i),0.);
	    
	    missingp4     = beam + target - meson;
	    
	    // Boost system //
	    boost_4vec  = CMVector(meson, target, beam);
	    boost_MMvec = CMVector(missingp4, target, beam);
	    boost_G3vec = CMVector(P, target, beam);
	    
	    cosine_CMS = boost_4vec.CosTheta();
	    Energy_CMS = boost_4vec.Energy();
	    
	    
	    // Missing Mass calculation (in MeV)
	    C1missingp4  = beam + target - C1;
	    C2missingp4  = beam + target - C2;
	    C3missingp4  = beam + target - C3;
	    
	    MissingM_C1  = C1missingp4.M();
	    MissingM_C2  = C2missingp4.M();
	    MissingM_C3  = C3missingp4.M();
	    MissingM     = boost_MMvec.M();
	    
	    
	    // Angular distribution of particles
	    phi1 = (boost_4vec).Phi()* TMath::RadToDeg();
	    phi2 = boost_MMvec.Phi()*TMath::RadToDeg();
	    phiP = P.Phi()*TMath::RadToDeg();
	    
	    theta1 = boost_4vec.Theta()* TMath::RadToDeg();
	    theta2 = boost_MMvec.Theta()*TMath::RadToDeg();
	    thetaP = boost_G3vec.Theta()*TMath::RadToDeg();
	    
	    //thetadiff  = (theta1 - thetaP); // Theta difference between pion candidate and 3rd gamma	      
	    //thetadiff2 = (theta2 - thetaP); // Theta difference between pion candidate and missing mass vector
	    //thetadiff  = boost_G3vec.TLorentzVector::Angle(boost_4vec.Vect())*TMath::RadToDeg(); // Theta difference between pion candidate and 3rd gamma	      
	    //thetadiff2 = boost_G3vec.TLorentzVector::Angle(boost_MMvec.Vect())*TMath::RadToDeg(); // Theta difference between pion candidate and missing mass vector
	    thetadiff  = P.TLorentzVector::Angle(meson.Vect())*TMath::RadToDeg(); // Theta difference between pion candidate and 3rd gamma	      
	    thetadiff2 = P.TLorentzVector::Angle(missingp4.Vect())*TMath::RadToDeg(); // Theta difference between pion candidate and missing mass vector
	    
	    //coplanarity  = phi1 - phiP; // Phi difference between pion candidate and 3rd gamma	      
	    //coplanarity2 = phi2 - phiP;	// Phi difference between pion candidate and missing mass vector
	    coplanarity  = (boost_G3vec.TLorentzVector::DeltaPhi(boost_4vec))*TMath::RadToDeg(); // Phi difference between pion candidate and 3rd gamma	      
	    coplanarity2 = (boost_G3vec.TLorentzVector::DeltaPhi(boost_MMvec))*TMath::RadToDeg();	// Phi difference between pion candidate and missing mass vector
	    
	    
	    //boost_4vec = CMVector(meson, target, -beam); // May need to unboost meson 4vec
	    //  if(coplanarity < -155 || coplanarity > 155){
	    //if(thetadiff2 > 0 && thetadiff2 < 20){
		tree->Fill();
		//}
		//}
	  }
	}
    }
  nEventsWritten++;
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
