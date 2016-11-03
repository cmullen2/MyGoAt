#include "GTreeMWPCHit.h"

GTreeMWPCHit::GTreeMWPCHit(GTreeManager *Manager, const TString& _Name)    :
    GTree(Manager,_Name),
    nChamberHitsin1(0)
{
    for(Int_t i=0; i<GTreeMWPCHit_MAX; i++)
    {/*
        clusterEnergy[i] = 0;
        theta[i] = 0;
        phi[i] = 0;
        time[i] = 0;
        clusterSize[i] = 0;
        centralCrystal[i] = -1;
        centralVeto[i] = -1;
        detectors[i] = 0;
        //Charged detector energies
        vetoEnergy[i] = 0;
        MWPC0Energy[i] = 0;
        MWPC1Energy[i] = 0;*/
    }
}

GTreeMWPCHit::~GTreeMWPCHit()
{

}

void    GTreeMWPCHit::SetBranchAdresses()
{ //It appears that ampersands are required here but it is not entirely clear why as they aren't in GTreeTrack
  //  inputTree->SetBranchAddress("MWPCHits",&MWPCHits);
    inputTree->SetBranchAddress("nChamberHitsin1",&nChamberHitsin1);
    inputTree->SetBranchAddress("Chamber1X", &Chamber1X);
    inputTree->SetBranchAddress("Chamber1Y", &Chamber1Y);
    inputTree->SetBranchAddress("Chamber1Z", &Chamber1Z);
    inputTree->SetBranchAddress("nChamberHitsin2",&nChamberHitsin2);
    inputTree->SetBranchAddress("Chamber2X", &Chamber2X);
    inputTree->SetBranchAddress("Chamber2Y", &Chamber2Y);
    inputTree->SetBranchAddress("Chamber2Z", &Chamber2Z);
    //inputTree->SetBranchAddress("Chamber1Hits", &Chamber1Hits);
	
    //inputTree->SetBranchAddress("Chamber1XPosition", &Chamber1XPosition);



}

void    GTreeMWPCHit::SetBranches()
{ //need to change the output branches and associated functions if want any chamber data out.
    

    //outputTree->Branch("MWPCHits",&MWPCHits,"MWPCHits/I");
   // outputTree->Branch("nChamberHitsin3",nChamberHitsin2,"nChamberHitsin3/I");
    outputTree->Branch("nChamberHitsin1",&nChamberHitsin1,"nChamberHitsin1/I");
    outputTree->Branch("Chamber1X",Chamber1X,"Chamber1X[nChamberHitsin1]/D"); 	//CAMChange
    outputTree->Branch("Chamber1Y",Chamber1Y,"Chamber1Y[nChamberHitsin1]/D");	//CAMChange
    outputTree->Branch("Chamber1Z",Chamber1Z,"Chamber1Z[nChamberHitsin1]/D");   //CAMChange
    outputTree->Branch("nChamberHitsin2",&nChamberHitsin2,"nChamberHitsin2/I");
    outputTree->Branch("Chamber2X",Chamber2X,"Chamber2X[nChamberHitsin2]/D");   //CAMChange
    outputTree->Branch("Chamber2Y",Chamber2Y,"Chamber2Y[nChamberHitsin2]/D");   //CAMChange
    outputTree->Branch("Chamber2Z",Chamber2Z,"Chamber2Z[nChamberHitsin2]/D");   //CAMChange

    //outputTree->Branch("Chamber1Hits",&Chamber1Hits,"Chamber1Hits/D");

    //outputTree->Branch("Chamber1XPosition",Chamber1XPosition,"Chamber1XPosition/D");



   


}


