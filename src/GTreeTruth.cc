#include "GTreeTruth.h"



GTreeTruth::GTreeTruth(GTreeManager *Manager, const TString& _Name)    :
    GTree(Manager, _Name),
    fNMC(0)
{

for(Int_t i=0; i<GTreeTruth_MAX; i++)
{
//	truthElab[i] 	= 0;
//	truthPlab[i] 	= 0;
//	dircos[i][3] 	= 0; 
}





}

GTreeTruth::~GTreeTruth()
{
}

void    GTreeTruth::SetBranchAdresses()
{
    inputTree->SetBranchAddress("fNMC", &fNMC);
    inputTree->SetBranchAddress("truthElab", &truthElab);
    inputTree->SetBranchAddress("truthPlab", &truthPlab);
    inputTree->SetBranchAddress("truthVertex", &truthVertex);
    inputTree->SetBranchAddress("truthBeam", &truthBeam);
    inputTree->SetBranchAddress("dircos", &dircos);

}

void    GTreeTruth::SetBranches()
{
    outputTree->Branch("fNMC", &fNMC, "fNMC/I");
    outputTree->Branch("truthElab", truthElab, "truthElab[fNMC]/D");
    outputTree->Branch("truthPlab", truthPlab, "truthPlab[fNMC]/D");
    outputTree->Branch("truthVertex", truthVertex, "truthVertex[3]/D");
    outputTree->Branch("truthBeam", truthBeam, "truthBeam[5]/D");
    outputTree->Branch("dircos", dircos, "dircos[fNMC][3]/F");

}


