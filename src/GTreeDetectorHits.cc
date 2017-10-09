#include "GTreeDetectorHits.h"



GTreeDetectorHits::GTreeDetectorHits(GTreeManager *Manager)    :
    GTree(Manager, TString("detectorHits")),
    nNaIHits(0),
    nPIDHits(0),
    nMWPCHits(0),
    nBaF2Hits(0),
    nVetoHits(0)
{

for(Int_t i=0; i<GTreeDetectors_MAX; i++)
{
	NaIEnergy[i] 	= 0;
	NaITime[i]   	= 0;
	PIDEnergy[i] 	= 0;
	PIDTime[i]   	= 0;
	BaF2Energy[i]	= 0;
	BaF2Time[i]	= 0;
	VetoEnergy[i]	= 0;
	VetoTime[i]	= 0;
}





}

GTreeDetectorHits::~GTreeDetectorHits()
{
}

void    GTreeDetectorHits::SetBranchAdresses()
{
    inputTree->SetBranchAddress("nNaIHits", &nNaIHits);
    inputTree->SetBranchAddress("NaIHits", NaIHits);
    inputTree->SetBranchAddress("NaICluster", NaICluster);
    inputTree->SetBranchAddress("NaIEnergy", &NaIEnergy);
    inputTree->SetBranchAddress("NaITime", &NaITime);
    inputTree->SetBranchAddress("nPIDHits", &nPIDHits);
    inputTree->SetBranchAddress("PIDHits", PIDHits);
    inputTree->SetBranchAddress("PIDEnergy", &PIDEnergy);
    inputTree->SetBranchAddress("PIDTime", &PIDTime);
    inputTree->SetBranchAddress("nMWPCHits", &nMWPCHits);
    inputTree->SetBranchAddress("MWPCHits", MWPCHits);
    inputTree->SetBranchAddress("nBaF2Hits", &nBaF2Hits);
    inputTree->SetBranchAddress("BaF2Hits", BaF2Hits);
    inputTree->SetBranchAddress("BaF2Cluster", BaF2Cluster);
    inputTree->SetBranchAddress("BaF2Energy", &BaF2Energy);
    inputTree->SetBranchAddress("BaF2Time", &BaF2Time);
    inputTree->SetBranchAddress("nVetoHits", &nVetoHits);
    inputTree->SetBranchAddress("VetoHits", VetoHits);
    inputTree->SetBranchAddress("VetoEnergy", &VetoEnergy);
    inputTree->SetBranchAddress("VetoTime", VetoTime);

}

void    GTreeDetectorHits::SetBranches()
{
    outputTree->Branch("nNaIHits", &nNaIHits, "nNaIHits/I");
    outputTree->Branch("NaIHits", NaIHits, "NaIHits[nNaIHits]/I");
    outputTree->Branch("NaICluster", NaICluster, "NaICluster[nNaIHits]/I");
    outputTree->Branch("NaIEnergy", NaIEnergy, "NaIEnergy[nNaIHits]/D");
    outputTree->Branch("NaITime", NaITime, "NaITime[nNaIHits]/D");
    outputTree->Branch("nPIDHits", &nPIDHits, "nPIDHits/I");
    outputTree->Branch("PIDHits", PIDHits, "PIDHits[nPIDHits]/I");
    outputTree->Branch("PIDEnergy", PIDEnergy, "PIDEnergy[nPIDHits]/D");
    outputTree->Branch("PIDTime", PIDTime, "PIDTime[nPIDHits]/D");
    outputTree->Branch("nMWPCHits", &nMWPCHits, "nMWPCHits/I");
    outputTree->Branch("MWPCHits", MWPCHits, "MWPCHits[nMWPCHits]/I");
    outputTree->Branch("nBaF2Hits", &nBaF2Hits, "nBaF2Hits/I");
    outputTree->Branch("BaF2Hits", BaF2Hits, "BaF2Hits[nBaF2Hits]/I");
    outputTree->Branch("BaF2Cluster", BaF2Cluster, "BaF2Cluster[nBaF2Hits]/I");
    outputTree->Branch("BaF2Energy", BaF2Energy, "BaF2Energy[nBaF2Hits]/D");
    outputTree->Branch("BaF2Time", BaF2Time, "BaF2Time[nBaF2Hits]/D");
    outputTree->Branch("nVetoHits", &nVetoHits, "nVetoHits/I");
    outputTree->Branch("VetoHits", VetoHits, "VetoHits[nVetoHits]/I");
    outputTree->Branch("VetoEnergy", VetoEnergy, "VetoEnergy[nVetoHits]/D");
    outputTree->Branch("VetoTime", VetoTime, "VetoTime[nVetoHits]/D");

}


