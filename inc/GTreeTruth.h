#ifndef __GTreeTruth_h__
#define __GTreeTruth_h__


#include "Rtypes.h"
#include "GTree.h"

#define GTreeTruth_MAX 4096


class GTreeParticle;
class GTreeMeson;


class  GTreeTruth    : public GTree
{
 private:
  Int_t		fNMC;
  Double_t	truthElab[GTreeTruth_MAX];
  Double_t	truthPlab[GTreeTruth_MAX];
  Double_t	truthVertex[GTreeTruth_MAX];
  Double_t	truthBeam[GTreeTruth_MAX];
  //  Float_t	dircos[GTreeTruth_MAX][GTreeTruth_MAX];  //dircos has 3 components for each particle
 Float_t  dircos[GTreeTruth_MAX];

 protected:
  virtual void    SetBranchAdresses();
  virtual void    SetBranches();

 public:
  GTreeTruth(GTreeManager *Manager, const TString& _Name);
  virtual ~GTreeTruth();

  virtual void            Clear()                         {fNMC=0;}

  Int_t		  GetfNMC()              	    const	{return fNMC;}
  const Double_t* GettruthElab()                    const       {return truthElab;}
  Double_t	  GettruthElab(const Int_t index)   const       {return truthElab[index];}  
  const Double_t* GettruthPlab()                    const       {return truthPlab;}
  Double_t	  GettruthPlab(const Int_t index)   const       {return truthPlab[index];}  
  const Double_t* GettruthVertex()                  const       {return truthVertex;}
  Double_t	  GettruthVertex(const Int_t index) const       {return truthVertex[index];}  
  const Double_t* GettruthBeam()                    const       {return truthBeam;}
  Double_t	  GettruthBeam(const Int_t index)   const       {return truthBeam[index];}  

  //  const Float_t*  Getdircos()                    const          {return dircos[3];}
  //  Float_t 	  Getdircos(const Int_t index, const Int_t index2)   const       {return dircos[index][index2];}  
 const Float_t* Getdircos()		const	{return dircos;}
 Float_t	Getdircos(const Int_t index)	const	{return dircos[index];}



friend class GTreeParticle;
friend class GTreeMeson;


}; 


#endif
