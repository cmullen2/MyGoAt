#ifndef __GTreeMWPCHit_h__
#define __GTreeMWPCHit_h__


#include <TLorentzVector.h>
#include "Rtypes.h"
#include "GTree.h"


//#define GTreeTrack_MAX 128 //Not sure about this line?
#define GTreeMWPCHit_MAX 256 //Not sure about this line?

class   GTreeParticle;
class   GTreeMeson;

class  GTreeMWPCHit : public GTree
{
public:
    enum
    {
        DETECTOR_NONE = 0,
        DETECTOR_NaI = 1,
        DETECTOR_PID = 2,
        DETECTOR_MWPC = 4,
        DETECTOR_BaF2 = 8,
        DETECTOR_PbWO4 = 16,
        DETECTOR_Veto = 32,
    };
private:
		//Not sure what all the _MAX stuff is about need to look to see if it has much of an effect elsewhere or if can disregard for MWPC.
		// May need it for the following params. It is needed or doubles otherwise seg faults due to memory issue
    Int_t	nChamberHitsin1;
    Double_t	Chamber1X[GTreeMWPCHit_MAX];
    Double_t 	Chamber1Y[GTreeMWPCHit_MAX];
    Double_t 	Chamber1Z[GTreeMWPCHit_MAX];
    Int_t       nChamberHitsin2;
    Double_t	Chamber2X[GTreeMWPCHit_MAX];
    Double_t	Chamber2Y[GTreeMWPCHit_MAX];
    Double_t	Chamber2Z[GTreeMWPCHit_MAX];



protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeMWPCHit(GTreeManager *Manager, const TString& _Name);
    virtual ~GTreeMWPCHit();

    virtual void    Clear()     {nChamberHitsin1 = 0;}


            Int_t           GetNMWPCHitsChrisChamber1()                        const	{return nChamberHitsin1;} //currently wrong, need one for each chamber.
           Int_t           GetNMWPCHitsChrisChamber2()                        const	{return nChamberHitsin2;} 
    const	    Double_t*	    GetMWPCChamber1X()			  const {return Chamber1X;}
		    Double_t	    GetMWPCChamber1X(const Int_t index)	  const {return Chamber1X[index];}
    const	    Double_t*	    GetMWPCChamber1Y()			  const {return Chamber1Y;}
		    Double_t	    GetMWPCChamber1Y(const Int_t index)	  const {return Chamber1Y[index];}
    const	    Double_t*	    GetMWPCChamber1Z()			  const {return Chamber1Z;}
		    Double_t	    GetMWPCChamber1Z(const Int_t index)	  const {return Chamber1Z[index];}
    const	    Double_t*	    GetMWPCChamber2X()			  const {return Chamber2X;}
		    Double_t	    GetMWPCChamber2X(const Int_t index)	  const {return Chamber2X[index];}
    const	    Double_t*	    GetMWPCChamber2Y()			  const {return Chamber2Y;}
		    Double_t	    GetMWPCChamber2Y(const Int_t index)	  const {return Chamber2Y[index];}
    const	    Double_t*	    GetMWPCChamber2Z()			  const {return Chamber2Z;}
		    Double_t	    GetMWPCChamber2Z(const Int_t index)	  const {return Chamber2Z[index];}


          


    friend  class GTreeParticle;
    friend  class GTreeMeson;
};

#endif
