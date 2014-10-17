#ifndef AmbuguityAndResolution_h
#define AmbuguityAndResolution_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

class AmbuguityAndResolution : public Selection {
 
 public:
  AmbuguityAndResolution(TString Name_, TString id_); 
  virtual ~AmbuguityAndResolution();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     PrimeVtx,
	     TauIsQuality,
	     TauPtCut,
	     TauIsIsoMVA2Medium,
	     TauIsIsoMVA2Tight,
	     TauIsIsoMVAMedium,
	     TauIsIsoMVATight,
	     TauIsIsoDBSumPtCorrMedium,
	     TauIsIsoDBSumPtCorrTight,
	     TauIsIsoCombThreeHitsMedium,
	     TauIsIsoCombThreeHitsTight,
	     TauIsNotElectron,
	     TauIsNotMuon,
	     MuonPt,
	     MuonisGlob,
	     MuonisPF,
	     MuonIso,
	     MuonNormChi2,
	     MuonHitsInMS,
	     MuonPixelHits,
	     MuonTrackLayer,
	     MuonNumOfStations,
	     MET,
	     charge,
	     TauDecay,
	     NMuon,
	     NTau,
	     NCuts};
  
 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
     virtual void  Finish();
 private:
  // Selection Variables

/*   std::vector<TH1D> NGoodVtx; */
/*   std::vector<TH1D> TauPt; */
/*   std::vector<TH1D> xmuon; */
/*   std::vector<TH1D> xmuonT; */
  std::vector<TH1D> ZMassPlus; 
  std::vector<TH1D> ZMassMins; 
  std::vector<TH1D> ZMassZero; 

  std::vector<TH1D> TauPtPlus; 
  std::vector<TH1D> TauPtMins; 
  std::vector<TH1D> TauPtZero; 

  std::vector<TH1D> PtBalancePlus; 
  std::vector<TH1D> PtBalanceMins; 
  std::vector<TH1D> PtBalanceZero; 


  std::vector<TH1D> ambuguity; 

  std::vector<TH1D> TauPtPlusResolution; 
  std::vector<TH1D> TauPtMinsResolution; 
  std::vector<TH1D> TauPtZeroResolution; 



  std::vector<TH2D> PtBalancePlusMins;
/*   std::vector<TH2D> PtBalancePlusAmbiguity; */
/*   std::vector<TH2D> PtBalanceMinsAmbiguity; */


};
#endif
