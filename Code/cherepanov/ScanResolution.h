#ifndef ScanResolution_h
#define ScanResolution_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "Tauola.h"
#include "LHAPDF/LHAPDF.h"
#include "tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"
#include<iostream>
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"




class ScanResolution : public Selection {

 public:
  ScanResolution(TString Name_, TString id_);
  virtual ~ScanResolution();

  virtual void  Configure();


  enum cuts {TriggerOk=0, 
	     PrimeVtx,
	     MuonisGlob,
	     TauIsQuality,
	     MuonPt,
	     TauPtCut,
	     MET,
	     MuonIso,
	     charge,
	     TauIsIso,
	     TauDecay,
	     NMuon,
	     NTau,
	     NCuts};



 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH1D> Tau3PiPt;
  std::vector<TH1D> Tau1theta;
  std::vector<TH1D> Tau2theta;  
  std::vector<TH1D> Tau1Deltatheta;
  std::vector<TH1D> Tau2Deltatheta;  
  std::vector<TH1D> DeltaTheta_1st;
  std::vector<TH1D> DeltaTheta_2nd;
  std::vector<TH1D> Tau1Pt;
  std::vector<TH1D> Tau2Pt;
  std::vector<TH1D> Energy_resolution1;
  std::vector<TH1D> Energy_resolution2;
  std::vector<TH1D> TruthResolution;
  std::vector<TH1D> TransverseEnergy_resolution;

  std::vector<TH1D> xmuon_plus;
  std::vector<TH1D> xmuonT_plus;
  std::vector<TH1D> xmuonTruthT_plus;
  std::vector<TH1D> xmuonTruth_plus;
  std::vector<TH1D> xmuonMC_plus;
  std::vector<TH1D> xmuonTMC_plus;

  std::vector<TH1D> res0;
  std::vector<TH1D> res5;
  std::vector<TH1D> res10;
  std::vector<TH1D> res15;
  std::vector<TH1D> res20;
  std::vector<TH1D> res30;
  std::vector<TH1D> res40;
  std::vector<TH1D> res50;
  std::vector<TH1D> res60;
  std::vector<TH1D> res70;
  std::vector<TH1D> res80;
  std::vector<TH1D> res90;
  std::vector<TH1D> res100;

  std::vector<TH1D> res110;
  std::vector<TH1D> res120;
  std::vector<TH1D> res150;

  std::vector<TH1D> xMC;
  std::vector<TH1D> xRec;
  std::vector<TH1D> xTruth;

  std::vector<TH1D> xp0;
  std::vector<TH1D> xp5;
  std::vector<TH1D> xp10;
  std::vector<TH1D> xp15;
  std::vector<TH1D> xp20;
  std::vector<TH1D> xp30;
  std::vector<TH1D> xp40;
  std::vector<TH1D> xp50;
  std::vector<TH1D> xp60;
  std::vector<TH1D> xp70;
  std::vector<TH1D> xp80;
  std::vector<TH1D> xp90;
  std::vector<TH1D> xp100;

  std::vector<TH1D> xp110;
  std::vector<TH1D> xp120;
  std::vector<TH1D> xp150;


  std::vector<TH1D> xm0;
  std::vector<TH1D> xm5;
  std::vector<TH1D> xm10;
  std::vector<TH1D> xm15;
  std::vector<TH1D> xm20;
  std::vector<TH1D> xm30;
  std::vector<TH1D> xm40;
  std::vector<TH1D> xm50;
  std::vector<TH1D> xm60;
  std::vector<TH1D> xm70;
  std::vector<TH1D> xm80;
  std::vector<TH1D> xm90;
  std::vector<TH1D> xm100;

  std::vector<TH1D> xm110;
  std::vector<TH1D> xm120;
  std::vector<TH1D> xm150;


  std::vector<TH2D> IsoVsMuonEnergy;
  std::vector<TH2D> IsoAbsVsMuonEnergy;
  std::vector<TH2D> IsoVsMuonTauEnergy;
  std::vector<TH1D> xmuon_minus;
  std::vector<TH1D> xmuonT_minus;
  std::vector<TH1D> xmuonTruthT_minus;
  std::vector<TH1D> xmuonTruth_minus;
  std::vector<TH1D> xmuonMC_minus;
  std::vector<TH1D> xmuonTMC_minus;   

  
  std::vector<TH1D> TruthDeltaTheta;


  std::vector<TH1D> xReco0;
  std::vector<TH1D> xReco5;
  std::vector<TH1D> xReco10;
  std::vector<TH1D> xReco15;
  std::vector<TH1D> xReco20;
  std::vector<TH1D> xReco30;
  std::vector<TH1D> xReco40;
  std::vector<TH1D> xReco50;
  std::vector<TH1D> xReco60;
  std::vector<TH1D> xReco70;
  std::vector<TH1D> xReco80;
  std::vector<TH1D> xReco90;
  std::vector<TH1D> xReco100;

  std::vector<TH1D> xReco110;
  std::vector<TH1D> xReco120;
  std::vector<TH1D> xReco150;


};
#endif
