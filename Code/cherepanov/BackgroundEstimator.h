#ifndef BackgroundEstimator_h
#define BackgroundEstimator_h

#include "Selection.h"
#include <vector>
#include "TString.h"


class BackgroundEstimator : public Selection {

 public:
  BackgroundEstimator(TString Name_, TString id_); 
  virtual ~BackgroundEstimator();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     PrimeVtx,
	     hasTau,
	     TauIsMediumIsolated,
	     TauHPSID,
	     hasMuon,
	     MuonIso,
	     MuonID,
	     PUJetID,
	     FlightLenghtSignificance,
	     MET,
	     charge,
	     Mass,
	     NCuts};


 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);

  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);

 private:
  // Selection Variables



  std::vector<TH1D> aqcd;
  std::vector<TH1D> bqcd;
  std::vector<TH1D> c1qcd;
  std::vector<TH1D> d1qcd;
  std::vector<TH1D> c2qcd;
  std::vector<TH1D> d2qcd;


  std::vector<TH1D> wj65os;
  std::vector<TH1D> wj65ss;

  std::vector<TH1D> wj25os;
  std::vector<TH1D> wj25ss;
  std::vector<TH1D> ChargeAll;

  std::vector<TH1D> TauIsolTightEfficiency;

  std::vector<TH2D> ChargeAllTauIso;
  std::vector<TH2D> ChargeAllTauHPS;
  std::vector<TH2D> ChargeAllMuonIso;
  std::vector<TH2D> ChargeAllMuonID;
  std::vector<TH2D> ChargeAllSign;
  std::vector<TH2D> ChargeAllMET;

  std::vector<TH1D> MissMass1;
  std::vector<TH1D> MissMass2;


  std::vector<TH1D> HPSTauDiscriminators;

 


};
#endif
