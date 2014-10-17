#ifndef WJetCorrectionFactor_h
#define WJetCorrectionFactor_h

#include "Selection.h"
#include <vector>
#include "TString.h"


class WJetCorrectionFactor : public Selection {

 public:
  WJetCorrectionFactor(TString Name_, TString id_); 
  virtual ~WJetCorrectionFactor();

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


 std::vector<TH1D> MetSSTauIso;
 std::vector<TH1D> MetOSTauIso;

 std::vector<TH1D> Met1;
 std::vector<TH1D> Met2;
 std::vector<TH1D> Met3;
 std::vector<TH1D> Met4;
 std::vector<TH1D> Met5;
 std::vector<TH1D> Met6;
 std::vector<TH1D> Met7;
 std::vector<TH1D> Met8;
 std::vector<TH1D> Met9;


 std::vector<TH1D> MetSSTauNonIso;
 std::vector<TH1D> MetOSTauNonIso;

 std::vector<TH1D> RatioTauNonIsoOS;
 std::vector<TH1D> RatioTauNonIsoSS;

 std::vector<TH1D> RatioTauIsoSS;
 std::vector<TH1D> RatioTauIsoOS;


 std::vector<TH1D> RatioTauMet2IsoSS;
 std::vector<TH1D> RatioTauMet2IsoOS;

 std::vector<TH1D> RatioTauMet3IsoSS;
 std::vector<TH1D> RatioTauMet3IsoOS;

 std::vector<TH1D> RatioTauMet4IsoSS;
 std::vector<TH1D> RatioTauMet4IsoOS;

 std::vector<TH1D> RatioTauMet5IsoSS;
 std::vector<TH1D> RatioTauMet5IsoOS;


 std::vector<TH1D> RatioTauMet6IsoSS;
 std::vector<TH1D> RatioTauMet6IsoOS;

 std::vector<TH1D> RatioTauMet7IsoSS;
 std::vector<TH1D> RatioTauMet7IsoOS;
 
 std::vector<TH1D> RatioTauMet8IsoSS;
 std::vector<TH1D> RatioTauMet8IsoOS;
 
 std::vector<TH1D> RatioTauMet9IsoSS;
 std::vector<TH1D> RatioTauMet9IsoOS;

 std::vector<TH1D> MetSSTauIsoSig2;
 std::vector<TH1D> MetOSTauIsoSig2;
 std::vector<TH1D> MetSSTauNonIsoSig2;
 std::vector<TH1D> MetOSTauNonIsoSig2;

 std::vector<TH1D> MetSSTauIsoSig1;
 std::vector<TH1D> MetOSTauIsoSig1;
 std::vector<TH1D> MetSSTauNonIsoSig1;
 std::vector<TH1D> MetOSTauNonIsoSig1;

 std::vector<TH1D> MetSSTauIsoSig0;
 std::vector<TH1D> MetOSTauIsoSig0;
 std::vector<TH1D> MetSSTauNonIsoSig0;
 std::vector<TH1D> MetOSTauNonIsoSig0;

 std::vector<TH2D>   FlightLengthMET;
 std::vector<TH2D>   SignigMET;


};
#endif
