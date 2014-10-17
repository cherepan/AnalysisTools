#ifndef Skimmer_h
#define Skimmer_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TRandom.h"

class Skimmer : public Selection {

 public:
  Skimmer(TString Name_, TString id_);
  virtual ~Skimmer();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     hasMuon,
	     hasTau,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  //  virtual void  Finish();
 private:
  // Selection Variables



  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;

  int channel;
  double jeteta,muoneta,TauTrackPtThreshold;

};
#endif
