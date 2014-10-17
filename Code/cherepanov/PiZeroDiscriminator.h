#ifndef PiZeroDiscriminator_h
#define PiZeroDiscriminator_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TRandom.h"

class PiZeroDiscriminator : public Selection {

 public:
  PiZeroDiscriminator(TString Name_, TString id_);
  virtual ~PiZeroDiscriminator();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     hasMuon,
	     hasTau,
	     hasValidPhoton,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  //  virtual void  Finish();
 private:
  // Selection Variables



  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;

  std::vector<TH1D> A1Mass;
  std::vector<TH1D> NPiZero;
  std::vector<TH1D> NGamma;
  std::vector<TH1D> NPhotons;
  std::vector<TH1D> PFTauNGamma;

  std::vector<TH1D> PiZeroMass;

  std::vector<TH1D> Pi0EnergyResolution;
  std::vector<TH1D> Pi0DirectionResolution;



  std::vector<TH1D> Pi0EnergyResolutionPair;
  std::vector<TH1D> Pi0DirectionResolutionPair;

  std::vector<TH1D> Pi0EnergyResolutionSingle;
  std::vector<TH1D> Pi0DirectionResolutionSingle;


  std::vector<TH1D> Pi0EnergyResolutionPairID;
  std::vector<TH1D> Pi0DirectionResolutionPairID;

  std::vector<TH1D> Pi0EnergyResolutionSingleID;
  std::vector<TH1D> Pi0DirectionResolutionSingleID;

  std::vector<TH1D> LeadingGammaEnergyResolution;
  std::vector<TH1D> LeadingGammaDirectionResolution;
  std::vector<TH2D> dRVsLeadingGammaEnergy;



  std::vector<TH1D> Pi0Energy;


  std::vector<TH1D> PairSingle;

  std::vector<TH1D> EnergySpectrSingle;
  std::vector<TH1D> EnergySpectrPair;
  std::vector<TH1D> LeadingGammaEnergy;

  std::vector<TH2D> dRVsPi0Energy;

  std::vector<TH2D> dRVsPi0EnergyPair;
  std::vector<TH2D> dRVsPi0EnergySingle;

  std::vector<TH2D> dRVsPi0EnergyPairCutLoose;
  std::vector<TH2D> dRVsPi0EnergySingleCutLoose;

  std::vector<TH2D> dRVsPi0EnergyPairCutTight;
  std::vector<TH2D> dRVsPi0EnergySingleCutTight;

  int channel;
  double jeteta,muoneta,TauTrackPtThreshold;

};
#endif
