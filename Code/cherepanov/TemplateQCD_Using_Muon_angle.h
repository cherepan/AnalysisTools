#ifndef TemplateQCD_Using_Muon_angle_h
#define TemplateQCD_Using_Muon_angle_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TauSpinerInterface.h"

#include "Tauola.h"
#include "LHAPDF/LHAPDF.h"
#include "tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"
#include<iostream>
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"


class TemplateQCD_Using_Muon_angle : public Selection {

 public:
  TemplateQCD_Using_Muon_angle(TString Name_, TString id_); 
  virtual ~TemplateQCD_Using_Muon_angle();

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
	     NCuts};
  
 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
 
 private:
  // Selection Variables

  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> TauPt;
  std::vector<TH1D> ZMass; 
  std::vector<TH1D> xmuon;
  std::vector<TH1D> xmuonT;



 


};
#endif
