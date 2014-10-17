#ifndef FitTest_h
#define FitTest_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TRandom.h"

class FitTest : public Selection {

 public:
  FitTest(TString Name_, TString id_);
  virtual ~FitTest();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     hasTag,
	     TagPtMin,
	     TagIso,
//	     numIsoTags,
	     TauPt,
	     TauEta,
	     TauIsIsolated,
	     TauFit,
	     deltaPhi,
	     Charge,
	     EventId,
	     //ZMassmax,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  virtual void  Finish();
 private:
  // Selection Variables


  TRandom *ran;
  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> TauFlightLength;
  std::vector<TH1D> TauFlightLengthTransverse;
  std::vector<TH1D> TauMomentum;
  std::vector<TH1D> TauMomentumTransverse;
  std::vector<TH1D> TauLife;
  std::vector<TH1D> TauLifeTransverse;
  std::vector<TH1D> ResTauFlightLength;
  std::vector<TH1D> ResTauFlightLengthTransverse;
  std::vector<TH1D> ResTrueTauMomentum;
  std::vector<TH1D> ResTauMomentumTransverse;
  std::vector<TH1D> ResTauMomentum;
  std::vector<TH1D> ResTrueTauMomentumTransverse;
  std::vector<TH1D> ResTrueTauLife;
  std::vector<TH1D> ResTrueTauLifeTransverse;
  std::vector<TH1D> ResTauLife;
  std::vector<TH1D> ResTauLifeTransverse;

  std::vector<TH1D> TauMomentumAvg;
  std::vector<TH1D> TauMomentumTransverseAvg;
  std::vector<TH1D> ResTauMomentumAvg;
  std::vector<TH1D> ResTauMomentumTransverseAvg;
  std::vector<TH1D> TauLifeAvg;
  std::vector<TH1D> TauLifeTransverseAvg;
  std::vector<TH1D> ResTauLifeAvg;
  std::vector<TH1D> ResTauLifeTransverseAvg;




  std::vector<TH1D> ZMass0;
  std::vector<TH1D> ZMass1;
  std::vector<TH1D> ZMass2;

  std::vector<TH1D> TauA1PtResol0;
  std::vector<TH1D> TauA1PtResol1;
  std::vector<TH1D> TauA1PtResol2;

  std::vector<TH1D> TauMuPtResol0;
  std::vector<TH1D> TauMuPtResol1;
  std::vector<TH1D> TauMuPtResol2;


  std::vector<TH1D> TauA1PhiResol0;
  std::vector<TH1D> TauA1PhiResol1;
  std::vector<TH1D> TauA1PhiResol2;

  std::vector<TH1D> TauMuPhiResol0;
  std::vector<TH1D> TauMuPhiResol1;
  std::vector<TH1D> TauMuPhiResol2;


  std::vector<TH1D> TauA1EtaResol0;
  std::vector<TH1D> TauA1EtaResol1;
  std::vector<TH1D> TauA1EtaResol2;

  std::vector<TH1D> TauMuEtaResol0;
  std::vector<TH1D> TauMuEtaResol1;
  std::vector<TH1D> TauMuEtaResol2;

  std::vector<TH1D> TauMuEnergyPull0;
  std::vector<TH1D> TauA1EnergyPull0;

  std::vector<TH1D> TauMuEnergyPull1;
  std::vector<TH1D> TauA1EnergyPull1;

  std::vector<TH1D> TauMuEnergyPull2;
  std::vector<TH1D> TauA1EnergyPull2;

  std::vector<TH1D> Probability0;
  std::vector<TH1D> Probability1;
  std::vector<TH1D> Probability2;

  std::vector<TH1D> ProbabilityOfCorrect;
  std::vector<TH2D> ProbabilityOfCorrectZPt;
  std::vector<TH2D> ProbabilityOfCorrectPtBalance;


  std::vector<TH1D> Charge_TPF_Plus;
  std::vector<TH1D> Charge_TPF_Minus;

  std::vector<TH1D> Charge_GEF_Plus;
  std::vector<TH1D> Charge_GEF_Minus;

  std::vector<TH1D> TruthA1Pt;
  std::vector<TH1D> TruthA1Phi;
  std::vector<TH1D> TruthA1Eta;

  std::vector<TH1D> TruthMuPt;
  std::vector<TH1D> TruthMuPhi;
  std::vector<TH1D> TruthMuEta;

  std::vector<TH1D> TruthA1PtAfter;
  std::vector<TH1D> TruthA1PhiAfter;
  std::vector<TH1D> TruthA1EtaAfter;

  std::vector<TH1D> TruthMuPtAfter;
  std::vector<TH1D> TruthMuPhiAfter;
  std::vector<TH1D> TruthMuEtaAfter;

  std::vector<TH1D> TruthA1PtAfterQC;
  std::vector<TH1D> TruthA1PhiAfterQC;
  std::vector<TH1D> TruthA1EtaAfterQC;

  std::vector<TH1D> TruthMuPtAfterQC;
  std::vector<TH1D> TruthMuPhiAfterQC;
  std::vector<TH1D> TruthMuEtaAfterQC;



  std::vector<TH1D> pvsvsignificance;
  std::vector<TH1D> pvMusvsignificance;



  std::vector<TH1D> Niter0;
  std::vector<TH1D> Niter1;
  std::vector<TH1D> Niter2;

  std::vector<TH2D> Chi2Dim;
  std::vector<TH2D> Chi2DimPlusMinus;


  //>>>>>>>.

  std::vector<TH1D> Niter;
  std::vector<TH1D> Probability;

  std::vector<TH1D> TauA1PhiResol;
  std::vector<TH1D> TauMuPhiResol;

  std::vector<TH1D> TauA1EtaResol;
  std::vector<TH1D> TauMuEtaResol;

  std::vector<TH1D> TauA1PtResol;
  std::vector<TH1D> TauMuPtResol;

  std::vector<TH1D> TauA1EResol;
  std::vector<TH1D> TauMuEResol;


  std::vector<TH2D> TauA1PhiResolVsProb;
  std::vector<TH2D> TauMuPhiResolVsProb;
	    
  std::vector<TH2D> TauA1EtaResolVsProb;
  std::vector<TH2D> TauMuEtaResolVsProb;
	    
  std::vector<TH2D> TauA1PtResolVsProb;
  std::vector<TH2D> TauMuPtResolVsProb;
	    
  std::vector<TH2D> TauA1EResolVsProb;
  std::vector<TH2D> TauMuEResolVsProb;
		


  std::vector<TH1D> TauA1PhiResolZero;
  std::vector<TH1D> TauMuPhiResolZero;

  std::vector<TH1D> TauA1EtaResolZero;
  std::vector<TH1D> TauMuEtaResolZero;

  std::vector<TH1D> TauA1PtResolZero;
  std::vector<TH1D> TauMuPtResolZero;

  std::vector<TH1D> TauA1EResolZero;
  std::vector<TH1D> TauMuEResolZero;


  std::vector<TH1D> TauA1PhiZero;
  std::vector<TH1D> TauA1EtaZero;
  std::vector<TH1D> TauA1PtZero;

  std::vector<TH1D> TauMuPhiZero;
  std::vector<TH1D> TauMuEtaZero;
  std::vector<TH1D> TauMuPtZero;

  std::vector<TH2D> AmbiguitySolverHvsSign;
  std::vector<TH1D> AmbiguitySolverH;






  std::vector<TH1D> TauA1PhiResolCorrectAmbiga;
  std::vector<TH1D> TauMuPhiResolCorrectAmbiga;

  std::vector<TH1D> TauA1EtaResolCorrectAmbiga;
  std::vector<TH1D> TauMuEtaResolCorrectAmbiga;

  std::vector<TH1D> TauA1PtResolCorrectAmbiga;
  std::vector<TH1D> TauMuPtResolCorrectAmbiga;

  std::vector<TH1D> TauA1EResolCorrectAmbiga;
  std::vector<TH1D> TauMuEResolCorrectAmbiga;


  std::vector<TH1D> Ambiguity;
  std::vector<TH1D> AmbiguityIsCorrect;


  std::vector<TH1D> RightlyTakenAmbiguity;
  std::vector<TH1D> WronglyTakenAmbiguity;
  std::vector<TH1D> CorrectAmbiguityTruth;
  std::vector<TH1D> AmbiguityEfficiency;


  std::vector<TH2D> TauMuEResolVsSignificance;

  std::vector<TH1D> TauMuEResolVsSignificanceSlices;

  std::vector<TH2D> TauMuEResolVsSignificance2;
  std::vector<TH2D> TauMuEResolVsChi2Probability;

  std::vector<TH2D> csum2Dim;
  std::vector<TH1D> csumZero;
  std::vector<TH1D> csumGF;
  std::vector<TH2D> csum2Prob;

  std::vector<TH1D> TauA1PtResolKFP;
  std::vector<TH1D> TauA1PhiResolKFP;
  std::vector<TH1D> TauA1EtaResolKFP;

  std::vector<TH1D> TauA1PtResolKFM;
  std::vector<TH1D> TauA1PhiResolKFM;
  std::vector<TH1D> TauA1EtaResolKFM;



  std::vector<TH1D> ZPtResol;
  std::vector<TH1D> ZPzResolRelative;
  std::vector<TH1D> ZEResolRelative;


  std::vector<TH1D> ZPtResolZPtCut;
  std::vector<TH1D> ZPzResolRelativeZPtCut;
  std::vector<TH1D> ZEResolRelativeZPtCut;



  std::vector<TH1D> ZMass;

  std::vector<TH1D> deltaPhiSignificance;

  std::vector<TH1D> idEff;
  std::vector<TH1D> idPassHPS;
  std::vector<TH1D> idPassKFMinus;
  std::vector<TH1D> idPassKFPlus;

  std::vector<TH1D> idPassGF;
  std::vector<TH1D> idPassGFZero;

  std::vector<TH1D> A1Mass;

  std::vector<TH1D> NiterQC;
  std::vector<TH1D> ProbabilityQC;

  std::vector<TH1D> TauA1PhiResolQC;
  std::vector<TH1D> TauMuPhiResolQC;

  std::vector<TH1D> TauA1EtaResolQC;
  std::vector<TH1D> TauMuEtaResolQC;

  std::vector<TH1D> TauA1PtResolQC;
  std::vector<TH1D> TauMuPtResolQC;

  std::vector<TH1D> TauA1EResolQC;
  std::vector<TH1D> TauMuEResolQC;

  std::vector<TH1D> TauA1PzResolQC;
  std::vector<TH1D> TauMuPzResolQC;




  std::vector<TH1D> ZPtTruth;
  std::vector<TH2D> ZPtTruthTauMuPtResolution;

  std::vector<TH1D> TauA1PhiResolZPtCut;
  std::vector<TH1D> TauMuPhiResolZPtCut;

  std::vector<TH1D> TauA1EtaResolZPtCut;
  std::vector<TH1D> TauMuEtaResolZPtCut;

  std::vector<TH1D> TauA1PtResolZPtCut;
  std::vector<TH1D> TauMuPtResolZPtCut;

  std::vector<TH1D> TauA1EResolZPtCut;
  std::vector<TH1D> TauMuEResolZPtCut;


  std::vector<TH2D> MissingEnergyCheckCorrectAmbiga;
  std::vector<TH2D> MissingEnergyCheckInCorrectAmbiga;


  //<<<<<<<.




  std::vector<TH2D> EtaTauMuEtaTauA1; 

  std::vector<TH2D> ResTauMuPtEtaTauMu;
  std::vector<TH2D> ResTauMuEEtaTauMu;
  std::vector<TH2D> ResTauMuPtPtTauMu;

  std::vector<TH2D> ResTauA1PtEtaTauA1;
  std::vector<TH2D> ResTauA1PtPtTauA1;
  std::vector<TH2D> ResTauA1EEtaTauA1;

  std::vector<TH2D> ResTauMuPtProb;
  std::vector<TH2D> ResTauMuEProb;
  std::vector<TH2D> ResTauA1PtProb;
  std::vector<TH2D> ResTauA1EProb;



  std::vector<TH1D> ResTauMuPtEtaTauMuProfile;
  std::vector<TH1D> ResTauMuEEtaTauMuProfile;
  std::vector<TH1D> ResTauMuPtPtTauMuProfile;

  std::vector<TH1D> ResTauA1PtEtaTauA1Profile;
  std::vector<TH1D> ResTauA1PtPtTauA1Profile;
  std::vector<TH1D> ResTauA1EEtaTauA1Profile;

  std::vector<TH1D> ResTauMuPtProbProfile;
  std::vector<TH1D> ResTauMuEProbProfile;
  std::vector<TH1D> ResTauA1PtProbProfile;
  std::vector<TH1D> ResTauA1EProbProfile;


  std::vector<TH1D> PhotonEnergyFraction1;
  std::vector<TH1D> PhotonEnergyFraction2;


  std::vector<TH1D> deltaAmbigaEnergyForWrongAmb;
  std::vector<TH1D> deltaAmbigaEnergyForCorreAmb;



  std::vector<TH1D> KFTauPhi_Minus;
  std::vector<TH1D> KFTauPhi_Plus;
  std::vector<TH1D> KFTauPhi_Ambiga;
  std::vector<TH1D> MeasuredPhi;
  std::vector<TH1D> KFTauPhi_Ambiga_DeltaPhiSignificanceCut;
  std::vector<TH2D> KFTauPhi_Ambiga_vsDeltaPhiSignificance;
  std::vector<TH2D> KFTauPhi_Ambiga_vsDeltaPhi;

  std::vector<TH1D> PhiRotationSignificancePhysicalTaus;
  std::vector<TH1D> PhiRotationSignificanceUnPhysicalTaus;



  int channel;
  double jeteta,muoneta,TauTrackPtThreshold;

};
#endif
