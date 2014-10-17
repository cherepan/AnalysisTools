#ifndef BasicSelection_h
#define BasicSelection_h

#include "TFile.h"
#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

class BasicSelection : public Selection {
 
 public:
  BasicSelection(TString Name_, TString id_); 
  virtual ~BasicSelection();

  virtual void  Configure();
  double tauPttest;
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
     virtual void  Finish();
     double TrigEff(Double_t m, double m0, double sigma, double alpha, double n, double norm);
     double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);

     double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);


 private:
  // Extra Variables

      std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSS;
      std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonTauIso;

      std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSScaled;
      std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonTauIsoScaled;

      std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits; 
      std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsScaled; 



/*       std::vector<TH1D> TauA1VisiblePtisTightIsolationDBSumPtCorr;  */
/*       std::vector<TH1D> TauA1VisiblePtisMediumIsolationDBSumPtCorr;  */
/*       std::vector<TH1D> TauA1VisiblePtisLooseIsolationDBSumPtCorr;  */
/*       std::vector<TH1D> TauA1VisiblePtisVLooseIsolationDBSumPtCorr;  */
/*       std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits;  */
/*       std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonMuonIso;  */
/*       std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits;  */
/*       std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByLooseIsolationMVA;  */
/*       std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByMediumIsolationMVA;  */
/*       std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByTightIsolationMVA;  */
/*       std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByLooseIsolationMVA2;  */
/*       std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByMediumIsolationMVA2;  */
/*       std::vector<TH1D> TauA1VisiblePtHPSPFTauDiscriminationByTightIsolationMVA2;  */





     std::vector<TH1D> TauA1VisiblePt;
     std::vector<TH1D> TauA1VisiblePtIsoCut;
     std::vector<TH1D> TauA1VisibleE;

     std::vector<TH1D> TauA1VisiblePx;
     std::vector<TH1D> TauA1VisiblePy;
     std::vector<TH1D> TauA1VisiblePz;

     std::vector<TH1D> TauA1VisibleP;
     //     std::vector<TH1D> TauA1VisiblePhi;
     std::vector<TH1D> TauA1VisibleTheta;
     std::vector<TH1D> TauA1VisibleM;



     std::vector<TH1D> TauA1VisiblePtA1LV;
     std::vector<TH1D> TauA1VisiblePtCorrected1;
     std::vector<TH1D> TauA1VisiblePtCorrected2;
     std::vector<TH1D> TauA1VisiblePtCorrected3;


     std::vector<TH1D> hTauJetdR;
     std::vector<TH1D> hPFJet_PUJetID_discr;
     std::vector<TH1D> hPFJet_PUJetID_looseWP;
     std::vector<TH1D> hPFJet_PUJetID_mediumWP;
     std::vector<TH1D> hPFJet_PUJetID_tightWP;


     std::vector<TH1D> TauMuVisiblePt;
     std::vector<TH1D> VisibleInvariantMass;

     std::vector<TH1D> TauA1VisibleEta;
     std::vector<TH1D> TauMuVisibleEta;

     std::vector<TH1D> TauA1VisiblePhi;
     std::vector<TH1D> TauMuVisiblePhi;


     std::vector<TH1D>   MET_hCorrT0pc_et;
     std::vector<TH1D>   MET_hCorrT0pcT1_et;
     std::vector<TH1D>   MET_hCorrT0rtTxy_et;


     std::vector<TH1D>   MET_hCorrT0rtT1Txy_et;
     std::vector<TH1D>   MET_hCorrT0pcTxy_et;
     std::vector<TH1D>   MET_hCorrT1_et;
     std::vector<TH1D>   MET_hCorrMVA_et;
     std::vector<TH1D>   MET_hCorrT0pcT1Txy_et;
     std::vector<TH1D>   MET_hCorrT0pcT1Txy_phi;



     std::vector<TH1D>   MET_hCorrT0pc_phi;
     std::vector<TH1D>   MET_hCorrT0pcT1_phi;
     std::vector<TH1D>   MET_hCorrT0rtTxy_phi;

     std::vector<TH1D>   MET_hCorrT0rtT1Txy_phi;
     std::vector<TH1D>   MET_hCorrT0pcTxy_phi;
     std::vector<TH1D>   MET_hCorrT1_phi;
     std::vector<TH1D>   MET_hCorrMVA_phi;

     std::vector<TH1D>   SSPionsPt;
     std::vector<TH1D>   SSPionsE;


     std::vector<TH1D>   OSPionsPt;
     std::vector<TH1D>   OSPionsE;

     std::vector<TH1D>   AllPionsPt;
     std::vector<TH1D>   AllPionsE;

     std::vector<TH1D>   SumPtOfjets;
     std::vector<TH1D>   LeadingJetPt;

     std::vector<TH1D>   LeadingJetdRMu;
     std::vector<TH1D>   LeadingJetdRA1;


     std::vector<TH1D>   EFitTauA1Pt;
     std::vector<TH1D>   EFitTauA1Phi;
     std::vector<TH1D>   EFitTauA1Eta;

     std::vector<TH1D>   EFitTauMuPt;
     std::vector<TH1D>   EFitTauMuPhi;
     std::vector<TH1D>   EFitTauMuEta;

     std::vector<TH1D>   EFitTauA1PtAmbPoint;
     std::vector<TH1D>   EFitTauA1PhiAmbPoint;
     std::vector<TH1D>   EFitTauA1EtaAmbPoint;

     std::vector<TH1D>   EFitTauMuPtAmbPoint;
     std::vector<TH1D>   EFitTauMuPhiAmbPoint;
     std::vector<TH1D>   EFitTauMuEtaAmbPoint;

     std::vector<TH1D>   PhiRotationSignificancePhysicalTaus;
     std::vector<TH1D>   PhiRotationSignificanceUnPhysicalTaus;


     std::vector<TH1D>   TauA1Pt;
     std::vector<TH1D>   TauA1E;

     std::vector<TH1D>   TauA1Phi;
     std::vector<TH1D>   TauA1Eta;

     std::vector<TH1D>   TauMuPt;
     std::vector<TH1D>   TauMuE;

     std::vector<TH1D>   TauMuPhi;
     std::vector<TH1D>   TauMuEta;

     std::vector<TH1D>   rtMetPhiMuPhi;
     std::vector<TH1D>   pcMetPhiMuPhi;



     std::vector<TH2D>   FlightLengthMET;
     std::vector<TH2D>   SignigMET;



/*   std::vector<TH1D> NGoodVtx; */
/*   std::vector<TH1D> TauPt; */
/*   std::vector<TH1D> xmuon; */
/*   std::vector<TH1D> xmuonT; */

/*      std::vector<TH1D> ZMass;  */
/*      std::vector<TH1D> ZMassAmbPoint;  */
/*      std::vector<TH1D> VisibleMass;  */
      

/*      std::vector<TH1D> TauEnergy;    */
/*      std::vector<TH1D> MuonEnergy;    */
/*      std::vector<TH1D> HPSPhi;    */


/*      std::vector<TH1D> MuonPhi;    */
/*      std::vector<TH1D> MuonEta;    */
/*      std::vector<TH1D> TauPhi;    */
/*      std::vector<TH1D> TauEta;    */
/*      std::vector<TH1D> FlightLenght; */
/*      std::vector<TH1D> deltaPhi;    */
/*      std::vector<TH1D> Pzeta;    */

/*      std::vector<TH1D> a1Mass;    */
/*      std::vector<TH1D> a1MassAlternat;    */
/*      std::vector<TH1D> Efrac;    */
/*      std::vector<TH1D> RecoTauValid;    */


/*      std::vector<TH1D> hGJ0; */
/*      std::vector<TH1D> hGJ1; */
/*      std::vector<TH1D> hGJ2; */

/*      std::vector<TH1D> rhGJ0; */
/*      std::vector<TH1D> rhGJ1; */
/*      std::vector<TH1D> rhGJ2; */
/*      std::vector<TH1D> DeltaTheta_1st; */
/*      std::vector<TH1D> DeltaTheta_2nd; */
/*      std::vector<TH1D> x; */
/*      std::vector<TH1D> xamb; */

/*      std::vector<TH1D> gamma; */
 
/*      std::vector<TH1D> taumupt; */
/*      std::vector<TH1D> tauptamb; */
/*      std::vector<TH1D> HPSTauPhi; */
/*      std::vector<TH1D> VisiblePhi; */

/*      std::vector<TH1D> TauPhiPlus; */
/*      std::vector<TH1D> TauPhiMinus; */
/*      std::vector<TH1D> TauPhiZero; */

/*      std::vector<TH1D> taumuphi; */
/*      std::vector<TH1D> taumutheta; */
/*      std::vector<TH1D> MuTauDeltaR; */
/*      std::vector<TH1D> EMiss; */

/*      std::vector<TH1D> MuonDistance; */


};
#endif
