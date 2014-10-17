#ifndef TauRecHelper_h
#define TauRecHelper_h


#include <vector>
#include "TLorentzVector.h"
#include "TComplex.h"

class TauRecHelper {

 public:
  TauRecHelper();
  ~TauRecHelper();

  void Initialize(TLorentzVector t, TLorentzVector mu);

  void SetParametersReco(TLorentzVector tau, TLorentzVector mu );
  void SetBoost(TLorentzVector TauA1, TLorentzVector TauMu );
  TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);
  std::vector<TLorentzVector> RecoEvent(TLorentzVector KFitZero,TLorentzVector KFitPlus, TLorentzVector KFitMinus, TLorentzVector Muon);




  double costheta(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3, TLorentzVector Z);
  double costheta1(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3, TLorentzVector Z);
  float cosbeta(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3);
  double cosbeta1(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3);
  std::vector<float> Sin2Cos2Gamma(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3);
  float cospsi(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3, TLorentzVector Z);





  double lambda(double x, double y, double z);
  float Scalar(TLorentzVector p1, TLorentzVector p2);


  //--------------------------- Hadronic current ---------------------
  float WA(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ);
  float WC(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ);
  float WD(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ);
  float WE(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ);
  TComplex  F(float s1, float s2,float QQ);
  float VV1(float SS1 ,float SS2, float SS3, float QQ);
  float VV2(float SS1 ,float SS2, float SS3, float QQ);
  float V1V2(float SS1 ,float SS2, float SS3, float QQ);
  float h0(float SS1 ,float SS2, float SS3, float QQ);
  float h(float SS1 ,float SS2, float SS3, float QQ);
  TComplex  BWa1(float QQ);
  TComplex  BWrho(float QQ);
  TComplex  BWrhoPrime(float QQ);
  float GammaA1(float QQ);
  float gForGammaA1(float QQ);
  float GammaRho(float QQ);
  float  GammaRhoPrime(float QQ);
  TComplex   Conjugate(TComplex a);






  bool isTauRecoValid();
  TLorentzVector Get1TauSolution(){return TauMu1_;}
  TLorentzVector Get2TauSolution(){return TauMu2_;}
  bool BothSidesAreOK(){return Flag_;}




 private:
 
  TLorentzVector KFitTau_;
  TLorentzVector RecoMuon_;
  TLorentzVector Tau1_;
  TLorentzVector Tau2_;

  TLorentzVector TauMu1_;
  TLorentzVector TauMu2_;

  bool Flag_;
};
#endif
