#include "TauRecHelper.h"
#include <iostream>


TauRecHelper::TauRecHelper(){
 
}
void 
TauRecHelper::SetParametersReco(TLorentzVector tau, TLorentzVector mu ){
 Initialize(tau,mu);
}
void 
TauRecHelper::SetBoost(TLorentzVector TauA1, TLorentzVector TauMu ){
  Tau1_ = TauA1;
  Tau2_ = TauMu;
}





TauRecHelper::~TauRecHelper(){
}



void 
TauRecHelper::Initialize(TLorentzVector t, TLorentzVector mu){
  RecoMuon_=mu;
  KFitTau_=t;
}


bool TauRecHelper::isTauRecoValid(){

      double zmass = 91.1876,mtau = 1.777;
      double ma1 = 1.23;

      TLorentzVector tau1,tau2,taumu1,taumu2;
      double e1 =0 ,e2 =0,pz1=0,pz2=0; 

      double A =  0.5*zmass*zmass - mtau*mtau - KFitTau_.Pt()* KFitTau_.Pt() ;
      double C = mtau*mtau + KFitTau_.Pt()* KFitTau_.Pt();
      double B =  A*A + C*KFitTau_.Pz()*KFitTau_.Pz();

      double a = floor(1*1000)/1000;
      double b = -floor(2*KFitTau_.E()*A/C*1000)/1000;
      double c = floor(B/C*1000)/1000;


      if(b*b - 4*a*c < 0){
	TauMu1_.SetXYZM(0,0,0,0);
	TauMu2_.SetXYZM(0,0,0,0);
	return false;}

      e1 = (-b + sqrt(b*b - 4*a*c))/2/a;
      e2 = (-b - sqrt(b*b - 4*a*c))/2/a;

      pz1 =  (e1*KFitTau_.E() + KFitTau_.Pt()* KFitTau_.Pt() - 0.5*zmass*zmass  + mtau*mtau)/(KFitTau_.Pz());
      pz2 =  (e2*KFitTau_.E() + KFitTau_.Pt()* KFitTau_.Pt() - 0.5*zmass*zmass  + mtau*mtau)/(KFitTau_.Pz());
    
      tau1.SetXYZM(-KFitTau_.Px(),-KFitTau_.Py(),pz1,mtau); 
      tau2.SetXYZM(-KFitTau_.Px(),-KFitTau_.Py(),pz2,mtau); 

      if(fabs(tau1.Theta() - RecoMuon_.Theta()) < fabs(tau2.Theta() - RecoMuon_.Theta())){taumu1 = tau1; taumu2 = tau2;}
      else{taumu1 = tau2; taumu2 = tau1;}


      TauMu1_ = taumu1;
      TauMu2_ = taumu2;
      return true;
}

std::vector<TLorentzVector> 
TauRecHelper::RecoEvent(TLorentzVector KFitZero,TLorentzVector KFitPlus, TLorentzVector KFitMinus, TLorentzVector Muon){
  std::vector<TLorentzVector> Output;// first a1 Tau, second muon Tau


  TLorentzVector muon    = Muon;
  TLorentzVector BalancedTau;
  if(fabs(KFitMinus.Pt() - muon.Pt()) <= fabs(KFitPlus.Pt() - muon.Pt()) ){ BalancedTau = KFitMinus; }
  else if(fabs(KFitMinus.Pt() - muon.Pt()) > fabs(KFitPlus.Pt() - muon.Pt())  ){ BalancedTau = KFitPlus; }
  
  TLorentzVector tau_3pi = BalancedTau;



  TLorentzVector tau1,tau2,taumu1,taumu2,Z;
  TLorentzVector z1,z2;
  double ma1 = 1.23;

  double zmass = 91.187,mtau = 1.777;
  double e1 =0 ,e2 =0,pz1=0,pz2=0; 
  double Tau3piPt=floor(tau_3pi.Pt()*100000)/100000;
  double Tau3piPz=floor(tau_3pi.Pz()*100000)/100000;
  double Tau3piE =/*floor(sqrt(mtau*mtau + Tau3piPz*Tau3piPz+Tau3piPt*Tau3piPt)*100000)/100000;//*/floor(tau_3pi.E()*100000)/100000;

  zmass = floor(zmass*1000)/1000;
  mtau = floor(mtau*1000)/1000;

  double A =  0.5*zmass*zmass - mtau*mtau - Tau3piPt* Tau3piPt ;
  //float A =  0.5*(zmass*zmass) + tau_3pi.Pt()* tau_3pi.Pt() ;
  //float B =  mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt() + A*A/tau_3pi.Pz()/tau_3pi.Pz();
  
  double C = mtau*mtau + Tau3piPt*Tau3piPt;
  double B =  A*A + C*Tau3piPz*Tau3piPz;


  A = floor(A*1000)/1000;
  B = floor(B*1000)/1000;
  C = floor(C*1000)/1000;
  

  double root = sqrt(A*A*Tau3piE*Tau3piE - C*(A*A + C*Tau3piPz*Tau3piPz));
  root = floor(root*1000)/1000;


  
  


  double a = floor(1*1000)/1000;
  double b = -floor(2*Tau3piE*A/C*1000)/1000;
  double c = floor(B/C*1000)/1000;


  e1 = (-b + sqrt(b*b - 4*a*c))/2/a;
  e2 = (-b - sqrt(b*b - 4*a*c))/2/a;

// alternative pz computation
  pz1 =  (e1*Tau3piE + Tau3piPt*Tau3piPt - 0.5*zmass*zmass  + mtau*mtau)/(Tau3piPz);
  pz2 =  (e2*Tau3piE + Tau3piPt*Tau3piPt - 0.5*zmass*zmass  + mtau*mtau)/(Tau3piPz);


  // alternative pz computation
  
  if(b*b - 4*a*c < 0){
    Flag_ = false;
  }
  if(b*b - 4*a*c >= 0){
    Flag_ = true;
  }

  
//     pz1 =  ((muon.Pz())/(fabs(muon.Pz())))*floor( sqrt(e1*e1 -  (mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt()))*1000)/1000;
//     pz2 =  ((muon.Pz())/(fabs(muon.Pz())))*floor( sqrt(e2*e2 -  (mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt()))*1000)/1000;



    
    tau1.SetXYZM(-tau_3pi.Px(),-tau_3pi.Py(),pz1,mtau); 
    tau2.SetXYZM(-tau_3pi.Px(),-tau_3pi.Py(),pz2,mtau); 

    z1 = tau1+tau_3pi;
    z2 = tau2+tau_3pi;
    
    if(fabs(tau1.Theta() - muon.Theta()) < fabs(tau2.Theta() - muon.Theta())){taumu1 = tau1; taumu2 = tau2;}
    else{taumu1 = tau2; taumu2 = tau1;}


    Output.push_back(tau_3pi);
    Output.push_back(taumu1);
    Output.push_back(taumu2);

  return Output;

}

double 
TauRecHelper::costheta(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3,TLorentzVector Z ){
  

  double zmass = Z.M();
  double mtau = 1.777;
  TLorentzVector a1 = p1+p2+p3;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();

  double x = 2*a1.E()/zmass;
  double ctheta = (2*x*mtau*mtau - mtau*mtau - QQ)/((mtau*mtau - QQ)*sqrt(1 - 4*mtau*mtau/zmass/zmass));
//   std::cout<<"p1 "<< p1.Px() << "  " <<p1.Py() << " "<< p1.Pz() << "  " <<p1.M()<<std::endl;
//   std::cout<<"p2 "<< p2.Px() << "  "<< p2.Py() << " "<< p2.Pz() << "  " <<p2.M()<<std::endl;
//   std::cout<<"p3 "<< p3.Px() << "  "<< p3.Py() << " "<< p3.Pz() << "  " <<p3.M()<<std::endl;
//   std::cout<<"a1 "<< a1.Px() << "  "<< a1.Py() << " "<< a1.Pz() << "  " <<a1.M()<<"  " <<a1.M()*a1.M() <<std::endl;
//   std::cout<<"QQ "<< QQ << " x "<< x <<std::endl;

//   std::cout<<"2*x*mtau*mtau/mtau*mtau - QQ "<< 2*x*mtau*mtau/(mtau*mtau - QQ)<< "  "<< x <<std::endl;
  return ctheta;
}


double 
TauRecHelper::costheta1(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3,TLorentzVector Z ){
  

  double zmass = Z.M();
  double mt = 1.777;
  TLorentzVector a1 = p1+p2+p3;
  double ma  = a1.M();
  double diffmass = mt*mt - ma*ma;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();

  double x = 2*a1.E()/zmass;
  double ctheta = 4*mt*mt*a1.E()/zmass/diffmass   - (mt*mt + ma*ma)/diffmass;
  return ctheta;

}


float 
TauRecHelper::cosbeta(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3){
  
  float mpi  = 0.139;
//   float E = p1.E() +  p2.E() +  p3.E(); 
//   float P = p1.P() +  p2.P() +  p3.P(); 
//   float QQ = E*E - P*P;


//   std::cout<<"  cosbeta --------------------- "<<std::endl;

  TLorentzVector a1 = p1+p2+p3;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();

  float B1 = (pow(p1.E()*a1.E()   - Scalar(p1,a1),2 ) - QQ*mpi*mpi)/QQ;
  float B2 = (pow(p2.E()*a1.E()   - Scalar(p2,a1),2 ) - QQ*mpi*mpi)/QQ;
  float B3 = (pow(p3.E()*a1.E()   - Scalar(p3,a1),2 ) - QQ*mpi*mpi)/QQ;

  float T = 0.5*sqrt(-lambda(B1,B2,B3));

  TLorentzVector p1Timesp2(p1.Py()*p2.Pz() - p1.Pz()*p2.Py(),p1.Pz()*p2.Px() - p1.Px()*p2.Pz(),p1.Px()*p2.Py() - p1.Py()*p2.Px(),1);
    
    
  float cbeta = Scalar(p3,p1Timesp2)/a1.P()/T;

//   std::cout<<"  T  "<< T<<std::endl;
//   std::cout<<"  B1  "<< B1<<std::endl;
//   std::cout<<"  B2  "<< B2<<std::endl;
//   std::cout<<"  B3  "<< B3<<std::endl;
//   std::cout<<"  QQ  "<< QQ<<std::endl;
//   std::cout<<"  Scalar(p3,p1Timesp2)  "<< Scalar(p3,p1Timesp2)<<std::endl;
//   std::cout<<"  cbeta  "<< cbeta  <<std::endl;


  return cbeta;

}


double 
TauRecHelper::cosbeta1(TLorentzVector p2,TLorentzVector p3, TLorentzVector p1){
  //std::cout<<"  cosbeta1  ================================= "<<std::endl;
  double mpi  = 0.139;
//   double E = p1.E() +  p2.E() +  p3.E(); 
//   double P = p1.P() +  p2.P() +  p3.P(); 
//   double QQ = E*E - P*P;

  TLorentzVector a1 = p1+p2+p3;
  //  float P = 
  TLorentzVector s12 = p1+p2;
  TLorentzVector s13 = p1+p3;
  TLorentzVector s23 = p2+p3;

  double QQ = a1.E()*a1.E() - a1.P()*a1.P();


  TLorentzVector p1Timesp2(p1.Py()*p2.Pz() - p1.Pz()*p2.Py(),p1.Pz()*p2.Px() - p1.Px()*p2.Pz(),p1.Px()*p2.Py() - p1.Py()*p2.Px(),1);
  float mm=a1.M()*a1.M();
  float mm12=s12.M()*s12.M();
  float mm13=s13.M()*s13.M();
  float mm23=s23.M()*s23.M();
  float mmpi=mpi*mpi;

  float l1  = lambda( mm, mm12 , mmpi);
  float l2  = lambda( mm, mm13 , mmpi);
  float l3  = lambda( mm, mm23 , mmpi);

 
  double cbeta = /*8*a1.M()*a1.M()**/Scalar(p3,p1Timesp2)*a1.P()/sqrt(-lambda(l1,l2,l3));

//   std::cout<<"  QQ  "<< QQ<<std::endl;
//   std::cout<<"  mm12  "<< mm12<<std::endl;
//   std::cout<<"  mm13  "<< mm13<<std::endl;
//   std::cout<<"  mm23  "<< mm23<<std::endl;
//   std::cout<<"  mm  "<< mm<<std::endl;

//   std::cout<<"  lambda1 = "<<l1<<std::endl;
//   std::cout<<"  lambda2 = "<<l2<<std::endl;
//   std::cout<<"  lambda3 = "<<l3<<std::endl;
//   std::cout<<"  lambda  = "<<-lambda(l1,l2,l3)<<std::endl;

  
//   std::cout<<"  Scalar(p3,p1Timesp2)*a1.P()  "<< Scalar(p3,p1Timesp2)<<std::endl;
//   std::cout<<"  a1.P()  "<< a1.P() <<std::endl;

//   std::cout<<"  a1.P()  "<< a1.P() << "   *  "<</*8*a1.M()*a1.M()**/Scalar(p3,p1Timesp2)/sqrt(-lambda(l1,l2,l3)) <<std::endl;



//   std::cout<<"  cbeta  "<< cbeta  <<std::endl;



  return cbeta;

}




double 
TauRecHelper::lambda(double x, double y, double z){
  
  return x*x +y*y +z*z - 2*x*y - 2*x*z - 2*z*y;

}




TLorentzVector 
TauRecHelper::Boost(TLorentzVector pB, TLorentzVector frame){
  
  


  float zmass = frame.M();
  float g  = sqrt(zmass*zmass + pow(frame.Pz(),2  ) +pow(frame.Px(),2  ) +pow(frame.Py(),2  )  )/zmass;
   // float g  = sqrt(zmass*zmass + pow(frame.Pz(),2  ) /*+pow(Frame1_.Px() + Frame2_.Px(),2  ) +pow(Frame1_.Py() + Frame2_.Py(),2  ) */ )/zmass;
  float b  = sqrt(1 - 1/g/g);


  int signx;
  int signy;
  int signz;


  float bx = frame.Px()/frame.E();
  float by = frame.Py()/frame.E();
  float bz = frame.Pz()/frame.E();


  float E = g*pB.E() - g*bx*pB.Px() - g*by*pB.Py() - g*bz*pB.Pz();
  float Px= -g*bx*pB.E() + (1+ (g-1)*bx*bx/b/b)*pB.Px() + ((g-1)*bx*by/b/b)*pB.Py() + ((g-1)*bx*bz/b/b)*pB.Pz();
  float Py= -g*by*pB.E() + ((g-1)*by*bx/b/b)*pB.Px() + (1 + (g-1)*by*by/b/b)*pB.Py() + ((g-1)*by*bz/b/b)*pB.Pz();
  float Pz= -g*bz*pB.E() + ((g-1)*bz*bx/b/b)*pB.Px() + ((g-1)*bz*by/b/b)*pB.Py()    + (1 + (g-1)*bz*bz/b/b)*pB.Pz();

//   float Pz = -sign*b*g*particleToBoost.E() + g*particleToBoost.Pz();
//   float E  = g*particleToBoost.E() - sign*b*g*particleToBoost.Pz();
//   //  std::cout<< "gamma:  "<< g << " beta:  "<<b<<" zmass " <<tau.M() <<std::endl;
  return TLorentzVector(Px,Py,Pz,E);
}


 


float 
TauRecHelper::Scalar(TLorentzVector p1, TLorentzVector p2){
  
  return p1.Px()*p2.Px() +  p1.Py()*p2.Py() +  p1.Pz()*p2.Pz();
}

std::vector<float> 
TauRecHelper::Sin2Cos2Gamma(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3){

  std::vector<float> sin2cos2;
  float mpi  = 0.139;
  TLorentzVector a1 = p1+p2+p3;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();

  float B1 = (pow(p1.E()*a1.E()   - Scalar(p1,a1),2 ) - QQ*mpi*mpi)/QQ;
  float B2 = (pow(p2.E()*a1.E()   - Scalar(p2,a1),2 ) - QQ*mpi*mpi)/QQ;
  float B3 = (pow(p3.E()*a1.E()   - Scalar(p3,a1),2 ) - QQ*mpi*mpi)/QQ;

  float T = 0.5*sqrt(-lambda(B1,B2,B3));

  float A1=(a1.E()*Scalar(a1,p1) - p1.E()*a1.P()*a1.P())/QQ;
  float A2=(a1.E()*Scalar(a1,p2) - p2.E()*a1.P()*a1.P())/QQ;
  float A3=(a1.E()*Scalar(a1,p3) - p3.E()*a1.P()*a1.P())/QQ;


  float cosgamma = A3/a1.P()/sqrt(B3)/sqrt(1 - cosbeta(p1,p2,p3)*cosbeta(p1,p2,p3));
  float singamma = -cosgamma*(B3*A1/A3 - 0.5*(B2 - B1 - B3))/T;

  sin2cos2.push_back(2*singamma*cosgamma);
  sin2cos2.push_back(2*cosgamma*cosgamma - 1);
  return sin2cos2;

}



float 
TauRecHelper::cospsi(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3,TLorentzVector Z ){
  float mtau =1.777;
  TLorentzVector a1 = p1+p2+p3;
  float QQ = a1.E()*a1.E() - a1.P()*a1.P();
  float cos = (costheta(p1,p2,p3,Z)*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(costheta(p1,p2,p3,Z)*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
  return cos;

}

 



//---------------------------------------  hadronic current ---------------------------
float TauRecHelper::WA(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ){

  float SS1 = (s2+s3).M()*(s2+s3).M();
  float SS2 = (s1+s3).M()*(s1+s3).M();
  float SS3 = (s1+s2).M()*(s1+s2).M();

 return  VV1(SS1,SS2,SS3,QQ)*F(SS1,SS2,QQ).Rho2() + VV2(SS1,SS2,SS3,QQ)*F(SS2,SS1,QQ).Rho2()  + 2*V1V2(SS1,SS2,SS3,QQ)*( F(SS2,SS1,QQ)*Conjugate(F(SS1,SS2,QQ)) ).Re();
}

float TauRecHelper::WC(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ){

  float SS1 = (s2+s3).M()*(s2+s3).M();
  float SS2 = (s1+s3).M()*(s1+s3).M();
  float SS3 = (s1+s2).M()*(s1+s2).M();

  return  (VV1(SS1,SS2,SS3,QQ) - 2*h(SS1,SS2,SS3,QQ))*F(SS1,SS2,QQ).Rho2() + (VV2(SS1,SS2,SS3,QQ) - 2*h(SS1,SS2,SS3,QQ))*F(SS2,SS1,QQ).Rho2() +   (2*V1V2(SS1,SS2,SS3,QQ) + 4*h(SS1,SS2,SS3,QQ))*( F(SS2,SS1,QQ)*Conjugate(F(SS1,SS2,QQ)) ).Re();
} 

float
TauRecHelper::WD(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ){
  float SS1 = (s2+s3).M()*(s2+s3).M();
  float SS2 = (s1+s3).M()*(s1+s3).M();
  float SS3 = (s1+s2).M()*(s1+s2).M();
  float mpi = 0.139;

  return  -sqrt(h(SS1,SS2,SS3,QQ))*(2*sqrt(VV1(SS1,SS2,SS3,QQ)  -h(SS1,SS2,SS3,QQ)) *F(SS1,SS2,QQ).Rho2() - 2*sqrt(VV2(SS1,SS2,SS3,QQ)  -h(SS1,SS2,SS3,QQ))*F(SS2,SS1,QQ).Rho2()  
			    + (QQ- mpi*mpi + SS3)*(SS1 - SS2 )*( F(SS2,SS1,QQ)*Conjugate(F(SS1,SS2,QQ)) ).Re()/QQ/sqrt(-h0(SS1,SS2,SS3,QQ)));

}

float
TauRecHelper::WE(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ){
  float SS1 = (s2+s3).M()*(s2+s3).M();
  float SS2 = (s1+s3).M()*(s1+s3).M();
  float SS3 = (s1+s2).M()*(s1+s2).M();
  return  -3*sqrt(-h(SS1,SS2,SS3,QQ)*h0(SS1,SS2,SS3,QQ))*( F(SS2,SS1,QQ)*Conjugate(F(SS1,SS2,QQ)) ).Im();

}


TComplex 
TauRecHelper::F(float s1, float s2,float QQ){


  float beta= -0.145;
  float fpi= 0.093;


  TComplex factor(0,-2*sqrt(2)/3/fpi/(1+beta));
  TComplex BreighWignerRhoPrime(beta*BWrhoPrime(s2).Re(), beta*BWrhoPrime(s2).Im());
  TComplex out = factor*BWa1(QQ)*(BWrho(s2) + BreighWignerRhoPrime);
  return out;
}


float
TauRecHelper::VV1(float SS1 ,float SS2, float SS3, float QQ){
  float mpi = 0.139;
  return  SS2 - 4*mpi*mpi + pow(SS3 - SS1,2)/4/QQ;
}

float
TauRecHelper::VV2(float SS1 ,float SS2, float SS3, float QQ){
  float mpi = 0.139;
  return  SS1 - 4*mpi*mpi + pow(SS3 - SS2,2)/4/QQ;
}

float
TauRecHelper::V1V2(float SS1 ,float SS2, float SS3, float QQ){
  float mpi = 0.139;
  return  (QQ/2 - SS3 - mpi*mpi/2) + (SS3 - SS1)*(SS3 - SS2)/4/QQ;
}


float
TauRecHelper::h0(float SS1 ,float SS2, float SS3, float QQ){
  float mpi = 0.139;
  return 4*mpi*mpi - pow(2*mpi*mpi - SS1 - SS2,2)/QQ;
}

float
TauRecHelper::h(float SS1 ,float SS2, float SS3, float QQ){
  float mpi = 0.139;
  return -(SS1*SS2*SS3 - mpi*mpi*pow(QQ - mpi*mpi,2))/h0(SS1,SS2,SS3,QQ)/QQ;
}



TComplex 
TauRecHelper::BWa1(float QQ){
  float m =  1.251;
  TComplex re,im;
  re = (m*m*pow(m*m - QQ,2))/(pow(m*m - QQ,2) - m*m*GammaA1(QQ)*GammaA1(QQ));
  im = m*m*m*GammaA1(QQ)/(pow(m*m - QQ,2) - m*m*GammaA1(QQ)*GammaA1(QQ));
  TComplex out(re,im);
  return out;
}



TComplex 
TauRecHelper::BWrho(float QQ){
  float m =  0.773;
  TComplex re,im;
  re = (m*m*pow(m*m - QQ,2))/(pow(m*m - QQ,2) - m*m*GammaRho(QQ)*GammaRho(QQ));
  im = m*m*m*GammaRho(QQ)/(pow(m*m - QQ,2) - m*m*GammaRho(QQ)*GammaRho(QQ));
  TComplex out(re,im);
  return out;
}


TComplex 
TauRecHelper::BWrhoPrime(float QQ){
  float m =  1.251;
  TComplex re,im;
  re = (m*m*pow(m*m - QQ,2))/(pow(m*m - QQ,2) - m*m*GammaRhoPrime(QQ)*GammaRhoPrime(QQ));
  im = m*m*m*GammaRhoPrime(QQ)/(pow(m*m - QQ,2) - m*m*GammaRhoPrime(QQ)*GammaRhoPrime(QQ));
  TComplex out(re,im);
  return out;
}


float
TauRecHelper::GammaA1(float QQ){
  float ma1 = 1.251;
  float ga1 = 0.599;
  float out = ga1*gForGammaA1(QQ)/gForGammaA1(ma1*ma1);
  return out;
}

float
TauRecHelper::gForGammaA1(float QQ){
  float mpi  = 0.139;
  float mrho = 0.773;
  float out;
  if(QQ > pow((mrho + mpi),2)){ out = QQ*(1.623 + 10.38/QQ - 9.34/QQ/QQ + 0.65/QQ/QQ/QQ);}
  else out = 4.1*pow((QQ - 9*mpi*mpi),3)*(1- 3.3*(QQ - 9*mpi*mpi) + 5.8*pow(QQ - 9*mpi*mpi,2));
  return out;
}


float
TauRecHelper::GammaRho(float QQ){
  float mpi  = 0.139;
  float mrho = 0.773;
  float grho = 0.599;
  float out  =grho*mrho*pow( sqrt(QQ - 4*mpi*mpi)/sqrt(mrho*mrho - mpi*mpi)   ,3)/sqrt(QQ);
  return out;
}


float
TauRecHelper::GammaRhoPrime(float QQ){

  float mrhoPrime = 1.370;
  float grhoPrime = 0.510;
  float out  =grhoPrime*QQ/mrhoPrime/mrhoPrime;
  return out;
}

TComplex 
TauRecHelper::Conjugate(TComplex a){
  return TComplex(a.Re(), -a.Im());
}


