
double ReturnUncert(double a, double b){

  return sqrt(pow(sqrt(a)/b,2) + pow(a*sqrt(b)/b/b,2) );

}

void BackgroundEstimator(){


  TFile *_file0 = TFile::Open("LOCAL_COMBINED_wjetcorrectionfactor_default.root");

   TH1F * MetSSTauIsoMCmu = (TH1F *)_file0->Get("wjetcorrectionfactor_default_MetSSTauIsoMC_WJ1");
   TH1F * MetOSTauIsoMCmu = (TH1F *)_file0->Get("wjetcorrectionfactor_default_MetOSTauIsoMC_WJ1");

   TH1F * MetSSTauNonIsoMCmu = (TH1F *)_file0->Get("wjetcorrectionfactor_default_MetSSTauNonIsoMC_WJ1");
   TH1F * MetOSTauNonIsoMCmu = (TH1F *)_file0->Get("wjetcorrectionfactor_default_MetOSTauNonIsoMC_WJ1");


   TH1F *RatioTauIsoSSMC  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauIsoSSMC_WJ1");
   TH1F *RatioTauIsoOSMC  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauIsoOSMC_WJ1");



   TH1F *RatioTauMet2IsoOS  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauMet2IsoOSMC_WJ1");
   TH1F *RatioTauMet3IsoOS  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauMet3IsoOSMC_WJ1");
   TH1F *RatioTauMet4IsoOS  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauMet4IsoOSMC_WJ1");
   TH1F *RatioTauMet5IsoOS  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauMet5IsoOSMC_WJ1");
   TH1F *RatioTauMet6IsoOS  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauMet6IsoOSMC_WJ1");
   TH1F *RatioTauMet7IsoOS  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauMet7IsoOSMC_WJ1");
   TH1F *RatioTauMet8IsoOS  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauMet8IsoOSMC_WJ1");
   TH1F *RatioTauMet8IsoSS  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauMet8IsoSSMC_WJ1");
   TH1F *RatioTauMet9IsoOS  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauMet9IsoOSMC_WJ1");
   TH1F *RatioTauMet9IsoSS  = (TH1F *)_file0->Get("wjetcorrectionfactor_default_RatioTauMet9IsoSSMC_WJ1");



//   TH1F * MetSSTauIsoMCmu = (TH1F *)_file0->Get("wjetcorrectionfactor_default_MetSSTauIsoSig2MC_WJ1");
//   TH1F * MetOSTauIsoMCmu = (TH1F *)_file0->Get("wjetcorrectionfactor_default_MetOSTauIsoSig2MC_WJ1");

//   TH1F * MetSSTauNonIsoMCmu = (TH1F *)_file0->Get("wjetcorrectionfactor_default_MetSSTauNonIsoMC_WJ1");
//   TH1F * MetOSTauNonIsoMCmu = (TH1F *)_file0->Get("wjetcorrectionfactor_default_MetOSTauNonIsoMC_WJ1");





  double N2565IsoWtoMu=0;   double dN2565IsoWtoMu=0;
  double N2565NonIsoWtoMu=0;  double dN2565NonIsoWtoMu=0;

  double N2565IsoWtoMuSS=0;   double dN2565IsoWtoMuSS=0;
  double N2565NonIsoWtoMuSS=0;  double dN2565NonIsoWtoMuSS=0;



  double NTauToMuIso=0; double dNTauToMuIso=0;
  double NTauToMuNonIso=0; double dNTauToMuNonIso=0;

  double NTauToMuIsoSS=0; double dNTauToMuIsoSS=0;
  double NTauToMuNonIsoSS=0; double dNTauToMuNonIsoSS=0;


  double N25Iso=0,N65Iso=0,N25NonIso=0,N65NonIso=0;
  double N25IsoSS=0,N65IsoSS=0,N25NonIsoSS=0,N65NonIsoSS=0;

  for(int i=0; i <MetOSTauIsoMCmu->GetNbinsX(); i++ ){

    if(i>0 && i <=30){

      N25Iso    +=MetOSTauIsoMCmu->GetBinContent(i);
      N25NonIso +=MetOSTauNonIsoMCmu->GetBinContent(i);

      N25IsoSS+=MetSSTauIsoMCmu->GetBinContent(i);
      N25NonIsoSS+=MetSSTauNonIsoMCmu->GetBinContent(i);
    }
    if(i>69 && i <100){

      N65Iso    +=MetOSTauIsoMCmu->GetBinContent(i);
      N65NonIso +=MetOSTauNonIsoMCmu->GetBinContent(i);

      N65IsoSS+=MetSSTauIsoMCmu->GetBinContent(i);
      N65NonIsoSS+=MetSSTauNonIsoMCmu->GetBinContent(i);

    }
  


  }


  printf("Met2: f2565 Iso  SS: %f  pm   %f \n",RatioTauMet2IsoOS->GetBinContent(1)/RatioTauMet2IsoOS->GetBinContent(2),ReturnUncert(RatioTauMet2IsoOS->GetBinContent(1),RatioTauMet2IsoOS->GetBinContent(2))); 
  printf("Met3: f2565 Iso  SS: %f  pm   %f \n",RatioTauMet3IsoOS->GetBinContent(1)/RatioTauMet3IsoOS->GetBinContent(2),ReturnUncert(RatioTauMet3IsoOS->GetBinContent(1),RatioTauMet3IsoOS->GetBinContent(2))); 
  printf("Met4: f2565 Iso  SS: %f  pm   %f \n",RatioTauMet4IsoOS->GetBinContent(1)/RatioTauMet4IsoOS->GetBinContent(2),ReturnUncert(RatioTauMet4IsoOS->GetBinContent(1),RatioTauMet4IsoOS->GetBinContent(2))); 
  printf("Met5: f2565 Iso  SS: %f  pm   %f \n",RatioTauMet5IsoOS->GetBinContent(1)/RatioTauMet5IsoOS->GetBinContent(2),ReturnUncert(RatioTauMet5IsoOS->GetBinContent(1),RatioTauMet5IsoOS->GetBinContent(2))); 
  printf("Met6: f2565 Iso  SS: %f  pm   %f \n",RatioTauMet6IsoOS->GetBinContent(1)/RatioTauMet6IsoOS->GetBinContent(2),ReturnUncert(RatioTauMet6IsoOS->GetBinContent(1),RatioTauMet6IsoOS->GetBinContent(2))); 
  printf("Met7: f2565 Iso  SS: %f  pm   %f \n",RatioTauMet7IsoOS->GetBinContent(1)/RatioTauMet7IsoOS->GetBinContent(2),ReturnUncert(RatioTauMet7IsoOS->GetBinContent(1),RatioTauMet7IsoOS->GetBinContent(2))); 
  printf("Met8: f2565 Iso  SS: %f  pm   %f \n",RatioTauMet8IsoOS->GetBinContent(1)/RatioTauMet8IsoOS->GetBinContent(2),ReturnUncert(RatioTauMet8IsoOS->GetBinContent(1),RatioTauMet8IsoOS->GetBinContent(2))); 
  printf("Met9: f2565 Iso  SS: %f  pm   %f \n",RatioTauMet9IsoOS->GetBinContent(1)/RatioTauMet9IsoOS->GetBinContent(2),ReturnUncert(RatioTauMet9IsoOS->GetBinContent(1),RatioTauMet9IsoOS->GetBinContent(2))); 


//   N25Iso = RatioTauIsoOSMC->GetBinContent(1);
//   N65Iso = RatioTauIsoOSMC->GetBinContent(2);

//   N25IsoSS = RatioTauIsoSSMC->GetBinContent(1);
//   N65IsoSS = RatioTauIsoSSMC->GetBinContent(2);

  N25Iso = RatioTauMet8IsoOS->GetBinContent(1);
  N65Iso = RatioTauMet8IsoOS->GetBinContent(2);

  N25IsoSS = RatioTauMet8IsoSS->GetBinContent(1);
  N65IsoSS = RatioTauMet8IsoSS->GetBinContent(2);



  printf("%f %f \n",N25Iso,N65Iso);

  N2565IsoWtoMu = N25Iso/N65Iso;dN2565IsoWtoMu = ReturnUncert(N25Iso,N65Iso);//N25Iso/N65Iso;  
  N2565NonIsoWtoMu = N25NonIso/N65NonIso;dN2565NonIsoWtoMu=ReturnUncert(N25NonIso,N65NonIso);//N25NonIso/N65NonIso;

  N2565IsoWtoMuSS = N25IsoSS/N65IsoSS;dN2565IsoWtoMuSS = ReturnUncert(N25IsoSS,N65IsoSS);//N25Iso/N65Iso;  
  N2565NonIsoWtoMuSS = N25NonIsoSS/N65NonIsoSS;dN2565NonIsoWtoMuSS=ReturnUncert(N25NonIsoSS,N65NonIsoSS);//N25NonIso/N65NonIso;



  printf("f2565 Iso  OS: %f  pm   %f ",N2565IsoWtoMu,dN2565IsoWtoMu);   printf("     f2565 nonIso  OS: %f  pm   %f  \n",N2565NonIsoWtoMu,dN2565NonIsoWtoMu);
  printf("f2565 Iso  SS: %f  pm   %f ",N2565IsoWtoMuSS,dN2565IsoWtoMuSS);   printf("    f2565 nonIso  SS: %f  pm   %f \n",N2565NonIsoWtoMuSS,dN2565NonIsoWtoMuSS);



  TFile *_file1 = TFile::Open("LOCAL_COMBINED_backgroundestimator_default.root");


  TH1F * aqcd = (TH1F *)_file1->Get("backgroundestimator_default_aqcdData");
  TH1F * bqcd = (TH1F *)_file1->Get("backgroundestimator_default_bqcdData");

  TH1F * bqcdTT = (TH1F *)_file1->Get("backgroundestimator_default_bqcdTT");
  TH1F * bqcdSig = (TH1F *)_file1->Get("backgroundestimator_default_bqcdSignal");
  TH1F * bqcdDY = (TH1F *)_file1->Get("backgroundestimator_default_bqcdMC_DY");
  TH1F * bqcdDYTT = (TH1F *)_file1->Get("backgroundestimator_default_bqcdMC_TT");
  TH1F * bqcdWW = (TH1F *)_file1->Get("backgroundestimator_default_bqcdEWK");


  TH1F * cqcd = (TH1F *)_file1->Get("backgroundestimator_default_c2qcdData");

  TH1F * cqcdTT = (TH1F *)_file1->Get("backgroundestimator_default_c2qcdTT");
  TH1F * cqcdSig = (TH1F *)_file1->Get("backgroundestimator_default_c2qcdSignal");
  TH1F * cqcdDY = (TH1F *)_file1->Get("backgroundestimator_default_c2qcdMC_DY");
  TH1F * cqcdDYTT = (TH1F *)_file1->Get("backgroundestimator_default_c2qcdMC_TT");
  TH1F * cqcdWW = (TH1F *)_file1->Get("backgroundestimator_default_c2qcdEWK");


  TH1F * dqcd = (TH1F *)_file1->Get("backgroundestimator_default_d2qcdData");

  TH1F * cqcdSignal = (TH1F *)_file1->Get("backgroundestimator_default_c2qcdSignal");
  TH1F * dqcdSignal = (TH1F *)_file1->Get("backgroundestimator_default_d2qcdSignal");


  TH1F * w65os = (TH1F *)_file1->Get("backgroundestimator_default_wj65osData");
  TH1F * w65ss = (TH1F *)_file1->Get("backgroundestimator_default_wj65ssData");


  TH1F * w65osTT = (TH1F *)_file1->Get("backgroundestimator_default_wj65osTT");
  TH1F * w65ssTT = (TH1F *)_file1->Get("backgroundestimator_default_wj65ssTT");

  TH1F * w65osDY = (TH1F *)_file1->Get("backgroundestimator_default_wj65osMC_DY");
  TH1F * w65ssDY = (TH1F *)_file1->Get("backgroundestimator_default_wj65ssMC_DY");


  TH1F * w65osSig = (TH1F *)_file1->Get("backgroundestimator_default_wj65osSignal");
  TH1F * w65ssSig = (TH1F *)_file1->Get("backgroundestimator_default_wj65ssSignal");

  TH1F * w65osDYTT = (TH1F *)_file1->Get("backgroundestimator_default_wj65osMC_TT");
  TH1F * w65ssDYTT = (TH1F *)_file1->Get("backgroundestimator_default_wj65ssMC_TT");


  TH1F * w65osWW = (TH1F *)_file1->Get("backgroundestimator_default_wj65osEWK");
  TH1F * w65ssWW = (TH1F *)_file1->Get("backgroundestimator_default_wj65ssEWK");

//    double WOSLow = (w65os->GetEntries() - w65osTT->GetEntries()*0.627  - w65osWW->GetEntries()*0.0579)* N2565IsoWtoMu;
//    double WSSLow = (w65ss->GetEntries() - w65ssTT->GetEntries()*0.627  - w65ssWW->GetEntries()*0.0579)* N2565IsoWtoMuSS;

//   double WOSLow = 2437.83;//(w65os->GetEntries() - w65osTT->GetEntries()*2.049 )* N2565IsoWtoMu;
//    // double WOSLow = (w65os->GetEntries() - w65osTT->GetEntries()*2.049 )* 0.25;
//   double WSSLow = 779.3;//(w65ss->GetEntries() - w65ssTT->GetEntries()*2.049 )* N2565IsoWtoMuSS;

  double WOSLow = (w65os->GetEntries() - w65osTT->GetEntries()*0.2049 -  w65osDY->GetEntries() * 2.28 - w65osDYTT->GetEntries() * 2.24  - w65osWW->GetEntries() * 0.03 -  w65osSig->GetEntries() * 2.24)*  N2565IsoWtoMu;
   // double WOSLow = (w65os->GetEntries() - w65osTT->GetEntries()*2.049 )* 0.25;
  double WSSLow = (w65ss->GetEntries() - w65ssTT->GetEntries()*0.2049 -  w65ssDY->GetEntries() * 2.28 - w65ssDYTT->GetEntries() * 2.24  - w65ssWW->GetEntries() * 0.03 -  w65ssSig->GetEntries() * 2.24 )* N2565IsoWtoMuSS;




//   double WOSLow = (w65os->GetEntries()  )* N2565IsoWtoMu;
//   double WSSLow = (w65ss->GetEntries()  )* N2565IsoWtoMuSS;


   double OS2SS = (cqcd->GetEntries() - backgroundestimator_default_c2qcdSignal->GetEntries()*2.33 - cqcdTT->GetEntries() * 2.24  - cqcdDY->GetEntries()* 2.28 - cqcdDYTT->GetEntries() *  2.24 - cqcdWW->GetEntries()*0.03)   /dqcd->GetEntries();

   double SSQCD = backgroundestimator_default_bqcdData->GetEntries() - WSSLow - bqcdTT->GetEntries() * 2.24  - bqcdSig->GetEntries()*2.24 - bqcdDY->GetEntries()* 2.28 - bqcdDYTT->GetEntries() *  2.24 - bqcdWW->GetEntries()*0.03;
   double OSQCD = SSQCD*OS2SS;


   printf("QCD events in B  = %f,  OS2SS = %f   pm  %f \n",SSQCD, OS2SS,ReturnUncert((cqcd->GetEntries() - backgroundestimator_default_c2qcdSignal->GetEntries()*2.33),dqcd->GetEntries()) );
   printf("OSQCD = %f   OSWLow = %f \n", OSQCD,WOSLow );


  printf("QCD Estimation  = %f   \n",OSQCD  );
  printf("WJ  Estimation  = %f   \n",WOSLow );

  printf("Total Selected  = %f   \n",aqcd->GetEntries() );

//   double ttbar = 22.00;
//   double ztautau = 28.77;
//   double murho = 0;//311.058;
//   double Fpi = 0;//1185.;
//   double dy  = 21.47;
//   double ww  = 12.26;

  double ttbar = 47.78;
  double ztautau = 88.28;
  double murho = 0;//311.058;
  double Fpi = 0;//1185.;
  double dy  = 24.65;
  double ww  = 11.62;

// ohne significance kut 
//   double ttbar = 42.963;
//   double ztautau = 76.8;
//   double murho = 0;//311.058;
//   double Fpi = 0;//1185.;
//   double dy  = 185.32;
//   double ww  = 20.144;


  printf("Signal Estimation = %f   \n",aqcd->GetEntries()  - OSQCD - WOSLow - ww - dy - ztautau  - ttbar - Fpi - murho);


 
  printf("Total selected %f \n", aqcd->GetEntries() );

}
  
//   KEY: TH1D     backgroundestimator_default_aqcdData;1  Data
//   KEY: TH1D     backgroundestimator_default_bqcdData;1  Data
//   KEY: TH1D     backgroundestimator_default_cqcdData;1  Data
//   KEY: TH1D     backgroundestimator_default_dqcdData;1  Data

//   KEY: TH1D     backgroundestimator_default_wj65osData;1        Data
//   KEY: TH1D     backgroundestimator_default_wj65ssData;1        Data
//   KEY: TH1D     backgroundestimator_default_wj25osData;1        Data
//   KEY: TH1D     backgroundestimator_default_wj25ssData;1        Data
