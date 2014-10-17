
#include "EFValidation.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauSolver.h"
//#include "TCanvas.h"
//#include "NTuple_Controller.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"

EFValidation::EFValidation(TString Name_, TString id_):
  Selection(Name_,id_)
  ,muoneta(2.1)
  ,jeteta(2.3)
  ,TauTrackPtThreshold(5.0)
		//  ,ran()
{
  verbose=false;
}

EFValidation::~EFValidation(){
  for(int j=0; j<Npassed.size(); j++){
    //std::cout << "EFValidation::~EFValidation Selection Summary before: "
	// << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 //<< Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  //std::cout << "EFValidation::~EFValidation()" << std::endl;
}

void  EFValidation::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==hasMuon)	      cut.at(hasMuon)=0;
    if(i==hasTau)	      cut.at(hasTau)=0;

  }
  TString hlabel;
  TString htitle;
  for(unsigned int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    Accumdist.push_back(std::vector<TH1D>());
    Nminus1dist.push_back(std::vector<TH1D>());


    TString c="_Cut_";
    if(i<10)c+="0";
    c+=i;
      if(i==TriggerOk){
      title.at(i)="$N_{\\mu} passing the Trigger$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#mu} passing the Trigger";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }  

    else if(i==hasMuon){
      title.at(i)="$Number of good muons$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of good muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }


    else if(i==hasTau){
      title.at(i)="$Number of good taus$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of good taus";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTau_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTau_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }

    //-----------
   /*else if(i==ZMassmax){
     title.at(i)="$M_{Z}<$";
     title.at(i)+=cut.at(ZMassmax);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{Z} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmax_",htitle,40,0,200,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmax_",htitle,40,0,200,hlabel,"Events"));
   }*/
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");


  // Setup Extra Histograms

  idAllTaus=HConfig.GetTH1D(Name+"_idAllTaus","idAllTaus",4,0.5,4.5,"idAllTaus","Events");
  idIdentifiedTau=HConfig.GetTH1D(Name+"_idIdentifiedTau","idIdentifiedTau",4,0.5,4.5," idIdentifiedTau","Events");
  idPassedEF=HConfig.GetTH1D(Name+"_idPassedEF","idPassedEF",4,0.5,4.5,"idPassedEF ","Events");


  TruthA1Pt=HConfig.GetTH1D(Name+"_TruthA1Pt","Truth A1 Pt",50,10,60,"A1 Pt, (GeV)","");
  TruthA1Phi=HConfig.GetTH1D(Name+"_TruthA1Phi","Truth A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  TruthA1Eta=HConfig.GetTH1D(Name+"_TruthA1Eta","Truth A1 Eta",50,-2,2,"A1 Eta","");
  
  TruthMuPt=HConfig.GetTH1D(Name+"_TruthMuPt","Truth Mu Pt",50,10,60,"#mu Pt, (GeV)","");
  TruthMuPhi=HConfig.GetTH1D(Name+"_TruthMuPhi","Truth Mu Phi",50,-3.14,3.14,"#mu Phi, (rad)","");
  TruthMuEta=HConfig.GetTH1D(Name+"_TruthMuEta","Truth Mu Eta",50,-2,2,"#mu Eta","");
  
  TruthTauMuPt=HConfig.GetTH1D(Name+"_TruthTauMuPt","Truth TauMu Pt",50,10,60,"#tau_{#mu} Pt, (GeV)","");
  TruthTauMuPhi=HConfig.GetTH1D(Name+"_TruthTauMuPhi","Truth TauMu Phi",50,-3.14,3.14,"#tau_{#mu} Phi, (rad)","");
  TruhTauMuEta=HConfig.GetTH1D(Name+"_TruthTauMuEta","Truth TauMu Eta",50,-2,2,"#tau_{#mu} Eta","");
  
  TruthTauA1Pt=HConfig.GetTH1D(Name+"_TruthTauA1Pt","Truth TauA1 Pt",50,10,60,"#tau_{a1} Pt, (GeV)","");
  TruthTauA1Phi=HConfig.GetTH1D(Name+"_TruthTauA1Phi","Truth TauA1 Phi",50,-3.14,3.14,"#tau_{a1} Phi, (rad)","");
  TruthTauA1Eta=HConfig.GetTH1D(Name+"_TruthTauA1Eta","Truth TauA1 Eta",50,-2,2,"#tau_{a1} Eta","");
  
  
  TauA1PtResolKFM=HConfig.GetTH1D(Name+"_TauA1PtResolKFM","Pt resolution of #tau_{a1} (TPF ambiguity minus)",50,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauA1PhiResolKFM=HConfig.GetTH1D(Name+"_TauA1PhiResolKFM","#phi resolution of #tau_{a1} (TPF ambiguity minus)",50,-1,1,"#delta #phi (mc - rec) ","");
  TauA1EtaResolKFM=HConfig.GetTH1D(Name+"_TauA1EtaResolKFM","#eta resolution of #tau_{a1} (TPF ambiguity minus)",50,-1,1,"#delta #eta (mc - rec) ","");
  

  TauA1PtResolKFP=HConfig.GetTH1D(Name+"_TauA1PtResolKFP","Pt resolution of #tau_{a1} (TPF ambiguity plus)",50,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauA1PhiResolKFP=HConfig.GetTH1D(Name+"_TauA1PhiResolKFP","#phi resolution of #tau_{a1} (TPF ambiguity plus)",50,-1,1,"#delta #phi (mc - rec) ","");
  TauA1EtaResolKFP=HConfig.GetTH1D(Name+"_TauA1EtaResolKFP","#eta resolution of #tau_{a1} (TPF ambiguity plus)",50,-1,1,"#delta #eta (mc - rec) ","");


  RecoA1Pt=HConfig.GetTH1D(Name+"_RecoA1Pt","Reco A1 Pt",50,10,60,"A1 Pt, (GeV)","");
  RecoA1Phi=HConfig.GetTH1D(Name+"_RecoA1Phi","Reco A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  RecoA1Eta=HConfig.GetTH1D(Name+"_RecoA1Eta","Reco A1 Eta",50,-2,2,"A1 Eta","");

  RecoMuPt=HConfig.GetTH1D(Name+"_RecoMuPt","Reco Mu Pt",50,10,60,"#mu Pt, (GeV)","");
  RecoMuPhi=HConfig.GetTH1D(Name+"_RecoMuPhi","Reco Mu Phi",50,-3.14,3.14,"#mu Phi, (rad)","");
  RecoMuEta=HConfig.GetTH1D(Name+"_RecoMuEta","Reco Mu Eta",50,-2,2,"#mu Eta","");


  EFitTauA1Pt=HConfig.GetTH1D(Name+"_FitTauA1Pt"," TauA1 Pt",50,10,60,"#tau_{a1} Pt, (GeV)","");
  EFitTauA1Phi=HConfig.GetTH1D(Name+"_FitTauA1Phi"," TauA1 Phi",50,-3.14,3.14,"#tau_{a1} Phi, (rad)","");
  EFitTauA1Eta=HConfig.GetTH1D(Name+"_FitTauA1Eta"," TauA1 Eta",50,-2,2,"#tau_{a1} Eta","");

  EFitTauMuPt=HConfig.GetTH1D(Name+"_FitTauMuPt"," TauMu Pt",50,10,60,"#tau_{#mu} Pt, (GeV)","");
  EFitTauMuPhi=HConfig.GetTH1D(Name+"_FitTauMuPhi"," TauMu Phi",50,-3.14,3.14,"#tau_{#mu} Phi, (rad)","");
  EFitTauMuEta=HConfig.GetTH1D(Name+"_FitTauMuEta"," TauMu Eta",50,-2,2,"#tau_{#mu} Eta","");



  EFitTauA1PtAmbPoint=HConfig.GetTH1D(Name+"_EFitTauA1PtAmbPoint"," TauA1 Pt",50,10,60,"#tau_{a1} Pt, (GeV)","");
  EFitTauA1PhiAmbPoint=HConfig.GetTH1D(Name+"_EFitTauA1PhiAmbPoint"," TauA1 Phi",50,-3.14,3.14,"#tau_{a1} Phi, (rad)","");
  EFitTauA1EtaAmbPoint=HConfig.GetTH1D(Name+"_EFitTauA1EtaAmbPoint"," TauA1 Eta",50,-2,2,"#tau_{a1} Eta","");

  EFitTauMuPtAmbPoint=HConfig.GetTH1D(Name+"_EFitTauMuPtAmbPoint"," TauMu Pt",50,10,60,"#tau_{#mu} Pt, (GeV)","");
  EFitTauMuPhiAmbPoint=HConfig.GetTH1D(Name+"_EFitTauMuPhiAmbPoint"," TauMu Phi",50,-3.14,3.14,"#tau_{#mu} Phi, (rad)","");
  EFitTauMuEtaAmbPoint=HConfig.GetTH1D(Name+"_EFitTauMuEtaAmbPoint"," TauMu Eta",50,-2,2,"#tau_{#mu} Eta","");



  A1Mass=HConfig.GetTH1D(Name+"_A1Mass","a1 mass",50,0.8,1.777,"M_{a1}, GeV","");
  
  TruthA1PtAfterFit=HConfig.GetTH1D(Name+"_TruthA1PtAfterFit","Truth A1 Pt",50,10,60,"A1 Pt, (GeV)","");
  TruthA1EAfterFit=HConfig.GetTH1D(Name+"_TruthA1EAfterFit","Truth A1 E",50,10,100,"A1 E, (GeV)","");
  TruthA1PhiAfterFit=HConfig.GetTH1D(Name+"_TruthA1PhiAfterFit","Truth A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  TruthA1EtaAfterFit=HConfig.GetTH1D(Name+"_TruthA1EtaAfterFit","Truth A1 Eta",50,-2,2,"A1 Eta","");
  
  TruthMuPtAfterFit=HConfig.GetTH1D(Name+"_TruthMuPtAfterFit","Truth Mu Pt",50,10,60,"#mu Pt, (GeV)","");
  TruthMuEAfterFit=HConfig.GetTH1D(Name+"_TruthMuEAfterFit","Truth Mu E",50,10,100,"#mu E, (GeV)","");
  TruthMuPhiAfterFit=HConfig.GetTH1D(Name+"_TruthMuPhiAfterFit","Truth A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  TruthMuEtaAfterFit=HConfig.GetTH1D(Name+"_TruthMuEtaAfterFit","Truth A1 Eta",50,-2,2,"A1 Eta","");
  
  EfficiencyOverA1Pt=HConfig.GetTH1D(Name+"_EfficiencyOverA1Pt","Truth A1 Pt",50,10,60,"A1 Pt, (GeV)","");
  EfficiencyOverA1Phi=HConfig.GetTH1D(Name+"_EfficiencyOverA1Phi","Truth A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  EfficiencyOverA1Eta=HConfig.GetTH1D(Name+"_EfficiencyOverA1Eta","Truth A1 Eta",50,-2,2,"A1 Eta","");


  EfficiencyOverMuPt=HConfig.GetTH1D(Name+"_EfficiencyOverMuPt","Truth Mu Pt",50,10,60,"#mu Pt, (GeV)","");
  EfficiencyOverMuPhi=HConfig.GetTH1D(Name+"_EfficiencyOverMuPhi","Truth A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  EfficiencyOverMuEta=HConfig.GetTH1D(Name+"_EfficiencyOverMuEta","Truth A1 Eta",50,-2,2,"A1 Eta","");


  
  TauA1PtResolution=HConfig.GetTH1D(Name+"_TauA1PtResolution","Pt resolution of #tau_{a1}",80,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauA1EResolution=HConfig.GetTH1D(Name+"_TauA1EResolution","E resolution of #tau_{a1}",80,-100,100,"#delta E (mc - rec), GeV ","");
  TauA1PhiResolution=HConfig.GetTH1D(Name+"_TauA1PhiResolution","#phi resolution of #tau_{a1}",50,-1,1,"#delta #phi (mc - rec) ","");
  TauA1EtaResolution=HConfig.GetTH1D(Name+"_TauA1EtaResolution","#eta resolution of #tau_{a1}",50,-1,1,"#delta #eta (mc - rec) ","");
  
  TauMuPtResolution=HConfig.GetTH1D(Name+"_TauMuPtResolution","Pt resolution of #tau_{#mu}",80,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauMuEResolution=HConfig.GetTH1D(Name+"_TauMuEResolution","E resolution of #tau_{#mu}",80,-100,100,"#delta E (mc - rec), GeV ","");
  TauMuPhiResolution=HConfig.GetTH1D(Name+"_TauMuPhiResolution","#phi resolution of #tau_{#mu}",50,-2,2,"#delta #phi (mc - rec) ","");
  TauMuEtaResolution=HConfig.GetTH1D(Name+"_TauMuEtaResolution","#eta resolution of #tau_{#mu}",50,-2,2,"#delta #eta (mc - rec) ","");

  TauA1PtResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauA1PtResolutionAmbPoint","Pt resolution of #tau_{a1}",80,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauA1PhiResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauA1PhiResolutionAmbPoint","#phi resolution of #tau_{a1}",50,-1,1,"#delta #phi (mc - rec) ","");
  TauA1EtaResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauA1EtaResolutionAmbPoint","#eta resolution of #tau_{a1}",50,-1,1,"#delta #eta (mc - rec) ","");

  TauMuPtResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauMuPtResolutionAmbPoint","Pt resolution of #tau_{#mu}",80,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauMuPhiResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauMuPhiResolutionAmbPoint","#phi resolution of #tau_{#mu}",50,-2,2,"#delta #phi (mc - rec) ","");
  TauMuEtaResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauMuEtaResolutionAmbPoint","#eta resolution of #tau_{#mu}",50,-2,2,"#delta #eta (mc - rec) ","");


  TauA1PhiResolVsProb=HConfig.GetTH2D(Name+"_TauA1PhiResolVsProb","Phi resolution of #tau_{a1} vs Probability  ",200,0,1,60,-1,1,"Probability ","#delta #phi (mc - rec)");
  TauMuPhiResolVsProb=HConfig.GetTH2D(Name+"_TauMuPhiResolVsProb","Phi resolution of #tau_{#mu} vs Probability ",200,0,1,60,-1,1,"Probability ","#delta #phi (mc - rec)");

  TauA1EtaResolVsProb=HConfig.GetTH2D(Name+"_TauA1EtaResolVsProb","Eta resolution of #tau_{a1} vs Probability  ",200,0,1,60,-1,1,"Probability ","#delta #eta (mc - rec)");
  TauMuEtaResolVsProb=HConfig.GetTH2D(Name+"_TauMuEtaResolVsProb","Eta resolution of #tau_{#mu} vs Probability  ",200,0,1,60,-1,1,"Probability ","#delta #eta (mc - rec)");

  TauA1PtResolVsProb=HConfig.GetTH2D(Name+"_TauA1PtResolVsProb","Pt resolution of #tau_{a1} vs Probability  ",200,0,1,80,-100,100,"Probability ","#delta pt (mc - rec)");
  TauMuPtResolVsProb=HConfig.GetTH2D(Name+"_TauMuPtResolVsProb","Pt resolution of #tau_{#mu} vs Probability  ",200,0,1,80,-100,100,"Probability ","#delta pt (mc - rec)");

  TauA1EResolVsProb=HConfig.GetTH2D(Name+"_TauA1EResolVsProb","E resolution of #tau_{a1} vs Probability  ",200,0,1,80,-100,100,"Probability ","#delta E (mc - rec)");
  TauMuEResolVsProb=HConfig.GetTH2D(Name+"_TauMuEResolVsProb","E resolution of #tau_{#mu} vs Probability  ",200,0,1,80,-100,100,"Probability ","#delta E (mc - rec)");


  TauA1PtResolVsPt=HConfig.GetTH2D(Name+"_TauA1PtResolVsPt","Pt resolution of #tau_{a1} vs a1 Pt",50,10,60,80,-100,100,"#tau_{a1} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsPt=HConfig.GetTH2D(Name+"_TauMuPtResolVsPt","Pt resolution of #tau_{#mu} vs #mu Pt",50,10,60,80,-100,100,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");

  TauA1PtResolVsEta=HConfig.GetTH2D(Name+"_TauA1PtResolVsEta","Pt resolution of #tau_{a1} vs a1 Eta",50,-2,2,80,-100,100,"a1 Eta ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsEta=HConfig.GetTH2D(Name+"_TauMuPtResolVsEta","Pt resolution of #tau_{#mu} vs  #mu Eta",50,-2,2,80,-100,100,"#mu Eta ","#delta Pt (mc - rec), GeV");
 
  TauA1EtaResolVsEta=HConfig.GetTH2D(Name+"_TauA1EtaResolVsEta","#eta resolution of #tau_{a1} vs a1 Eta",50,-2,2,50,-1,1,"a1 Eta ","#delta #eta (mc - rec)");
  TauMuEtaResolVsEta=HConfig.GetTH2D(Name+"_TauMuEtaResolVsEta","#eta resolution of #tau_{#mu} vs  #mu Eta",50,-2,2,50,-1,1,"#mu Eta ","#delta #eta (mc - rec)");

  
  TauA1PtResolVsZPt=HConfig.GetTH2D(Name+"_TauA1PtResolVsZPt","Pt resolution of #tau_{a1} vs Z Pt",50,0,20,80,-100,100,"#tau_{a1} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsZPt=HConfig.GetTH2D(Name+"_TauMuPtResolVsZPt","Pt resolution of #tau_{#mu} vs Z  Pt",50,0,20,80,-100,100,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");
 
  TauA1EtaResolVsZPt=HConfig.GetTH2D(Name+"_TauA1EtaResolVsZPt","Pt resolution of #tau_{a1} vs Z Pt",50,0,20,80,-100,100,"a1 Eta ","#delta #eta (mc - rec)");
  TauMuEtaResolVsZPt=HConfig.GetTH2D(Name+"_TauMuEtaResolVsZPt","Pt resolution of #tau_{#mu} vs  Z Pt",50,0,20,80,-100,100,"#mu Eta ","#delta #eta (mc - rec)");

  TauMuPtResolVsZPtMean=HConfig.GetTH1D(Name+"_TauMuPtResolVsZPtMean","Pt resolution of #tau_{#mu} vs Z  Pt",50,0,20,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsZPtRMS=HConfig.GetTH1D(Name+"_TauMuPtResolVsZPtRMS","Pt resolution of #tau_{#mu} vs Z  Pt",50,0,20,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauA1PtResolVsEtaRMS=HConfig.GetTH1D(Name+"_TauA1PtResolVsEtaRMS","Pt resolution of #tau_{a1} vs a1 Eta",50,-2,2,"a1 Eta ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsPtRMS=HConfig.GetTH1D(Name+"_TauMuPtResolVsPtRMS","Pt resolution of #tau_{#mu} vs #mu Pt",50,10,60,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");


  TauA1PtResolVsPtRMS=HConfig.GetTH1D(Name+"_TauA1PtResolVsPtRMS","Pt resolution of #tau_{a1} vs a1 Pt",50,10,60,"#tau_{a1} Pt, GeV ","#delta Pt (mc - rec), GeV");


  TauA1PhiResolVsPVSVSignificance=HConfig.GetTH2D(Name+"_TauA1PhiResolVsPVSVSignificance","Phi resolution of #tau_{a1} vs PVSV significance ",50,0,20,60,-1,1,"#sigma (PV-SV) ","#delta #phi (mc - rec)");
  TauA1EtaResolVsPVSVSignificance=HConfig.GetTH2D(Name+"_TauA1EtaResolVsPVSVSignificance","Eta resolution of #tau_{a1} vs PVSV significance ",50,0,20,60,-1,1,"#sigma (PV-SV) ","#delta #eta (mc - rec)");
  TauA1PtResolVsPVSVSignificance=HConfig.GetTH2D(Name+"_TauA1PtResolVsPVSVSignificance","Pt resolution of #tau_{a1} vs PVSV significance ",50,0,20,80,-100,100,"#sigma (PV-SV) ","#delta #pt (mc - rec), GeV");
  TauA1EResolVsPVSVSignificance=HConfig.GetTH2D(Name+"_TauA1EResolVsPVSVSignificance","E resolution of #tau_{a1} vs PVSV significance ",50,0,20,80,-100,100,"#sigma (PV-SV) ","#delta E (mc - rec), GeV");
  
  PVSVSignificanceCutEfficiency=HConfig.GetTH1D(Name+"_PVSVSignificanceCutEfficiency","PVSVSignificanceCutEfficiency",50,0,20,"#sigma (PV-SV)","#epsilon of a cut");
  TauA1PtResolVsPVSVSignificanceRMS=HConfig.GetTH1D(Name+"_TauA1PtResolVsPVSVSignificanceRMS","Pt resolution of #tau_{a1} vs PVSV significance ",50,0,20,"#sigma (PV-SV) ","#delta #pt (mc - rec), GeV");


  TauA1PtResolVsPtAmbPoint=HConfig.GetTH2D(Name+"_TauA1PtResolVsPtAmbPoint","Pt resolution of #tau_{a1} vs a1 Pt",50,10,60,80,-100,100,"#tau_{a1} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsPtAmbPoint=HConfig.GetTH2D(Name+"_TauMuPtResolVsPtAmbPoint","Pt resolution of #tau_{#mu} vs #mu Pt",50,10,60,80,-100,100,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauA1PtResolVsEtaAmbPoint=HConfig.GetTH2D(Name+"_TauA1PtResolVsEtaAmbPoint","Pt resolution of #tau_{a1} vs a1 Eta",50,-2,2,80,-100,100,"a1 Eta ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsEtaAmbPoint=HConfig.GetTH2D(Name+"_TauMuPtResolVsEtaAmbPoint","Pt resolution of #tau_{#mu} vs  #mu Eta",50,-2,2,80,-100,100,"#mu Eta ","#delta Pt (mc - rec), GeV");
  TauA1EtaResolVsEtaAmbPoint=HConfig.GetTH2D(Name+"_TauA1EtaResolVsEtaAmbPoint","Pt resolution of #tau_{#mu} vs  #mu Eta",50,-2,2,50,-2,2,"#mu Eta ","#delta Pt (mc - rec), GeV");
  TauMuEtaResolVsEtaAmbPoint=HConfig.GetTH2D(Name+"_TauMuEtaResolVsEtaAmbPoint","Pt resolution of #tau_{#mu} vs  #mu Eta",50,-2,2,50,-2,2,"#mu Eta ","#delta Pt (mc - rec), GeV");

  TauA1PtResolVsZPtAmbPoint=HConfig.GetTH2D(Name+"_TauA1PtResolVsZPtAmbPoint","Pt resolution of #tau_{a1} vs Z Pt",50,0,20,80,-100,100,"#tau_{a1} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsZPtAmbPoint=HConfig.GetTH2D(Name+"_TauMuPtResolVsZPtAmbPoint","Pt resolution of #tau_{#mu} vs Z  Pt",50,0,20,80,-100,100,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");

  TauA1EtaResolVsZPtAmbPoint=HConfig.GetTH2D(Name+"_TauA1EtaResolVsZPtAmbPoint","Pt resolution of #tau_{a1} vs Z Pt",50,0,20,80,-100,100,"a1 Eta ","#delta #eta (mc - rec)");
  TauMuEtaResolVsZPtAmbPoint=HConfig.GetTH2D(Name+"_TauMuEtaResolVsZPtAmbPoint","Pt resolution of #tau_{#mu} vs  Z Pt",50,0,20,80,-100,100,"#mu Eta ","#delta #eta (mc - rec)");


  TauA1PhiResolVsPhiRotSignificance=HConfig.GetTH2D(Name+"_TauA1PhiResolVsPhiRotSignificance","Phi resolution of #tau_{a1} vs PVSV significance ",50,0,10,60,-1,1,"#sigma (phi rotation) ","#delta #phi (mc - rec)");
  TauA1EtaResolVsPhiRotSignificance=HConfig.GetTH2D(Name+"_TauA1EtaResolVsPhiRotSignificance","Eta resolution of #tau_{a1} vs PVSV significance ",50,0,10,60,-1,1,"#sigma (phi rotation) ","#delta #eta (mc - rec)");
  TauA1PtResolVsPhiRotSignificance=HConfig.GetTH2D(Name+"_TauA1PtResolVsPhiRotSignificance","Pt resolution of #tau_{a1} vs PVSV significance ",50,0,10,80,-100,100,"#sigma (phi rotation) ","#delta #pt (mc - rec), GeV");
  TauA1EResolVsPhiRotSignificance=HConfig.GetTH2D(Name+"_TauA1EResolVsPhiRotSignificance","E resolution of #tau_{a1} vs PVSV significance ",50,0,10,80,-100,100,"#sigma (phi rotation) ","#delta E (mc - rec), GeV");
  
  PhiRotSignificanceCutEfficiency=HConfig.GetTH1D(Name+"_PhiRotSignificanceCutEfficiency","PhiRotSignificanceCutEfficiency",50,0,10,"#sigma (phi rotation) ","#epsilon of a cut");


  TauA1PhiResolVsPVSVSignificanceAmbPoint=HConfig.GetTH2D(Name+"_TauA1PhiResolVsPVSVSignificanceAmbPoint","Phi resolution of #tau_{a1} vs PVSV significance ",50,0,10,60,-1,1,"#sigma (PV-SV) ","#delta #phi (mc - rec)");
  TauA1EtaResolVsPVSVSignificanceAmbPoint=HConfig.GetTH2D(Name+"_TauA1EtaResolVsPVSVSignificanceAmbPoint","Eta resolution of #tau_{a1} vs PVSV significance ",50,0,10,60,-1,1,"#sigma (PV-SV) ","#delta #phi (mc - rec)");
  TauA1PtResolVsPVSVSignificanceAmbPoint=HConfig.GetTH2D(Name+"_TauA1PtResolVsPVSVSignificanceAmbPoint","Pt resolution of #tau_{a1} vs PVSV significanc ",50,0,10,80,-100,100,"#sigma (PV-SV)","#delta #pt (mc - rec), GeV");
  TauA1EResolVsPVSVSignificanceAmbPoint=HConfig.GetTH2D(Name+"_TauA1EResolVsPVSVSignificanceAmbPoint","E resolution of #tau_{a1} vs PVSV significance ",50,0,10,80,-100,100,"#sigma (PV-SV) ","#delta E (mc - rec), GeV");


  
  ProbabilityOfCorrect=HConfig.GetTH1D(Name+"_ProbabilityOfCorrect","ProbabilityOfCorrect",300,0,1,"Probability","Events");

  ProbabilityOfCorrectndf2=HConfig.GetTH1D(Name+"_ProbabilityOfCorrectndf2","ProbabilityOfCorrectndf2",300,0,1,"Probability ndf 2","Events");
  ProbabilityOfCorrectndf3=HConfig.GetTH1D(Name+"_ProbabilityOfCorrectndf3","ProbabilityOfCorrectndf3",300,0,1,"Probability ndf 3","Events");
  Chi2OfCorrect=HConfig.GetTH1D(Name+"_Chi2OfCorrect","Chi2OfCorrect",300,0,50,"Probability","Events");


  ProbabilityOfAmbiguityPoint=HConfig.GetTH1D(Name+"_ProbabilityOfAmbiguityPoint","ProbabilityOfAmbiguityPoint",300,0,1,"Probability","Events");
  Chi2Dim=HConfig.GetTH2D(Name+"_Chi2Dim","Chi2Dim",100,0,1,100,0,1,"Chi2Dim","Events");
  csum2Dim=HConfig.GetTH2D(Name+"_csum2Dim","csum2Dim",200,0,1,200,0,1,"csum, GeV","Events");
  CorrectAmbiguityTruth=HConfig.GetTH1D(Name+"_CorrectAmbiguityTruth","CorrectAmbiguityTruth",2,0.5,2.5,"Correct Ambiguity Truth ","Events");

  MeasuredPhi=HConfig.GetTH1D(Name+"_MeasuredPhi","Measured TauA1 phi",40,-3.14,3.14,"  #phi, (rad) ","");

  pvsvsignificance=HConfig.GetTH1D(Name+"_pvsvsignificance","pvsv significance",40,0,10,"  #sigma (PV-SV) ","");
  PhiRotationSignificanceUnPhysicalTaus=HConfig.GetTH1D(Name+"_PhiRotationSignificanceUnPhysicalTaus","Significance of #phi angle rotation",20,0,10,"#sigma(#phi)","");


  csumProb2Dim=HConfig.GetTH2D(Name+"_csumProb2Dim","csumProb2Dim",200,0,0.5,200,0,0.2,"csum, GeV","probability");
  ProbabilityOfCorrectVsZPt=HConfig.GetTH2D(Name+"_ProbabilityOfCorrectVsZPt","ProbabilityOfCorrectVsZPt",400,0,0.5,100,0,50,"probability","Z pt, GeV");
  ProbabilityOfCorrectVsChi2=HConfig.GetTH2D(Name+"_ProbabilityOfCorrectVsChi2","ProbabilityOfCorrectVsChi2",400,0,1,100,0,50,"probability","#chi^{2}");


  ChannelEfficiency=HConfig.GetTH1D(Name+"_ChannelEfficiency","ChannelEfficiency",4,0.5,4.5,"  0 - signal, 1 - 4#pi, 2 - K#pi#pi, 3 - KK#pi ","");


  TauMuPxPull=HConfig.GetTH1D(Name+"_TauMuPxPull","TauMuPxPull",50,-100,100,"TauMu Px, pull ","");
  TauMuPyPull=HConfig.GetTH1D(Name+"_TauMuPyPull","TauMuPyPull",50,-100,100,"TauMu Py, pull ","");
  TauMuPzPull=HConfig.GetTH1D(Name+"_TauMuPzPull","TauMuPzPull",50,-100,100,"TauMu Pz, pull ","");

  TauA1PxPull=HConfig.GetTH1D(Name+"_TauA1PxPull","TauA1PxPull",200,-50,50,"TauA1 Px, pull ","");
  TauA1PyPull=HConfig.GetTH1D(Name+"_TauA1PyPull","TauA1PyPull",200,-50,50,"TauA1 Py, pull ","");
  TauA1PzPull=HConfig.GetTH1D(Name+"_TauA1PzPull","TauA1PzPull",200,-50,50,"TauA1 Pz, pull ","");

  TauA1PMinusPull=HConfig.GetTH1D(Name+"_TauA1PMinusPull","TauA1PMinusPull",200,-50,50,"TauA1 P Minus, pull ","");
  TauA1PPlusPull=HConfig.GetTH1D(Name+"_TauA1PPlusPull","TauA1PPlusPull",200,-50,50,"TauA1 P Plus, pull ","");
  TauA1PAmbigaPull=HConfig.GetTH1D(Name+"_TauA1PAmbigaPull","TauA1PAmbigaPull",200,-50,50,"TauA1 P Ambiga, pull ","");

  TauA1PMinusPPropagatorPull=HConfig.GetTH1D(Name+"_TauA1PMinusPPropagatorPull","TauA1PMinusPPropagatorPull",200,-50,50,"TauA1 P Plus, pull ","");
  TauA1PPlusPPropagatorPull=HConfig.GetTH1D(Name+"_TauA1PPlusPPropagatorPull","TauA1PPlusPPropagatorPull",200,-50,50,"TauA1 P Plus, pull ","");

  TestPropagation=HConfig.GetTH1D(Name+"_TestPropagation","TestPropagation",200,-50,50,"TauA1 P Ambiga, pull ","");

  TauA1SFPxPull=HConfig.GetTH1D(Name+"_TauA1SFPxPull","TauA1SFPxPull",50,-100,100,"TauA1 Px, pull ","");
  TauA1SFPyPull=HConfig.GetTH1D(Name+"_TauA1SFPyPull","TauA1SFPyPull",50,-100,100,"TauA1 Py, pull ","");
  TauA1SFPzPull=HConfig.GetTH1D(Name+"_TauA1SFPzPull","TauA1SFPzPull",50,-100,100,"TauA1 Pz, pull ","");

//   TauA1PPlusPull=HConfig.GetTH1D(Name+"_TauA1PPlusPull","TauA1PPlusPull",50,-100,100,"TauA1 P Plus, pull ","");
//   TauA1PMinusPull=HConfig.GetTH1D(Name+"_TauA1PMinusPull","TauA1PMinusPull",50,-100,100,"TauA1 P Minus, pull ","");



  SVxPull=HConfig.GetTH1D(Name+"_SVxPull","SVxPull",50,-5,5,"SVx, pull ","");
  SVyPull=HConfig.GetTH1D(Name+"_SVyPull","SVyPull",50,-5,5,"SVy, pull ","");
  SVzPull=HConfig.GetTH1D(Name+"_SVzPull","SVzPull",50,-5,5,"SVz, pull ","");

  PVxPull=HConfig.GetTH1D(Name+"_PVxPull","PVxPull",50,-5,5,"PVx, pull ","");
  PVyPull=HConfig.GetTH1D(Name+"_PVyPull","PVyPull",50,-5,5,"PVy, pull ","");
  PVzPull=HConfig.GetTH1D(Name+"_PVzPull","PVzPull",50,-5,5,"PVz, pull ","");



  A1PxPull=HConfig.GetTH1D(Name+"_A1PxPull","A1PxPull",50,-5,5,"A1 Px, pull ","");
  A1PyPull=HConfig.GetTH1D(Name+"_A1PyPull","A1PyPull",50,-5,5,"A1 Py, pull ","");
  A1PzPull=HConfig.GetTH1D(Name+"_A1PzPull","A1PzPull",50,-5,5,"A1 Pz, pull ","");
  A1PPull=HConfig.GetTH1D(Name+"_A1PPull","A1PPull",50,-5,5,"A1 P, pull ","");


  GJReco=HConfig.GetTH1D(Name+"_GJReco","GJReco",100,-1,1,"GJ, reco ","");  
  dGJReco=HConfig.GetTH1D(Name+"_dGJReco","dGJReco",100,-1,1,"GJ, reco ","");  
  GJRecoZoom=HConfig.GetTH1D(Name+"_GJRecoZoom","GJRecoZoom",100,0.8,1,"GJ, reco ","");  
  GJTruth=HConfig.GetTH1D(Name+"_GJTruth","GJTruth",100,0.8,1,"GJ, Truth ","");  

  GJResol=HConfig.GetTH1D(Name+"_GJResol","GJResol",100,-1,1,"GJ, resol ","");  
  GJPull=HConfig.GetTH1D(Name+"_GJPull","GJPull",100,-0.5,0.5,"GJ, pull ","");  
  GJPull2=HConfig.GetTH1D(Name+"_GJPull2","GJPull2",100,-0.5,0.5,"GJ2, pull ","");

  pullpx1px2=HConfig.GetTH1D(Name+"_pullpx1px2","pullpx1px2",100,-5,5,"px1px2, pull ","");
  pullsvxpvx=HConfig.GetTH1D(Name+"_pullsvxpvx","pullsvxpvx",100,-5,5,"pullsvxpvx, pull ","");


  pullpx1=HConfig.GetTH1D(Name+"_pullpx1","pullpx1",100,-5,5,"px1, pull ","");
  pullpx2=HConfig.GetTH1D(Name+"_pullpx2","pullpx2",100,-5,5,"px2, pull ","");


  pulltaudirx=HConfig.GetTH1D(Name+"_pulltaudirx","pulltaudirx",100,-5,5,"x tau dir, pull ","");
  pulltaudiry=HConfig.GetTH1D(Name+"_pulltaudiry","pulltaudiry",100,-5,5,"y tau dir, pull ","");
  pulltaudirz=HConfig.GetTH1D(Name+"_pulltaudirz","pulltaudirz",100,-5,5,"z tau dir, pull ","");

  NormedXpulltaudir=HConfig.GetTH1D(Name+"_NormedXpulltaudir","NormedXpulltaudir",100,-5,5,"xx tau dir, pull ","");
  NormedYpulltaudir=HConfig.GetTH1D(Name+"_NormedYpulltaudir","NormedYpulltaudir",100,-5,5,"yy tau dir, pull ","");
  NormedZpulltaudir=HConfig.GetTH1D(Name+"_NormedZpulltaudir","NormedZpulltaudir",100,-5,5,"zz tau dir, pull ","");

  pulltaudirmag=HConfig.GetTH1D(Name+"_pulltaudirmag","pulltaudirmag",100,-5,5,"mag tau dir, pull ","");


  pullscalarproductz=HConfig.GetTH1D(Name+"_pullscalarproductz","pullscalarproductz",100,-5,5,", pull ","");

  pullsvxaddpvx=HConfig.GetTH1D(Name+"_pullsvxaddpvx","pullsvxaddpvx",100,-5,5,"pullsvxaddpvx, pull ","");


  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    //std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
  //  ran = new TRandom();
}



void  EFValidation::Store_ExtraDist(){


 Extradist1d.push_back(&TruthA1Pt);
 Extradist1d.push_back(&TruthA1Phi);
 Extradist1d.push_back(&TruthA1Eta);
 
 Extradist1d.push_back(&TruthMuPt);
 Extradist1d.push_back(&TruthMuPhi);
 Extradist1d.push_back(&TruthMuEta);

 Extradist1d.push_back(&TruthTauMuPt);
 Extradist1d.push_back(&TruthTauMuPhi);
 Extradist1d.push_back(&TruhTauMuEta);
 
 Extradist1d.push_back(&TruthTauA1Pt);
 Extradist1d.push_back(&TruthTauA1Phi);
 Extradist1d.push_back(&TruthTauA1Eta);
 
 
 Extradist1d.push_back(&TauA1PtResolKFM);
 Extradist1d.push_back(&TauA1PhiResolKFM);
 Extradist1d.push_back(&TauA1EtaResolKFM);
 
 
 Extradist1d.push_back(&TauA1PtResolKFP);
 Extradist1d.push_back(&TauA1PhiResolKFP);
 Extradist1d.push_back(&TauA1EtaResolKFP);
 
 
 Extradist1d.push_back(&RecoA1Pt);
 Extradist1d.push_back(&RecoA1Phi);
 Extradist1d.push_back(&RecoA1Eta);
 
 Extradist1d.push_back(&RecoMuPt);
 Extradist1d.push_back(&RecoMuPhi);
 Extradist1d.push_back(&RecoMuEta);
 
 Extradist1d.push_back(&A1Mass);
 
 Extradist1d.push_back(&TruthA1PtAfterFit);
 Extradist1d.push_back(&TruthA1EAfterFit);
 Extradist1d.push_back(&TruthA1PhiAfterFit);
 Extradist1d.push_back(&TruthA1EtaAfterFit);
 
 Extradist1d.push_back(&TruthMuPtAfterFit);
 Extradist1d.push_back(&TruthMuEAfterFit);
 Extradist1d.push_back(&TruthMuPhiAfterFit);
 Extradist1d.push_back(&TruthMuEtaAfterFit);
 
 
 Extradist1d.push_back(&EFitTauA1Pt);
 Extradist1d.push_back(&EFitTauA1Phi);
 Extradist1d.push_back(&EFitTauA1Eta);

 Extradist1d.push_back(&EFitTauMuPt);
 Extradist1d.push_back(&EFitTauMuPhi);
 Extradist1d.push_back(&EFitTauMuEta);



 Extradist1d.push_back(&EFitTauA1PtAmbPoint);
 Extradist1d.push_back(&EFitTauA1PhiAmbPoint);
 Extradist1d.push_back(&EFitTauA1EtaAmbPoint);

 Extradist1d.push_back(&EFitTauMuPtAmbPoint);
 Extradist1d.push_back(&EFitTauMuPhiAmbPoint);
 Extradist1d.push_back(&EFitTauMuEtaAmbPoint);

 
 Extradist1d.push_back(&TauA1PtResolution);
 Extradist1d.push_back(&TauA1EResolution);
 Extradist1d.push_back(&TauA1PhiResolution);
 Extradist1d.push_back(&TauA1EtaResolution);
 
 Extradist1d.push_back(&TauMuPtResolution);
 Extradist1d.push_back(&TauMuEResolution);
 Extradist1d.push_back(&TauMuPhiResolution);
 Extradist1d.push_back(&TauMuEtaResolution);
 
 
 Extradist2d.push_back(&TauA1PhiResolVsProb);
 Extradist2d.push_back(&TauMuPhiResolVsProb);
 Extradist2d.push_back(&TauA1EtaResolVsProb);
 Extradist2d.push_back(&TauMuEtaResolVsProb);
 Extradist2d.push_back(&TauA1PtResolVsProb);
 Extradist2d.push_back(&TauMuPtResolVsProb);
 Extradist2d.push_back(&TauA1EResolVsProb);
 Extradist2d.push_back(&TauMuEResolVsProb);
      
 Extradist2d.push_back(&TauA1PtResolVsPt);
 Extradist2d.push_back(&TauMuPtResolVsPt);
 Extradist2d.push_back(&TauA1PtResolVsEta);
 Extradist2d.push_back(&TauMuPtResolVsEta);
 Extradist2d.push_back(&TauA1EtaResolVsEta);
 Extradist2d.push_back(&TauMuEtaResolVsEta);
 
 
 Extradist2d.push_back(&TauA1PtResolVsZPt);
 Extradist2d.push_back(&TauMuPtResolVsZPt);
 
 Extradist2d.push_back(&TauA1EtaResolVsZPt);
 Extradist2d.push_back(&TauMuEtaResolVsZPt);
 
 

 
 Extradist2d.push_back(&TauA1PhiResolVsPVSVSignificance);
 Extradist2d.push_back(&TauA1EtaResolVsPVSVSignificance);
 Extradist2d.push_back(&TauA1PtResolVsPVSVSignificance);
 Extradist2d.push_back(&TauA1EResolVsPVSVSignificance);
 
 Extradist2d.push_back(&TauA1PtResolVsPtAmbPoint);
 Extradist2d.push_back(&TauMuPtResolVsPtAmbPoint);
 Extradist2d.push_back(&TauA1PtResolVsEtaAmbPoint);
 Extradist2d.push_back(&TauMuPtResolVsEtaAmbPoint);
 Extradist2d.push_back(&TauA1EtaResolVsEtaAmbPoint);
 Extradist2d.push_back(&TauMuEtaResolVsEtaAmbPoint);
 
 Extradist2d.push_back(&TauA1PtResolVsZPtAmbPoint);
 Extradist2d.push_back(&TauMuPtResolVsZPtAmbPoint);
 
 Extradist2d.push_back(&TauA1EtaResolVsZPtAmbPoint);
 Extradist2d.push_back(&TauMuEtaResolVsZPtAmbPoint);
 
 Extradist1d.push_back(&TauA1PtResolutionAmbPoint);
 Extradist1d.push_back(&TauA1PhiResolutionAmbPoint);
 Extradist1d.push_back(&TauA1EtaResolutionAmbPoint);
 
 Extradist1d.push_back(&TauMuPtResolutionAmbPoint);
 Extradist1d.push_back(&TauMuPhiResolutionAmbPoint);
 Extradist1d.push_back(&TauMuEtaResolutionAmbPoint);
 
 Extradist2d.push_back(&TauA1PhiResolVsPhiRotSignificance);
 Extradist2d.push_back(&TauA1EtaResolVsPhiRotSignificance);
 Extradist2d.push_back(&TauA1PtResolVsPhiRotSignificance);
 Extradist2d.push_back(&TauA1EResolVsPhiRotSignificance);
 
 
 Extradist2d.push_back(&TauA1PhiResolVsPVSVSignificanceAmbPoint);
 Extradist2d.push_back(&TauA1EtaResolVsPVSVSignificanceAmbPoint);
 Extradist2d.push_back(&TauA1PtResolVsPVSVSignificanceAmbPoint);
 Extradist2d.push_back(&TauA1EResolVsPVSVSignificanceAmbPoint);
 
 Extradist2d.push_back(&Chi2Dim);
 Extradist2d.push_back(&csum2Dim);
 Extradist1d.push_back(&ProbabilityOfCorrect);


 Extradist1d.push_back(&idAllTaus);
 Extradist1d.push_back(&idIdentifiedTau);
 Extradist1d.push_back(&idPassedEF);
 Extradist1d.push_back(&CorrectAmbiguityTruth);
 Extradist1d.push_back(&MeasuredPhi);
 Extradist1d.push_back(&pvsvsignificance);

 Extradist1d.push_back(&PhiRotationSignificanceUnPhysicalTaus);

 Extradist2d.push_back(&csumProb2Dim);
 Extradist2d.push_back(&ProbabilityOfCorrectVsZPt);
 Extradist1d.push_back(&ProbabilityOfAmbiguityPoint);



 Extradist1d.push_back(&TauA1PtResolVsPVSVSignificanceRMS);
 Extradist1d.push_back(&TauA1PtResolVsPtRMS);
 Extradist1d.push_back(&TauMuPtResolVsZPtRMS);
 Extradist1d.push_back(&TauA1PtResolVsEtaRMS);
 Extradist1d.push_back(&TauMuPtResolVsPtRMS);
 Extradist1d.push_back(&TauMuPtResolVsZPtMean);

 Extradist1d.push_back(&EfficiencyOverA1Pt);
 Extradist1d.push_back(&EfficiencyOverA1Phi);
 Extradist1d.push_back(&EfficiencyOverA1Eta);


 Extradist1d.push_back(&EfficiencyOverMuPt);
 Extradist1d.push_back(&EfficiencyOverMuPhi);
 Extradist1d.push_back(&EfficiencyOverMuEta);

 Extradist1d.push_back(&ProbabilityOfCorrectndf2);
 Extradist1d.push_back(&ProbabilityOfCorrectndf3);
 Extradist1d.push_back(&Chi2OfCorrect);


 Extradist2d.push_back(&ProbabilityOfCorrectVsChi2);
 Extradist1d.push_back(&PhiRotSignificanceCutEfficiency);
 Extradist1d.push_back(&PVSVSignificanceCutEfficiency);


 Extradist1d.push_back(&ChannelEfficiency);


 Extradist1d.push_back(&TauMuPxPull);
 Extradist1d.push_back(&TauMuPyPull);
 Extradist1d.push_back(&TauMuPzPull);

 Extradist1d.push_back(&TauA1PxPull);
 Extradist1d.push_back(&TauA1PyPull);
 Extradist1d.push_back(&TauA1PzPull);
 Extradist1d.push_back(&TauA1PMinusPull);
 Extradist1d.push_back(&TauA1PPlusPull);
 Extradist1d.push_back(&TauA1PAmbigaPull);


 Extradist1d.push_back(&TauA1PMinusPPropagatorPull);
 Extradist1d.push_back(&TauA1PPlusPPropagatorPull);
 Extradist1d.push_back(&A1PxPull);
 Extradist1d.push_back(&A1PyPull);
 Extradist1d.push_back(&A1PzPull);
 Extradist1d.push_back(&A1PPull);

 Extradist1d.push_back(&TestPropagation);

 Extradist1d.push_back(&TauA1SFPxPull);
 Extradist1d.push_back(&TauA1SFPyPull);
 Extradist1d.push_back(&TauA1SFPzPull);

//  Extradist1d.push_back(&TauA1PMinusPull);
//  Extradist1d.push_back(&TauA1PPlusPull);



 Extradist1d.push_back(&SVxPull);
 Extradist1d.push_back(&SVyPull);
 Extradist1d.push_back(&SVzPull);


 Extradist1d.push_back(&PVxPull);
 Extradist1d.push_back(&PVyPull);
 Extradist1d.push_back(&PVzPull);


 Extradist1d.push_back(&GJReco); 
 Extradist1d.push_back(&dGJReco); 
 Extradist1d.push_back(&GJTruth);


 Extradist1d.push_back(&GJResol);
 Extradist1d.push_back(&GJPull);
 Extradist1d.push_back(&GJPull2);

 Extradist1d.push_back(&pullpx1px2);
 Extradist1d.push_back(&pullsvxpvx);



 Extradist1d.push_back(&pullpx1);
 Extradist1d.push_back(&pullpx2);


 Extradist1d.push_back(&pullsvxaddpvx);

 Extradist1d.push_back(&pulltaudirx);
 Extradist1d.push_back(&pulltaudiry);
 Extradist1d.push_back(&pulltaudirz);


 Extradist1d.push_back(&NormedXpulltaudir);
 Extradist1d.push_back(&NormedYpulltaudir);
 Extradist1d.push_back(&NormedZpulltaudir);


 Extradist1d.push_back(&pulltaudirmag);

 Extradist1d.push_back(&pullscalarproductz);

}

void  EFValidation::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
    std::cout << "EFValidation Ntp->GetMCID(): " << Ntp->GetMCID() << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){/* std::cout << "failed to find id" <<std::endl;*/ return;}

  if(verbose)std::cout << "void  EFValidation::doEvent() Mu" << std::endl;


//   // check out triggers
//   for(unsigned int itr = 0; itr <  Ntp->NHLTTriggers(); itr ++){
//     std::cout << "Trigger names  " <<Ntp->HTLTriggerName(itr) << std::endl;
//   }

  int idToFill(0);
  if(id ==998)idToFill = 1;
  if(id ==10230833)idToFill = 2;
  if(id ==10231433)idToFill = 3;
  if(id ==10231833)idToFill = 4;
  idAllTaus.at(t).Fill(idToFill,1);


  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("eta2p1_LooseIsoPFTau"))value.at(TriggerOk)+=1;
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);

//   for(unsigned int iTrig =0; iTrig <Ntp-> NHLTTriggers(); iTrig ++){

//     std::cout<<" Triggers  "<< Ntp->HTLTriggerName(iTrig)<<std::endl;
//     }
  // Apply Selection
  std::vector<unsigned int> mu_idx_good, mu_idx_pt, mu_idx_iso;
  unsigned mu_idx(999);
  //double mu_pt(0);
  for(unsigned int i=0;i<Ntp->NMuons();i++){
    //std::cout << "Checking if muon number = " << i << " is good..." << std::endl;
    if(Ntp->isGoodMuon(i) && fabs(Ntp->Muons_p4(i).Eta())<2.4 && Ntp->Muons_p4(i).Pt() > 17){
      mu_idx_good.push_back(i);
    }  
  }



  value.at(hasMuon)=mu_idx_good.size();
  pass.at(hasMuon)=(value.at(hasMuon)>cut.at(hasMuon));
  
  std::vector<unsigned int> tau_idx_pt, tau_idx_eta, tau_idx_iso, tau_idx_good ;
  unsigned int tau_idx(999);
  for(unsigned int i=0;i<Ntp->NPFTaus();i++){
    if(Ntp->PFTau_p4(i).Pt()> 17){
      if(fabs(Ntp->PFTau_p4(i).Eta())<2.7){
	if(Ntp->PFTau_hpsDecayMode(i)==10 && Ntp->PFTau_isHPSByDecayModeFinding(i)  ){ 
	//	if(Ntp->PFTau_hpsDecayMode(i)==10 && Ntp->PFTau_isHPSByDecayModeFinding(i)  ){ 
	  tau_idx_good.push_back(i);
	}
      }
    }
  }



  value.at(hasTau)=tau_idx_good.size();
  pass.at(hasTau)=(value.at(hasTau)>cut.at(hasTau));

  std::cout<<"ntaus  "<<tau_idx_good.size()<<std::endl;
  std::cout<<"nmuons  "<<mu_idx_good.size()<<std::endl;

   tau_idx =   tau_idx_good.at(0);
   mu_idx=   mu_idx_good.at(0);



  double wobs(1),w(1),w1(1);
   if(!Ntp->isData()){
     w1*=Ntp->EvtWeight3D();
   }
   else{w1=1;}

  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
    if(verbose)std::cout << "MC type: " << Ntp->GetMCID() <<std::endl;
    //std::cout << "tau_idx: " << tau_idx << "; mu_idx: " << mu_idx << std::endl;
    //std::cout << "tau_idx_iso.size: " << tau_idx_iso.size() << " mu_idx_iso.size: " << mu_idx_iso.size() << std::endl;
    int isInPhysicalRegion(0);
    double significan(0);
    if (tau_idx!=999 && mu_idx!=999){
      //>>>>>>>>> search for a matched jet 
      unsigned int matchedJet(0);
      for(unsigned int ijet =0; ijet < Ntp->NPFJets(); ijet++){
	
	double delta = 999.;
	if(Ntp->PFTau_p4(tau_idx).DeltaR(Ntp->PFJet_p4(ijet)) < delta ){
	  delta = Ntp->PFTau_p4(tau_idx).DeltaR(Ntp->PFJet_p4(ijet));
	  matchedJet=ijet;
	}
      }
      //<<<<<<<<< search for a matched jet


      //>>>>>>>>>>>>>>>>>>>>>>>Fill histos for efficiency 
  
      if(Ntp->CheckDecayID(2,5)){
	std::cout<<" 1 "<<std::endl;
	TLorentzVector TruthA1 = Ntp->GetTruthTauProductLV(5, 20213);
	TLorentzVector TruthMu = Ntp->GetTruthTauProductLV(2, 13);
	
	TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
    
	TruthA1Pt.at(t).Fill(TruthA1.Pt(),1);
	TruthA1Phi.at(t).Fill(TruthA1.Phi(),1);
	TruthA1Eta.at(t).Fill(TruthA1.Eta(),1);

	TruthMuPt.at(t).Fill(TruthMu.Pt(),1);
	TruthMuPhi.at(t).Fill(TruthMu.Phi(),1);
	TruthMuEta.at(t).Fill(TruthMu.Eta(),1);

	TruthTauMuPt.at(t).Fill(TruthMu.Pt(),1);
	TruthTauMuPhi.at(t).Fill(TruthMu.Phi(),1);
	TruhTauMuEta.at(t).Fill(TruthMu.Eta(),1);
	
	TruthTauA1Pt.at(t).Fill(TruthA1.Pt(),1);
	TruthTauA1Phi.at(t).Fill(TruthA1.Phi(),1);
	TruthTauA1Eta.at(t).Fill(TruthA1.Eta(),1);

	double angle = TruthTauA1.Angle(TruthA1.Vect());
	double rootAmb=sqrt((TruthA1.M()*TruthA1.M()+TruthA1.P()*TruthA1.P())*(pow(TruthA1.M()*TruthA1.M()-TruthTauA1.M()*TruthTauA1.M(),2)-4*TruthTauA1.M()*TruthTauA1.M()*TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle)));
	double PtruthMinus  = ((TruthA1.M()*TruthA1.M()+TruthTauA1.M()*TruthTauA1.M() )*TruthA1.P()*cos(angle)  -rootAmb )/(TruthA1.M()*TruthA1.M() + TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle))/2;
	double PtruthPlus  = ((TruthA1.M()*TruthA1.M()+TruthTauA1.M()*TruthTauA1.M() )*TruthA1.P()*cos(angle)  +rootAmb )/(TruthA1.M()*TruthA1.M() + TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle))/2;
	
	int TruthAmibiga =0;
	
	if( fabs(TruthTauA1.P()  -PtruthMinus ) < fabs(TruthTauA1.P()  - PtruthPlus))TruthAmibiga =1;
	if( fabs(TruthTauA1.P()  -PtruthMinus ) > fabs(TruthTauA1.P()  - PtruthPlus))TruthAmibiga =2;
	CorrectAmbiguityTruth.at(t).Fill(TruthAmibiga,1.);
      }
  

      //<<<<<<<<<<<<<<<<<<<<<<<Fill histos for efficiency 
      
      int idToFill(0);
      if(id ==998)idToFill = 1;
      if(id ==10230833)idToFill = 2;
      if(id ==10231433)idToFill = 3;
      if(id ==10231833)idToFill = 4;
      if(Ntp->PFTau_hpsDecayMode(tau_idx) == 10)idIdentifiedTau.at(t).Fill(idToFill,1);

      std::vector<bool> Tau_FitOk; 
      std::vector<double> Tau_FitChi2;
      std::vector<TLorentzVector> TauA1_TPF;


      std::vector<bool> EventFit_Ok;
      std::vector<double> SignificanceOfTauRotation;

      std::vector<double> DeltaPhiSign;

      std::vector<double> deltaphiRotation;

      std::vector<TLorentzVector> TauA1_EF;
      std::vector<TLorentzVector> TauMu_EF;
      std::vector<TLorentzVector> TauMuNetrino_EF;
      std::vector<TLorentzVector> TauA1Neutrino_TPF;

      std::vector<TLorentzVector> Z_EF;

      std::vector<double> Chi2EventFit;
      std::vector<double> Chi2ProbabilityEventFit;
      std::vector<double> Chi2ProbabilityEventFitndf2;
      std::vector<double> Chi2ProbabilityEventFitndf3;


      std::vector<double> NiterationsEventFit;
      std::vector<double> csum_GF;

      double flightsignificanceA1Side;
      double flightsignificanceMuSide;
      std::vector<LorentzVectorParticle> theTauAmbiguityVector;

      std::vector<LorentzVectorParticle> theTauA1TLV;
      std::vector<LorentzVectorParticle> theTauMuTLV;

      TVector3 pv;
      TMatrixTSym<double> PVcov;
      TVector3 sv;
      TMatrixTSym<double> SVcov;
    
      TLorentzVector A1Reco = Ntp->PFTau_a1_lvp(tau_idx).LV();
      TLorentzVector MuReco = Ntp->Muons_p4(mu_idx);

      RecoA1Pt.at(t).Fill(A1Reco.Pt(),1);
      RecoA1Phi.at(t).Fill(A1Reco.Phi(),1);
      RecoA1Eta.at(t).Fill(A1Reco.Eta(),1);

      RecoMuPt.at(t).Fill(MuReco.Pt(),1);
      RecoMuPhi.at(t).Fill(MuReco.Phi(),1);
      RecoMuEta.at(t).Fill(MuReco.Eta(),1);

      A1Mass.at(t).Fill(A1Reco.M(),1);

      //////////////////////////////
      //loop over ambiguity points
      MeasuredPhi.at(t).Fill(Ntp->MeasuredTauDirection(tau_idx).Phi(),1);

      for(unsigned int j=0;j< MultiProngTauSolver::NAmbiguity;j++){

	LorentzVectorParticle theTau;
	LorentzVectorParticle theZ;
	std::vector<LorentzVectorParticle> daughter;
	std::vector<LorentzVectorParticle> theZdaughter;

	LorentzVectorParticle  LVPEF_TauMu;
	LorentzVectorParticle  LVPEF_TauA1;
	       
	TLorentzVector NeutrinoA1(0,0,0,0);
	TLorentzVector NeutrinoMu(0,0,0,0);
	       
	TLorentzVector TauA1EventFit(0,0,0,0);
	TLorentzVector TauA1ThreeProngFit(0,0,0,0);
	TLorentzVector TauMuEventFit(0,0,0,0);
	TLorentzVector Z_sol(0,0,0,0);
	double deltaphi(0);
	       
	double LC_Eventchi2(0);
	double LC_Eventchi2Probability(0);
	double LC_Eventchi2Probabilityndf2(0);
	double LC_Eventchi2Probabilityndf3(0);


	double LC_chi2(0);
	double phisign(0);
	int NiterationsEF(0);
	double csum(0);
	double SEA1(0);
	double SEMU(0);

	bool  ThreeProngFitSuccess =false;
	bool  EventFitSuccess =false;
	ThreeProngFitSuccess=Ntp->ThreeProngTauFit(tau_idx,j,theTau,daughter,LC_chi2,phisign);
	MeasuredPhi.at(t).Fill(Ntp->MeasuredTauDirection(tau_idx).Phi(),1);
	if(ThreeProngFitSuccess){

	  //>>>>>>>>>>> Fill Resol plots dor TPF
	  if(Ntp->CheckDecayID(2,5)){

	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);

	    TLorentzVector TruthA1 = Ntp->GetTruthTauProductLV(5, 20213);
	    TLorentzVector TruthMu = Ntp->GetTruthTauProductLV(2, 13);

	    double angle = TruthTauA1.Angle(TruthA1.Vect());
	    double rootAmb=sqrt((TruthA1.M()*TruthA1.M()+TruthA1.P()*TruthA1.P())*(pow(TruthA1.M()*TruthA1.M()-TruthTauA1.M()*TruthTauA1.M(),2)-4*TruthTauA1.M()*TruthTauA1.M()*TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle)));
	    double PtruthMinus  = ((TruthA1.M()*TruthA1.M()+TruthTauA1.M()*TruthTauA1.M() )*TruthA1.P()*cos(angle)  -rootAmb )/(TruthA1.M()*TruthA1.M() + TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle))/2;
	    double PtruthPlus  = ((TruthA1.M()*TruthA1.M()+TruthTauA1.M()*TruthTauA1.M() )*TruthA1.P()*cos(angle)  +rootAmb )/(TruthA1.M()*TruthA1.M() + TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle))/2;
	    //  PlotSmallParameter.at(t).Fill(4*1.777*1.777*TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle)/(1.777*1.777 - TruthA1.M()*TruthA1.M())/(1.777*1.777 - TruthA1.M()*TruthA1.M()),1);
	    int TruthAmibiga =0;
     
	    std::cout<<"---- small param "<<4*TruthTauA1.M()*TruthTauA1.M()*TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle)/pow(TruthA1.M()*TruthA1.M()-TruthTauA1.M()*TruthTauA1.M(),2) <<std::endl;

	    if( fabs(TruthTauA1.P()  -PtruthMinus ) < fabs(TruthTauA1.P()  - PtruthPlus))TruthAmibiga =1;
	    if( fabs(TruthTauA1.P()  -PtruthMinus ) > fabs(TruthTauA1.P()  - PtruthPlus))TruthAmibiga =2;




	    if(j==1){
	      TauA1PtResolKFM.at(t).Fill(TruthTauA1.Pt() - theTau.LV().Pt(),1);
	      TauA1PhiResolKFM.at(t).Fill(TruthTauA1.Phi() - theTau.LV().Phi(),1);
	      TauA1EtaResolKFM.at(t).Fill(TruthTauA1.Eta() - theTau.LV().Eta(),1);

	      if(TruthAmibiga==1){
		double deltaPTau = sqrt(pow(theTau.LV().Px(),2) * theTau.Covariance(3,3)  +pow(theTau.LV().Py(),2) * theTau.Covariance(4,4) + pow(theTau.LV().Pz(),2) * theTau.Covariance(5,5) )/theTau.LV().P();
		TauA1PMinusPull.at(t).Fill((PtruthMinus - theTau.LV().P() )/deltaPTau,1);



	      }
	    }
	    if(j==2){
	      TauA1PtResolKFP.at(t).Fill(TruthTauA1.Pt() -  theTau.LV().Pt(),1);
	      TauA1PhiResolKFP.at(t).Fill(TruthTauA1.Phi() - theTau.LV().Phi(),1);
	      TauA1EtaResolKFP.at(t).Fill(TruthTauA1.Eta() - theTau.LV().Eta(),1);

	      if(TruthAmibiga==2){
		double deltaPTau = sqrt(pow(theTau.LV().Px(),2) * theTau.Covariance(3,3)  +pow(theTau.LV().Py(),2) * theTau.Covariance(4,4) + pow(theTau.LV().Pz(),2) * theTau.Covariance(5,5) )/theTau.LV().P();
		TauA1PPlusPull.at(t).Fill((PtruthPlus - theTau.LV().P() )/deltaPTau,1);


	      }
	    }





	    TVector3 pv=Ntp->PFTau_TIP_primaryVertex_pos(tau_idx);
	    TMatrixTSym<double> pvcov=Ntp->PFTau_TIP_primaryVertex_cov(tau_idx);

	    TVector3 TruthTauA1SV1 = Ntp->GetTruthTauProductVertex(5,211);
	   
	    SVxPull.at(t).Fill((TruthTauA1SV1.X() - theTau.Vertex().X() )/sqrt(theTau.Covariance(0,0)*1.26*1.26),1);
	    SVyPull.at(t).Fill((TruthTauA1SV1.Y() - theTau.Vertex().Y() )/sqrt(theTau.Covariance(1,1)),1);
	    SVzPull.at(t).Fill((TruthTauA1SV1.Z() - theTau.Vertex().Z() )/sqrt(theTau.Covariance(2,2)),1);

	    if(sqrt(theTau.Covariance(3,3)) > 0  ) TauA1PxPull.at(t).Fill((TruthTauA1.Px() - theTau.LV().Px() )/sqrt(theTau.Covariance(3,3)),1);
	    if(sqrt(theTau.Covariance(4,4)) > 0  ) TauA1PyPull.at(t).Fill((TruthTauA1.Py() - theTau.LV().Py() )/sqrt(theTau.Covariance(4,4)),1);
	    if(sqrt(theTau.Covariance(5,5)) > 0  ) TauA1PzPull.at(t).Fill((TruthTauA1.Pz() - theTau.LV().Pz() )/sqrt(theTau.Covariance(5,5)),1);


	    TLorentzVector TruthA1LV = Ntp->GetTruthTauProductLV(5, 20213);
	    TLorentzVector TruthMuLV = Ntp->GetTruthTauProductLV(2, 13);

	    if(sqrt(daughter.at(0).Covariance(3,3)) > 0  ) A1PxPull.at(t).Fill((TruthA1LV.Px() - daughter.at(0).LV().Px() )/sqrt(daughter.at(0).Covariance(3,3)),1);
	    if(sqrt(daughter.at(0).Covariance(4,4)) > 0  ) A1PyPull.at(t).Fill((TruthA1LV.Py() - daughter.at(0).LV().Py() )/sqrt(daughter.at(0).Covariance(4,4)),1);
	    if(sqrt(daughter.at(0).Covariance(5,5)) > 0  ) A1PzPull.at(t).Fill((TruthA1LV.Pz() - daughter.at(0).LV().Pz() )/sqrt(daughter.at(0).Covariance(5,5)),1);
	    TVector3 PVTruth = Ntp->MCSignalParticle_Poca(0);
	  

	    PVxPull.at(t).Fill((PVTruth.X() - pv.X() )/sqrt(pvcov(0,0)),1);
	    PVyPull.at(t).Fill((PVTruth.Y() - pv.Y() )/sqrt(pvcov(1,1)),1);
	    PVzPull.at(t).Fill((PVTruth.Z() - pv.Z() )/sqrt(pvcov(2,2)),1);

	    TVector3 TauDir = theTau.Vertex()-pv;
	    TVector3 TauDirTruth = TruthTauA1SV1-PVTruth;

	    TVector3 TauDirTruthNormed(TauDirTruth.X()/TauDirTruth.Mag(),TauDirTruth.Y()/TauDirTruth.Mag(),TauDirTruth.Z()/TauDirTruth.Mag());


	    double svxpvx = pv.X()*theTau.Vertex().X();
	    double dsvxpvx = sqrt( pv.X()*pv.x()* theTau.Covariance(0,0)/1.26/1.26 +  theTau.Vertex().Y()*theTau.Vertex().Y()*pvcov(0,0));

	    double svxaddpvx = pv.X() + theTau.Vertex().X();
	    double dsvxaddpvx = sqrt(  theTau.Covariance(0,0)*1.26*1.26 +  pvcov(0,0));


	    pullsvxpvx.at(t).Fill( (svxpvx - PVTruth.X()*TruthTauA1SV1.X())/dsvxpvx,1);
	    pullsvxaddpvx.at(t).Fill( (svxaddpvx - (PVTruth.X() + TruthTauA1SV1.X()))/dsvxaddpvx,1);


	    double x = TauDir.X();
	    double y = TauDir.Y();
	    double z = TauDir.Z();
	    double R = TauDir.Mag();

	    double dx  = (theTau.Covariance(0,0)  + pvcov(0,0));
	    double dy  = (theTau.Covariance(1,1)  + pvcov(1,1));
	    double dz  = (theTau.Covariance(2,2)  + pvcov(2,2));


	    double dR = 1.26*sqrt(x*x*dx + y*y*dy + z*z*dz)/R;



	    double ddx = sqrt(dx/R/R + x*x*dR*dR/R/R/R/R);//sqrt( pow((1 - x*x/R/R)/R,2)*dx  + pow(x*y/R/R/R,2)*dy +  pow(x*z/R/R/R,2)*dz   );
	    double ddy = sqrt( pow((1 - y*y/R/R)/R,2)*dy  + pow(x*y/R/R/R,2)*dx +  pow(y*z/R/R/R,2)*dz   );
	    double ddz = sqrt( pow((1 - z*z/R/R)/R,2)*dz  + pow(x*z/R/R/R,2)*dx +  pow(y*z/R/R/R,2)*dy   );
 



// 	    pulltaudirx.at(t).Fill((x - TauDirTruth.X())/sqrt(dx),1);
// 	    pulltaudiry.at(t).Fill((y - TauDirTruth.Y())/sqrt(dy),1);
// 	    pulltaudirz.at(t).Fill((z - TauDirTruth.Z())/sqrt(dz),1);


	    double dRdx = x/R;
	    double dRdy = y/R;
	    double dRdz = z/R;

	    TVector3 TauDirNormed(x/R,y/R,z/R);
	    TVector3 TauDirError(	  (pow((1/R - x*dRdx/R/R),2)*dx), (pow((1/R - y*dRdy/R/R),2)*dy),(pow((1/R - x*dRdz/R/R),2)*dz) );


	    std::cout<<" x - TauDirTruth.X())/sqrt(dx) "<< (x - TauDirTruth.X())/sqrt(dx)  <<" R*(x/R - TauDirTruth.X()/TauDirTruth.Mag())/sqrt(dx)  "<< (x - TauDirTruth.X()*R/TauDirTruth.Mag())/sqrt(dx) << " (x/R - TauDirTruth.X()/TauDirTruth.Mag())/ddx  "<< (x/R - TauDirTruth.X()/TauDirTruth.Mag())/ddx <<std::endl;

	    std::cout<<" R   "<<  R   <<" TauDirTruth.Mag()  "<< TauDirTruth.Mag() <<  " ratio   " <<R/TauDirTruth.Mag()<<std::endl;
	    std::cout<<" x - TauDirTruth.X()   "<<  x - TauDirTruth.X()   <<" x - TauDirTruth.X()*ratio  "<< (x - TauDirTruth.X()*R/TauDirTruth.Mag()) <<std::endl;




	    pulltaudirx.at(t).Fill((x - TauDirTruth.X())/sqrt(dx),1);
	    pulltaudiry.at(t).Fill((y - TauDirTruth.Y())/sqrt(dy),1);
	    pulltaudirz.at(t).Fill((z - TauDirTruth.Z())/sqrt(dz),1);

	    NormedXpulltaudir.at(t).Fill((x - R*TauDirTruth.X()/TauDirTruth.Mag())/sqrt(dx),1);
	    NormedYpulltaudir.at(t).Fill((y/R - TauDirTruth.Y()/R)/ddy,1);
	    NormedZpulltaudir.at(t).Fill((z/R - TauDirTruth.Z()/R)/ddz,1);


	    pulltaudirmag.at(t).Fill((R - TauDirTruth.Mag())/dR,1);


	    pullscalarproductz.at(t).Fill((x*daughter.at(0).LV().Px() + y*daughter.at(0).LV().Py() + z*daughter.at(0).LV().Pz()  - TauDirTruth.X()*TruthA1LV.Px()-TauDirTruth.Y()*TruthA1LV.Py()-TauDirTruth.Z()*TruthA1LV.Pz())/sqrt(
																       daughter.at(0).LV().Px()*daughter.at(0).LV().Px()*dx + x*x*daughter.at(0).Covariance(3,3) +
																       daughter.at(0).LV().Py()*daughter.at(0).LV().Py()*dy + y*y*daughter.at(0).Covariance(4,4) +
																       daughter.at(0).LV().Pz()*daughter.at(0).LV().Pz()*dz + z*z*daughter.at(0).Covariance(5,5)),1);

	    TLorentzVector RefitA1 = daughter.at(0).LV();


	    double px1 =  RefitA1.Px();
	    double px2 = theTau.LV().Px();

	    double dpx1 =  daughter.at(0).Covariance(3,3);
	    double dpx2 = theTau.Covariance(3,3);

	    double py1 = RefitA1.Py();
	    double py2 = theTau.LV().Py();

	    double dpy1 = daughter.at(0).Covariance(4,4);
	    double dpy2 =theTau.Covariance(4,4);

	    double pz1 = RefitA1.Pz();
	    double pz2 = theTau.LV().Pz();

	    double dpz1 = daughter.at(0).Covariance(5,5);
	    double dpz2 = theTau.Covariance(5,5);

	    double p1 = RefitA1.P();           
	    double p2 = theTau.LV().P();

	    double dp1dpx = px1/p1;
	    double dp1dpy = py1/p1;
	    double dp1dpz = pz1/p1;
	    
	    double dp2dpx = px2/p2;
	    double dp2dpy = py2/p2;
	    double dp2dpz = pz2/p2;

	    TVector3 A1Dir(px1/p1,py1/p1,pz1/p1);
	    TVector3 A1DirError(  (pow((1/p1 - px1*dp1dpx/p1/p1),2)*dpx1),  (pow((1/p1 - py1*dp1dpy/p1/p1),2)*dpy1),  (pow((1/p1 - pz1*dp1dpz/p1/p1),2)*dpz1));
	    TVector3 TruthA1DirNormed(TruthA1LV.Px()/TruthA1LV.P(),TruthA1LV.Py()/TruthA1LV.P(),TruthA1LV.Pz()/TruthA1LV.P());


	    double dcosGJ =  sqrt(pow(px2/p1/p2 - px1*px2*dp1dpx/p1/p1/p2,2)*dpx1*dpx1  + pow(py2/p1/p2 - py1*py2*dp1dpy/p1/p1/p2,2)*dpy1*dpy1   + pow(pz2/p1/p2 - pz1*pz2*dp1dpz/p1/p1/p2,2)*dpz1*dpz1+
				  pow(px1/p1/p2 - px1*px2*dp2dpx/p1/p1/p2,2)*dpx2*dpx2  + pow(py1/p1/p2 - py1*py2*dp2dpy/p1/p1/p2,2)*dpy2*dpy2   + pow(pz1/p1/p2 - pz1*pz2*dp2dpz/p1/p1/p2,2)*dpz2*dpz2);  	
    
	    
	    double acosTruth = (TruthA1LV.Px()*TruthTauA1.Px() + TruthA1LV.Py()*TruthTauA1.Py() + TruthA1LV.Pz()*TruthTauA1.Pz())/TruthA1LV.P()/TruthTauA1.P();
	    double acosReco  = (RefitA1.Px()*theTau.LV().Px() + RefitA1.Py()*theTau.LV().Py() + RefitA1.Pz()*theTau.LV().Pz())/RefitA1.P()/theTau.LV().P();
	  

	    double cosReco2 =  TauDirNormed*A1Dir ;
	    
	    double dcosReco2 = sqrt(A1Dir.X()*A1Dir.X()*TauDirError.X() + TauDirNormed.X()*TauDirNormed.X()*A1DirError.X() +    
				    A1Dir.Y()*A1Dir.Y()*TauDirError.Y() + TauDirNormed.Y()*TauDirNormed.Y()*A1DirError.Y() +    
				    A1Dir.Z()*A1Dir.Z()*TauDirError.Z() + TauDirNormed.Z()*TauDirNormed.Z()*A1DirError.Z());


	    double px1px2 =  TauDirNormed.X()*A1Dir.X();
	    double dpx1px2 = sqrt(A1Dir.X()*A1Dir.X()* TauDirError.X() + TauDirNormed.X()*TauDirNormed.X()*A1DirError.X());


	    pullpx1.at(t).Fill((px1  - TruthA1LV.Px())/dpx1);
	    pullpx2.at(t).Fill((px2  - TruthTauA1.Px())/dpx2);
	    pullpx1px2.at(t).Fill( (px1px2 -TruthA1DirNormed.X()*TauDirTruthNormed.X() )/dpx1px2,1 );
            
	    GJReco.at(t).Fill(cosReco2,1);
	    dGJReco.at(t).Fill(dcosReco2,1);
	    GJRecoZoom.at(t).Fill(cosReco2,1);
	    GJTruth.at(t).Fill(acosTruth,1);

	    

	    GJPull.at(t).Fill((acosReco-acosTruth)/dcosGJ,1);
	    GJPull2.at(t).Fill((cosReco2-acosTruth)/dcosReco2,1);

	    //------------  check propagation of A1 momenta uncertainty

	    

              double angleRec = TauDir.Angle(A1Dir);
	      GJResol.at(t).Fill(angleRec - angle,1);
	      double deltaPa1 = sqrt(pow(daughter.at(0).LV().Px(),2) * daughter.at(0).Covariance(3,3)  +pow(daughter.at(0).LV().Py(),2) * daughter.at(0).Covariance(4,4) + pow(daughter.at(0).LV().Pz(),2) * daughter.at(0).Covariance(5,5) )/daughter.at(0).LV().P();
	      if(j ==0 ) A1PPull.at(t).Fill((TruthA1LV.P() - daughter.at(0).LV().P() )/deltaPa1,1);
	      TLorentzVector a1lv=daughter.at(0).LV();
	      double m1 = 1.777;
	      double m2 = a1lv.M();
	      double sin = TMath::Sin(angleRec);
	      double root = sqrt(  (m2*m2 +  a1lv.P()*a1lv.P()) *( pow(m2*m2-m1*m1,2) - 4*m1*m1*a1lv.P()*a1lv.P()*sin*sin ) );
	      double D = 0.5/(m2*m2 +a1lv.P()*a1lv.P()*sin*sin );

	      double drootdp = (a1lv.P()*pow(m2*m2 - m1*m1,2) - 8*m1*m1 *a1lv.P()*a1lv.P()*a1lv.P()*sin*sin - m1*m1*m2*m2*a1lv.P()*sin*sin)/root;
	      double dDdp = -a1lv.P()*sin*sin/(m2*m2 + a1lv.P()*a1lv.P()*sin*sin);

	      double TauPPlus = D*((m1*m1 + m2*m2 )*a1lv.P()*sqrt(1 - sin*sin) + root);
	      double TauPMinus = D*((m1*m1 + m2*m2 )*a1lv.P()*sqrt(1 - sin*sin) - root);

	      double deltaTauPplus = ((deltaPa1)*(1 - sin*sin)*(m1*m1+m2*m2) *(D  + a1lv.P()*dDdp) + (root * dDdp + D*drootdp));
	      double deltaTauPminus = ((deltaPa1)*(1 - sin*sin)*(m1*m1+m2*m2) *(D  + a1lv.P()*dDdp) - (root * dDdp + D*drootdp));


	      std::cout<<" D "<< D<< "  D2 " << m2*m2 +a1lv.P()*a1lv.P()*sin*sin  <<"  root "<<root<<"  sin " << sin <<"  drootdp  " << drootdp <<" dDdp  " << dDdp << "  angle  "<<angle << "  rootAmb  " <<rootAmb <<std::endl;
	      std::cout<<" TruthA1.P() "<< TruthA1.P() <<"  a1lv.P() "<<a1lv.P() <<"  TruthA1.M() " << TruthA1.M() <<"  drootdp  " << drootdp <<" dDdp  " << dDdp << "  angle  "<<angle << "  rootAmb  " <<rootAmb <<std::endl;

	      std::cout<<"TruthAmibiga  "<<TruthAmibiga <<" PtruthMinus  "<<PtruthMinus<<" PtruthPlus  " << PtruthPlus <<"   TruthTauA1.P() " <<TruthTauA1.P()<<"   " <<TauPPlus << "  "<< TauPMinus <<std::endl;

	      if(TruthAmibiga == 2 && j==2){
			TauA1PPlusPPropagatorPull.at(t).Fill((PtruthPlus - TauPPlus )/deltaTauPplus,1);
	      }
	      
	      if(TruthAmibiga == 1 && j==1){
			TauA1PMinusPPropagatorPull.at(t).Fill((PtruthMinus - TauPMinus )/deltaTauPminus,1);
	      }
	    
	        //  --- ambiguity point 

	      double A  = 2*m1*m1*(m1*m1 + m2*m2)/(4*m1*m1*m2*m2 + pow(m2*m2 - m1*m1,2));
	      double B  = pow( (m2*m2 - m1*m1)/2/m1,2);

	      double ATruth  = 2*m1*m1*(m1*m1 + TruthA1LV.M()*TruthA1LV.M())/(4*m1*m1*TruthA1LV.M()*TruthA1LV.M() + pow(TruthA1LV.M()*TruthA1LV.M() - m1*m1,2));
	      double BTruth  = pow( (TruthA1LV.M()*TruthA1LV.M() - m1*m1)/2/m1,2);


	      double AmbiguityRoot   = sqrt(a1lv.P()*a1lv.P() - B);
	      double TauPAmbiguity = A*AmbiguityRoot ;
	      double dAmbiguityRootdp = 4*m1*m1*a1lv.P()/AmbiguityRoot;
	      double deltaTauPAmbiguityPoint = A*a1lv.P()*deltaPa1/sqrt(a1lv.P()*a1lv.P() - B);

	      double TruthAmibuityTauP = ATruth*sqrt(TruthA1LV.P()*TruthA1LV.P()-BTruth);
	      std::cout<<" A  "<< A <<" B  "<< B <<"  TauPAmbiguity  " << TauPAmbiguity  <<" TruthTauA1.P()  " <<TruthTauA1.P() << " AmbiguityRoot  "<<AmbiguityRoot <<std::endl;
	      std::cout<< "EFValidation tau cov:   " << theTau.Covariance(3,3) << "  "<<theTau.Covariance(4,4) << "  "<<theTau.Covariance(5,5) <<std::endl;
	      if(j ==0 ){
		double deltaPTau = sqrt(pow(theTau.LV().Px(),2) * theTau.Covariance(3,3)  +pow(theTau.LV().Py(),2) * theTau.Covariance(4,4) + pow(theTau.LV().Pz(),2) * theTau.Covariance(5,5) )/theTau.LV().P();

		TauA1PAmbigaPull.at(t).Fill((TruthAmibuityTauP - theTau.LV().P() )/deltaPTau,1);
	

	      }
	      if(j==0 ) TestPropagation.at(t).Fill((sqrt(5* TruthA1LV.P()*TruthA1LV.P()+1) - sqrt(5*daughter.at(0).LV().P()*daughter.at(0).LV().P()+1) )/(5*daughter.at(0).LV().P()*deltaPa1/sqrt(5*daughter.at(0).LV().P()*daughter.at(0).LV().P()+1)),1);


	        //  --- ambiguity point 
	    
	    //------------  check propagation of A1 momenta uncertainty




	      //-----check fit uncert 



	  }
	  
	  NeutrinoA1=daughter.at(1).LV();
	  TauA1ThreeProngFit=theTau.LV();



	  EventFitSuccess = Ntp->EventFit(tau_idx,mu_idx,theTau,theZ,theZdaughter,LC_Eventchi2,NiterationsEF,csum);
	  if(EventFitSuccess){
	    LC_Eventchi2Probability = TMath::Prob(LC_Eventchi2,1);
	    LC_Eventchi2Probabilityndf2 = TMath::Prob(LC_Eventchi2,2);
	    LC_Eventchi2Probabilityndf3 = TMath::Prob(LC_Eventchi2,3);

	    TauA1EventFit =theZdaughter.at(0).LV();
	    TauMuEventFit =theZdaughter.at(1).LV();

	    LVPEF_TauMu = theZdaughter.at(0);
	    LVPEF_TauA1 = theZdaughter.at(1);
	      

	    NeutrinoMu=TauMuEventFit-Ntp->Muons_p4(mu_idx);
	    Z_sol=theZ.LV();
	  }
	}

	SignificanceOfTauRotation.push_back(phisign);
	TauA1Neutrino_TPF.push_back(NeutrinoA1);
	TauMuNetrino_EF.push_back(NeutrinoMu);

	Tau_FitOk.push_back(ThreeProngFitSuccess);
	Tau_FitChi2.push_back(LC_chi2);
  
	Chi2ProbabilityEventFit.push_back(LC_Eventchi2Probability);
	Chi2ProbabilityEventFitndf2.push_back(LC_Eventchi2Probabilityndf2);
	Chi2ProbabilityEventFitndf3.push_back(LC_Eventchi2Probabilityndf3);

	EventFit_Ok.push_back(EventFitSuccess);
  	Z_EF.push_back(Z_sol);
 	Chi2EventFit.push_back(LC_Eventchi2);
	theTauAmbiguityVector.push_back(theTau);

 	theTauA1TLV.push_back(LVPEF_TauA1);
 	theTauMuTLV.push_back(LVPEF_TauMu);


 	TauA1_EF.push_back(TauA1EventFit);
 	TauMu_EF.push_back(TauMuEventFit);
	NiterationsEventFit.push_back(NiterationsEF);
	csum_GF.push_back(csum);
	DeltaPhiSign.push_back(phisign);   
	deltaphiRotation.push_back(deltaphi);
      }
      
      PhiRotationSignificanceUnPhysicalTaus.at(t).Fill(SignificanceOfTauRotation.at(0));
    


      int AmbiguitySolution =0;
      bool AmbPoint(false);
      if(Ntp->AmbiguitySolver(Tau_FitOk,EventFit_Ok,Chi2ProbabilityEventFit,AmbiguitySolution, AmbPoint)){

 	double flightsignificance=Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(tau_idx),Ntp->PFTau_TIP_primaryVertex_cov(tau_idx),theTauAmbiguityVector.at(AmbiguitySolution).Vertex(),theTauAmbiguityVector.at(AmbiguitySolution).VertexCov());
 	pvsvsignificance.at(t).Fill(flightsignificance,1);


	std::cout<<" AmbPoint "<<AmbPoint   <<"  "<<AmbiguitySolution  <<"  "<<Chi2ProbabilityEventFit.at(AmbiguitySolution)<< std::endl; 
	int idToFill(0);
	if(id ==998)idToFill = 1;
	if(id ==10230833)idToFill = 2;
	if(id ==10231433)idToFill = 3;
	if(id ==10231833)idToFill = 4;
	idPassedEF.at(t).Fill(idToFill,1);
      


	TLorentzVector EventFitTauA1 = TauA1_EF.at(AmbiguitySolution);
	TLorentzVector EventFitTauMu = TauMu_EF.at(AmbiguitySolution);



	TLorentzVector EventFitZ = EventFitTauA1+EventFitTauMu;

	if(!AmbPoint){
	  EFitTauA1Pt.at(t).Fill(EventFitTauA1.Pt(),1);
	  EFitTauA1Phi.at(t).Fill(EventFitTauA1.Phi(),1);
	  EFitTauA1Eta.at(t).Fill(EventFitTauA1.Eta(),1);

	  EFitTauMuPt.at(t).Fill(EventFitTauMu.Pt(),1);
	  EFitTauMuPhi.at(t).Fill(EventFitTauMu.Phi(),1);
	  EFitTauMuEta.at(t).Fill(EventFitTauMu.Eta(),1);
	}


	if(AmbPoint){
	  EFitTauA1PtAmbPoint.at(t).Fill(EventFitTauA1.Pt(),1);
	  EFitTauA1PhiAmbPoint.at(t).Fill(EventFitTauA1.Phi(),1);
	  EFitTauA1EtaAmbPoint.at(t).Fill(EventFitTauA1.Eta(),1);

	  EFitTauMuPtAmbPoint.at(t).Fill(EventFitTauMu.Pt(),1);
	  EFitTauMuPhiAmbPoint.at(t).Fill(EventFitTauMu.Phi(),1);
	  EFitTauMuEtaAmbPoint.at(t).Fill(EventFitTauMu.Eta(),1);
	}

	  if(Ntp->CheckDecayID(2,5)){
	    TLorentzVector TruthA1 = Ntp->GetTruthTauProductLV(5, 20213);
	    TLorentzVector TruthMu = Ntp->GetTruthTauProductLV(2, 13);
	
	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	    

	    TVector3 TruthTauA1SV = Ntp->GetTruthTauProductVertex(5,211);
	    std::cout<<" TruthTauA1SV  "<<TruthTauA1SV.X()<<" TruthTauA1SV  "<<TruthTauA1SV.Y()<<" TruthTauA1SV  "<<TruthTauA1SV.Z()<<std::endl;


	    std::cout<<"(TruthTauA1SV.X() - theTauMuTLV.at(AmbiguitySolution).Vertex().X()  " << TruthTauA1SV.X() - theTauMuTLV.at(AmbiguitySolution).Vertex().X() <<std::endl;
	    std::cout<<" sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(0,0) " << theTauA1TLV.at(AmbiguitySolution).Covariance(0,0)<<std::endl;


	    TruthA1PtAfterFit.at(t).Fill(TruthA1.Pt(),1);
	    TruthA1EAfterFit.at(t).Fill(TruthA1.E(),1);
	    TruthA1PhiAfterFit.at(t).Fill(TruthA1.Phi(),1);
	    TruthA1EtaAfterFit.at(t).Fill(TruthA1.Eta(),1);

	    TruthMuPtAfterFit.at(t).Fill(TruthMu.Pt(),1);
	    TruthMuEAfterFit.at(t).Fill(TruthMu.E(),1);
	    TruthMuPhiAfterFit.at(t).Fill(TruthMu.Phi(),1);
	    TruthMuEtaAfterFit.at(t).Fill(TruthMu.Eta(),1);



	    if(!AmbPoint){
	      TauA1PtResolution.at(t).Fill(TruthTauA1.Pt()-EventFitTauA1.Pt(),1);
	      TauA1EResolution.at(t).Fill(TruthTauA1.E()-EventFitTauA1.E(),1);
	      TauA1PhiResolution.at(t).Fill(TruthTauA1.Phi()-EventFitTauA1.Phi(),1);
	      TauA1EtaResolution.at(t).Fill(TruthTauA1.Eta()-EventFitTauA1.Eta(),1);
	      
	      TauMuPtResolution.at(t).Fill(TruthTauMu.Pt()-EventFitTauMu.Pt(),1);
	      TauMuEResolution.at(t).Fill(TruthTauMu.E()-EventFitTauMu.E(),1);
	      TauMuPhiResolution.at(t).Fill(TruthTauMu.Phi()-EventFitTauMu.Phi(),1);
	      TauMuEtaResolution.at(t).Fill(TruthTauMu.Eta()-EventFitTauMu.Eta(),1);


	      TauA1PhiResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauA1.Phi() - EventFitTauA1.Phi(),1);
	      TauMuPhiResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauMu.Phi() -EventFitTauMu.Phi(),1);
	      TauA1EtaResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauMuEtaResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauMu.Eta() -EventFitTauMu.Eta(),1);
	      TauA1PtResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);
	      TauA1EResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauA1.E() - EventFitTauA1.E(),1);
	      TauMuEResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauMu.E() - EventFitTauMu.E(),1);
	    

	      TauA1PtResolVsPt.at(t).Fill(A1Reco.Pt(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsPt.at(t).Fill(MuReco.Pt(),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);
	      TauA1PtResolVsEta.at(t).Fill(A1Reco.Eta(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsEta.at(t).Fill(MuReco.Eta(),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);
	      TauA1EtaResolVsEta.at(t).Fill(A1Reco.Eta(),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauMuEtaResolVsEta.at(t).Fill(MuReco.Eta(),TruthTauMu.Eta() - EventFitTauMu.Eta(),1);


	      TauA1PtResolVsZPt.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsZPt.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauMu.Pt()-EventFitTauMu.Pt(),1);

	      TauA1EtaResolVsZPt.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauMuEtaResolVsZPt.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauMu.Eta() - EventFitTauMu.Eta(),1);


 	      std::cout<<"theTauMuTLV  "<<theTauMuTLV.size()<<std::endl;
 	      std::cout<<"theTauA1TLV  "<<theTauA1TLV.size()<<std::endl;


	      std::cout<<" (TruthTauMu.Px() - theTauMuTLV.at(AmbiguitySolution).LV().Px()  "<< TruthTauMu.Px() - theTauMuTLV.at(AmbiguitySolution).LV().Px() <<"   " <<sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(3,3))<<std::endl;
	      std::cout<<" (TruthTauMu.Px() - theTauMuTLV.at(AmbiguitySolution).LV().Px()  "<< TruthTauMu.Px() - theTauMuTLV.at(AmbiguitySolution).LV().Px() <<"   " <<theTauMuTLV.at(AmbiguitySolution).Covariance(3,3)<<std::endl;

 	      std::cout<<" EventFitTauA1.Px()   "<< EventFitTauA1.Px() <<"   " <<" EventFitTauA1.Py()   "<<EventFitTauA1 .Py() <<"   " <<" EventFitTauA1.Pz()   "<< EventFitTauA1.Pz() <<"   " <<std::endl;
 	      std::cout<<" EventFitTauMu.Px()   "<< EventFitTauMu.Px() <<"   " <<" EventFitTauMu.Py()   "<< EventFitTauMu.Py() <<"   " <<" EventFitTauMu.Pz()   "<< EventFitTauMu.Pz() <<"   " <<std::endl;
 	      std::cout<<" EventFit Zmass    "<< (EventFitTauMu +  EventFitTauA1).M()<< " pt1 " << (EventFitTauMu +  EventFitTauA1).Pt() <<" pt2 "  <<EventFitTauMu.Pt()  -  EventFitTauA1.Pt() <<std::endl;


 	      TauA1SFPxPull.at(t).Fill((TruthTauA1.Px() - theTauAmbiguityVector.at(AmbiguitySolution).LV().Px() )/sqrt(theTauAmbiguityVector.at(AmbiguitySolution).Covariance(3,3)),1);
 	      TauA1SFPyPull.at(t).Fill((TruthTauA1.Py() - theTauAmbiguityVector.at(AmbiguitySolution).LV().Py() )/sqrt(theTauAmbiguityVector.at(AmbiguitySolution).Covariance(4,4)),1);
 	      TauA1SFPzPull.at(t).Fill((TruthTauA1.Pz() - theTauAmbiguityVector.at(AmbiguitySolution).LV().Pz() )/sqrt(theTauAmbiguityVector.at(AmbiguitySolution).Covariance(5,5)),1);


 	      TauMuPxPull.at(t).Fill((TruthTauMu.Px() - theTauMuTLV.at(AmbiguitySolution).LV().Px() )/sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(3,3)),1);
 	      TauMuPyPull.at(t).Fill((TruthTauMu.Py() - theTauMuTLV.at(AmbiguitySolution).LV().Py() )/sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(4,4)),1);
 	      TauMuPzPull.at(t).Fill((TruthTauMu.Pz() - theTauMuTLV.at(AmbiguitySolution).LV().Pz() )/sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(5,5)),1);

	      

	      std::cout<<"sqrt(theTauA1TLV.at(AmbiguitySolution).Covariance(3,3))  " <<sqrt(theTauA1TLV.at(AmbiguitySolution).Covariance(3,3)) <<std::endl;
	      std::cout<<"(TruthTauA1.Px() - theTauA1TLV.at(AmbiguitySolution).LV().Px()  " <<TruthTauA1.Px() - theTauA1TLV.at(AmbiguitySolution).LV().Px()<<std::endl;

	      TauA1PhiResolVsPVSVSignificance.at(t).Fill(flightsignificance,TruthTauA1.Phi() - EventFitTauA1.Phi(),1);
	      TauA1EtaResolVsPVSVSignificance.at(t).Fill(flightsignificance,TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauA1PtResolVsPVSVSignificance.at(t).Fill(flightsignificance,TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauA1EResolVsPVSVSignificance.at(t).Fill(flightsignificance,TruthTauA1.E() - EventFitTauA1.E(),1);
	    }
	    if(AmbPoint){

	      ProbabilityOfAmbiguityPoint.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),1);

	      TauA1PtResolVsPtAmbPoint.at(t).Fill(A1Reco.Pt(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsPtAmbPoint.at(t).Fill(MuReco.Pt(),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);
	      TauA1PtResolVsEtaAmbPoint.at(t).Fill(A1Reco.Eta(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsEtaAmbPoint.at(t).Fill(MuReco.Eta(),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);
	      TauA1EtaResolVsEtaAmbPoint.at(t).Fill(A1Reco.Eta(),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauMuEtaResolVsEtaAmbPoint.at(t).Fill(MuReco.Eta(),TruthTauMu.Eta() - EventFitTauMu.Eta(),1);

	      TauA1PtResolVsZPtAmbPoint.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsZPtAmbPoint.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);

	      TauA1EtaResolVsZPtAmbPoint.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauMuEtaResolVsZPtAmbPoint.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauMu.Eta() - EventFitTauMu.Eta(),1);

	      TauA1PtResolutionAmbPoint.at(t).Fill(TruthTauA1.Pt()-EventFitTauA1.Pt() ,1);
	      TauA1PhiResolutionAmbPoint.at(t).Fill(TruthTauA1.Phi()-EventFitTauA1.Phi(),1);
	      TauA1EtaResolutionAmbPoint.at(t).Fill(TruthTauA1.Eta()-EventFitTauA1.Eta(),1);
	      
	      TauMuPtResolutionAmbPoint.at(t).Fill(TruthTauMu.Pt()-EventFitTauMu.Pt(),1);
	      TauMuPhiResolutionAmbPoint.at(t).Fill(TruthTauMu.Phi()-EventFitTauMu.Phi(),1);
	      TauMuEtaResolutionAmbPoint.at(t).Fill(TruthTauMu.Eta()-EventFitTauMu.Eta(),1);

	      TauA1PhiResolVsPhiRotSignificance.at(t).Fill(SignificanceOfTauRotation.at(AmbiguitySolution),TruthTauA1.Phi() - EventFitTauA1.Phi(),1);
	      TauA1EtaResolVsPhiRotSignificance.at(t).Fill(SignificanceOfTauRotation.at(AmbiguitySolution),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauA1PtResolVsPhiRotSignificance.at(t).Fill(SignificanceOfTauRotation.at(AmbiguitySolution),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauA1EResolVsPhiRotSignificance.at(t).Fill(SignificanceOfTauRotation.at(AmbiguitySolution),TruthTauA1.E() - EventFitTauA1.E(),1);


	      TauA1PhiResolVsPVSVSignificanceAmbPoint.at(t).Fill(flightsignificance,TruthTauA1.Phi() - EventFitTauA1.Phi(),1);
	      TauA1EtaResolVsPVSVSignificanceAmbPoint.at(t).Fill(flightsignificance,TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauA1PtResolVsPVSVSignificanceAmbPoint.at(t).Fill(flightsignificance,TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauA1EResolVsPVSVSignificanceAmbPoint.at(t).Fill(flightsignificance,TruthTauA1.E() - EventFitTauA1.E(),1);

	    }
	  }
      }
      

      if(Tau_FitOk.at(1) && Tau_FitOk.at(2)){
	if(EventFit_Ok.at(1) && EventFit_Ok.at(2)){
	  if(Ntp->CheckDecayID(2,5)){
	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	    

	    int CorrectAmbiguity; int IncorrectAmbiguity;
	    if(fabs(TauA1_EF.at(1).E() - TruthTauA1.E()) < fabs(TauA1_EF.at(2).E() - TruthTauA1.E())){  CorrectAmbiguity=1; IncorrectAmbiguity =2;}
	    if(fabs(TauA1_EF.at(2).E() - TruthTauA1.E()) < fabs(TauA1_EF.at(1).E() - TruthTauA1.E())){  CorrectAmbiguity=2;IncorrectAmbiguity =1; }
	    //----------------
	    Chi2Dim.at(t).Fill(TMath::Prob(Chi2EventFit.at(CorrectAmbiguity),1),TMath::Prob(Chi2EventFit.at(IncorrectAmbiguity),1),1);
	    csum2Dim.at(t).Fill(csum_GF.at(CorrectAmbiguity),csum_GF.at(IncorrectAmbiguity),1);
	    csumProb2Dim.at(t).Fill(csum_GF.at(CorrectAmbiguity),Chi2ProbabilityEventFit.at(CorrectAmbiguity),1);
	    if((TruthTauA1 + TruthTauMu).Pt() < 3)	    ProbabilityOfCorrect.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),1);



	    ProbabilityOfCorrectVsChi2.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),Chi2EventFit.at(CorrectAmbiguity),1);
	    ProbabilityOfCorrectndf2.at(t).Fill(Chi2ProbabilityEventFitndf2.at(CorrectAmbiguity),1);
	    ProbabilityOfCorrectndf3.at(t).Fill(Chi2ProbabilityEventFitndf3.at(CorrectAmbiguity),1);
	    Chi2OfCorrect.at(t).Fill(Chi2EventFit.at(CorrectAmbiguity),1);
	    std::cout<<" chiw EFVaildation code --->   "<< Chi2EventFit.at(CorrectAmbiguity) << std::endl;
	    ProbabilityOfCorrectVsZPt.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),(TruthTauA1 + TruthTauMu).Pt(),1);
	  }
	}
      }


    }
  }
}




void  EFValidation::Finish(){
  unsigned int t=1;


 
  for(int iBin = 1; iBin < TauA1PtResolVsPtRMS.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauA1PtResolVsPt.at(t).GetNbinsX();
    TauA1PtResolVsPtRMS.at(t).SetBinContent(iBin,TauA1PtResolVsPt.at(t).ProjectionY(" ",iBin, NMax)->GetRMS());
    TauA1PtResolVsPtRMS.at(t).SetBinError(iBin,TauA1PtResolVsPt.at(t).ProjectionY(" ",iBin, NMax)->GetRMSError());
  }


  for(int iBin = 1; iBin < TauMuPtResolVsPtRMS.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauMuPtResolVsPt.at(t).GetNbinsX();
    TauMuPtResolVsPtRMS.at(t).SetBinContent(iBin,TauMuPtResolVsPt.at(t).ProjectionY(" ",iBin, NMax)->GetRMS());
    TauMuPtResolVsPtRMS.at(t).SetBinError(iBin,TauMuPtResolVsPt.at(t).ProjectionY(" ",iBin, NMax)->GetRMSError());
  }

  for(int iBin = 1; iBin < TauA1PtResolVsEtaRMS.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauA1PtResolVsEta.at(t).GetNbinsX();
    TauA1PtResolVsEtaRMS.at(t).SetBinContent(iBin,TauA1PtResolVsEta.at(t).ProjectionY(" ",iBin, iBin )->GetRMS());
    TauA1PtResolVsEtaRMS.at(t).SetBinError(iBin,TauA1PtResolVsEta.at(t).ProjectionY(" ",iBin, iBin)->GetRMSError());
  }

  for(int iBin = 1; iBin < TauMuPtResolVsZPtRMS.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauMuPtResolVsZPt.at(t).GetNbinsX();
    TauMuPtResolVsZPtRMS.at(t).SetBinContent(iBin,TauMuPtResolVsZPt.at(t).ProjectionY(" ",1, iBin)->GetRMS());
    TauMuPtResolVsZPtRMS.at(t).SetBinError(iBin,TauMuPtResolVsZPt.at(t).ProjectionY(" ",1, iBin)->GetRMSError());
  }

  for(int iBin = 1; iBin < TauMuPtResolVsZPtMean.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauMuPtResolVsZPt.at(t).GetNbinsX();
    TauMuPtResolVsZPtMean.at(t).SetBinContent(iBin,TauMuPtResolVsZPt.at(t).ProjectionY(" ",1, iBin)->GetMean());
    TauMuPtResolVsZPtMean.at(t).SetBinError(iBin,TauMuPtResolVsZPt.at(t).ProjectionY(" ",1, iBin)->GetMeanError());
  }

  for(int iBin = 1; iBin < TauA1PtResolVsPVSVSignificanceRMS.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauA1PtResolVsPVSVSignificance.at(t).GetNbinsX();
    TauA1PtResolVsPVSVSignificanceRMS.at(t).SetBinContent(iBin,TauA1PtResolVsPVSVSignificance.at(t).ProjectionY(" ",iBin, NMax)->GetRMS());
    TauA1PtResolVsPVSVSignificanceRMS.at(t).SetBinError(iBin,TauA1PtResolVsPVSVSignificance.at(t).ProjectionY(" ",iBin, NMax)->GetRMSError());
  }
  
  for(int iBin = 1; iBin < TauA1PtResolVsPhiRotSignificance.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauA1PtResolVsPhiRotSignificance.at(t).GetNbinsX();
    
    double  nom = TauA1PtResolVsPhiRotSignificance.at(t).ProjectionY(" ",0, iBin)->GetEntries();
    double  denom = TauA1PtResolVsPhiRotSignificance.at(t).ProjectionY(" ",0, NMax)->GetEntries();

    
    if(denom!=0){     PhiRotSignificanceCutEfficiency.at(t).SetBinContent(iBin,nom/denom);}
    else{PhiRotSignificanceCutEfficiency.at(t).SetBinContent(iBin,0);}
    
    
  }
  
  for(int iBin = 1; iBin < TauA1PtResolVsPVSVSignificance.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauA1PtResolVsPVSVSignificance.at(t).GetNbinsX();
    double  nom = TauA1PtResolVsPVSVSignificance.at(t).ProjectionY(" ",iBin,NMax)->GetEntries();
    double  denom = TauA1PtResolVsPVSVSignificance.at(t).ProjectionY(" ",0, NMax)->GetEntries();
    
    if(denom!=0){     PVSVSignificanceCutEfficiency.at(t).SetBinContent(iBin,nom/denom);}
    else{PVSVSignificanceCutEfficiency.at(t).SetBinContent(iBin,0);}
  }



  ChannelEfficiency.at(1).Divide(&idPassedEF.at(1), &idAllTaus.at(1));
  ChannelEfficiency.at(2).Divide(&idPassedEF.at(2), &idAllTaus.at(2));
  ChannelEfficiency.at(3).Divide(&idPassedEF.at(3), &idAllTaus.at(3));
  ChannelEfficiency.at(4).Divide(&idPassedEF.at(4), &idAllTaus.at(4));
  EfficiencyOverA1Pt.at(t).Divide(&TruthA1PtAfterFit.at(t),&TruthA1Pt.at(t));
  EfficiencyOverA1Phi.at(t).Divide(&TruthA1PhiAfterFit.at(t),&TruthA1Phi.at(t));
  EfficiencyOverA1Eta.at(t).Divide(&TruthA1EtaAfterFit.at(t),&TruthA1Eta.at(t));


  EfficiencyOverMuPt.at(t).Divide(&TruthMuPtAfterFit.at(t),&TruthMuPt.at(t));
  EfficiencyOverMuPhi.at(t).Divide(&TruthMuPhiAfterFit.at(t),&TruthMuPhi.at(t));
  EfficiencyOverMuEta.at(t).Divide(&TruthMuEtaAfterFit.at(t),&TruthMuEta.at(t));


  


  
  Selection::Finish();
}



