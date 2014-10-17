#include "TemplatesData.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"
#include "TauRecHelper.h"

TemplatesData::TemplatesData(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

TemplatesData::~TemplatesData(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "TemplatesData::~TemplatesData Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "TemplatesData::~TemplatesData()" << std::endl;
}

void  TemplatesData::Configure(){
  // Setup Cut Values

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)                   cut.at(TriggerOk)=1;
    if(i==PrimeVtx)                    cut.at(PrimeVtx)=1;
    if(i==hasTau)                      cut.at(hasTau)=1;
    if(i==TauIsMediumIsolated)         cut.at(TauIsMediumIsolated)=1;
    if(i==TauHPSID)                    cut.at(TauHPSID)=10;
    if(i==hasMuon)                     cut.at(hasMuon)= 1;
    if(i==MuonIso)                     cut.at(MuonIso)=0.12;
    if(i==MuonID)                      cut.at(MuonID)=1;
    if(i==PUJetID)                     cut.at(PUJetID)=0;
    if(i==FlightLenghtSignificance)    cut.at(FlightLenghtSignificance)=2.2; // 0.05  for 0.3 cone
    if(i==MET)                         cut.at(MET)=30; // 0.05  for 0.3 cone
    if(i==charge)                      cut.at(charge)=-1; // 0.05  for 0.3 cone




 

  }

  TString hlabel;
  TString htitle;
  for(unsigned int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
      if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>$";
      title.at(i)+=cut.at(PrimeVtx);
      title.at(i)+=")";
      htitle=title.at(i); 
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
   else if(i==hasTau){
      title.at(i)="At least 1 good tau";
      hlabel="hasTau";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTau_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTau_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }

   else if(i==TauIsMediumIsolated){
      title.at(i)=" TauIsMediumIsolated";
      hlabel="TauIsMediumIsolated";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsMediumIsolated_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsMediumIsolated_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
   else if(i==TauHPSID){
      title.at(i)="TauHPSID";
      hlabel="TauHPSID";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauHPSID_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauHPSID_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
   else if(i==hasMuon){
      title.at(i)="At least 1 good mu";
      hlabel="hasMuon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasMuon_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasMuon_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }



   else if(i==MuonIso){
      title.at(i)="MuonIsIsolated";
      hlabel="MuonIsolation";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIso_",htitle,25,0,0.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIso_",htitle,25,0,0.2,hlabel,"Events"));
    }
   else if(i==MuonID){
      title.at(i)="MuonID";
      hlabel="MuonID";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
   else if(i==PUJetID){
      title.at(i)="PUJetID";
      hlabel="PUJetID";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PUJetID_",htitle,50,-1,1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PUJetID_",htitle,50,-1,1,hlabel,"Events"));
    }


   else if(i==FlightLenghtSignificance){
      title.at(i)="$ Sigma_{PV-SV} > $";
      title.at(i)+=cut.at(FlightLenghtSignificance);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#sigma_{PV-SV}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_FlightLenghtSignificance_",htitle,20,0,10,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_FlightLenghtSignificance_",htitle,20,0,10,hlabel,"Events"));
    }

   else if(i==MET){
      title.at(i)="$M_{T} < $";
      title.at(i)+=cut.at(MET);
      title.at(i)+=" GeV ";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{T}, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MET_",htitle,15,0,100,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MET_",htitle,15,0,100,hlabel,"Events"));
    }  
   else if(i==charge){
      title.at(i)="$Pair charge = $";
      title.at(i)+=cut.at(charge);
      hlabel="charge";    
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,3,-1.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,3,-1.5,1.5,hlabel,"Events"));
    } 






  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  TauA1VisiblePt=HConfig.GetTH1D(Name+"_TauA1VisiblePt","TauA1VisiblePt",15,15,50,"Pt_{a1}, GeV","");
  TauA1VisiblePtIsoCut=HConfig.GetTH1D(Name+"_TauA1VisiblePtIsoCut","TauA1VisiblePtIsoCut",15,15,50,"Pt_{a1}, GeV","");

  TauA1VisibleE=HConfig.GetTH1D(Name+"_TauA1VisibleE","TauA1VisibleE",15,15,80,"E_{a1}, GeV","");

  TauA1VisiblePtCorrected1=HConfig.GetTH1D(Name+"_TauA1VisiblePtCorrected1","TauA1VisiblePtCorrected1",15,15,50,"Pt_{a1}, GeV","");
  TauA1VisiblePtCorrected2=HConfig.GetTH1D(Name+"_TauA1VisiblePtCorrected2","TauA1VisiblePtCorrected2",15,15,50,"Pt_{a1}, GeV","");
  TauA1VisiblePtCorrected3=HConfig.GetTH1D(Name+"_TauA1VisiblePtCorrected3","TauA1VisiblePtCorrected3",15,15,50,"Pt_{a1}, GeV","");


  TauA1VisiblePtA1LV=HConfig.GetTH1D(Name+"_TauA1VisiblePtA1LV","TauA1VisiblePtA1LV",15,15,50,"Pt_{a1}, GeV","");
  TauMuVisiblePt=HConfig.GetTH1D(Name+"_TauMuVisiblePt","TauMuVisiblePt",15,15,50,"Pt_{#mu}, GeV","");

  TauA1VisibleEta=HConfig.GetTH1D(Name+"_TauA1VisibleEta","",10,-2.5,2.5,"#eta_{a1}","");
  TauMuVisibleEta=HConfig.GetTH1D(Name+"_TauMuVisibleEta","",10,-2.5,2.5,"#eta_{#mu}","");

  TauA1VisiblePhi=HConfig.GetTH1D(Name+"_TauA1VisiblePhi","",10,-3.15,3.15,"#phi_{a1}","");
  TauMuVisiblePhi=HConfig.GetTH1D(Name+"_TauMuVisiblePhi","",10,-3.15,3.15,"#phi_{#mu}","");


  MET_hCorrT0pc_et=HConfig.GetTH1D(Name+"_MET_hCorrT0pc_et","",20,0,80," met  ,GeV","");
  MET_hCorrT0pcT1_et=HConfig.GetTH1D(Name+"_MET_hCorrT0pcT1_et","",20,0,80," met  ,GeV","");
  MET_hCorrT0rtTxy_et=HConfig.GetTH1D(Name+"_MET_hCorrT0rtTxy_et","",20,0,80," met  ,GeV","");
  
  MET_hCorrT0rtT1Txy_et=HConfig.GetTH1D(Name+"_MET_hCorrT0rtT1Txy_et","",20,0,80," met  ,GeV","");
  MET_hCorrT0pcTxy_et=HConfig.GetTH1D(Name+"_MET_hCorrT0pcTxy_et","",20,0,80," met  ,GeV","");
  MET_hCorrT1_et=HConfig.GetTH1D(Name+"_MET_hCorrT1_et","",20,0,80," met  ,GeV","");
  MET_hCorrMVA_et=HConfig.GetTH1D(Name+"_MET_hCorrMVA_et","",20,0,80," met  ,GeV","");



  MET_hCorrT0pc_phi=HConfig.GetTH1D(Name+"_MET_hCorrT0pc_phi","",10,-3.14,3.14,"met #phi, rad  ","");
  MET_hCorrT0pcT1_phi=HConfig.GetTH1D(Name+"_MET_hCorrT0pcT1_phi","",10,-3.14,3.14,"met #phi, rad   ","");
  MET_hCorrT0rtTxy_phi=HConfig.GetTH1D(Name+"_MET_hCorrT0rtTxy_phi","",10,-3.14,3.14,"met #phi, rad   ","");

  MET_hCorrT0rtT1Txy_phi=HConfig.GetTH1D(Name+"_MET_hCorrT0rtT1Txy_phi","",10,-3.14,3.14,"met #phi, rad   ","");
  MET_hCorrT0pcTxy_phi=HConfig.GetTH1D(Name+"_MET_hCorrT0pcTxy_phi","",10,-3.14,3.14,"met #phi, rad   ","");
  MET_hCorrT1_phi=HConfig.GetTH1D(Name+"_MET_hCorrT1_phi","",10,-3.14,3.14,"met #phi, rad   ","");
  MET_hCorrMVA_phi=HConfig.GetTH1D(Name+"_MET_hCorrMVA_phi","",10,-3.14,3.14,"met #phi, rad   ","");


  MET_hCorrT0pcT1Txy_et=HConfig.GetTH1D(Name+"_MET_hCorrT0pcT1Txy_et","",20,0,80," met  ,GeV","");
  MET_hCorrT0pcT1Txy_phi=HConfig.GetTH1D(Name+"_MET_hCorrT0pcT1Txy_phi","",10,-3.14,3.14,"met #phi, rad   ","");

  VisibleInvariantMass=HConfig.GetTH1D(Name+"_VisibleInvariantMass","VisibleInvariantMass",15,20,120,"M_{a1-#mu}, GeV","");


  hTauJetdR=HConfig.GetTH1D(Name+"_hTauJetdR","",20,0,0.2,"#DeltaR (tau-pfjet)","");
  hPFJet_PUJetID_discr=HConfig.GetTH1D(Name+"_hPFJet_PUJetID_discr","",50,-1,1,"pujetId  mva output","");
  hPFJet_PUJetID_looseWP=HConfig.GetTH1D(Name+"_hPFJet_PUJetID_looseWP","",2,-0.5,1.5,"PUJetIDLoose","");
  hPFJet_PUJetID_mediumWP=HConfig.GetTH1D(Name+"_hPFJet_PUJetID_mediumWP","",2,-0.5,1.5,"PUJetIDMedium","");
  hPFJet_PUJetID_tightWP=HConfig.GetTH1D(Name+"_hPFJet_PUJetID_tightWP","",2,-0.05,1.5,"PUJetIDTight","");


  SSPionsPt=HConfig.GetTH1D(Name+"_SSPionsPt","",20,0,80,"SS pions Pt, GeV","");
  SSPionsE=HConfig.GetTH1D(Name+"_SSPionsE","",20,0,50,"SS pions E, GeV","");
  OSPionsPt=HConfig.GetTH1D(Name+"_OSPionsPt","",20,0,50,"OS pion Pt, GeV","");
  OSPionsE=HConfig.GetTH1D(Name+"_OSPionsE","",20,0,80,"OS pion E, GeV","");


  AllPionsPt=HConfig.GetTH1D(Name+"_AllPionsPt","",20,0,50,"All pions Pt, GeV","");
  AllPionsE=HConfig.GetTH1D(Name+"_AllPionsE","",20,0,80,"All pions E, GeV","");


  TauA1VisiblePx=HConfig.GetTH1D(Name+"_TauA1VisiblePx","",15,0,60,"Px_{a1}, GeV","");
  TauA1VisiblePy=HConfig.GetTH1D(Name+"_TauA1VisiblePy","",15,0,60,"Py_{a1}, GeV","");
  TauA1VisiblePz=HConfig.GetTH1D(Name+"_TauA1VisiblePz","",15,0,60,"Pz_{a1}, GeV","");

  TauA1VisibleP=HConfig.GetTH1D(Name+"_TauA1VisibleP","",15,15,60,"P, GeV, GeV","");
  //  TauA1VisiblePhi=HConfig.GetTH1D(Name+"_TauA1VisiblePhi","",30,-3.14,3.14,"TauA1VisiblePhi","");
  TauA1VisibleTheta=HConfig.GetTH1D(Name+"_TauA1VisibleTheta","",15,0,3.14,"#theta_{a1}","");
  TauA1VisibleM=HConfig.GetTH1D(Name+"_TauA1VisibleM","",15,0.5,1.777," A1 mass, GeV","");

  SumPtOfjets=HConfig.GetTH1D(Name+"_SumPtOfjets","",20,0,50," pT of All Jets, GeV ","");
  LeadingJetPt=HConfig.GetTH1D(Name+"_LeadingJetPt","",20,0,50," ","");
  
  LeadingJetdRMu=HConfig.GetTH1D(Name+"_LeadingJetdRMu","",20,0,1," ","");
  LeadingJetdRA1=HConfig.GetTH1D(Name+"_LeadingJetdRA1","",20,0,1," ","");


  EFitTauA1Pt=HConfig.GetTH1D(Name+"_EFitTauA1Pt"," TauA1 Pt",10,10,60,"#tau_{a1} Pt, (GeV)","");
  EFitTauA1Phi=HConfig.GetTH1D(Name+"_EFitTauA1Phi"," TauA1 Phi",10,-3.14,3.14,"#tau_{a1} Phi, (rad)","");
  EFitTauA1Eta=HConfig.GetTH1D(Name+"_EFitTauA1Eta"," TauA1 Eta",10,-2,2,"#tau_{a1} Eta","");

  EFitTauMuPt=HConfig.GetTH1D(Name+"_EFitTauMuPt"," TauMu Pt",10,10,60,"#tau_{#mu} Pt, (GeV)","");
  EFitTauMuPhi=HConfig.GetTH1D(Name+"_EFitTauMuPhi"," TauMu Phi",10,-3.14,3.14,"#tau_{#mu} Phi, (rad)","");
  EFitTauMuEta=HConfig.GetTH1D(Name+"_EFitTauMuEta"," TauMu Eta",10,-2,2,"#tau_{#mu} Eta","");



  EFitTauA1PtAmbPoint=HConfig.GetTH1D(Name+"_EFitTauA1PtAmbPoint"," TauA1 Pt",20,10,60,"#tau_{a1} Pt, (GeV)","");
  EFitTauA1PhiAmbPoint=HConfig.GetTH1D(Name+"_EFitTauA1PhiAmbPoint"," TauA1 Phi",20,-3.14,3.14,"#tau_{a1} Phi, (rad)","");
  EFitTauA1EtaAmbPoint=HConfig.GetTH1D(Name+"_EFitTauA1EtaAmbPoint"," TauA1 Eta",20,-2,2,"#tau_{a1} Eta","");

  EFitTauMuPtAmbPoint=HConfig.GetTH1D(Name+"_EFitTauMuPtAmbPoint"," TauMu Pt",20,10,60,"#tau_{#mu} Pt, (GeV)","");
  EFitTauMuPhiAmbPoint=HConfig.GetTH1D(Name+"_EFitTauMuPhiAmbPoint"," TauMu Phi",20,-3.14,3.14,"#tau_{#mu} Phi, (rad)","");
  EFitTauMuEtaAmbPoint=HConfig.GetTH1D(Name+"_EFitTauMuEtaAmbPoint"," TauMu Eta",20,-2,2,"#tau_{#mu} Eta","");

  TauA1Pt=HConfig.GetTH1D(Name+"_TauA1Pt"," TauA1 Pt",20,10,60,"#tau_{a1} Pt, (GeV)","");
  TauA1Phi=HConfig.GetTH1D(Name+"_TauA1Phi"," TauA1 Phi",20,-3.14,3.14,"#tau_{a1} Phi, (rad)","");
  TauA1Eta=HConfig.GetTH1D(Name+"_TauA1Eta"," TauA1 Eta",20,-2,2,"#tau_{a1} Eta","");

  TauMuPt=HConfig.GetTH1D(Name+"_TauMuPt"," TauMu Pt",20,10,60,"#tau_{#mu} Pt, (GeV)","");
  TauMuPhi=HConfig.GetTH1D(Name+"_TauMuPhi"," TauMu Phi",20,-3.14,3.14,"#tau_{#mu} Phi, (rad)","");
  TauMuEta=HConfig.GetTH1D(Name+"_TauMuEta"," TauMu Eta",20,-2,2,"#tau_{#mu} Eta","");


  RecOmegaA1p=HConfig.GetTH1D(Name+"_RecOmegaA1p","RecOmegaA1p",20,-1,1,"#omega_{a1}","");
  RecOmegaA1m=HConfig.GetTH1D(Name+"_RecOmegaA1m","RecOmegaA1m",20,-1,1,"#omega_{a1}","");
  RecOmegaA1z=HConfig.GetTH1D(Name+"_RecOmegaA1z","RecOmegaA1z",20,-1,1,"#omega_{a1}","");

  RecoCosTheta_CosBeta_ZFrame_Plus=HConfig.GetTH2D(Name+"_RecoCosTheta_CosBeta_ZFrame_Plus","RecoCosTheta_CosBeta_ZFrame_Plus",20,-1,1,20,-1,1,"cos#theta","cos#beta");
  RecoCosTheta_CosBeta_ZFrame_Minus=HConfig.GetTH2D(Name+"_RecoCosTheta_CosBeta_ZFrame_Minus","RecoCosTheta_CosBeta_ZFrame_Minus",20,-1,1,20,-1,1,"cos#theta","cos#beta");
  RecoCosTheta_CosBeta_ZFrame_Zero=HConfig.GetTH2D(Name+"_RecoCosTheta_CosBeta_ZFrame_Zero","RecoCosTheta_CosBeta_ZFrame_Zero",20,-1,1,20,-1,1,"cos#theta","cos#beta");

  
  RecOmegaMup=HConfig.GetTH1D(Name+"_RecOmegaMup","RecOmegaMup",20,-1,2,"#omega_{#mu}","");
  RecOmegaMum=HConfig.GetTH1D(Name+"_RecOmegaMum","RecOmegaMum",20,-1,2,"#omega_{#mu}","");
  RecOmegaMuz=HConfig.GetTH1D(Name+"_RecOmegaMuz","RecOmegaMuz",20,-1,2,"#omega_{#mu}","");
  
  
  RecOmegaCombinedp=HConfig.GetTH1D(Name+"_RecOmegaCombinedp","RecOmegaCombinedp",20,-1,2,"#omega_{#mu}","");
  RecOmegaCombinedm=HConfig.GetTH1D(Name+"_RecOmegaCombinedm","RecOmegaCombinedm",20,-1,2,"#omega_{#mu}","");
  RecOmegaCombinedz=HConfig.GetTH1D(Name+"_RecOmegaCombinedz","RecOmegaCombinedz",20,-1,2,"#omega_{#mu}","");


//   TauA1VisiblePtisTightIsolationDBSumPtCorr=HConfig.GetTH1D(Name+"_TauA1VisiblePtisTightIsolationDBSumPtCorr","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtisMediumIsolationDBSumPtCorr=HConfig.GetTH1D(Name+"_TauA1VisiblePtisMediumIsolationDBSumPtCorr","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtisLooseIsolationDBSumPtCorr=HConfig.GetTH1D(Name+"_TauA1VisiblePtisLooseIsolationDBSumPtCorr","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtisVLooseIsolationDBSumPtCorr=HConfig.GetTH1D(Name+"_TauA1VisiblePtisVLooseIsolationDBSumPtCorr","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtHPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtHPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtHPSPFTauDiscriminationByLooseIsolationMVA=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByLooseIsolationMVA","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtHPSPFTauDiscriminationByMediumIsolationMVA=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByMediumIsolationMVA","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtHPSPFTauDiscriminationByTightIsolationMVA=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByTightIsolationMVA","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtHPSPFTauDiscriminationByLooseIsolationMVA2=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByLooseIsolationMVA2","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtHPSPFTauDiscriminationByMediumIsolationMVA2=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByMediumIsolationMVA2","",15,15,50,"Pt_{a1}, GeV","");
//   TauA1VisiblePtHPSPFTauDiscriminationByTightIsolationMVA2=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByTightIsolationMVA2","",15,15,50,"Pt_{a1}, GeV","");
 
  TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSS=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSS","",15,15,50,"Pt_{a1}, GeV","");
  TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonTauIso=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonTauIso","",15,15,50,"Pt_{a1}, GeV","");

  TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSScaled=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSScaled","",15,15,50,"Pt_{a1}, GeV","");
  TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonTauIsoScaled=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonTauIsoScaled","",15,15,50,"Pt_{a1}, GeV","");


  PhiRotationSignificanceUnPhysicalTaus=HConfig.GetTH1D(Name+"_PhiRotationSignificanceUnPhysicalTaus","Significance of #phi angle rotation",20,0,10,"#sigma(#phi)","");
  PhiRotationSignificancePhysicalTaus=HConfig.GetTH1D(Name+"_PhiRotationSignificancePhysicalTaus","Significance of #phi angle rotation",20,0,10,"#sigma(#phi)","");

 TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits","",15,15,50,"Pt_{a1}, GeV","");
  TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsScaled=HConfig.GetTH1D(Name+"_TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsScaled","",15,15,50,"Pt_{a1}, GeV","");
//   ZMass =HConfig.GetTH1D(Name+"_ZMass","ZMass",20,40,160,"M_{#tau-#mu}, GeV","");
//   ZMassAmbPoint =HConfig.GetTH1D(Name+"_ZMassAmbPoint","ZMassAmbPoint",20,40,160,"M_{#tau-#mu}, GeV","");
//   VisibleMass =HConfig.GetTH1D(Name+"_VisibleMass","VisibleMass",20,20,160,"M_{#tau-#mu}, GeV","");
//   tauptamb =HConfig.GetTH1D(Name+"_tauptamb","tauptamb",30,15,60,"pT_{#tau_{#mu}}","");
//   taumupt =HConfig.GetTH1D(Name+"_taumupt","taumupt",30,15,60,"pT_{#tau_{#mu}}","");
//   taumuphi=HConfig.GetTH1D(Name+"_taumuphi","taumuphi",20,-3.16,3.16,"#phi_{#tau_{#mu}}","");
//   taumutheta =HConfig.GetTH1D(Name+"_taumutheta","taumutheta",20,-2.1,2.11,"#eta_{#tau_{#mu}}","");




//   MuonPhi =HConfig.GetTH1D(Name+"_MuonPhi","MuonPhi",20,-3.16,3.16,"#phi_{#mu}","");
//   TauPhi =HConfig.GetTH1D(Name+"_TauPhi","TauPhi",20,-3.16,3.16,"#phi_{#tau}","");

//   MuonEta =HConfig.GetTH1D(Name+"_MuonEta","MuonEta",20,-2.1,2.1,"#eta_{#mu}","");
//   TauEta =HConfig.GetTH1D(Name+"_TauEta","TauEta",20,-2.1,2.11,"#eta_{#tau}","");

//   deltaPhi =HConfig.GetTH1D(Name+"_deltaPhi","deltaPhi",30,0,6.28,"#Delta#phi","");

//   FlightLenght =HConfig.GetTH1D(Name+"_FlightLenght","FlightLenght",20,0,2,"L, cm","");
//   Pzeta =HConfig.GetTH1D(Name+"_Pzeta","Pzeta",50,-150,100,"p_{#zeta}^{mis} - 0.85Pp_{#zeta}^{vis} ","");
//   //  Jak1VsJak2 =HConfig.GetTH2D(Name+"_Jak1VsJak2","Jak1VsJak2",20,-0.5,20.5,20,-0.5,20.5,"Jak1VsJak2","");

//   hGJ0 =HConfig.GetTH1D(Name+"_hGJ0","hGJ0",25,0,0.05,"hGJ0","");
//   hGJ1 =HConfig.GetTH1D(Name+"_hGJ1","hGJ1",25,0,0.05,"hGJ1","");
//   hGJ2 =HConfig.GetTH1D(Name+"_hGJ2","hGJ2",25,0,0.05,"hGJ2","");

//   rhGJ0 =HConfig.GetTH1D(Name+"_rhGJ0","rhGJ0",25,0,2,"rhGJ0","");
//   rhGJ1 =HConfig.GetTH1D(Name+"_rhGJ1","rhGJ1",25,0,2,"rhGJ1","");
//   rhGJ2 =HConfig.GetTH1D(Name+"_rhGJ2","rhGJ2",25,0,2,"rhGJ2","");

//   a1Mass =HConfig.GetTH1D(Name+"_a1Mass","a1Mass",25,0,2,"M_{a1}, GeV","");
//   a1MassAlternat =HConfig.GetTH1D(Name+"_a1MassAlternat","a1MassAlternat",25,0,2,"M_{a1}, GeV","");
//   Efrac =HConfig.GetTH1D(Name+"_Efrac","Efrac",25,0,1.5,"Et^{vis}/Et^{full}","");
//   EMiss =HConfig.GetTH1D(Name+"_EMiss","EMiss",25,0,100,"E_{T}, GeV","");
//   MuonDistance =HConfig.GetTH1D(Name+"_MuonDistance","MuonDistance",25,0,0.5,"L_{#mu}, cm","");
//   DeltaTheta_1st=HConfig.GetTH1D(Name+"_DeltaTheta_1st","DeltaTheta_1st",20,-3.16,3.16,"#Delta #theta_{1}, rad","");
//   DeltaTheta_2nd=HConfig.GetTH1D(Name+"_DeltaTheta_2nd","DeltaTheta_2nd",20,-3.16,3.16,"#Delta #theta_{1}, rad","");
//   MuTauDeltaR=HConfig.GetTH1D(Name+"_MuTauDeltaR","MuTauDeltaR",30,0,1,"MuTauDeltaR","");
//   TauPhiPlus=HConfig.GetTH1D(Name+"_TauPhiPlus","TauPhiPlus",20,-3.16,3.16,"TauPhiPlus","");
//   TauPhiMinus=HConfig.GetTH1D(Name+"_TauPhiMinus","TauPhiMinus",20,-3.16,3.16,"TauPhiMinus","");
//   TauPhiZero=HConfig.GetTH1D(Name+"_TauPhiZero","TauPhiZero",20,-3.16,3.16,"TauPhiZero","");
//   x=HConfig.GetTH1D(Name+"_x","x",20,0,1.5,"x","");
//   xamb=HConfig.GetTH1D(Name+"_xamb","xamb",20,0,1.5,"xamb","");

//   RecoTauValid=HConfig.GetTH1D(Name+"_RecoTauValid","RecoTauValid",2,-0.5,1.5,"RecoTauValid","");

//   TauEnergy=HConfig.GetTH1D(Name+"_TauEnergy","TauEnergy",30,20,150,"E_{#tau}, GeV","");
//   MuonEnergy=HConfig.GetTH1D(Name+"_MuonEnergy","MuonEnergy",30,20,150,"E_{#mu}, GeV","");

//   gamma = HConfig.GetTH1D(Name+"_gamma","gamma",30,-1.5,1.5,"gamma","");
//   HPSPhi=HConfig.GetTH1D(Name+"_HPSPhi","HPSPhi",20,-3.16,3.16,"HPSPhi","");
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);

} 

 
// TString name,TString title,int nbinsx, double minx, double maxx, 
// 				       int nbinsy, double miny, double maxy, TString xaxis, TString yaxis)

void  TemplatesData::Store_ExtraDist(){

  
  Extradist1d.push_back(&TauA1VisiblePt);
  Extradist1d.push_back(&TauA1VisibleE);
  Extradist1d.push_back(&TauMuVisiblePt);
  Extradist1d.push_back(&VisibleInvariantMass);


  Extradist1d.push_back(&TauA1VisibleEta);
  Extradist1d.push_back(&TauMuVisibleEta);

  //  Extradist1d.push_back(&TauA1VisiblePhi);
  Extradist1d.push_back(&TauMuVisiblePhi);


//   Extradist1d.push_back(&MET_hCorrT0pc_et);
//   Extradist1d.push_back(&MET_hCorrT0pcT1_et);
//   Extradist1d.push_back(&MET_hCorrT0rtTxy_et);

  Extradist1d.push_back(&MET_hCorrT0rtT1Txy_et);
//   Extradist1d.push_back(&MET_hCorrT0pcTxy_et);
//   Extradist1d.push_back(&MET_hCorrT1_et);
//   Extradist1d.push_back(&MET_hCorrMVA_et);


  
//   Extradist1d.push_back(&MET_hCorrT0pc_phi);
//   Extradist1d.push_back(&MET_hCorrT0pcT1_phi);
//   Extradist1d.push_back(&MET_hCorrT0rtTxy_phi);

  Extradist1d.push_back(&MET_hCorrT0rtT1Txy_phi);
//   Extradist1d.push_back(&MET_hCorrT0pcTxy_phi);
//   Extradist1d.push_back(&MET_hCorrT1_phi);
//   Extradist1d.push_back(&MET_hCorrMVA_phi);

//   Extradist1d.push_back(&MET_hCorrT0pcT1Txy_et);
//   Extradist1d.push_back(&MET_hCorrT0pcT1Txy_phi);
  Extradist1d.push_back(&TauA1VisiblePtCorrected1);
  Extradist1d.push_back(&TauA1VisiblePtCorrected2);
  Extradist1d.push_back(&TauA1VisiblePtCorrected3);

  Extradist1d.push_back(&hTauJetdR);
  Extradist1d.push_back(&hPFJet_PUJetID_discr);
  Extradist1d.push_back(&hPFJet_PUJetID_looseWP);
  Extradist1d.push_back(&hPFJet_PUJetID_mediumWP);
  Extradist1d.push_back(&hPFJet_PUJetID_tightWP);

  Extradist1d.push_back(&SSPionsPt);
  Extradist1d.push_back(&SSPionsE);


  Extradist1d.push_back(&OSPionsPt);
  Extradist1d.push_back(&OSPionsE);

  Extradist1d.push_back(&AllPionsPt);
  Extradist1d.push_back(&AllPionsE);

  Extradist1d.push_back(&TauA1VisiblePtA1LV);

  Extradist1d.push_back(&TauA1VisiblePx);
  Extradist1d.push_back(&TauA1VisiblePy);
  Extradist1d.push_back(&TauA1VisiblePz);

  Extradist1d.push_back(&TauA1VisibleP);
  Extradist1d.push_back(&TauA1VisiblePhi);
  Extradist1d.push_back(&TauA1VisibleTheta);
  Extradist1d.push_back(&TauA1VisibleM);

  Extradist1d.push_back(&SumPtOfjets);
  Extradist1d.push_back(&LeadingJetPt);

//   Extradist1d.push_back(&LeadingJetdRMu);
//   Extradist1d.push_back(&LeadingJetdRA1);

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

  //  Extradist1d.push_back(&PhiRotationSignificancePhysicalTaus);
  Extradist1d.push_back(&PhiRotationSignificanceUnPhysicalTaus);

  Extradist1d.push_back(&TauA1Pt);
  Extradist1d.push_back(&TauA1Phi);
  Extradist1d.push_back(&TauA1Eta);
	
  Extradist1d.push_back(&TauMuPt);
  Extradist1d.push_back(&TauMuPhi);
  Extradist1d.push_back(&TauMuEta);
  Extradist1d.push_back(&TauA1VisiblePtIsoCut);

//   Extradist1d.push_back(&TauA1VisiblePtisTightIsolationDBSumPtCorr);
//   Extradist1d.push_back(&TauA1VisiblePtisMediumIsolationDBSumPtCorr);
//   Extradist1d.push_back(&TauA1VisiblePtisLooseIsolationDBSumPtCorr);
//   Extradist1d.push_back(&TauA1VisiblePtisVLooseIsolationDBSumPtCorr);
//   Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits);
  Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits);
  Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsScaled);
  //  Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonMuonIso);
//   Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits);
//   Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByLooseIsolationMVA);
//   Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByMediumIsolationMVA);
//   Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByTightIsolationMVA);
//   Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByLooseIsolationMVA2);
//   Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByMediumIsolationMVA2);
//   Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByTightIsolationMVA2);


    Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSS);
    Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonTauIso);

    Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSScaled);
    Extradist1d.push_back(&TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonTauIsoScaled);

    Extradist1d.push_back(&RecOmegaA1p);
    Extradist1d.push_back(&RecOmegaA1m);
    Extradist1d.push_back(&RecOmegaA1z);

    Extradist2d.push_back(&RecoCosTheta_CosBeta_ZFrame_Plus);
    Extradist2d.push_back(&RecoCosTheta_CosBeta_ZFrame_Minus);
    Extradist2d.push_back(&RecoCosTheta_CosBeta_ZFrame_Zero);


    Extradist1d.push_back(&RecOmegaMup);
    Extradist1d.push_back(&RecOmegaMum);
    Extradist1d.push_back(&RecOmegaMuz);

    Extradist1d.push_back(&RecOmegaCombinedp);
    Extradist1d.push_back(&RecOmegaCombinedm);
    Extradist1d.push_back(&RecOmegaCombinedz);


  //  Extradist1d.push_back(&NGoodVtx);
  //  Extradist1d.push_back(&TauPt);
  //   Extradist1d.push_back(&ZMass);
  //   Extradist1d.push_back(&ZMassAmbPoint);
  
  //   Extradist1d.push_back(&MuonPhi);
  //   Extradist1d.push_back(&TauPhi);
  
  //   Extradist1d.push_back(&MuonEta);
  //   Extradist1d.push_back(&TauEta);  
  
  //   Extradist1d.push_back(&deltaPhi);
  //   Extradist1d.push_back(&FlightLenght);
  //   Extradist1d.push_back(&Pzeta);
  //   Extradist1d.push_back(&hGJ0);
  //   Extradist1d.push_back(&hGJ1);
  //   Extradist1d.push_back(&hGJ2);
  
  //   Extradist1d.push_back(&rhGJ0);
  //   Extradist1d.push_back(&rhGJ1);
  //   Extradist1d.push_back(&rhGJ2);
  
  //   Extradist1d.push_back(&RecoTauValid);
  
  //   Extradist1d.push_back(&a1Mass);
  //   Extradist1d.push_back(&a1MassAlternat);
  //   Extradist1d.push_back(&Efrac);
  //   Extradist1d.push_back(&DeltaTheta_1st);
  //   Extradist1d.push_back(&DeltaTheta_2nd);
  
  
  //   Extradist1d.push_back(&taumupt);
  //   Extradist1d.push_back(&tauptamb);
  //   Extradist1d.push_back(&taumuphi);
  //   Extradist1d.push_back(&taumutheta);
  //   Extradist1d.push_back(&MuonDistance);
  //   Extradist1d.push_back(&EMiss);
  //   Extradist1d.push_back(&MuTauDeltaR);
  //   Extradist1d.push_back(&VisibleMass);

  //   Extradist1d.push_back(&TauEnergy);
  //   Extradist1d.push_back(&MuonEnergy);
  //   Extradist1d.push_back(&HPSPhi);
  //   Extradist1d.push_back(&x);
  //   Extradist1d.push_back(&xamb);

//   Extradist1d.push_back(&gamma);

}

void  TemplatesData::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id " << id <<std::endl; return;}
  

  // Apply Selection
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }


  unsigned int TauCandidate;
  unsigned int MuonCandidate;
  
  
  ///////////////////////////////////
  // Fill vector of object candidate
  std::vector<unsigned int> mu_idx_good;
  unsigned mu_idx(999);
  for(unsigned int i=0;i<Ntp->NMuons();i++){
    if(Ntp->isGoodMuon(i) && fabs(Ntp->Muons_p4(i).Eta())<2.4 && Ntp->Muons_p4(i).Pt() > 20){
      mu_idx_good.push_back(i);
    }  
  }


  std::vector<unsigned int> tau_idx_good;
  unsigned int tau_idx(999);
  for(unsigned int i=0;i<Ntp->NPFTaus();i++){
    if(Ntp->PFTau_p4(i).Pt()> 20){
      if(fabs(Ntp->PFTau_p4(i).Eta())<2.7){
	if(Ntp->PFTau_hpsDecayMode(i)==10 && Ntp->PFTau_isHPSByDecayModeFinding(i) && Ntp->PFTau_isHPSAgainstMuonTight(i) && Ntp->PFTau_isHPSAgainstElectronsTight(i)  ){
	  if(  Ntp->PFTau_TIP_hassecondaryVertex(i)){
	  tau_idx_good.push_back(i);
	  }
	}
      }
    }
  }


  ////////////////////////////////////////////////////
  // Find the leading pt ones if more than one found 

  double tempmupt(0);
  for(unsigned int i =0; i<mu_idx_good.size(); i++){
    if(Ntp->Muons_p4(mu_idx_good.at(i)).Pt() > tempmupt){
      tempmupt = Ntp->Muons_p4(mu_idx_good.at(i)).Pt();
      MuonCandidate = mu_idx_good.at(i);
    }
  }
    
  double temptaupt(0);
  for(unsigned int i =0; i<tau_idx_good.size(); i++){
    if(Ntp->PFTau_p4(tau_idx_good.at(i)).Pt() >temptaupt ){
      temptaupt = Ntp->PFTau_p4(tau_idx_good.at(i)).Pt();
      TauCandidate= tau_idx_good.at(i);
    }
  }

  //check muon is from the same vertex 
  double diffPV(999.);
  int indexpv(-1);
  
  value.at(PrimeVtx)=nGoodVtx;
  
  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("eta2p1_LooseIsoPFTau"))value.at(TriggerOk)+=1;
  
  value.at(hasTau) =  tau_idx_good.size();
  value.at(hasMuon) = mu_idx_good.size();


  if((value.at(hasTau)==cut.at(hasTau)) && (value.at(hasMuon) ==cut.at(hasMuon))){
    
    unsigned int PFJetIndex(0);
    double TauJetdR(999.);
    for(unsigned int ijet = 0; ijet < Ntp->NPFJets(); ijet++){
      if(Ntp->PFTau_p4(TauCandidate).DeltaR(Ntp->PFJet_p4(ijet))  <TauJetdR ){
	TauJetdR=Ntp->PFTau_p4(TauCandidate).DeltaR(Ntp->PFJet_p4(ijet));
	PFJetIndex = ijet;
      }
    }

    //  value.at(TauIsMediumIsolated) =   Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate);
    value.at(TauIsMediumIsolated) =   Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate);
    value.at(TauHPSID) = Ntp->PFTau_hpsDecayMode(TauCandidate);
    
    value.at(MuonIso) =  ((Ntp->Muon_sumChargedHadronPt04(MuonCandidate) + max(0.,Ntp->Muon_sumNeutralHadronEt04(MuonCandidate) + Ntp->Muon_sumPhotonEt04(MuonCandidate)  - 0.5*Ntp->Muon_sumPUPt04(MuonCandidate)))/ Ntp->Muons_p4(MuonCandidate).Pt());
    value.at(MuonID) =  Ntp->isMuonID(MuonCandidate); 


    value.at(FlightLenghtSignificance)=Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate),Ntp->PFTau_TIP_primaryVertex_cov(TauCandidate),Ntp->PFTau_a1_lvp(TauCandidate).Vertex(),Ntp->PFTau_a1_lvp(TauCandidate).VertexCov());
    
    value.at(MET) =   sqrt(2*Ntp->Muons_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muons_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) );
    value.at(charge) =Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate);
  
    value.at(PUJetID) = Ntp->PFJet_PUJetID_discr(PFJetIndex);


  }

  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  pass.at(hasTau)=(value.at(hasTau)==cut.at(hasTau));
  pass.at(TauIsMediumIsolated)=(value.at(TauIsMediumIsolated) == cut.at(TauIsMediumIsolated));
  pass.at(TauHPSID)=(value.at(TauHPSID) ==cut.at(TauHPSID));

  pass.at(hasMuon)=(value.at(hasMuon) ==cut.at(hasMuon));
  pass.at(MuonIso)=(value.at(MuonIso) <= cut.at(MuonIso));
  pass.at(MuonID) = (value.at(MuonID) == cut.at(MuonID));
  pass.at(PUJetID) = (value.at(PUJetID) > cut.at(PUJetID));


  pass.at(FlightLenghtSignificance)=(value.at(FlightLenghtSignificance)>=cut.at(FlightLenghtSignificance));
  pass.at(MET)=(value.at(MET)<=cut.at(MET));
  pass.at(charge)=(value.at(charge) == cut.at(charge));


  
  double wobs=1;
  double w;
  
    







  if(!Ntp->isData()){w = Ntp->EvtWeight3D();}
  else{w=1;}
  
    
  if(!pass.at(charge) and   !pass.at(MuonIso)  ){
    //if( value.at(MuonIso) > 0.25){
    //      std::cout<< value.at(MuonIso) <<std::endl;
      if(Ntp->isData() && value.at(MuonIso) > 0.2 && value.at(MuonIso) < 0.5){
      
	if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id for QCD: "<< DataMCType::QCD <<std::endl; return;}
	pass.at(charge)=true;
	pass.at(MuonIso) = true;
	//  }
    }
  }
  

//   //  std::cout<<"before id=: "<<id<<std::endl;
//    if(pass.at(charge) and  !pass.at(MuonIso)){
//      if(id==20){
//        // if(checkLooseIsolation)
//        {
//  	if(!HConfig.GetHisto(id!=20,DataMCType::W_lnu,t)){ std::cout << "failed to find id for WJ: "<< DataMCType::W_lnu <<std::endl; return;}
//  	//std::cout<<"after id=: "<<id<<std::endl;
//  	//	 std::cout<<"=== t:"<<t<<std::endl;
//  	pass.at(charge)=true;
//  	pass.at(MuonIso) = true;
//  	//	pass.at(FlightLenghtSignificance)= true;
	
//        }
//      }
//    }else if(pass.at(charge) and  pass.at(MuonIso)){
//      if(id==20){
//        //if(checkLooseIsolation)
//        {
//  	if(!HConfig.GetHisto(id!=20,DataMCType::W_lnu,t)){ std::cout << "failed to find id for WJ: "<< DataMCType::W_lnu <<std::endl; return;}
//  	pass.at(charge)=true;
//  	pass.at(MuonIso) = true;
//  	//	pass.at(FlightLenghtSignificance)= true;
//        }
//      }
  
//    }
  
//   if(pass.at(charge) and  !pass.at(TauIsMediumIsolated)){
//     if(id==998){
//          if(!HConfig.GetHisto(id!=998,998,t)){  return;}
// 	 pass.at(charge)=true;
// 	 pass.at(TauIsMediumIsolated) = true;
// 	 //	 pass.at(FlightLenghtSignificance)= true;
     
//     }
//   }
//   else if(pass.at(charge) and  pass.at(TauIsMediumIsolated)){
//     if(id==998){
//       if(!HConfig.GetHisto(id!=998,998,t)){  return;}
//       pass.at(charge)=true;
//       pass.at(TauIsMediumIsolated) = true;
//       //   pass.at(FlightLenghtSignificance)= true;
      
//     }

//   }
  
    
  if( pass.at(TriggerOk) &&
      pass.at(PrimeVtx) &&
      pass.at(hasTau) &&
      pass.at(TauHPSID)&&
      pass.at(MuonIso)&&
      pass.at(MuonID) &&
      pass.at(PUJetID)&&
      pass.at(FlightLenghtSignificance)&&
      pass.at(MET)&&
      pass.at(charge) ){

    TLorentzVector TauScaled;
    if(!Ntp->isData()){   TauScaled = Ntp->PFTau_p4(TauCandidate)*(1.012 + 0.001 * TMath::Min(TMath::Max( Ntp->PFTau_p4(TauCandidate).Pt()-32.,0.),18.)) ;}
    else TauScaled = Ntp->PFTau_p4(TauCandidate);


    TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
    TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsScaled.at(t).Fill(TauScaled.Pt(),w);




   }


  if( pass.at(TriggerOk) &&
      pass.at(PrimeVtx) &&
      pass.at(hasTau) &&
      pass.at(TauHPSID)&&
      pass.at(MuonIso)&&
      pass.at(MuonID) &&
      pass.at(PUJetID)&&
      pass.at(FlightLenghtSignificance)&&
      pass.at(MET)&&
      !pass.at(charge) ){

    TLorentzVector TauScaled;
    if(!Ntp->isData()){   TauScaled = Ntp->PFTau_p4(TauCandidate)*(1.012 + 0.001 * TMath::Min(TMath::Max( Ntp->PFTau_p4(TauCandidate).Pt()-32.,0.),18.)) ;}
    else TauScaled = Ntp->PFTau_p4(TauCandidate);


      
      
      TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSS.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
      TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSScaled.at(t).Fill(TauScaled.Pt(),w);

  }
  if( pass.at(TriggerOk) &&
      pass.at(PrimeVtx) &&
      pass.at(hasTau) &&
      pass.at(TauHPSID)&&
      pass.at(TauIsMediumIsolated)&&
      pass.at(MuonIso)&&
      pass.at(MuonID) &&
      pass.at(PUJetID)&&
      pass.at(FlightLenghtSignificance)&&
      pass.at(MET)&&
      !pass.at(charge) ){
      TLorentzVector TauScaled = Ntp->PFTau_p4(TauCandidate)*(1.012 + 0.001 * TMath::Min(TMath::Max( Ntp->PFTau_p4(TauCandidate).Pt()-32.,0.),18.)) ;
      
      
      TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonTauIso.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
      TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsSSNonTauIsoScaled.at(t).Fill(TauScaled.Pt(),w);

  }





  bool status=AnalysisCuts(t,w,wobs); 
  
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){

//   std::cout<<" debug 1  " <<std::endl;

//   for(int ipv =0; ipv < Ntp->NVtx(); ipv++){
//     if(sqrt(pow(Ntp->Vtx(ipv).X() - Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate).X(),2)+ 
// 	    pow(Ntp->Vtx(ipv).Y() - Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate).Y(),2)+ 
// 	    pow(Ntp->Vtx(ipv).Z() - Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate).Z(),2)) < diffPV ){
      
//       diffPV = sqrt(pow(Ntp->Vtx(ipv).X() - Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate).X(),2)+ 
// 		    pow(Ntp->Vtx(ipv).Y() - Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate).Y(),2)+ 
// 		    pow(Ntp->Vtx(ipv).Z() - Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate).Z(),2));
//       indexpv = ipv;
//     }
    
//   }
//   double deltaRMuTrack(999.);
//   std::cout<<"--- indexpv  "<<indexpv <<"  Ntp->NPVTracks(indexpv) "<< Ntp->NPVTracks(indexpv) <<std::endl;

//   for(int ipvtrack =0; ipvtrack < Ntp->NPVTracks(indexpv);ipvtrack++){
    
//     if(Ntp->Vtx_TracksP4(indexpv,ipvtrack).DeltaR(Ntp->Muons_p4(MuonCandidate)) < deltaRMuTrack){
//       deltaRMuTrack=Ntp->Vtx_TracksP4(indexpv,ipvtrack).DeltaR(Ntp->Muons_p4(MuonCandidate));

//     }
//   }

//   std::cout<<"--- deltaRMuTrack  "<<deltaRMuTrack <<" diffPV  "<< diffPV <<std::endl;






    unsigned int nGoodVtx=0;
      for(unsigned int i=0;i<Ntp->NVtx();i++){
	if(Ntp->isVtxGood(i))nGoodVtx++;
      }

    TLorentzVector TauScaled;
    if(!Ntp->isData()){   TauScaled = Ntp->PFTau_p4(TauCandidate)*(1.012 + 0.001 * TMath::Min(TMath::Max( Ntp->PFTau_p4(TauCandidate).Pt()-32.,0.),18.)) ;}
    else TauScaled = Ntp->PFTau_p4(TauCandidate);



 //      if(Ntp->PFTau_isTightIsolationDBSumPtCorr(TauCandidate)) TauA1VisiblePtisTightIsolationDBSumPtCorr.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_isMediumIsolationDBSumPtCorr(TauCandidate)) TauA1VisiblePtisMediumIsolationDBSumPtCorr.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_isLooseIsolationDBSumPtCorr(TauCandidate)) TauA1VisiblePtisLooseIsolationDBSumPtCorr.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_isVLooseIsolationDBSumPtCorr(TauCandidate)) TauA1VisiblePtisVLooseIsolationDBSumPtCorr.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits(TauCandidate)) TauA1VisiblePtHPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate)) TauA1VisiblePtHPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits(TauCandidate)) TauA1VisiblePtHPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA(TauCandidate)) TauA1VisiblePtHPSPFTauDiscriminationByLooseIsolationMVA.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA(TauCandidate)) TauA1VisiblePtHPSPFTauDiscriminationByMediumIsolationMVA.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA(TauCandidate)) TauA1VisiblePtHPSPFTauDiscriminationByTightIsolationMVA.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2(TauCandidate)) TauA1VisiblePtHPSPFTauDiscriminationByLooseIsolationMVA2.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2(TauCandidate)) TauA1VisiblePtHPSPFTauDiscriminationByMediumIsolationMVA2.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
//       if(Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA2(TauCandidate)) TauA1VisiblePtHPSPFTauDiscriminationByTightIsolationMVA2.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);




      if(Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate))      TauA1VisiblePtIsoCut.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
     
      //    NGoodVtx.at(t).Fill(nGoodVtx,w);;
      //      std::cout<<"  id "<<id<<std::endl;


//        std::cout<<" TauJetdR "<<TauJetdR <<std::endl;
//        std::cout<<" Ntp->PFJet_PUJetID_discr(PFJetIndex) "<<Ntp->PFJet_PUJetID_discr(PFJetIndex) <<std::endl;
      unsigned int PFJetIndex(0);
      double TauJetdR(999.);
      for(unsigned int ijet = 0; ijet < Ntp->NPFJets(); ijet++){
	if(Ntp->PFTau_p4(TauCandidate).DeltaR(Ntp->PFJet_p4(ijet))  <TauJetdR ){
	  TauJetdR=Ntp->PFTau_p4(TauCandidate).DeltaR(Ntp->PFJet_p4(ijet));
	  PFJetIndex = ijet;
	}
      }
      unsigned int PFJetMuonIndex(0);
      double MuJetdR(999.);
      for(unsigned int ijet = 0; ijet < Ntp->NPFJets(); ijet++){
	//	std::cout<<" jet pt  "<< Ntp->PFJet_p4(ijet).Pt() << std::endl;
	if(Ntp->Muons_p4(MuonCandidate).DeltaR(Ntp->PFJet_p4(ijet))  < MuJetdR){
	  MuJetdR=Ntp->Muons_p4(MuonCandidate).DeltaR(Ntp->PFJet_p4(ijet));
	  PFJetMuonIndex= ijet;
	}
      }

      //////////////////////////////////////////////////
      //  find all the jets that belong to signal vertex
      std::vector<int> PFJetsIndices; 
      for( unsigned int ijet = 0; ijet < Ntp->NPFJets(); ijet++){

	if(ijet !=PFJetIndex && ijet !=PFJetMuonIndex) {

	  if(Ntp->PFJet_PUJetID_discr(ijet) > 0.4){

	    PFJetsIndices.push_back(ijet);

	  }
	}
      }
// 	std::cout<<"=== mu pt "<< Ntp->Muons_p4(MuonCandidate).Pt() << std::endl;
// 	std::cout<<"=== a1 pt "<< Ntp->PFTau_p4(TauCandidate).Pt() << std::endl;

//       std::cout<<" PFJetsIndices.size()  " <<PFJetsIndices.size() <<std::endl;
      TLorentzVector SumAllJets(0,0,0,0);
      double MaxJetPt(0);
      int LeadJetIndex;

      for(unsigned int ibadjet =0; ibadjet <  PFJetsIndices.size(); ibadjet++){

	//std::cout<<"=== jet pt  reduced  "<< Ntp->PFJet_p4(PFJetsIndices.at(ibadjet)).Pt() << std::endl;
	SumAllJets+=Ntp->PFJet_p4(PFJetsIndices.at(ibadjet));
	if(Ntp->PFJet_p4(PFJetsIndices.at(ibadjet)).Pt() > MaxJetPt){
	  LeadJetIndex = PFJetsIndices.at(ibadjet);
	  MaxJetPt = Ntp->PFJet_p4(PFJetsIndices.at(ibadjet)).Pt();

	}
      }

      SumPtOfjets.at(t).Fill(SumAllJets.Pt(),w);
      if(PFJetsIndices.size()!=0)  LeadingJetPt.at(t).Fill(Ntp->PFJet_p4(LeadJetIndex).Pt(),w);
      
      if(PFJetsIndices.size()!=0)  LeadingJetdRMu.at(t).Fill(Ntp->PFJet_p4(LeadJetIndex).DeltaR(Ntp->Muons_p4(MuonCandidate)),w);
      if(PFJetsIndices.size()!=0)  LeadingJetdRA1.at(t).Fill(Ntp->PFJet_p4(LeadJetIndex).DeltaR(Ntp->PFTau_p4(TauCandidate)),w);
	

      TLorentzVector SSPion1,SSPion2, OSPion;
      unsigned int ssindex2(0);
      
      if(Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,1) == 1 && Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,2) == -1){
	SSPion1 = Ntp->PFTau_PionsP4(TauCandidate,0);
	SSPion2 = Ntp->PFTau_PionsP4(TauCandidate,1);
	OSPion = Ntp->PFTau_PionsP4(TauCandidate,2);

      }
      if(Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,1) == -1 && Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,2) == 1){
	SSPion1 = Ntp->PFTau_PionsP4(TauCandidate,0);
	SSPion2 = Ntp->PFTau_PionsP4(TauCandidate,2);
	OSPion = Ntp->PFTau_PionsP4(TauCandidate,1);
	
      }
      if(Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,1) == -1 && Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,2) == -1){
	SSPion1 = Ntp->PFTau_PionsP4(TauCandidate,1);
	SSPion2 = Ntp->PFTau_PionsP4(TauCandidate,2);
	OSPion = Ntp->PFTau_PionsP4(TauCandidate,0);
      }
      
      SSPionsPt.at(t).Fill((SSPion1 + SSPion2).Pt(),w);
      SSPionsE.at(t).Fill((SSPion1 + SSPion2).E(),w);
      
      
      OSPionsPt.at(t).Fill(OSPion.Pt(),w);
      OSPionsE.at(t).Fill(OSPion.E(),w);
      
      
      AllPionsPt.at(t).Fill((SSPion1 + SSPion2+OSPion).Pt(),w);
      AllPionsE.at(t).Fill((SSPion1 + SSPion2+OSPion).E(),w);
      

      
      hTauJetdR.at(t).Fill(TauJetdR,w);
      hPFJet_PUJetID_discr.at(t).Fill(Ntp->PFJet_PUJetID_discr(PFJetIndex),w);
      hPFJet_PUJetID_looseWP.at(t).Fill(Ntp->PFJet_PUJetID_looseWP(PFJetIndex),w);
      hPFJet_PUJetID_mediumWP.at(t).Fill(Ntp->PFJet_PUJetID_mediumWP(PFJetIndex),w);
      hPFJet_PUJetID_tightWP.at(t).Fill(Ntp->PFJet_PUJetID_tightWP(PFJetIndex),w);
  

      ////////////////////////////////////
     //compute trigger efficiency weights 
      double wTau(1.);

      if(id==998 or id == 10230833){
	if(fabs(Ntp->PFTau_p4(TauCandidate).Eta()) < 1.5){
	  wTau = TrigEff(Ntp->PFTau_p4(TauCandidate).Pt(),18.604910,    0.276042,    0.137039,    2.698437,    0.940721)/TrigEff(Ntp->PFTau_p4(TauCandidate).Pt(),18.532997,    1.027880,    2.262950,    1.003322,    5.297292);
	}
	if(fabs(Ntp->PFTau_p4(TauCandidate).Eta()) > 1.5){
	  wTau = TrigEff(Ntp->PFTau_p4(TauCandidate).Pt(),18.701715,    0.216523,    0.148111,    2.245081,    0.895320)/TrigEff(Ntp->PFTau_p4(TauCandidate).Pt(),18.212782,    0.338119,    0.122828,    12.577926,    0.893975);
	}
      }
      double wMuon(1.);
      if(id==998 or id == 10230833){
	if(Ntp->Muons_p4(MuonCandidate).Eta() < -1.2){
	  wMuon = TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),5.9977, 	7.64004e-05 ,	6.4951e-08, 	1.57403 ,	0.865325)/TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),16.0051, 	2.45144e-05, 	4.3335e-09, 	1.66134, 	0.87045);
	}
	if(Ntp->Muons_p4(MuonCandidate).Eta() >  -1.2 && Ntp->Muons_p4(MuonCandidate).Eta() < -0.8 ){
	  wMuon = TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),17.3974, 	0.804001, 	1.47145, 	1.24295, 	0.928198)/TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),17.3135, 	0.747636, 	1.21803, 	1.40611, 	0.934983);
	}
	if(Ntp->Muons_p4(MuonCandidate).Eta() >  -0.8 && Ntp->Muons_p4(MuonCandidate).Eta() <0 ){
	  wMuon = TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),16.4307,	0.226312, 	0.265553, 	1.55756, 	0.974462)/TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),15.9556, 	0.0236127, 	0.00589832, 	1.75409, 	0.981338);
	}
	if(Ntp->Muons_p4(MuonCandidate).Eta() >  0 && Ntp->Muons_p4(MuonCandidate).Eta() <0.8 ){
	  wMuon = TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),17.313, 	        0.662731, 	1.3412,  	1.05778, 	1.26624)/TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),15.9289, 	0.0271317, 	0.00448573, 	1.92101, 	0.978625);
	}
	
	if(Ntp->Muons_p4(MuonCandidate).Eta() >  0.8 && Ntp->Muons_p4(MuonCandidate).Eta() <1.2 ){
	  wMuon = TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),16.9966, 	0.550532, 	0.807863, 	1.55402, 	0.885134)/TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),16.5678, 	0.328333, 	0.354533, 	1.67085, 	0.916992);
	}
	
	if(Ntp->Muons_p4(MuonCandidate).Eta() >  1.2){
	  wMuon = TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),15.9962, 	0.000106195, 	4.95058e-08, 	1.9991, 	0.851294)/TrigEff(Ntp->Muons_p4(MuonCandidate).Pt(),15.997, 	7.90069e-05, 	4.40036e-08 ,	1.66272, 	0.884502);
	}
	
      }
      // IsoMu eta < -1.2 - 2012ABCD 	15.9977 	7.64004e-05 	6.4951e-08 	1.57403 	0.865325
      // IsoMu17 eta < -1.2 - 53X MC 	16.0051 	2.45144e-05 	4.3335e-09 	1.66134 	0.87045
      

      // IsoMu -1.2 < eta < -0.8 - 2012ABCD 	17.3974 	0.804001 	1.47145 	1.24295 	0.928198
      // IsoMu17 -1.2 < eta < -0.8 - 53X MC 	17.3135 	0.747636 	1.21803 	1.40611 	0.934983

      
      // IsoMu -0.8 < eta < 0 - 2012ABCD 	16.4307 	0.226312 	0.265553 	1.55756 	0.974462
      // IsoMu17 -0.8 < eta < 0 - 53X MC 	15.9556 	0.0236127 	0.00589832 	1.75409 	0.981338


      // IsoMu 0 < eta < 0.8 - 2012ABCD 	17.313 	        0.662731 	1.3412  	1.05778 	1.26624
      // IsoMu17 0 < eta < 0.8 - 53X MC 	15.9289 	0.0271317 	0.00448573 	1.92101 	0.978625


      // IsoMu 0.8 < eta < 1.2 - 2012ABCD 	16.9966 	0.550532 	0.807863 	1.55402 	0.885134
      // IsoMu17 0.8 < eta < 1.2 - 53X MC 	16.5678 	0.328333 	0.354533 	1.67085 	0.916992


      // IsoMu eta > 1.2 - 2012ABCD 	15.9962 	0.000106195 	4.95058e-08 	1.9991 	0.851294
      // IsoMu17 eta > 1.2 - 53X MC 	15.997 	7.90069e-05 	4.40036e-08 	1.66272 	0.884502


      // ----  energy scale  p4' = p4 * (1.012 + 0.001 * TMath::Min(TMath::Max(pT-32.,0.),18.) 
      // Barrel Data /afs/cern.ch/user/b/benitezj/www/Summer13Studies/TauTrigger/muTauABCD_June30/TauTriggerEfficiency_muTau_Barrel_Summer13ReReco_Data.root : 
      // 18.604910,    0.276042,    0.137039,    2.698437,    0.940721

      // EndCap Data /afs/cern.ch/user/b/benitezj/www/Summer13Studies/TauTrigger/muTauABCD_June30/TauTriggerEfficiency_muTau_EndCap_Summer13ReReco_Data.root : 
      // 18.701715,    0.216523,    0.148111,    2.245081,    0.895320



      // Barrle MC   /afs/cern.ch/user/b/benitezj/www/Summer13Studies/TauTrigger/muTauABCD_June30/TauTriggerEfficiency_muTau_Barrel_Summer13ReReco_MC.root : 
      // 18.532997,    1.027880,    2.262950,    1.003322,    5.297292

      // EndCap MC/afs/cern.ch/user/b/benitezj/www/Summer13Studies/TauTrigger/muTauABCD_June30/TauTriggerEfficiency_muTau_EndCap_Summer13ReReco_MC.root : 
      // 18.212782,    0.338119,    0.122828,    12.577926,    0.893975

//       std::cout<<" Tau Weight "<<wTau <<std::endl;
//       std::cout<<" Muon  Weight "<<wMuon <<std::endl;
//       std::cout<<" tau times muon "<<wTau*wMuon << "   Tau Pt " << Ntp->PFTau_p4(TauCandidate).Pt()<<std::endl;

//       std::cout<<"--------------- "<<std::endl;

//       std::cout<<" P - PCorrected "<<Ntp->PFTau_p4(TauCandidate).P() - TauScaled.P() <<std::endl;


      TauA1VisiblePt.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w);
      TauMuVisiblePt.at(t).Fill(Ntp->Muons_p4(MuonCandidate).Pt(),w);
      TauA1VisibleE.at(t).Fill(Ntp->PFTau_p4(TauCandidate).E(),w);

      TauA1VisiblePtCorrected1.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pt(),w/wTau/wMuon);
      TauA1VisiblePtCorrected2.at(t).Fill(TauScaled.Pt(),w);
      TauA1VisiblePtCorrected3.at(t).Fill(TauScaled.Pt(),w/wTau/wMuon);



      TauA1VisiblePx.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Px(),w);
      TauA1VisiblePy.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Py(),w);
      TauA1VisiblePz.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Pz(),w);

      TauA1VisibleP.at(t).Fill(Ntp->PFTau_p4(TauCandidate).P(),w);
      //      TauA1VisiblePhi.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Phi(),w);
      TauA1VisibleTheta.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Theta(),w);
      TauA1VisibleM.at(t).Fill(Ntp->PFTau_p4(TauCandidate).M(),w);



      TauA1VisiblePtA1LV.at(t).Fill(Ntp->PFTau_a1_lvp(TauCandidate).LV().Pt(),w);

      TauA1VisibleEta.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Eta(),w);
      TauMuVisibleEta.at(t).Fill(Ntp->Muons_p4(MuonCandidate).Eta(),w);

      TauA1VisiblePhi.at(t).Fill(Ntp->PFTau_p4(TauCandidate).Phi(),w);
      TauMuVisiblePhi.at(t).Fill(Ntp->Muons_p4(MuonCandidate).Phi(),w);



      VisibleInvariantMass.at(t).Fill((Ntp->Muons_p4(MuonCandidate) + Ntp->PFTau_p4(TauCandidate)).M(),w);
      
      MET_hCorrT0pc_et.at(t).Fill(Ntp->MET_CorrT0pc_et(),w);
      MET_hCorrT0pcT1_et.at(t).Fill(Ntp->MET_CorrT0pcT1_et(),w);
      MET_hCorrT0rtTxy_et.at(t).Fill(Ntp->MET_CorrT0rtTxy_et(),w);
      MET_hCorrT0pcT1Txy_et.at(t).Fill(Ntp->MET_CorrT0pcT1Txy_et(),w);


      MET_hCorrT0rtT1Txy_et.at(t).Fill(Ntp->MET_CorrT0rtT1Txy_et(),w);
      MET_hCorrT0pcTxy_et.at(t).Fill(Ntp->MET_CorrT0pcTxy_et(),w);
      MET_hCorrT1_et.at(t).Fill(Ntp->MET_CorrT1_et(),w);
      MET_hCorrMVA_et.at(t).Fill(Ntp->MET_CorrMVA_et(),w);



      
      MET_hCorrT0pc_phi.at(t).Fill(Ntp->MET_CorrT0pc_phi(),w);
      MET_hCorrT0pcT1_phi.at(t).Fill(Ntp->MET_CorrT0pcT1_phi(),w);
      MET_hCorrT0rtTxy_phi.at(t).Fill(Ntp->MET_CorrT0rtTxy_phi(),w);
      MET_hCorrT0pcT1Txy_phi.at(t).Fill(Ntp->MET_CorrT0pcT1Txy_phi(),w);


      MET_hCorrT0rtT1Txy_phi.at(t).Fill(Ntp->MET_CorrT0rtT1Txy_phi(),w);
      MET_hCorrT0pcTxy_phi.at(t).Fill(Ntp->MET_CorrT0pcTxy_phi(),w);
      MET_hCorrT1_phi.at(t).Fill(Ntp->MET_CorrT1_phi(),w);
      MET_hCorrMVA_phi.at(t).Fill(Ntp->MET_CorrMVA_phi(),w);
      


      ////////////////////////////////
      // Run kinematic fit ///////////

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
      bool ptbalance(true);
      for(unsigned int ibadjet =0; ibadjet <  PFJetsIndices.size(); ibadjet++){

	if(Ntp->PFJet_p4(PFJetsIndices.at(ibadjet)).Pt() > 20){
	  ptbalance = false;
	}
     }


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
	ThreeProngFitSuccess=Ntp->ThreeProngTauFit(TauCandidate,j,theTau,daughter,LC_chi2,phisign);
	if(ThreeProngFitSuccess){

	  NeutrinoA1=daughter.at(1).LV();
	  TauA1ThreeProngFit=theTau.LV();
	  
	  EventFitSuccess = Ntp->EventFit(TauCandidate,MuonCandidate,theTau,theZ,theZdaughter,LC_Eventchi2,NiterationsEF,csum);
	  if(EventFitSuccess){
	    LC_Eventchi2Probability = TMath::Prob(LC_Eventchi2,1);
	    LC_Eventchi2Probabilityndf2 = TMath::Prob(LC_Eventchi2,2);
	    LC_Eventchi2Probabilityndf3 = TMath::Prob(LC_Eventchi2,3);


	    TauA1EventFit =theZdaughter.at(0).LV();
	    TauMuEventFit =theZdaughter.at(1).LV();

	    LVPEF_TauMu = theZdaughter.at(0);
	    LVPEF_TauA1 = theZdaughter.at(1);
	      

	    NeutrinoMu=TauMuEventFit-Ntp->Muons_p4(MuonCandidate);
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
	TauRecHelper Helper;
	TLorentzVector EventFitTauA1 = TauA1_EF.at(AmbiguitySolution);
	TLorentzVector EventFitTauMu = TauMu_EF.at(AmbiguitySolution);
	
	
	TauA1Pt.at(t).Fill(EventFitTauA1.Pt(),w);
	TauA1Phi.at(t).Fill(EventFitTauA1.Phi(),w);
	TauA1Eta.at(t).Fill(EventFitTauA1.Eta(),w);
	
	TauMuPt.at(t).Fill(EventFitTauMu.Pt(),w);
	TauMuPhi.at(t).Fill(EventFitTauMu.Phi(),w);
	TauMuEta.at(t).Fill(EventFitTauMu.Eta(),w);

// 	ZPt.at(t).Fill(EventFitTauMu.Pt(),1);
// 	ZPhi.at(t).Fill(EventFitTauMu.Phi(),1);
// 	ZEta.at(t).Fill(EventFitTauMu.Eta(),1);


	//	if(ptbalance)
	  {
	  if(!AmbPoint){
	    EFitTauA1Pt.at(t).Fill(EventFitTauA1.Pt(),w);
	    EFitTauA1Phi.at(t).Fill(EventFitTauA1.Phi(),w);
	    EFitTauA1Eta.at(t).Fill(EventFitTauA1.Eta(),w);
	    
	    EFitTauMuPt.at(t).Fill(EventFitTauMu.Pt(),w);
	    EFitTauMuPhi.at(t).Fill(EventFitTauMu.Phi(),w);
	    EFitTauMuEta.at(t).Fill(EventFitTauMu.Eta(),w);
	    PhiRotationSignificancePhysicalTaus.at(t).Fill(SignificanceOfTauRotation.at(AmbiguitySolution),w);
	    
	  }

	  if(AmbPoint){
	    EFitTauA1PtAmbPoint.at(t).Fill(EventFitTauA1.Pt(),w);
	    EFitTauA1PhiAmbPoint.at(t).Fill(EventFitTauA1.Phi(),w);
	    EFitTauA1EtaAmbPoint.at(t).Fill(EventFitTauA1.Eta(),w);
	    
	    EFitTauMuPtAmbPoint.at(t).Fill(EventFitTauMu.Pt(),w);
	    EFitTauMuPhiAmbPoint.at(t).Fill(EventFitTauMu.Phi(),w);
	    EFitTauMuEtaAmbPoint.at(t).Fill(EventFitTauMu.Eta(),w);
	  }
	  }

      TLorentzVector SSPion1,SSPion2, OSPion;
      unsigned int ssindex2(0);
      
      if(Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,1) == 1 && Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,2) == -1){
	SSPion1 = Ntp->PFTau_PionsP4(TauCandidate,0);
	SSPion2 = Ntp->PFTau_PionsP4(TauCandidate,1);
	OSPion = Ntp->PFTau_PionsP4(TauCandidate,2);

      }
      if(Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,1) == -1 && Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,2) == 1){
	SSPion1 = Ntp->PFTau_PionsP4(TauCandidate,0);
	SSPion2 = Ntp->PFTau_PionsP4(TauCandidate,2);
	OSPion = Ntp->PFTau_PionsP4(TauCandidate,1);
	
      }
      if(Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,1) == -1 && Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,2) == -1){
	SSPion1 = Ntp->PFTau_PionsP4(TauCandidate,1);
	SSPion2 = Ntp->PFTau_PionsP4(TauCandidate,2);
	OSPion = Ntp->PFTau_PionsP4(TauCandidate,0);
      }


	TLorentzVector A1LabFrame = Ntp->PFTau_p4(TauCandidate);
	TLorentzVector EventFitZ = EventFitTauA1+EventFitTauMu;
	



	//	TLorentzVector TruthZ_LabFrame  = TruthTauA1_LabFrame + TruthTauMu_LabFrame;
	TLorentzVector OSPion_ZFrame = Helper.Boost(OSPion,EventFitZ);
	TLorentzVector SSPion1_ZFrame = Helper.Boost(SSPion1,EventFitZ);
	TLorentzVector SSPion2_ZFrame = Helper.Boost(SSPion2,EventFitZ);
	TLorentzVector RecoA1_ZFrame = Helper.Boost(A1LabFrame,EventFitZ);


	TLorentzVector LV12  = OSPion_ZFrame  + SSPion1_ZFrame;
	TLorentzVector LV23  = SSPion1_ZFrame + SSPion2_ZFrame;
	TLorentzVector LV13  = OSPion_ZFrame  + SSPion2_ZFrame;

	RecoCosTheta_CosBeta_ZFrame_Plus.at(t).Fill(Helper.costheta1(OSPion_ZFrame,SSPion1_ZFrame,SSPion2_ZFrame,EventFitZ),Helper.cosbeta(SSPion1_ZFrame,SSPion2_ZFrame,OSPion_ZFrame),w);
	RecoCosTheta_CosBeta_ZFrame_Minus.at(t).Fill(Helper.costheta1(OSPion_ZFrame,SSPion1_ZFrame,SSPion2_ZFrame,EventFitZ),Helper.cosbeta(SSPion1_ZFrame,SSPion2_ZFrame,OSPion_ZFrame),w);
	RecoCosTheta_CosBeta_ZFrame_Zero.at(t).Fill(Helper.costheta1(OSPion_ZFrame,SSPion1_ZFrame,SSPion2_ZFrame,EventFitZ),Helper.cosbeta(SSPion1_ZFrame,SSPion2_ZFrame,OSPion_ZFrame),w);


	TLorentzVector pi1 = SSPion1_ZFrame;
	TLorentzVector pi2 = SSPion2_ZFrame;
	TLorentzVector pi3 = OSPion_ZFrame;
	TLorentzVector a = pi1+pi2+pi3;
	float mtau = 1.777;
	
	float	  cospsi  = Helper.cospsi(pi1,pi2,pi3,EventFitZ);
	float	  sinpsi  = sqrt(1 - cospsi*cospsi);
	float	  sin2psi = 2*sinpsi*cospsi;
	
	float	  sin2gamma = Helper.Sin2Cos2Gamma(pi1,pi2,pi3).at(0);
	float	  cos2gamma = Helper.Sin2Cos2Gamma(pi1,pi2,pi3).at(1);
	
	
	float	  cosbeta=Helper.cosbeta(pi1,pi2,pi3);
	float	  sinbeta = sqrt(1  - cosbeta*cosbeta);
	
	  
	float	  costheta=Helper.costheta1(pi1,pi2,pi3,EventFitZ);
	float	  sintheta = sqrt(1 - costheta*costheta);

	double RR  = mtau*mtau/a.M()/a.M();
	float U = 0.5*(3*cospsi*cospsi - 1)*(1 - RR);
	float V = 0.5*(3*cospsi*cospsi - 1)*(1 + RR)*costheta + 0.5*3*sin2psi*sintheta*sqrt(RR);


	float Wa =Helper.WA(pi1,pi2,pi3,a.M()*a.M());
	float Wc =Helper.WC(pi1,pi2,pi3,a.M()*a.M());
	float Wd =Helper.WD(pi1,pi2,pi3,a.M()*a.M());
	float We =Helper.WE(pi1,pi2,pi3,a.M()*a.M());

	float fa1 = (2  + RR + 0.5*(3*cosbeta*cosbeta - 1)*U)*Wa/3 - 0.5*sinbeta*sinbeta*cos2gamma*U*Wc + 0.5*sinbeta*sinbeta*sin2gamma*U*Wd + cospsi*cosbeta*We;
	float ga1 = (costheta*(RR -2) - 0.5*(3*cosbeta*cosbeta - 1)*V)*Wa/3 + 0.5*sinbeta*sinbeta*cos2gamma*V*Wc - 0.5*sinbeta*sinbeta*sin2gamma*V*Wd -cosbeta*(costheta*cospsi + sintheta*sinpsi*sqrt(RR))*We;

// 	std::cout<<" V: "<<V <<"  U: "<<U<<std::endl;
// 	// 	  std::cout<<"cospsi: "<<cospsi << "  sinpsi: " << sinpsi << "  sin2psi:  "<< sin2psi<< "  sin2gamma: " <<sin2gamma <<"cosbeta: "<<cosbeta << "  sinbeta: " << sinbeta << "  costheta:  "<< costheta<< "  sintheta: " <<sintheta<<std::endl;
// 	// 	  std::cout<<" costheta1: "<<Helper.costheta1(pi1,pi2,pi3,EventFitZ)<<std::endl;
// 	std::cout<<"WA: "<<Wa << "  WC: " << Wc << "  WD:  "<< Wd<< "  WE: " <<We <<std::endl;
// 	std::cout<<"Unitarity: "<<" Wa2 "<<Wa << "  =  " << sqrt(Wc*Wc  + Wd*Wd + We*We) <<std::endl;
	
// 	std::cout<<"f "<<fa1 << "  ga1: " << ga1<<std::endl;
	
	float om = ga1/fa1;
	float tempx;
	if(EventFitTauMu.E()!=0){ tempx = Ntp->Muons_p4(MuonCandidate).E()/EventFitTauMu.E();}
	else{ tempx = -1;}
	float x = 2*tempx-1;
	RecOmegaA1p.at(t).Fill(ga1/fa1,w);
	RecOmegaA1m.at(t).Fill(ga1/fa1,w);
	RecOmegaA1z.at(t).Fill(ga1/fa1,w);

	RecOmegaMup.at(t).Fill(x,w);
	RecOmegaMum.at(t).Fill(x,w);
	RecOmegaMuz.at(t).Fill(x,w);

	RecOmegaCombinedp.at(t).Fill(x+om/(1+x*om),w);
	RecOmegaCombinedm.at(t).Fill(x+om/(1+x*om),w);
	RecOmegaCombinedz.at(t).Fill(x+om/(1+x*om),w);





      }
  }
}




void  TemplatesData::Finish(){
  unsigned int t;
//   //----- selection with tauisolation on  
//    if(Nminus0.at(0).at(0).Integral()!=0){ if(HConfig.GetHisto(false,DataMCType::QCD,t))ScaleAllHistOfType(t,1789.46/Nminus0.at(0).at(t).Integral());}
//    if(Nminus0.at(0).at(2).Integral()!=0){ if(HConfig.GetHisto(false,20,t)){ScaleAllHistOfType(t,752.57/Nminus0.at(0).at(t).Integral());}}
//    if(Nminus0.at(0).at(8).Integral()!=0){ if(HConfig.GetHisto(false,998,t)){ScaleAllHistOfType(t,7149.98/Nminus0.at(0).at(t).Integral());}}
  

//muon pt cut 19  //----- selection with tauisolation on  
//    if(Nminus0.at(0).at(0).Integral()!=0){ if(HConfig.GetHisto(false,DataMCType::QCD,t))ScaleAllHistOfType(t,1555.67/Nminus0.at(0).at(t).Integral());}
//    if(Nminus0.at(0).at(2).Integral()!=0){ if(HConfig.GetHisto(false,20,t)){ScaleAllHistOfType(t,708.637/Nminus0.at(0).at(t).Integral());}}
//    if(Nminus0.at(0).at(8).Integral()!=0){ if(HConfig.GetHisto(false,998,t)){ScaleAllHistOfType(t,6579.22/Nminus0.at(0).at(t).Integral());}}



  //----- selection with tauisolation off  
//   if(Nminus0.at(0).at(0).Integral()!=0){ if(HConfig.GetHisto(false,DataMCType::QCD,t))ScaleAllHistOfType(t,40888.45/Nminus0.at(0).at(t).Integral());}
//   if(Nminus0.at(0).at(2).Integral()!=0){ if(HConfig.GetHisto(false,20,t)){ScaleAllHistOfType(t,17067.8/Nminus0.at(0).at(t).Integral());}}
//   if(Nminus0.at(0).at(8).Integral()!=0){ if(HConfig.GetHisto(false,998,t)){ScaleAllHistOfType(t,12969.86/Nminus0.at(0).at(t).Integral());}}
  
  Selection::Finish();
}



double TemplatesData::TrigEff(Double_t m, double m0, double sigma, double alpha, double n, double norm){
  
  const double sqrtPiOver2 = 1.2533141373;
  const double sqrt2 = 1.4142135624;
  double sig = fabs((double) sigma);
  double t = (m - m0)/sig;
  if(alpha < 0)    t = -t;
  double absAlpha = fabs(alpha/sig);
  double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
  double b = absAlpha - n/absAlpha;
  double ApproxErf;
  double arg = absAlpha / sqrt2;
  if (arg > 5.) ApproxErf = 1;
  else if (arg < -5.) ApproxErf = -1;
  else ApproxErf = TMath::Erf(arg);
  double leftArea = (1 + ApproxErf) * sqrtPiOver2;
  double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
  double area = leftArea + rightArea;
  if( t <= absAlpha ){
    arg = t / sqrt2;
    if(arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
  }
  else{
    return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
  }
 
}

//       TLorentzVector pion1 = Ntp->KFTau_Initial_pions(FirstLeadingTau,0);
//       TLorentzVector pion2 = Ntp->KFTau_Initial_pions(FirstLeadingTau,1);
//       TLorentzVector pion3 = Ntp->KFTau_Initial_pions(FirstLeadingTau,2);
//       TLorentzVector AlternatA1 = pion1+ pion2 + pion3;

//       TLorentzVector tau_3pi0 = Ntp->KFTau_TauFit_p4(FirstLeadingTau,0);
//       TLorentzVector tau_3pi1 = Ntp->KFTau_TauFit_p4(FirstLeadingTau,1);
//       TLorentzVector tau_3pi2 = Ntp->KFTau_TauFit_p4(FirstLeadingTau,2);

//       TLorentzVector tau_3piVis0 = Ntp->KFTau_TauVis_p4(FirstLeadingTau,0);
//       TLorentzVector tau_3piVis1 = Ntp->KFTau_TauVis_p4(FirstLeadingTau,1);
//       TLorentzVector tau_3piVis2 = Ntp->KFTau_TauVis_p4(FirstLeadingTau,2);
//       TLorentzVector muon    = Ntp->Muons_p4(FirstLeadingMuon);
//       TLorentzVector PtbalancedTau;
//       int BalancedAmbiguity;
//       bool balance = false;
//       if(fabs(tau_3pi1.Pt() - muon.Pt()) <= fabs(tau_3pi2.Pt() - muon.Pt()) ){ PtbalancedTau = tau_3pi1;balance = true;BalancedAmbiguity = 1;}
//       else if(fabs(tau_3pi1.Pt() - muon.Pt()) > fabs(tau_3pi2.Pt() - muon.Pt())  ){ PtbalancedTau = tau_3pi2; balance = true;BalancedAmbiguity = 2;}

//       std::cout<<" +  "<<tau_3pi2.Pt() <<std::endl;
//       std::cout<<" -  "<<tau_3pi1.Pt() <<std::endl;
//       std::cout<<" 0  "<<tau_3pi0.Pt() <<std::endl;


//       TLorentzVector a1Fit= Ntp->KFTau_TauVis_p4(FirstLeadingTau,BalancedAmbiguity);

  
//       TLorentzVector tau_3pi = PtbalancedTau;//Ntp->KFTau_TauFit_p4(FirstLeadingTau,0);
//       double EnergyFraction;
//       if(tau_3pi.Et()!=0)  EnergyFraction = a1Fit.Et()/tau_3pi.Et();
      
//       double DeltaPhi2 = fabs(tau_3pi.Phi() - muon.Phi())/2;



//       double axisPhi = tau_3pi.Phi() - DeltaPhi2;

//       if(tau_3pi.Phi() - DeltaPhi2 ==  muon.Phi() + DeltaPhi2){axisPhi=tau_3pi.Phi() - DeltaPhi2;}
//       else if(tau_3pi.Phi() + DeltaPhi2 == muon.Phi() - DeltaPhi2){axisPhi = tau_3pi.Phi() + DeltaPhi2;}
//       a1MassAlternat .at(t).Fill(AlternatA1.M(),1); 
//       a1Mass.at(t).Fill(a1Fit.M(),1); 
//       Efrac.at(t).Fill(EnergyFraction,w); 



//       double GJ0    = acos(  (tau_3pi0.Px()*tau_3piVis0.Px()  + tau_3pi0.Py()*tau_3piVis0.Py() + tau_3pi0.Pz()*tau_3piVis0.Pz())/ tau_3pi0.P()/tau_3piVis0.P() );
//       double GJ1    = acos(  (tau_3pi1.Px()*tau_3piVis1.Px()  + tau_3pi1.Py()*tau_3piVis1.Py() + tau_3pi1.Pz()*tau_3piVis1.Pz())/ tau_3pi1.P()/tau_3piVis1.P() );
//       double GJ2    = acos(  (tau_3pi2.Px()*tau_3piVis2.Px()  + tau_3pi2.Py()*tau_3piVis2.Py() + tau_3pi2.Pz()*tau_3piVis2.Pz())/ tau_3pi2.P()/tau_3piVis2.P() );
      
//       double MaxGJ0 =  asin( (1.777*1.777 - tau_3piVis0.M())/2/1.777/tau_3piVis0.P());
//       double MaxGJ1 =  asin( (1.777*1.777 - tau_3piVis1.M())/2/1.777/tau_3piVis1.P());


//       double MaxGJ2 =  asin( (1.777*1.777 - tau_3piVis2.M())/2/1.777/tau_3piVis2.P());
 
//       hGJ0.at(t).Fill( GJ0,1); 
//       hGJ1.at(t).Fill( GJ1,1); 
//       hGJ2.at(t).Fill( GJ2,1); 
 
//       rhGJ0.at(t).Fill( GJ0/MaxGJ0, 1); 
//       rhGJ1.at(t).Fill( GJ1/MaxGJ1 ,1); 
//       rhGJ2.at(t).Fill( GJ2/MaxGJ2 ,1); 

//       TLorentzVector tau,Z,ZVis,ZAmb;
//       double zmass = 91.1876,mtau = 1.777;
      

//       double energy;
//       double theta = muon.Theta();
      
//       if(sin(theta)!=0){
// 	energy = sqrt(mtau*mtau*sin(theta)*sin(theta) + tau_3pi.Pt()*tau_3pi.Pt())/sin(theta);
//       }else energy = -1;
      
//       tau.SetXYZM(-tau_3pi.Px(), -tau_3pi.Py(),sqrt(energy*energy - mtau*mtau)*cos(theta),mtau);
      
//       TauEnergy.at(t).Fill(tau_3pi.E(),w); 
//       MuonEnergy.at(t).Fill(muon.E(),w); 

//       Z  = tau_3pi + muon;
//       ZAmb  = tau_3pi0 + muon;
//       ZVis  = Ntp->PFTau_p4(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau)) + muon;

//       double psi = Z.Pt()*cos(axisPhi - Z.Phi());
//       double psiMVS = Ntp->MET_CorrMVA_et()*cos(axisPhi - Ntp->MET_CorrMVA_phi());
//       Pzeta.at(t).Fill(psiMVS - 0.85*psi,w); 

//       //  deltaPhi.at(t).Fill(abs(1  - abs(muon.Phi() - tau_3pi.Phi())/3.1415));
//       if(Ntp->KFTau_IsAmbiguityValid(FirstLeadingTau,BalancedAmbiguity)==1)        ZMass.at(t).Fill(Z.M(),w); 
//       if(Ntp->KFTau_IsAmbiguityValid(FirstLeadingTau,0)==1)      ZMassAmbPoint.at(t).Fill(ZAmb.M(),w); 
//       VisibleMass.at(t).Fill(ZVis.M(),w); 
//       HPSPhi.at(t).Fill(Ntp->PFTau_p4(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau)).Phi(),w); 
//       tauPttest = Z.M();
      
//       FlightLenght.at(t).Fill( (Ntp->KFTau_Fit_PrimaryVertex(FirstLeadingTau) - Ntp->KFTau_Fit_SecondaryVertex(FirstLeadingTau)).Mag(),w); 

//       MuonDistance.at(t).Fill( (Ntp->KFTau_Fit_PrimaryVertex(FirstLeadingTau) - Ntp->Muon_Poca(FirstLeadingMuon)).Mag(),w); 

//       MuonPhi.at(t).Fill(muon.Phi(),w); 
//       MuonEta.at(t).Fill(muon.Eta(),w); 
      
//       TauPhi.at(t).Fill(tau_3pi.Phi(),w); 



//       TauPhiPlus.at(t).Fill(Ntp->KFTau_TauFit_p4(FirstLeadingTau,1).Phi(),w);
//       TauPhiMinus.at(t).Fill(Ntp->KFTau_TauFit_p4(FirstLeadingTau,2).Phi(),w);
//       TauPhiZero.at(t).Fill(Ntp->KFTau_TauFit_p4(FirstLeadingTau,0).Phi(),w);



//       TauEta.at(t).Fill(tau_3pi.Eta(),w); 
//       MuTauDeltaR.at(t).Fill(tau_3pi.DeltaR(muon),w); 
//       deltaPhi.at(t).Fill(fabs(tau_3pi.Phi() - muon.Phi() ),w); 
//       //---------------------------------   Reco atumu---------------------------------------------------------------------------------


//       TLorentzVector tauCand1,tauCand2,taumu1(0,0,0,0),taumu2(0,0,0,0),ZFull;
//       TLorentzVector TauMu11(0,0,0,0),TauMu12(0,0,0,0),Tau3pi,BestTau,WorstTau;
//       TLorentzVector z1,z2;
// //     double ma1 = 1.23;
// //     double zmass = 91.187,mtau = 1.777;

//       TauRecHelper RcTH0(tau_3pi0, muon);
//       TauRecHelper RcTH1(PtbalancedTau, muon);

//       std::cout<<" bal  "<<PtbalancedTau.Pt() <<std::endl;


//       if(RcTH1.isTauRecoValid()){
//       taumu1 = RcTH1.Get1TauSolution();
//       taumu2 = RcTH1.Get2TauSolution();
//       taumupt.at(t).Fill(taumu1.Pt(),w); 
//       taumuphi.at(t).Fill(taumu1.Phi(),w); 
//       taumutheta.at(t).Fill(taumu1.Eta(),w); 
//       x.at(t).Fill(muon.E()/taumu1.E(),w);


//       }


//       if(RcTH0.isTauRecoValid()){
//       TauMu11 = RcTH0.Get1TauSolution();
//       TauMu12 = RcTH0.Get2TauSolution();
//       xamb.at(t).Fill(muon.E()/TauMu11.E(),w);
//       }

//       tauptamb.at(t).Fill(tau_3pi0.Pt(),w);
//       std::cout<<"tau reco balanced Pt: "<< taumu1.Pt()<<" ambpoint taupt: "<<TauMu11.Pt()<<std::endl;

//       z1 = taumu1+PtbalancedTau;
//       z2 = taumu2+PtbalancedTau;

 
//       DeltaTheta_1st.at(t).Fill(taumu1.Theta() - muon.Theta(),w); 
//       DeltaTheta_2nd.at(t).Fill(taumu2.Theta() - muon.Theta(),w); 
//       EMiss.at(t).Fill(Ntp->MET_CorrMVA_et(),1); 
//       //------------------------------------------------------------------------------------------------------------------

//       TLorentzVector PlusPiFor(0,0,0,0), MinusPiFor(0,0,0,0), a1LVFor(0,0,0,0);
//       int nplusFor =0,nminusFor=0;
//       double gammaFor=0;

//       for(unsigned int ip =0; ip < Ntp->KFTau_NDaughter(FirstLeadingTau); ip++){


// 	//	std::cout<<" partPDGID "<<Ntp->KFTau_Daughter_pdgid(FirstLeadingTau,ip)<<"  pt " << Ntp->KFTau_Daughter_p4(FirstLeadingTau,ip).Pt()<<std::endl;
// 	if(Ntp->KFTau_Daughter_pdgid(FirstLeadingTau,ip) == 211){PlusPiFor+=Ntp->KFTau_Daughter_p4(FirstLeadingTau,ip);nplusFor++;a1LVFor+=PlusPiFor; }
// 	if(Ntp->KFTau_Daughter_pdgid(FirstLeadingTau,ip) == -211){MinusPiFor+=Ntp->KFTau_Daughter_p4(FirstLeadingTau,ip);nminusFor++;a1LVFor+=MinusPiFor;}
// 	if(nplusFor+nminusFor==3) break;

	
//       }
//       //    std::cout<<"nplusFor "<<nplusFor<<" nminusFor "<<nminusFor<<" PlusPiFor.Pt()  "<< PlusPiFor.Pt()<< "  a1LVFor.Pt()  "<<a1LVFor.Pt()<<std::endl;
//       if(nplusFor+nminusFor==3 && nplusFor==1) gammaFor=2*PlusPiFor.Pt()/a1LVFor.Pt()-1;
//       if(nplusFor+nminusFor==3 && nminusFor==1) gammaFor=2*MinusPiFor.Pt()/a1LVFor.Pt()-1; 
   
//       gamma.at(t).Fill(gammaFor,w);




//       RecoTauValid.at(t).Fill(RcTH0.isTauRecoValid());
//       //   if(taumu1.E()!=0)xmuon.at(t).Fill(muon.E()/taumu1.E(),w);  
      
//       //     TauPt.at(t).Fill(tau.Pt(),w); 
//       //     xmuon.at(t).Fill(muon.E()/tau.E(),w); 
//       //     xmuonT.at(t).Fill(muon.Et()/tau.Et(),w); 
//       //       std::cout<<"Tau vertex x,y,z  "<<Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).X()<<"  "<<Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).Y()<<"  "<<Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).Z()<<std::endl;
//       //       std::cout<<"muon track index  "<<Ntp->Muon_Track_idx(FirstLeadingMuon)<<std::endl;
//       //       bool matched = false;
//       //       for(int iVert = 0; iVert < Ntp->NVtx(); iVert++){
//       // // 	std::cout<<"Vertex position x,y,z  "<< Ntp->Vtx(iVert).X()<<"  "<<Ntp->Vtx(iVert).Y()<<"  "<<Ntp->Vtx(iVert).Z()<<std::endl; 
//       // // 	std::cout<<"N tracks per vertex  "<<Ntp->Vtx_Track_idx(iVert).size()<<std::endl;
      
//       // 	for(int iTrack = 0; iTrack <Ntp->Vtx_Track_idx(iVert).size(); iTrack++){
//       // 	  //	std::cout<<"track index  "<<Ntp->Vtx_Track_idx(iVert).at(iTrack)<<std::endl;
      
//       // 	  //	if(Ntp->Vtx_Track_idx(iVert).at(iTrack) == Ntp->Muon_Track_idx(FirstLeadingMuon)){std::cout<<"Vertex position x,y,z  "<< Ntp->Vtx(iVert).X()<<"  "<<Ntp->Vtx(iVert).Y()<<"  "<<Ntp->Vtx(iVert).Z()<<std::endl;matched = true;}
      
//       // 	}
      
      
//       //       }
//       // //       if(matched)std::cout<<" mathced  "<<std::endl;
//       //       else std::cout<<"not  mathced  "<<std::endl;
      
      
      
//       //       if(id == 33){
//       //           for(int iz =0; iz<Ntp->NMCSignalParticles(); iz++){
//       // 	    if(Ntp->MCSignalParticle_Tauidx(iz).size()!=0){
//       // 	      Jak1VsJak2.at(t).Fill(Ntp->MCTau_JAK(0),Ntp->MCTau_JAK(1));
//       // 	      std::cout<<"Jak1 = "<< Ntp->MCTau_JAK(0)<< " Jak2 =  "<< Ntp->MCTau_JAK(1)<<std::endl;
//       // 	    }
//       // 	  }
//       //     }
      
//     }
//   }
// }


