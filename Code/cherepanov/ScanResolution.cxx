#include "ScanResolution.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"
#include "Tauola.h"
#include "LHAPDF/LHAPDF.h"
#include "tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"
#include<iostream>
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"





ScanResolution::ScanResolution(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

ScanResolution::~ScanResolution(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ScanResolution::~ScanResolution Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ScanResolution::~ScanResolution()" << std::endl;
}

void  ScanResolution::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==MuonisGlob)         cut.at(MuonisGlob)=1;
    if(i==TauIsQuality)       cut.at(TauIsQuality)=1;
    if(i==MuonPt)             cut.at(MuonPt)=16; //18
    if(i==TauPtCut)           cut.at(TauPtCut)=25;
    if(i==MET)                cut.at(MET)=40;
    if(i==MuonIso)            cut.at(MuonIso)=2.7;
    //    if(i==MuonIso)            cut.at(MuonIso)=0.1;
    if(i==TauIsIso)           cut.at(TauIsIso)=1;
    if(i==charge)             cut.at(charge)=-1;
    if(i==TauDecay)           cut.at(TauDecay)=10;
    if(i==NMuon)              cut.at(NMuon)=1;
    if(i==NTau)               cut.at(NTau)=1;



 

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
    else if(i==MuonisGlob){
      title.at(i)="Muon is global ";
      hlabel="MuonisGlob";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonisGlob_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonisGlob_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==TauIsQuality){
      title.at(i)="TauIsQuality ";
      hlabel="TauIsQuality ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsQuality_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsQuality_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
   else if(i==MuonPt){
      title.at(i)="$MuonPt > $";
      title.at(i)+=cut.at(MuonPt);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#mu_{pT}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonPt_",htitle,20,0,50,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonPt_",htitle,20,0,50,hlabel,"Events"));
    }
   else if(i==TauPtCut){
      title.at(i)="$TauPtCut > $";
      title.at(i)+=cut.at(TauPtCut);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#tau_{pT}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauPtCut_",htitle,20,0,50,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauPtCut_",htitle,20,0,50,hlabel,"Events"));
    }
   else if(i==MET){
      title.at(i)="$MET < $";
      title.at(i)+=cut.at(MET);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="MET";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MET_",htitle,10,0,60,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MET_",htitle,10,0,60,hlabel,"Events"));
    }  
   else if(i==MuonIso){
      title.at(i)="$MuonIso < $";
      title.at(i)+=cut.at(MuonIso);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#MuonIso";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIso_",htitle,30,0,10,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIso_",htitle,30,0,10,hlabel,"Events"));

//       Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIso_",htitle,30,0,1,hlabel,"Events"));
//       Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIso_",htitle,30,0,1,hlabel,"Events"));
    }
   else if(i==charge){
      title.at(i)="$opposite charge$";
      title.at(i)+=cut.at(charge);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="charge";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,3,-1.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,3,-1.5,1.5,hlabel,"Events"));
    } 

   else if(i==TauIsIso){
      title.at(i)="$TauIsIso == $";
      title.at(i)+=cut.at(TauIsIso);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsIso";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIso_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIso_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
   else if(i==TauDecay){
      title.at(i)="$TauDecay == $";
      title.at(i)+=cut.at(TauDecay);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauDecay";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauDecay_",htitle,11,-0.5,10.5,hlabel,""));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauDecay_",htitle,11,-0.5,10.5,hlabel,""));
    } 
   else if(i==NMuon){
      title.at(i)="$NMuon == $";
      title.at(i)+=cut.at(NMuon);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="NMuon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuon_",htitle,5,-0.5,4.5,hlabel,""));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuon_",htitle,5,-0.5,4.5,hlabel,""));
    } 
   else if(i==NTau){
      title.at(i)="$NTau == $";
      title.at(i)+=cut.at(NTau);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="NTau";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTau_",htitle,5,-0.5,4.5,hlabel,""));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTau_",htitle,5,-0.5,4.5,hlabel,""));
    } 



  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Accumulative Cuts Passed","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  Tau3PiPt=HConfig.GetTH1D(Name+"_Tau3PiPt","Tau3PiPt",30,0,100,"Tau3PiPt","Events");
  Tau1Pt =HConfig.GetTH1D(Name+"_Tau1Pt","Tau1Pt",20,0,80,"Tau1Pt","Events");
  Tau2Pt =HConfig.GetTH1D(Name+"_Tau2Pt","Tau2Pt",20,0,80,"Tau2Pt","Events");
  Tau1theta=HConfig.GetTH1D(Name+"_Tau1theta","Tau1theta",20,0,3.14,"Tau1theta","Events");
  Tau2theta=HConfig.GetTH1D(Name+"_Tau2theta","Tau2theta",20,0,3.14,"Tau2theta","Events");  
  Tau1Deltatheta=HConfig.GetTH1D(Name+"_Tau1Deltatheta","Tau1Deltatheta",20,0,3.14,"Tau1Deltatheta","Events");
  Tau2Deltatheta=HConfig.GetTH1D(Name+"_Tau2Deltatheta","Tau2Deltatheta",20,0,3.14,"Tau2Deltatheta","Events");  
  DeltaTheta_1st=HConfig.GetTH1D(Name+"_DeltaTheta_1st","DeltaTheta_1st",20,0,3.14,"DeltaTheta_1st","Events");
  DeltaTheta_2nd=HConfig.GetTH1D(Name+"_DeltaTheta_2nd","DeltaTheta_2nd",20,0,3.14,"DeltaTheta_2nd","Events");
  Energy_resolution1=HConfig.GetTH1D(Name+"_Energy_resolution1","Energy_resolution1",30,-100,100,"Energy_resolution1","Events");
  Energy_resolution2=HConfig.GetTH1D(Name+"_Energy_resolution2","Energy_resolution2",30,-100,100,"Energy_resolution2","Events");
  TransverseEnergy_resolution=HConfig.GetTH1D(Name+"_TransverseEnergy_resolution","TransverseEnergy_resolution",30,100,100,"TransverseEnergy_resolution","Events");


  res0=HConfig.GetTH1D(Name+"_res0","res0",30,-100,100,"res0","");
  res5=HConfig.GetTH1D(Name+"_res5","res5",30,-100,100,"res5","");
  res10=HConfig.GetTH1D(Name+"_res10","res10",30,-100,100,"res10","");
  res15=HConfig.GetTH1D(Name+"_res15","res15",30,-100,100,"res15","");
  res20=HConfig.GetTH1D(Name+"_res20","res20",30,-100,100,"res20","");
  res30=HConfig.GetTH1D(Name+"_res30","res30",30,-100,100,"res30","");
  res40=HConfig.GetTH1D(Name+"_res40","res40",30,-100,100,"res40","");
  res50=HConfig.GetTH1D(Name+"_res50","res50",30,-100,100,"res50","");
  res60=HConfig.GetTH1D(Name+"_res60","res60",30,-100,100,"res60","");
  res70=HConfig.GetTH1D(Name+"_res70","res70",30,-100,100,"res70","");
  res80=HConfig.GetTH1D(Name+"_res80","res80",30,-100,100,"res80","");
  res90=HConfig.GetTH1D(Name+"_res90","res90",30,-100,100,"res90","");
  res100=HConfig.GetTH1D(Name+"_res100","res100",30,-100,100,"res100","");

  res110=HConfig.GetTH1D(Name+"_res110","res110",30,-100,100,"res110","");
  res120=HConfig.GetTH1D(Name+"_res120","res120",30,-100,100,"res120","");
  res150=HConfig.GetTH1D(Name+"_res150","res150",30,-100,100,"res150","");


  xMC=HConfig.GetTH1D(Name+"_xMC","xMC",14,0.,1.4,"xMC","");
  xRec=HConfig.GetTH1D(Name+"_xRec","xRec",14,0.,1.4,"xRec","");
  xTruth=HConfig.GetTH1D(Name+"_xTruth","xTruth",14,0.,1.4,"xTruth","");

  xp0=HConfig.GetTH1D(Name+"_xp0","xp0",14,0.,1.4,"xp0","");
  xp5=HConfig.GetTH1D(Name+"_xp5","xp5",14,0.,1.4,"xp5","");
  xp10=HConfig.GetTH1D(Name+"_xp10","xp10",14,0.,1.4,"xp10","");
  xp15=HConfig.GetTH1D(Name+"_xp15","xp15",14,0.,1.4,"xp15","");
  xp20=HConfig.GetTH1D(Name+"_xp20","xp20",14,0.,1.4,"xp20","");
  xp30=HConfig.GetTH1D(Name+"_xp30","xp30",14,0.,1.4,"xp30","");
  xp40=HConfig.GetTH1D(Name+"_xp40","xp40",14,0.,1.4,"xp40","");
  xp50=HConfig.GetTH1D(Name+"_xp50","xp50",14,0.,1.4,"xp50","");
  xp60=HConfig.GetTH1D(Name+"_xp60","xp60",14,0.,1.4,"xp60","");
  xp70=HConfig.GetTH1D(Name+"_xp70","xp70",14,0.,1.4,"xp70","");
  xp80=HConfig.GetTH1D(Name+"_xp80","xp80",14,0.,1.4,"xp80","");
  xp90=HConfig.GetTH1D(Name+"_xp90","xp90",14,0.,1.4,"xp90","");
  xp100=HConfig.GetTH1D(Name+"_xp100","xp100",14,0.,1.4,"xp100","");
  xp110=HConfig.GetTH1D(Name+"_xp110","xp110",14,0.,1.4,"xp110","");
  xp120=HConfig.GetTH1D(Name+"_xp120","xp120",14,0.,1.4,"xp120","");
  xp150=HConfig.GetTH1D(Name+"_xp150","xp150",14,0.,1.4,"xp150","");

  xm0=HConfig.GetTH1D(Name+"_xm0","xm0",14,0.,1.4,"xm0","");
  xm5=HConfig.GetTH1D(Name+"_xm5","xm5",14,0.,1.4,"xm5","");
  xm10=HConfig.GetTH1D(Name+"_xm10","xm10",14,0.,1.4,"xm10","");
  xm15=HConfig.GetTH1D(Name+"_xm15","xm15",14,0.,1.4,"xm15","");
  xm20=HConfig.GetTH1D(Name+"_xm20","xm20",14,0.,1.4,"xm20","");
  xm30=HConfig.GetTH1D(Name+"_xm30","xm30",14,0.,1.4,"xm30","");
  xm40=HConfig.GetTH1D(Name+"_xm40","xm40",14,0.,1.4,"xm40","");
  xm50=HConfig.GetTH1D(Name+"_xm50","xm50",14,0.,1.4,"xm50","");
  xm60=HConfig.GetTH1D(Name+"_xm60","xm60",14,0.,1.4,"xm60","");
  xm70=HConfig.GetTH1D(Name+"_xm70","xm70",14,0.,1.4,"xm70","");
  xm80=HConfig.GetTH1D(Name+"_xm80","xm80",14,0.,1.4,"xm80","");
  xm90=HConfig.GetTH1D(Name+"_xm90","xm90",14,0.,1.4,"xm90","");
  xm100=HConfig.GetTH1D(Name+"_xm100","xm100",14,0.,1.4,"xm100","");
  xm110=HConfig.GetTH1D(Name+"_xm110","xm110",14,0.,1.4,"xm110","");
  xm120=HConfig.GetTH1D(Name+"_xm120","xm120",14,0.,1.4,"xm120","");
  xm150=HConfig.GetTH1D(Name+"_xm150","xm150",14,0.,1.4,"xm150","");

  xReco0=HConfig.GetTH1D(Name+"_xReco0","xReco0",14,0.,1.4,"xReco0","");
  xReco5=HConfig.GetTH1D(Name+"_xReco5","xReco5",14,0.,1.4,"xReco5","");
  xReco10=HConfig.GetTH1D(Name+"_xReco10","xReco10",14,0.,1.4,"xReco10","");
  xReco15=HConfig.GetTH1D(Name+"_xReco15","xReco15",14,0.,1.4,"xReco15","");
  xReco20=HConfig.GetTH1D(Name+"_xReco20","xReco20",14,0.,1.4,"xReco20","");
  xReco30=HConfig.GetTH1D(Name+"_xReco30","xReco30",14,0.,1.4,"xReco30","");
  xReco40=HConfig.GetTH1D(Name+"_xReco40","xReco40",14,0.,1.4,"xReco40","");
  xReco50=HConfig.GetTH1D(Name+"_xReco50","xReco50",14,0.,1.4,"xReco50","");
  xReco60=HConfig.GetTH1D(Name+"_xReco60","xReco60",14,0.,1.4,"xReco60","");
  xReco70=HConfig.GetTH1D(Name+"_xReco70","xReco70",14,0.,1.4,"xReco70","");
  xReco80=HConfig.GetTH1D(Name+"_xReco80","xReco80",14,0.,1.4,"xReco80","");
  xReco90=HConfig.GetTH1D(Name+"_xReco90","xReco90",14,0.,1.4,"xReco90","");
  xReco100=HConfig.GetTH1D(Name+"_xReco100","xReco100",14,0.,1.4,"xReco100","");
  xReco110=HConfig.GetTH1D(Name+"_xReco110","xReco110",14,0.,1.4,"xReco110","");
  xReco120=HConfig.GetTH1D(Name+"_xReco120","xReco120",14,0.,1.4,"xReco120","");
  xReco150=HConfig.GetTH1D(Name+"_xReco150","xReco150",14,0.,1.4,"xReco150","");




  xmuon_plus=HConfig.GetTH1D(Name+"_xmuon_plus","xmuon_plus",14,0.,1.4,"xmuon_plus","Events");
  xmuon_minus=HConfig.GetTH1D(Name+"_xmuon_minus","xmuon_minus",14,0.,1.4,"xmuon_minus","Events");

  xmuonT_plus=HConfig.GetTH1D(Name+"_xmuonT_plus","xmuonT_plus",15,0,2,"xmuonT_plus","Events");
  xmuonT_minus=HConfig.GetTH1D(Name+"_xmuonT_minus","xmuonT_minus",15,0,2,"xmuonT_minus","Events");

  xmuonTruthT_plus=HConfig.GetTH1D(Name+"_xmuonTruthT_plus","xmuonTruthT_plus",15,0,2,"xmuonTruthT_plus","Events");
  xmuonTruthT_minus=HConfig.GetTH1D(Name+"_xmuonTruthT_minus","xmuonTruthT_minus",15,0,2,"xmuonTruthT_minus","Events");

  xmuonTruth_plus=HConfig.GetTH1D(Name+"_xmuonTruth_plus","xmuonTruth_plus",14,0.,1.4,"xmuonTruth_plus","Events");
  xmuonTruth_minus=HConfig.GetTH1D(Name+"_xmuonTruth_minus","xmuonTruth_minus",14,0.,1.4,"xmuonTruth_minus","Events");

  xmuonMC_plus=HConfig.GetTH1D(Name+"_xmuonMC_plus","xmuonMC_plus",14,0.,1.4,"xmuonMC_plus","Events");
  xmuonMC_minus=HConfig.GetTH1D(Name+"_xmuonMC_minus","xmuonMC_minus",14,0.,1.4,"xmuonMC_minus","Events");

  xmuonTMC_plus=HConfig.GetTH1D(Name+"_xmuonTMC_plus","xmuonTMC_plus",15,0,2,"xmuonTMC_plus","Events");
  xmuonTMC_minus=HConfig.GetTH1D(Name+"_xmuonTMC_minus","xmuonTMC_minus",15,0,2,"xmuonTMC_minus","Events");

  IsoVsMuonEnergy=HConfig.GetTH2D(Name+"_IsoVsMuonEnergy","IsoVsMuonEnergy",55,15,60,30,0,1,"IsoVsMuonEnergy","");
  IsoVsMuonTauEnergy=HConfig.GetTH2D(Name+"_IsoVsMuonTauEnergy","IsoVsMuonTauEnergy",55,15,60,30,0,1,"IsoVsMuonTauEnergy","");

  IsoAbsVsMuonEnergy=HConfig.GetTH2D(Name+"_IsoAbsVsMuonEnergy","IsoAbsVsMuonEnergy",30,0,30,30,15,45,"IsoAbsVsMuonEnergy","");

  TruthResolution=HConfig.GetTH1D(Name+"_TruthResolution","TruthResolution",30,-100,100,"TruthResolution","Events");
  TruthDeltaTheta=HConfig.GetTH1D(Name+"_TruthDeltaTheta","TruthDeltaTheta",20,0,3.14,"TruthDeltaTheta","Events");





  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}


// TString name,TString title,int nbinsx, double minx, double maxx, 
// 				       int nbinsy, double miny, double maxy, TString xaxis, TString yaxis)

void  ScanResolution::Store_ExtraDist(){

 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist1d.push_back(&Tau3PiPt);
 Extradist1d.push_back(&Tau1Pt);
 Extradist1d.push_back(&Tau2Pt);
 Extradist1d.push_back(&Tau1theta);	    
 Extradist1d.push_back(&Tau2theta);      
 Extradist1d.push_back(&Tau1Deltatheta); 
 Extradist1d.push_back(&Tau2Deltatheta); 
 Extradist1d.push_back(&DeltaTheta_1st);
 Extradist1d.push_back(&DeltaTheta_2nd);
 Extradist1d.push_back(&Energy_resolution1);
 Extradist1d.push_back(&Energy_resolution2);
 Extradist1d.push_back(&TransverseEnergy_resolution);

 Extradist1d.push_back(&xmuon_plus);	      
 Extradist1d.push_back(&xmuonT_plus);      
 Extradist1d.push_back(&xmuonTruthT_plus); 
 Extradist1d.push_back(&xmuonTruth_plus);  
 Extradist1d.push_back(&xmuonMC_plus);     
 Extradist1d.push_back(&xmuonTMC_plus);    
				      
 Extradist1d.push_back(&xmuon_minus);      
 Extradist1d.push_back(&xmuonT_minus);     
 Extradist1d.push_back(&xmuonTruthT_minus);
 Extradist1d.push_back(&xmuonTruth_minus); 
 Extradist1d.push_back(&xmuonMC_minus);    
 Extradist1d.push_back(&xmuonTMC_minus);   

 Extradist1d.push_back(&TruthResolution);
 Extradist1d.push_back(&TruthDeltaTheta);

 Extradist2d.push_back(&IsoVsMuonEnergy);
 Extradist2d.push_back(&IsoVsMuonTauEnergy);
 Extradist2d.push_back(&IsoAbsVsMuonEnergy);

 Extradist1d.push_back(&xp0);  
 Extradist1d.push_back(&xp5);  
 Extradist1d.push_back(&xp10); 
 Extradist1d.push_back(&xp15); 
 Extradist1d.push_back(&xp20); 
 Extradist1d.push_back(&xp30); 
 Extradist1d.push_back(&xp40); 
 Extradist1d.push_back(&xp50); 
 Extradist1d.push_back(&xp60); 
 Extradist1d.push_back(&xp70); 
 Extradist1d.push_back(&xp80); 
 Extradist1d.push_back(&xp90); 
 Extradist1d.push_back(&xp100);
 Extradist1d.push_back(&xp110);
 Extradist1d.push_back(&xp120);
 Extradist1d.push_back(&xp150);

 Extradist1d.push_back(&xm0);  
 Extradist1d.push_back(&xm5);  
 Extradist1d.push_back(&xm10); 
 Extradist1d.push_back(&xm15); 
 Extradist1d.push_back(&xm20); 
 Extradist1d.push_back(&xm30); 
 Extradist1d.push_back(&xm40); 
 Extradist1d.push_back(&xm50); 
 Extradist1d.push_back(&xm60); 
 Extradist1d.push_back(&xm70); 
 Extradist1d.push_back(&xm80); 
 Extradist1d.push_back(&xm90); 
 Extradist1d.push_back(&xm100);
 Extradist1d.push_back(&xm110);
 Extradist1d.push_back(&xm120);
 Extradist1d.push_back(&xm150);

 Extradist1d.push_back(&res0);  
 Extradist1d.push_back(&res5);  
 Extradist1d.push_back(&res10); 
 Extradist1d.push_back(&res15); 
 Extradist1d.push_back(&res20); 
 Extradist1d.push_back(&res30); 
 Extradist1d.push_back(&res40); 
 Extradist1d.push_back(&res50); 
 Extradist1d.push_back(&res60); 
 Extradist1d.push_back(&res70); 
 Extradist1d.push_back(&res80); 
 Extradist1d.push_back(&res90); 
 Extradist1d.push_back(&res100);
 Extradist1d.push_back(&res110);
 Extradist1d.push_back(&res120);
 Extradist1d.push_back(&res150);


 Extradist1d.push_back(&xReco0);
 Extradist1d.push_back(&xReco5);
 Extradist1d.push_back(&xReco10);
 Extradist1d.push_back(&xReco15);
 Extradist1d.push_back(&xReco20);
 Extradist1d.push_back(&xReco30);
 Extradist1d.push_back(&xReco40);
 Extradist1d.push_back(&xReco50);
 Extradist1d.push_back(&xReco60);
 Extradist1d.push_back(&xReco70);
 Extradist1d.push_back(&xReco80);
 Extradist1d.push_back(&xReco90);
 Extradist1d.push_back(&xReco100);
 Extradist1d.push_back(&xReco110);
 Extradist1d.push_back(&xReco120);
 Extradist1d.push_back(&xReco150);


 Extradist1d.push_back(&xMC);
 Extradist1d.push_back(&xRec);
 Extradist1d.push_back(&xTruth);





}

void  ScanResolution::doEvent(){
 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  
  // Apply Selection
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
 

  unsigned int HighestPtMuonIndex=0;
  unsigned int HighestPtTauIndex=0;
  unsigned int SecondPtMuonIndex =0;
  float muonPt=0;
  float tauPt=0;
  unsigned int iTau=0;
  unsigned int iMuon=0;
  unsigned int nTau=0,nMuon=0;

  for(iTau=0;iTau < Ntp->NKFTau();iTau++){
    if(Ntp->KFTau_discriminatorByQC(iTau))nTau++;
    if(Ntp->KFTau_TauFit_p4(iTau).Pt() > tauPt){
      tauPt = Ntp->KFTau_TauFit_p4(iTau).Pt();
      HighestPtTauIndex = iTau; 
    }
  }
  
  if(Ntp->NMuons()!=0){
    for(iMuon=0; iMuon<Ntp->NMuons(); iMuon++){
      if(Ntp->Muon_isGlobalMuon(iMuon))nMuon++;
      if(Ntp->Muons_p4(iMuon).Pt() > muonPt){
	muonPt = Ntp->Muons_p4(iMuon).Pt();
	SecondPtMuonIndex = HighestPtMuonIndex;
	HighestPtMuonIndex = iMuon;
      }
    }
  }
  
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  value.at(TriggerOk)=1;
  pass.at(TriggerOk)=true;

  if(Ntp->NKFTau()!=0 && Ntp->NMuons()!=0){

  value.at(MuonisGlob) = Ntp->Muon_isGlobalMuon(HighestPtMuonIndex);
  pass.at(MuonisGlob)  = (value.at(MuonisGlob)==cut.at(MuonisGlob));

  value.at(TauIsQuality) =  Ntp->KFTau_discriminatorByQC(HighestPtTauIndex);
  //  value.at(TauIsQuality) =  Ntp->KFTau_discriminatorByKFit(HighestPtTauIndex);//Ntp->KFTau_discriminatorByQC(HighestPtTauIndex);
  pass.at(TauIsQuality)  =  (value.at(TauIsQuality)==cut.at(TauIsQuality));

  value.at(TauPtCut) = Ntp->KFTau_TauFit_p4(HighestPtTauIndex).Pt();
  pass.at(TauPtCut)=true;//(value.at(TauPtCut)>=cut.at(TauPtCut));

  value.at(MuonPt) = Ntp->Muons_p4(HighestPtMuonIndex).Pt();
  pass.at(MuonPt)=true;//(value.at(MuonPt) >=cut.at(MuonPt));

  value.at(MET) =   sqrt(2*Ntp->Muons_p4(HighestPtMuonIndex).Pt()*Ntp->MET_et()*(1  - cos(Ntp->Muons_p4(HighestPtMuonIndex).Phi() - Ntp->MET_phi())) );
  pass.at(MET)=(value.at(MET)<=cut.at(MET));

  value.at(TauIsIso) =   Ntp->PFTau_isMediumIsolation(Ntp->KFTau_MatchedHPS_idx(HighestPtTauIndex));
  pass.at(TauIsIso)=true;//(value.at(TauIsIso) ==cut.at(TauIsIso));


  //  value.at(MuonIso) = (Ntp->Muon_emEt05(HighestPtMuonIndex) + Ntp->Muon_hadEt05(HighestPtMuonIndex) + Ntp->Muon_sumPt05(HighestPtMuonIndex))/Ntp->Muons_p4(HighestPtMuonIndex).Pt();
  value.at(MuonIso) = (Ntp->Muon_emEt05(HighestPtMuonIndex) + Ntp->Muon_hadEt05(HighestPtMuonIndex) + Ntp->Muon_sumPt05(HighestPtMuonIndex));
  pass.at(MuonIso)=true;//(value.at(MuonIso)<=cut.at(MuonIso));

  value.at(charge) =Ntp->KFTau_Fit_charge(HighestPtTauIndex)*Ntp->Muon_Charge(HighestPtMuonIndex);
  pass.at(charge)=(Ntp->KFTau_Fit_charge(HighestPtTauIndex)*Ntp->Muon_Charge(HighestPtMuonIndex) == -1);
    

  value.at(TauDecay)=Ntp->PFTau_hpsDecayMode(Ntp->KFTau_MatchedHPS_idx(HighestPtTauIndex));
  pass.at(TauDecay)=true;//(Ntp->PFTau_hpsDecayMode(Ntp->KFTau_MatchedHPS_idx(HighestPtTauIndex)) == cut.at(TauDecay));
  
  value.at(NTau)=nTau;
  pass.at(NTau)=(nTau == cut.at(NTau));
    
  value.at(NMuon)=nMuon;
  pass.at(NMuon)= (nMuon == cut.at(NMuon));
  




  }

  double wobs=1;
  double w;

  double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
  double UnSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::UnSpin);
  double FlipSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin);


  if(!Ntp->isData()){w = Ntp->EvtWeight3D();}
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs); 
 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){

    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);;
    if(Ntp->NKFTau()!=0 && Ntp->NMuons()!=0){

    TLorentzVector tau_3pi = Ntp->KFTau_TauFit_p4(HighestPtTauIndex);
    TLorentzVector muon    = Ntp->Muons_p4(HighestPtMuonIndex);

    TLorentzVector tau1,tau2,taumu1,taumu2,Z;
    double e1,e2,pz1,pz2;

    TLorentzVector z1,z2;

    double zmass = 91.1876,mtau = 1.777;

    double A =  0.5*(zmass*zmass - 2*mtau*mtau) -tau_3pi.Pt()* tau_3pi.Pt() ;
    double B =  mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt() + A*A/tau_3pi.Pz()/tau_3pi.Pz();
    double root = sqrt(A*A*tau_3pi.E()*tau_3pi.E() - B*tau_3pi.Pz()*tau_3pi.Pz()*(mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt()));



    e1 = (A*tau_3pi.E() + root)/(mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt());
    e2 = (A*tau_3pi.E() - root)/(mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt());

    pz1 = sqrt(e1*e1 -  (mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt()));
    pz2 = sqrt(e2*e2 -  (mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt()));
    
    tau1.SetXYZM(-tau_3pi.Px(),-tau_3pi.Py(),pz1,mtau); 
    tau2.SetXYZM(-tau_3pi.Px(),-tau_3pi.Py(),pz2,mtau); 

    z1 = tau1+tau_3pi;
    z2 = tau2+tau_3pi;
    Tau3PiPt.at(t).Fill(tau_3pi.Pt(),w);
    Tau1Pt.at(t).Fill(tau1.Pt(),w);
    Tau2Pt.at(t).Fill(tau2.Pt(),w);
    Tau1theta.at(t).Fill(tau1.Theta(),w);	    
    Tau2theta.at(t).Fill(tau2.Theta(),w);      
    Tau1Deltatheta.at(t).Fill(tau1.Theta() - muon.Theta(),w); 
    Tau2Deltatheta.at(t).Fill(tau2.Theta() - muon.Theta(),w); 
    
    if(fabs(tau1.Theta() - muon.Theta()) < fabs(tau2.Theta() - muon.Theta())){taumu1 = tau1; taumu2 = tau2;}
    else{taumu1 = tau2; taumu2 = tau1;}
    DeltaTheta_1st.at(t).Fill(taumu1.Theta() - muon.Theta(),w); 
    DeltaTheta_2nd.at(t).Fill(taumu2.Theta() - muon.Theta(),w); 
    //    if(taumu1.E()!=0 && taumu1.Et()!=0){


    IsoVsMuonEnergy.at(t).Fill((Ntp->Muon_emEt05(HighestPtMuonIndex) + Ntp->Muon_hadEt05(HighestPtMuonIndex) + Ntp->Muon_sumPt05(HighestPtMuonIndex))/Ntp->Muons_p4(HighestPtMuonIndex).Pt(),muon.E(),w); 
    IsoVsMuonTauEnergy.at(t).Fill((Ntp->Muon_emEt05(HighestPtMuonIndex) + Ntp->Muon_hadEt05(HighestPtMuonIndex) + Ntp->Muon_sumPt05(HighestPtMuonIndex))/Ntp->Muons_p4(HighestPtMuonIndex).Pt(),taumu1.E(),w); 
    IsoAbsVsMuonEnergy.at(t).Fill(Ntp->Muon_emEt05(HighestPtMuonIndex) + Ntp->Muon_hadEt05(HighestPtMuonIndex) + Ntp->Muon_sumPt05(HighestPtMuonIndex),muon.E(),w); 

    xRec.at(t).Fill(muon.E()/taumu1.E(),w); 

    xmuon_plus.at(t).Fill(muon.E()/taumu1.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
    xmuonT_plus.at(t).Fill(muon.Et()/taumu1.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 


    xmuon_minus.at(t).Fill(muon.E()/taumu1.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
    xmuonT_minus.at(t).Fill(muon.Et()/taumu1.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 


    
    //------------------------------- comment out for work with data

    for(int iz =0; iz<Ntp->NMCSignalParticles(); iz++){
      TLorentzVector TruthTauMu;
      TLorentzVector TruthTauPi;
      
      TLorentzVector TruthMu;
      TLorentzVector RecoiledTauMu;
      TLorentzVector RecoiledZ;
      int TauMuIndex =0;
      
      double e1sim,e2sim,pz1sim,pz2sim;
      TLorentzVector z1sim,z2sim;
      TLorentzVector tau1sim,tau2sim,taumu1sim,taumu2sim,Zsim;
      double Erec,Emc;
      bool signal =false;
      double x[16] = {0,0.05,0.10,0.15,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.0,1.1,1.2,1.5};
      double E[16];

      if(Ntp->MCTau_JAK(0) == 2 and (Ntp->MCTau_JAK(1) ==5 or Ntp->MCTau_JAK(1) ==8 or Ntp->MCTau_JAK(1) ==14 or Ntp->MCTau_JAK(1) ==18) ){TruthTauMu = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(0)); TruthTauPi = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(1)); signal =true;TauMuIndex = 0;}
      else if((Ntp->MCTau_JAK(0) ==5 or Ntp->MCTau_JAK(0) ==8 or Ntp->MCTau_JAK(0) ==14 or Ntp->MCTau_JAK(0) ==18) and Ntp->MCTau_JAK(1) ==2){TruthTauMu = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(1)); TruthTauPi = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(0)); signal =true;TauMuIndex = 1;}
      
      
      if(signal){
	for(int iProd1 =0; iProd1 < Ntp->NMCTauDecayProducts(TauMuIndex); iProd1++ ){
	 	  if(abs( Ntp->MCTauandProd_pdgid(TauMuIndex,iProd1))==13){
		    TruthMu = Ntp->MCTauandProd_p4(TauMuIndex,iProd1);
		  }
	}
	    
      
 
	     double A1 =  0.5*(zmass*zmass - 2*mtau*mtau) -TruthTauPi.Pt()*TruthTauPi.Pt() ;
	     double B1 =  mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt() + A1*A1/TruthTauPi.Pz()/TruthTauPi.Pz();
	     double root1 = sqrt(A1*A1*TruthTauPi.E()*TruthTauPi.E() - B1*TruthTauPi.Pz()*TruthTauPi.Pz()*(mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt()));
	     
	     e1sim = (A1*TruthTauPi.E() + root1)/(mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt());
	     e2sim = (A1*TruthTauPi.E() - root1)/(mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt());
	    
	    pz1sim = sqrt(e1sim*e1sim -  (mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt()));
	    pz2sim = sqrt(e2sim*e2sim -  (mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt()));
	 
	 
	    tau1sim.SetXYZM(-TruthTauPi.Px(),-TruthTauPi.Py(),pz1sim,mtau); 
	    tau2sim.SetXYZM(-TruthTauPi.Px(),-TruthTauPi.Py(),pz2sim,mtau); 

	    z1sim = tau1sim+TruthTauPi;
	    z2sim = tau2sim+TruthTauPi;
	    if(fabs(tau1sim.Theta() - TruthMu.Theta()) < fabs(tau2sim.Theta() - TruthMu.Theta())){taumu1sim = tau1sim; taumu2sim = tau2sim;}
	    else{taumu1sim = tau2sim; taumu2sim = tau1sim;}


	    xmuonTruthT_plus.at(t).Fill(TruthMu.Et()/TruthTauMu.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	    xmuonTruth_plus.at(t).Fill(TruthMu.E()/taumu1sim.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
  
	    xmuonTruthT_minus.at(t).Fill(TruthMu.Et()/TruthTauMu.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	    xmuonTruth_minus.at(t).Fill(TruthMu.E()/taumu1sim.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	    xTruth.at(t).Fill(TruthMu.E()/taumu1sim.E(),w); 

	    xMC.at(t).Fill(TruthMu.E()/TruthTauMu.E(),w); 
	    xmuonMC_plus.at(t).Fill(TruthMu.E()/TruthTauMu.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	    xmuonTMC_plus.at(t).Fill(TruthMu.Et()/TruthTauMu.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
  
	    xmuonMC_minus.at(t).Fill(TruthMu.E()/TruthTauMu.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	    xmuonTMC_minus.at(t).Fill(TruthMu.Et()/TruthTauMu.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 

	    TruthResolution.at(t).Fill(taumu1sim.E() - TruthTauMu.E(),w); 
	    TruthDeltaTheta.at(t).Fill(taumu1sim.Theta() - TruthMu.Theta(),w); 


      }
      
      if(signal){
	res0.at(t).Fill(x[0]*(taumu1.E() -TruthTauMu.E() ),w); 
	res5.at(t).Fill(x[1]*(taumu1.E() -TruthTauMu.E() ),w); 
	res10.at(t).Fill(x[2]*(taumu1.E() -TruthTauMu.E() ),w); 
	res15.at(t).Fill(x[3]*(taumu1.E() -TruthTauMu.E() ),w); 
	res20.at(t).Fill(x[4]*(taumu1.E() -TruthTauMu.E() ),w); 
	res30.at(t).Fill(x[5]*(taumu1.E() -TruthTauMu.E() ),w); 
	res40.at(t).Fill(x[6]*(taumu1.E() -TruthTauMu.E() ),w); 
	res50.at(t).Fill(x[7]*(taumu1.E() -TruthTauMu.E() ),w); 
	res60.at(t).Fill(x[8]*(taumu1.E() -TruthTauMu.E() ),w); 
	res70.at(t).Fill(x[9]*(taumu1.E() -TruthTauMu.E() ),w); 
	res80.at(t).Fill(x[10]*(taumu1.E() -TruthTauMu.E() ),w); 
	res90.at(t).Fill(x[11]*(taumu1.E() -TruthTauMu.E() ),w); 
	res100.at(t).Fill(x[12]*(taumu1.E() -TruthTauMu.E() ),w); 

	res110.at(t).Fill(x[13]*(taumu1.E() -TruthTauMu.E() ),w); 
	res120.at(t).Fill(x[14]*(taumu1.E() -TruthTauMu.E() ),w); 
	res150.at(t).Fill(x[15]*(taumu1.E() -TruthTauMu.E() ),w); 

	E[0]  = TruthTauMu.E() + x[0]*(taumu1.E() - TruthTauMu.E());
	E[1]  = TruthTauMu.E() + x[1]*(taumu1.E() - TruthTauMu.E());
	E[2]  = TruthTauMu.E() + x[2]*(taumu1.E() - TruthTauMu.E());
	E[3]  = TruthTauMu.E() + x[3]*(taumu1.E() - TruthTauMu.E());
	E[4]  = TruthTauMu.E() + x[4]*(taumu1.E() - TruthTauMu.E());
	E[5]  = TruthTauMu.E() + x[5]*(taumu1.E() - TruthTauMu.E());
	E[6]  = TruthTauMu.E() + x[6]*(taumu1.E() - TruthTauMu.E());
	E[7]  = TruthTauMu.E() + x[7]*(taumu1.E() - TruthTauMu.E());
	E[8]  = TruthTauMu.E() + x[8]*(taumu1.E() - TruthTauMu.E());
	E[9]  = TruthTauMu.E() + x[9]*(taumu1.E() - TruthTauMu.E());
	E[10]  = TruthTauMu.E() + x[10]*(taumu1.E() - TruthTauMu.E());
	E[11]  = TruthTauMu.E() + x[11]*(taumu1.E() - TruthTauMu.E());
	E[12]  = TruthTauMu.E() + x[12]*(taumu1.E() - TruthTauMu.E());
	E[13]  = TruthTauMu.E() + x[13]*(taumu1.E() - TruthTauMu.E());
	E[14]  = TruthTauMu.E() + x[14]*(taumu1.E() - TruthTauMu.E());
	E[15]  = TruthTauMu.E() + x[15]*(taumu1.E() - TruthTauMu.E());


	xp0.at(t).Fill(muon.E()/E[0],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp5.at(t).Fill(muon.E()/E[1],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp10.at(t).Fill(muon.E()/E[2],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp15.at(t).Fill(muon.E()/E[3],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp20.at(t).Fill(muon.E()/E[4],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp30.at(t).Fill(muon.E()/E[5],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp40.at(t).Fill(muon.E()/E[6],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp50.at(t).Fill(muon.E()/E[7],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp60.at(t).Fill(muon.E()/E[8],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp70.at(t).Fill(muon.E()/E[9],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp80.at(t).Fill(muon.E()/E[10],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp90.at(t).Fill(muon.E()/E[11],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp100.at(t).Fill(muon.E()/E[12],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp110.at(t).Fill(muon.E()/E[13],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp120.at(t).Fill(muon.E()/E[14],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	xp150.at(t).Fill(muon.E()/E[15],w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 



	xm0.at(t).Fill(muon.E()/E[0],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm5.at(t).Fill(muon.E()/E[1],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm10.at(t).Fill(muon.E()/E[2],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm15.at(t).Fill(muon.E()/E[3],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm20.at(t).Fill(muon.E()/E[4],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm30.at(t).Fill(muon.E()/E[5],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm40.at(t).Fill(muon.E()/E[6],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm50.at(t).Fill(muon.E()/E[7],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm60.at(t).Fill(muon.E()/E[8],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm70.at(t).Fill(muon.E()/E[9],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm80.at(t).Fill(muon.E()/E[10],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm90.at(t).Fill(muon.E()/E[11],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm100.at(t).Fill(muon.E()/E[12],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm110.at(t).Fill(muon.E()/E[13],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm120.at(t).Fill(muon.E()/E[14],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	xm150.at(t).Fill(muon.E()/E[15],w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 

	xReco0.at(t).Fill(muon.E()/E[0],w); 
	xReco5.at(t).Fill(muon.E()/E[1],w); 
	xReco10.at(t).Fill(muon.E()/E[2],w); 
	xReco15.at(t).Fill(muon.E()/E[3],w); 
	xReco20.at(t).Fill(muon.E()/E[4],w); 
	xReco30.at(t).Fill(muon.E()/E[5],w); 
	xReco40.at(t).Fill(muon.E()/E[6],w); 
	xReco50.at(t).Fill(muon.E()/E[7],w); 
	xReco60.at(t).Fill(muon.E()/E[8],w); 
	xReco70.at(t).Fill(muon.E()/E[9],w); 
	xReco80.at(t).Fill(muon.E()/E[10],w); 
	xReco90.at(t).Fill(muon.E()/E[11],w); 
	xReco100.at(t).Fill(muon.E()/E[12],w); 
	xReco110.at(t).Fill(muon.E()/E[13],w); 
	xReco120.at(t).Fill(muon.E()/E[14],w); 
	xReco150.at(t).Fill(muon.E()/E[15],w); 



	Energy_resolution1.at(t).Fill(taumu1.E() -TruthTauMu.E(),w); 
	Energy_resolution2.at(t).Fill(taumu2.E() -TruthTauMu.E(),w); 
	TransverseEnergy_resolution.at(t).Fill(taumu1.Et() -TruthTauMu.Et(),w); 
      }
    
      
    }

    }
    
    //      //------------------------------- comment out for work with data
    
  }
  
}
 




 

