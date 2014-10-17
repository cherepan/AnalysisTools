#include "AmbuguityAndResolution.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"

AmbuguityAndResolution::AmbuguityAndResolution(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

AmbuguityAndResolution::~AmbuguityAndResolution(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "AmbuguityAndResolution::~AmbuguityAndResolution Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "AmbuguityAndResolution::~AmbuguityAndResolution()" << std::endl;
}

void  AmbuguityAndResolution::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==TauIsQuality)       cut.at(TauIsQuality)=1;
    if(i==TauPtCut)           cut.at(TauPtCut)=20;
    if(i==TauIsIsoMVA2Medium)              cut.at(TauIsIsoMVA2Medium)=1;
    if(i==TauIsIsoMVA2Tight)               cut.at(TauIsIsoMVA2Tight)=1;
    if(i==TauIsIsoMVAMedium)               cut.at(TauIsIsoMVAMedium)=1;
    if(i==TauIsIsoMVATight)                cut.at(TauIsIsoMVATight)=1;
    if(i==TauIsIsoDBSumPtCorrMedium)       cut.at(TauIsIsoDBSumPtCorrMedium)=1;
    if(i==TauIsIsoDBSumPtCorrTight)        cut.at(TauIsIsoDBSumPtCorrTight)=1;
    if(i==TauIsIsoCombThreeHitsMedium)     cut.at(TauIsIsoCombThreeHitsMedium)=1;
    if(i==TauIsIsoCombThreeHitsTight)      cut.at(TauIsIsoCombThreeHitsTight)=1;
    if(i==TauIsNotElectron)   cut.at(TauIsNotElectron)=1;
    if(i==TauIsNotMuon)       cut.at(TauIsNotMuon)=1;
    if(i==MuonPt)             cut.at(MuonPt)=18; //18
    if(i==MuonisGlob)         cut.at(MuonisGlob)=1;
    if(i==MuonisPF)           cut.at(MuonisPF)=1;
    if(i==MuonIso)            cut.at(MuonIso)=0.12; // 0.05  for 0.3 cone
    if(i==MuonNormChi2)       cut.at(MuonNormChi2)=10;   // <10
    if(i==MuonHitsInMS)       cut.at(MuonHitsInMS)=0;    // >0
    if(i==MuonPixelHits)      cut.at(MuonPixelHits)=0;   // >0
    if(i==MuonTrackLayer)     cut.at(MuonTrackLayer)=5;  // >5
    if(i==MuonNumOfStations)  cut.at(MuonNumOfStations)=1;  // >1
    if(i==MET)                cut.at(MET)=30;             // < 40
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
    else if(i==TauIsQuality){
      title.at(i)="TauIsQuality ";
      hlabel="TauIsQuality ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsQuality_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsQuality_",htitle,2,-0.5,1.5,hlabel,"Events"));
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
   else if(i==TauIsIsoMVA2Medium){
      title.at(i)="TauIsIsoMVA2Medium == ";
      title.at(i)+=cut.at(TauIsIsoMVA2Medium);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsIsoMVA2Medium";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsoMVA2Medium_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsoMVA2Medium_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
   else if(i==TauIsIsoMVA2Tight){
      title.at(i)=" TauIsIsoMVA2Tight== ";
      title.at(i)+=cut.at(TauIsIsoMVA2Tight);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsIsoMVA2Tight";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsoMVA2Tight_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsoMVA2Tight_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 

   else if(i==TauIsIsoMVAMedium){
      title.at(i)="TauIsIsoMVAMedium == ";
      title.at(i)+=cut.at(TauIsIsoMVAMedium);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsIsoMVAMedium";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsoMVAMedium_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsoMVAMedium_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
   else if(i==TauIsIsoMVATight){
      title.at(i)=" TauIsIsoMVATight== ";
      title.at(i)+=cut.at(TauIsIsoMVATight);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsIsoMVATight";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsoMVATight_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsoMVATight_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 

   else if(i==TauIsIsoDBSumPtCorrMedium){
      title.at(i)="TauIsIsoDBSumPtCorrMedium == ";
      title.at(i)+=cut.at(TauIsIsoDBSumPtCorrMedium);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsIsoDBSumPtCorrMedium";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsoDBSumPtCorrMedium_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsoDBSumPtCorrMedium_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
   else if(i==TauIsIsoDBSumPtCorrTight){
      title.at(i)="TauIsIsoDBSumPtCorrTight == ";
      title.at(i)+=cut.at(TauIsIsoDBSumPtCorrTight);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsIsoDBSumPtCorrTight";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsoDBSumPtCorrTight_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsoDBSumPtCorrTight_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
   else if(i==TauIsIsoCombThreeHitsMedium){
      title.at(i)="TauIsIsoCombThreeHitsMedium == ";
      title.at(i)+=cut.at(TauIsIsoCombThreeHitsMedium);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsIsoCombThreeHitsMedium";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsoCombThreeHitsMedium_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsoCombThreeHitsMedium_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
   else if(i==TauIsIsoCombThreeHitsTight){
      title.at(i)="TauIsIsoCombThreeHitsTight == ";
      title.at(i)+=cut.at(TauIsIsoCombThreeHitsTight);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsIsoCombThreeHitsTight";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsoCombThreeHitsTight_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsoCombThreeHitsTight_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 


   else if(i==TauIsNotElectron){
      title.at(i)="$TauIsNotElectron ==1 $";
      title.at(i)+=cut.at(TauIsNotElectron);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsNotElectron";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsNotElectron_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsNotElectron_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
   else if(i==TauIsNotMuon){
      title.at(i)="$TauIsNotMuon ==1 $";
      title.at(i)+=cut.at(TauIsNotMuon);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsNotMuon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsNotMuon_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsNotMuon_",htitle,2,-0.5,1.5,hlabel,"Events"));
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

   else if(i==MuonisGlob){
      title.at(i)="Muon is global ";
      hlabel="MuonisGlob";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonisGlob_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonisGlob_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
   else if(i==MuonisPF){
      title.at(i)="MuonisPF ";
      hlabel="MuonisPF";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonisPF_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonisPF_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
   else if(i==MuonIso){
      title.at(i)="$MuonIso < $";
      title.at(i)+=cut.at(MuonIso);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#MuonIso";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIso_",htitle,60,0,0.3,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIso_",htitle,60,0,0.3,hlabel,"Events"));
    }

   else if(i==MuonNormChi2){
      title.at(i)="$MuonNormChi2 < $";
      title.at(i)+=cut.at(MuonNormChi2);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="MuonNormChi2";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonNormChi2_",htitle,20,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonNormChi2_",htitle,20,0,20,hlabel,"Events"));
    }  
   else if(i==MuonHitsInMS){
      title.at(i)="$MuonHitsInMS < $";
      title.at(i)+=cut.at(MuonHitsInMS);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="MuonHitsInMS";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonHitsInMS_",htitle,10,0,10,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonHitsInMS_",htitle,10,0,10,hlabel,"Events"));
    }  

   else if(i==MuonPixelHits){
      title.at(i)="$MuonPixelHits < $";
      title.at(i)+=cut.at(MuonPixelHits);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="MuonPixelHits";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonPixelHits_",htitle,10,0,10,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonPixelHits_",htitle,10,0,10,hlabel,"Events"));
    }  

   else if(i==MuonTrackLayer){
      title.at(i)="$MuonTrackLayer < $";
      title.at(i)+=cut.at(MuonTrackLayer);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="MuonTrackLayer";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonTrackLayer_",htitle,10,0,10,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonTrackLayer_",htitle,10,0,10,hlabel,"Events"));
    }  
   else if(i==MuonNumOfStations){
      title.at(i)="$MuonNumOfStations < $";
      title.at(i)+=cut.at(MuonNumOfStations);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="MuonNumOfStations";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonNumOfStations_",htitle,10,0,10,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonNumOfStations_",htitle,10,0,10,hlabel,"Events"));
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
//   NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Vertex","Events");
//   TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",20,0,80,"TauPt","Events");
//   xmuon=HConfig.GetTH1D(Name+"_xmuon","xmuon",15,0.2,1.4,"xmuon","Events");
//   xmuonT=HConfig.GetTH1D(Name+"_xmuonT","xmuonT",15,0,2,"xmuonT","Events");


ZMassPlus=HConfig.GetTH1D(Name+"_ZMassPlus","ZMassPlus",30,40,160,"M_{#tau-#mu}, GeV","");
ZMassMins=HConfig.GetTH1D(Name+"_ZMassMins","ZMassMins",30,40,160,"M_{#tau-#mu}, GeV","");
ZMassZero=HConfig.GetTH1D(Name+"_ZMassZero","ZMassZero",30,40,160,"M_{#tau-#mu}, GeV","");

TauPtPlus=HConfig.GetTH1D(Name+"_TauPtPlus","TauPtPlus",30,15,50,"TauPtPlus, GeV","");
TauPtMins=HConfig.GetTH1D(Name+"_TauPtMins","TauPtMins",30,15,50,"TauPtMins, GeV","");
TauPtZero=HConfig.GetTH1D(Name+"_TauPtZero","TauPtZero",30,15,50,"TauPtZero, GeV","");

PtBalancePlus=HConfig.GetTH1D(Name+"_PtBalancePlus","PtBalancePlus",30,-40,40,"PtBalancePlus, GeV","");
PtBalanceMins=HConfig.GetTH1D(Name+"_PtBalanceMins","PtBalanceMins",30,-40,40,"PtBalanceMins, GeV","");
PtBalanceZero=HConfig.GetTH1D(Name+"_PtBalanceZero","PtBalanceZero",30,-40,40,"PtBalanceZero, GeV","");


ambuguity=HConfig.GetTH1D(Name+"_ambuguity","ambuguity",100,0,100," GeV","");

TauPtPlusResolution=HConfig.GetTH1D(Name+"_TauPtPlusResolution","TauPtPlusResolution",30,-40,40,"TauPtPlusResolution, GeV","");
TauPtMinsResolution=HConfig.GetTH1D(Name+"_TauPtMinsResolution","TauPtMinsResolution",30,-40,40,"TauPtMinsResolution, GeV","");
TauPtZeroResolution=HConfig.GetTH1D(Name+"_TauPtZeroResolution","TauPtZeroResolution",30,-40,40,"TauPtZeroResolution, GeV","");



PtBalancePlusMins=HConfig.GetTH2D(Name+"_PtBalancePlusMins","PtBalancePlusMins",30,-40,40,30,-40,40,"PtBalancePlusMins","");
// PtBalancePlusAmbiguity=HConfig.GetTH2D(Name+"_PtBalancePlusAmbiguity","PtBalancePlusAmbiguity",30,-40,40,100,0,100,"PtBalancePlusAmbiguity","");
// PtBalanceMinsAmbiguity=HConfig.GetTH2D(Name+"_PtBalanceMinsAmbiguity","PtBalanceMinsAmbiguity",30,-40,40,100,0,100,"PtBalanceMinsAmbiguity","");



  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
} 

 
// TString name,TString title,int nbinsx, double minx, double maxx, 
// 				       int nbinsy, double miny, double maxy, TString xaxis, TString yaxis)

void  AmbuguityAndResolution::Store_ExtraDist(){


//  Extradist1d.push_back(&NGoodVtx);
//  Extradist1d.push_back(&TauPt);
//  Extradist1d.push_back(&xmuon);
//  Extradist1d.push_back(&xmuonT);
  Extradist1d.push_back(&ZMassPlus); 		   
  Extradist1d.push_back(&ZMassMins); 		   
  Extradist1d.push_back(&ZMassZero); 		   
			   
  Extradist1d.push_back(&TauPtPlus); 		   
  Extradist1d.push_back(&TauPtMins); 		   
  Extradist1d.push_back(&TauPtZero); 		   
			   
  Extradist1d.push_back(&PtBalancePlus); 	   
  Extradist1d.push_back(&PtBalanceMins); 	   
  Extradist1d.push_back(&PtBalanceZero); 	   
			   
			   
  Extradist1d.push_back(&ambuguity); 		   
			   
  Extradist1d.push_back(&TauPtPlusResolution);   
  Extradist1d.push_back(&TauPtMinsResolution);   
  Extradist1d.push_back(&TauPtZeroResolution);   
			   
			   
		   
  Extradist2d.push_back(&PtBalancePlusMins);	   
//   Extradist2d.push_back(&PtBalancePlusAmbiguity);
//   Extradist2d.push_back(&PtBalanceMinsAmbiguity);





 
//  Extradist1d.push_back(&TruthRatio);
//  Extradist1d.push_back(&TruthTransverseRatio);

}

void  AmbuguityAndResolution::doEvent(){

  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  
  // Apply Selection
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
 


  unsigned int FirstLeadingTau=0;
  unsigned int SecondLeadingTau=0;
  unsigned int FirstLeadingMuon=0;
  unsigned int SecondLeadingMuon =0;
  unsigned int nTau=0,nMuon=0;

  float muonPt=0;
  float muonPt2=0;
  float tauPt=0;
  float tauPt2=0;

  int ambiguity = 0;



  for(unsigned int iTau=0;iTau < Ntp->NKFTau();iTau++){
    if(Ntp->KFTau_discriminatorByQC(iTau,ambiguity))nTau++;
    if(Ntp->KFTau_TauFit_p4(iTau,ambiguity).Pt() > tauPt){
      tauPt = Ntp->KFTau_TauFit_p4(iTau,ambiguity).Pt();
      SecondLeadingTau = FirstLeadingTau;
      FirstLeadingTau = iTau; 
    }
    if(Ntp->KFTau_TauFit_p4(iTau,ambiguity).Pt()  > tauPt2 && iTau!=FirstLeadingTau){
      tauPt2= Ntp->KFTau_TauFit_p4(iTau,ambiguity).Pt();
      SecondLeadingTau = iTau;
    }
  }
  
  unsigned int MatchedVertex;
  for(int iVert = 0; iVert < Ntp->NVtx(); iVert++){
    if(sqrt(pow(Ntp->Vtx(iVert).X() - Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).X(),2) 
	    + pow(Ntp->Vtx(iVert).Y() - Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).Y(),2)  
	    + pow(Ntp->Vtx(iVert).Z() - Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).Z(),2)) < 0.001 ){MatchedVertex = iVert;}
    
  }


//   std::cout<<"Tau     Vertex: "<<Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).X()<<"  "<<Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).Y()<<"  "<<Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).Z()<<std::endl;
//   std::cout<<"Matched Vertex: "<<Ntp->Vtx(MatchedVertex).X()<<"  "<<Ntp->Vtx(MatchedVertex).Y()<<"  "<<Ntp->Vtx(MatchedVertex).Z()<<std::endl;



  std::vector<int> MuonsFromTauVertex; 
  if(Ntp->NMuons()!=0){
    for(int iTrack = 0; iTrack <Ntp->Vtx_Track_idx(MatchedVertex).size(); iTrack++){
      for(unsigned int iMuon=0; iMuon<Ntp->NMuons(); iMuon++){
	if(Ntp->Muon_Track_idx(iMuon) == Ntp->Vtx_Track_idx(MatchedVertex).at(iTrack))MuonsFromTauVertex.push_back(iMuon);
	  }
    }
  }
  
  if(MuonsFromTauVertex.size()!=0){
    for(unsigned int iMuon=0; iMuon!=MuonsFromTauVertex.size(); iMuon++){
      if(Ntp->Muon_isGlobalMuon(iMuon))nMuon++;
      if(Ntp->Muons_p4(iMuon).Pt() > muonPt){
	muonPt = Ntp->Muons_p4(iMuon).Pt();
	SecondLeadingMuon = FirstLeadingMuon;
	FirstLeadingMuon = iMuon;
      }
      if(Ntp->Muons_p4(iMuon).Pt() > muonPt2 && iMuon != FirstLeadingMuon){
	muonPt2 =  Ntp->Muons_p4(iMuon).Pt();
	SecondLeadingMuon = iMuon;
      }
	
    }
    //  std::cout<<"NMuons: "<<MuonsFromTauVertex.size()<<"  MuonPt   "<< Ntp->Muons_p4(FirstLeadingMuon).Pt() <<std::endl;
    //  if(MuonsFromTauVertex.size() > 1)   std::cout<<"SecondLeadingMuonPt: "<< Ntp->Muons_p4(SecondLeadingMuon).Pt() <<std::endl;
  }
//   if(Ntp->NKFTau()!=0)  std::cout<<"1 TauPt : "<<Ntp->KFTau_TauFit_p4(FirstLeadingTau,ambiguity).Pt()<<std::endl;
//   if(Ntp->NKFTau()>1)  std::cout<<"2 TauPt : "<<Ntp->KFTau_TauFit_p4(SecondLeadingTau,ambiguity).Pt()<<std::endl;




  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  value.at(TriggerOk)=1;
  pass.at(TriggerOk)=true; 

  if(Ntp->NKFTau()!=0 && MuonsFromTauVertex.size()!=0 ){



  value.at(TauIsQuality) =  Ntp->KFTau_discriminatorByQC(FirstLeadingTau,ambiguity);
  //  value.at(TauIsQuality) =  Ntp->KFTau_discriminatorByKFit(FirstLeadingTau);//Ntp->KFTau_discriminatorByQC(FirstLeadingTau);
  pass.at(TauIsQuality)  =  (value.at(TauIsQuality)==cut.at(TauIsQuality));

  value.at(TauPtCut) = Ntp->KFTau_TauFit_p4(FirstLeadingTau,ambiguity).Pt();
  pass.at(TauPtCut)=(value.at(TauPtCut)>=cut.at(TauPtCut));



  value.at(TauIsIsoMVA2Medium) =   Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
  pass.at(TauIsIsoMVA2Medium)=true;//(value.at(TauIsIsoMVA2Medium) == cut.at(TauIsIsoMVA2Medium));



//   value.at(TauIsIsoMVA2Tight) =   Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA2(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
   pass.at(TauIsIsoMVA2Tight)=true;//(value.at(TauIsIsoMVA2Tight) == cut.at(TauIsIsoMVA2Tight));



  value.at(TauIsIsoMVAMedium) =   Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
  pass.at(TauIsIsoMVAMedium)=true;//(value.at(TauIsIsoMVAMedium) == cut.at(TauIsIsoMVAMedium));


//   value.at(TauIsIsoMVATight) =   Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
   pass.at(TauIsIsoMVATight)=true;//(value.at(TauIsIsoMVATight) == cut.at(TauIsIsoMVATight));



  value.at(TauIsIsoDBSumPtCorrMedium) =   Ntp->PFTau_isMediumIsolationDBSumPtCorr(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
  pass.at(TauIsIsoDBSumPtCorrMedium)=true;//(value.at(TauIsIsoDBSumPtCorrMedium) == cut.at(TauIsIsoDBSumPtCorrMedium));

  value.at(TauIsIsoDBSumPtCorrTight) =   Ntp->PFTau_isTightIsolationDBSumPtCorr(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
  pass.at(TauIsIsoDBSumPtCorrTight)=(value.at(TauIsIsoDBSumPtCorrTight) == cut.at(TauIsIsoDBSumPtCorrTight));



  value.at(TauIsIsoCombThreeHitsMedium) =   Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
  pass.at(TauIsIsoCombThreeHitsMedium)=true;//(value.at(TauIsIsoCombThreeHitsMedium) == cut.at(TauIsIsoCombThreeHitsMedium));

  value.at(TauIsIsoCombThreeHitsTight) =   Ntp->PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
  pass.at(TauIsIsoCombThreeHitsTight)=true;//(value.at(TauIsIsoCombThreeHitsTight) == cut.at(TauIsIsoCombThreeHitsTight));



  value.at(TauIsNotElectron)=Ntp->PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
  pass.at(TauIsNotElectron)=(value.at(TauIsNotElectron)==cut.at(TauIsNotElectron));
										     
  value.at(TauIsNotMuon)=Ntp->PFTau_isHPSAgainstMuonMedium2(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
  pass.at(TauIsNotMuon)=(value.at(TauIsNotMuon)==cut.at(TauIsNotMuon));


  value.at(MuonPt) = Ntp->Muons_p4(FirstLeadingMuon).Pt();
  pass.at(MuonPt)=(value.at(MuonPt) >=cut.at(MuonPt));
 
  value.at(MuonisGlob) = Ntp->Muon_isGlobalMuon(FirstLeadingMuon);
  pass.at(MuonisGlob)  = (value.at(MuonisGlob)==cut.at(MuonisGlob));

  value.at(MuonisPF) = Ntp->Muon_isPFMuon(FirstLeadingMuon);
  pass.at(MuonisPF)  = (value.at(MuonisPF)==cut.at(MuonisPF));

  value.at(MuonIso) = (  (Ntp->Muon_sumChargedHadronPt04(FirstLeadingMuon) + max(0.,Ntp->Muon_sumNeutralHadronEt04(FirstLeadingMuon) + Ntp->Muon_sumPhotonEt04(FirstLeadingMuon)  - 0.5*Ntp->Muon_sumPUPt04(FirstLeadingMuon)))/ Ntp->Muons_p4(FirstLeadingMuon).Pt());
  pass.at(MuonIso)=(value.at(MuonIso)<=cut.at(MuonIso));

  value.at(MuonNormChi2) = Ntp->Muon_normChi2(FirstLeadingMuon);
  pass.at(MuonNormChi2)  = (value.at(MuonNormChi2)<cut.at(MuonNormChi2));

  value.at(MuonHitsInMS) = Ntp->Muon_hitPattern_numberOfValidMuonHits(FirstLeadingMuon);
  pass.at(MuonHitsInMS)  = (value.at(MuonHitsInMS)>cut.at(MuonHitsInMS));

 
//   value.at(MuonPixelHits) = Ntp->Muon_numberofValidPixelHits(FirstLeadingMuon);
  pass.at(MuonPixelHits)  = true;//(value.at(MuonPixelHits)>cut.at(MuonPixelHits));

//   value.at(MuonTrackLayer) = Ntp->Muon_trackerLayersWithMeasurement(FirstLeadingMuon);
   pass.at(MuonTrackLayer)  = true;//(value.at(MuonTrackLayer)>cut.at(MuonTrackLayer));

   value.at(MuonNumOfStations) = Ntp->Muon_numberOfMatchedStations(FirstLeadingMuon);
   pass.at(MuonNumOfStations)  = (value.at(MuonNumOfStations)>cut.at(MuonNumOfStations));

  value.at(MET) =   sqrt(2*Ntp->Muons_p4(FirstLeadingMuon).Pt()*Ntp->MET_et()*(1  - cos(Ntp->Muons_p4(FirstLeadingMuon).Phi() - Ntp->MET_phi())) );
  pass.at(MET)=(value.at(MET)<=cut.at(MET));





  value.at(charge) =Ntp->KFTau_Fit_charge(FirstLeadingTau)*Ntp->Muon_Charge(FirstLeadingMuon);
  pass.at(charge)=(Ntp->KFTau_Fit_charge(FirstLeadingTau)*Ntp->Muon_Charge(FirstLeadingMuon) == -1);


  value.at(TauDecay)=Ntp->PFTau_hpsDecayMode(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau));
  pass.at(TauDecay)=true;//(Ntp->PFTau_hpsDecayMode(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau)) == cut.at(TauDecay));
  
  value.at(NTau)=nTau;
  pass.at(NTau)=true;//(nTau == cut.at(NTau));
    
  value.at(NMuon)=nMuon;
  pass.at(NMuon)= (nMuon == cut.at(NMuon));
  


  }

  double wobs=1;
  double w;



  if(!pass.at(charge) and  !pass.at(MuonIso)){
    if(Ntp->isData()){
         if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id for QCD: "<< DataMCType::QCD <<std::endl; return;}
      pass.at(charge)=true;
      pass.at(MuonIso) = true;
     
    }
  }
  



  if(!Ntp->isData()){w = Ntp->EvtWeight3D();}
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs); 
 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
  


    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
         if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    //    NGoodVtx.at(t).Fill(nGoodVtx,w);;
    if(Ntp->NKFTau()!=0 && Ntp->NMuons()!=0){


   
      TLorentzVector tau_3piZero = Ntp->KFTau_TauFit_p4(FirstLeadingTau,1);
      TLorentzVector tau_3piPlus = Ntp->KFTau_TauFit_p4(FirstLeadingTau,2);
      TLorentzVector tau_3piMins = Ntp->KFTau_TauFit_p4(FirstLeadingTau,0);

      TLorentzVector tau_a1Initial = Ntp->KFTau_a1Initial_p4(FirstLeadingTau);

      TLorentzVector muon    = Ntp->Muons_p4(FirstLeadingMuon);
      std::cout<<"ambig Mins "<< tau_3piMins.Pt()<<std::endl;
      std::cout<<"ambig Zero "<< tau_3piZero.Pt()<<std::endl;
      std::cout<<"ambig Plus "<< tau_3piPlus.Pt()<<std::endl;
      TLorentzVector tau,ZPlus, Zmins, ZZero;
      double zmass = 91.1876,mtau = 1.777;
      double ma1 = 1.23;
      
      double energy;
      double theta = muon.Theta();

      ZPlus  = tau_3piPlus + muon;
      Zmins  = tau_3piMins + muon;
      ZZero  = tau_3piZero + muon;
      


      //  deltaPhi.at(t).Fill(abs(1  - abs(muon.Phi() - tau_3pi.Phi())/3.1415));

      ZMassPlus.at(t).Fill(ZPlus.M(),w); 
      ZMassMins.at(t).Fill(Zmins.M(),w); 
      ZMassZero.at(t).Fill(ZZero.M(),w); 

      TauPtPlus.at(t).Fill(tau_3piPlus.Pt(),w); 
      TauPtMins.at(t).Fill(tau_3piMins.Pt(),w); 
      TauPtZero.at(t).Fill(tau_3piZero.Pt(),w); 


      PtBalancePlus.at(t).Fill(tau_3piPlus.Pt()-muon.Pt(),w); 
      PtBalanceMins.at(t).Fill(tau_3piMins.Pt()-muon.Pt(),w); 
      PtBalanceZero.at(t).Fill(tau_3piZero.Pt()-muon.Pt(),w); 


      double     ambig= (mtau*mtau - ma1*ma1)*tau_a1Initial.E()/ma1/ma1;
      ambuguity.at(t).Fill(ambig,w); 


      PtBalancePlusMins.at(t).Fill(tau_3piPlus.Pt()-muon.Pt(),tau_3piMins.Pt()-muon.Pt(),w); 
//       PtBalancePlusAmbiguity.at(t).Fill(tau_3piPlus.Pt()-muon.Pt(),ambig,w); 
//       PtBalanceMinsAmbiguity.at(t).Fill(tau_3piMins.Pt()-muon.Pt(),ambig,w); 


    //------------------------------- comment out for work with data
    
          for(int iz =0; iz<Ntp->NMCSignalParticles(); iz++){
 
   // int   iz = 0;
      TLorentzVector TruthTauMu;
      TLorentzVector TruthTauPi;
      
      TLorentzVector TruthMu;
      TLorentzVector RecoiledTauMu;
      TLorentzVector RecoiledZ;
      int TauMuIndex =0;
      int TauA1Index =0;
      
      double e1sim,e2sim,pz1sim,pz2sim;
      TLorentzVector z1sim,z2sim;
      TLorentzVector tau1sim,tau2sim,taumu1sim,taumu2sim,Zsim;
      if(Ntp->MCSignalParticle_Tauidx(iz).size()!=0){
      bool signal =false;
      if(Ntp->MCTau_JAK(0) == 2 and (Ntp->MCTau_JAK(1) ==5 or Ntp->MCTau_JAK(1) ==8 or Ntp->MCTau_JAK(1) ==14 or Ntp->MCTau_JAK(1) ==18) ){TruthTauMu = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(0)); TruthTauPi = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(1)); signal =true;TauMuIndex = 0;TauA1Index=1;   }
      else if((Ntp->MCTau_JAK(0) ==5 or Ntp->MCTau_JAK(0) ==8 or Ntp->MCTau_JAK(0) ==14 or Ntp->MCTau_JAK(0) ==18) and Ntp->MCTau_JAK(1) ==2){TruthTauMu = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(1)); TruthTauPi = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(0)); signal =true;TauMuIndex = 1;TauA1Index=0;     }
              
      
      if(signal){
      
	for(int iProd1 =0; iProd1 < Ntp->NMCTauDecayProducts(TauMuIndex); iProd1++ ){
	
	 	  if(abs( Ntp->MCTauandProd_pdgid(TauMuIndex,iProd1))==13){
		    TruthMu = Ntp->MCTauandProd_p4(TauMuIndex,iProd1);
		  }

 

	}
	//------- matching ------
      


	TauPtPlusResolution.at(t).Fill(tau_3piPlus.Pt()-TruthTauPi.Pt(),w); 
	TauPtMinsResolution.at(t).Fill(tau_3piMins.Pt()-TruthTauPi.Pt(),w); 
	TauPtZeroResolution.at(t).Fill(tau_3piZero.Pt()-TruthTauPi.Pt(),w); 



      }


      }

 }














      //     TauPt.at(t).Fill(tau.Pt(),w); 
      //     xmuon.at(t).Fill(muon.E()/tau.E(),w); 
      //     xmuonT.at(t).Fill(muon.Et()/tau.Et(),w); 
//       std::cout<<"Tau vertex x,y,z  "<<Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).X()<<"  "<<Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).Y()<<"  "<<Ntp->KFTau_Fit_InitialPrimaryVertex(FirstLeadingTau).Z()<<std::endl;
//       std::cout<<"muon track index  "<<Ntp->Muon_Track_idx(FirstLeadingMuon)<<std::endl;
// 	bool matched = false;
//       for(int iVert = 0; iVert < Ntp->NVtx(); iVert++){
// // 	std::cout<<"Vertex position x,y,z  "<< Ntp->Vtx(iVert).X()<<"  "<<Ntp->Vtx(iVert).Y()<<"  "<<Ntp->Vtx(iVert).Z()<<std::endl;
// // 	std::cout<<"N tracks per vertex  "<<Ntp->Vtx_Track_idx(iVert).size()<<std::endl;

// 	for(int iTrack = 0; iTrack <Ntp->Vtx_Track_idx(iVert).size(); iTrack++){
// 	  //	std::cout<<"track index  "<<Ntp->Vtx_Track_idx(iVert).at(iTrack)<<std::endl;
	
// 	if(Ntp->Vtx_Track_idx(iVert).at(iTrack) == Ntp->Muon_Track_idx(FirstLeadingMuon)){std::cout<<"Vertex position x,y,z  "<< Ntp->Vtx(iVert).X()<<"  "<<Ntp->Vtx(iVert).Y()<<"  "<<Ntp->Vtx(iVert).Z()<<std::endl;matched = true;}
	
// 	}


//       }
//       if(matched)std::cout<<" mathced  "<<std::endl;
//       else std::cout<<"not  mathced  "<<std::endl;










    }

  }
}


   void  AmbuguityAndResolution::Finish(){
     unsigned int t;
     if(Nminus0.at(0).at(t).Integral()!=0) if(HConfig.GetHisto(false,DataMCType::QCD,t))ScaleAllHistOfType(t,1129./Nminus0.at(0).at(t).Integral());
     Selection::Finish();
   }



 
