#include "WJetCorrectionFactor.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"
#include <algorithm>
WJetCorrectionFactor::WJetCorrectionFactor(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

WJetCorrectionFactor::~WJetCorrectionFactor(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "WJetCorrectionFactor::~WJetCorrectionFactor Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "WJetCorrectionFactor::~WJetCorrectionFactor()" << std::endl;
}

void  WJetCorrectionFactor::Configure(){
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
    if(i==Mass)                        cut.at(Mass)=-1;




 

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
      hlabel="$sigma_{PV-SV}";
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

   else if(i==Mass){
      title.at(i)="Pair Mass";
      title.at(i)+=cut.at(Mass);
      hlabel="Mass";    
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mass_",htitle,3,-1.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mass_",htitle,3,-1.5,1.5,hlabel,"Events"));
    } 





  }

  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  MetSSTauIso=HConfig.GetTH1D(Name+"_MetSSTauIso","MetSSTauIso",100,0.,100,"MetSSTauIso","");     
  MetOSTauIso=HConfig.GetTH1D(Name+"_MetOSTauIso","MetOSTauIso",100,0.,100,"MetOSTauIso","");     

  MetSSTauNonIso=HConfig.GetTH1D(Name+"_MetSSTauNonIso","MetSSTauNonIso",100,0.,100,"MetSSTauNonIso","");     
  MetOSTauNonIso=HConfig.GetTH1D(Name+"_MetOSTauNonIso","MetOSTauNonIso",100,0.,100,"MetOSTauNonIso","");    


  MetSSTauIsoSig2=HConfig.GetTH1D(Name+"_MetSSTauIsoSig2","MetSSTauIsoSig2",100,0.,100,"MetSSTauIsoSig2","");     
  MetOSTauIsoSig2=HConfig.GetTH1D(Name+"_MetOSTauIsoSig2","MetOSTauIsoSig2",100,0.,100,"MetOSTauIsoSig2","");     
  MetSSTauNonIsoSig2=HConfig.GetTH1D(Name+"_MetSSTauNonIsoSig2","MetSSTauNonIsoSig2",100,0.,100,"MetSSTauNonIsoSig2","");     
  MetOSTauNonIsoSig2=HConfig.GetTH1D(Name+"_MetOSTauNonIsoSig2","MetOSTauNonIsoSig2",100,0.,100,"MetOSTauNonIsoSig2","");     
  
  MetSSTauIsoSig1=HConfig.GetTH1D(Name+"_MetSSTauIsoSig1","MetSSTauIsoSig1",100,0.,100,"MetSSTauIsoSig1","");     
  MetOSTauIsoSig1=HConfig.GetTH1D(Name+"_MetOSTauIsoSig1","MetOSTauIsoSig1",100,0.,100,"MetOSTauIsoSig1","");     
  MetSSTauNonIsoSig1=HConfig.GetTH1D(Name+"_MetSSTauNonIsoSig1","MetSSTauNonIsoSig1",100,0.,100,"MetSSTauNonIsoSig1","");     
  MetOSTauNonIsoSig1=HConfig.GetTH1D(Name+"_MetOSTauNonIsoSig1","MetOSTauNonIsoSig1",100,0.,100,"MetOSTauNonIsoSig1","");     

  MetSSTauIsoSig0=HConfig.GetTH1D(Name+"_MetSSTauIsoSig0","MetSSTauIsoSig0",100,0.,100,"MetSSTauIsoSig0","");     
  MetOSTauIsoSig0=HConfig.GetTH1D(Name+"_MetOSTauIsoSig0","MetOSTauIsoSig0",100,0.,100,"MetOSTauIsoSig0","");     
  MetSSTauNonIsoSig0=HConfig.GetTH1D(Name+"_MetSSTauNonIsoSig0","MetSSTauNonIsoSig0",100,0.,100,"MetSSTauNonIsoSig0","");     
  MetOSTauNonIsoSig0=HConfig.GetTH1D(Name+"_MetOSTauNonIsoSig0","MetOSTauNonIsoSig0",100,0.,100,"MetOSTauNonIsoSig0","");    
 
  RatioTauNonIsoOS=HConfig.GetTH1D(Name+"_RatioTauNonIsoOS","RatioTauNonIsoOS",3,0.5,3.5,"","");     
  RatioTauNonIsoSS=HConfig.GetTH1D(Name+"_RatioTauNonIsoSS","RatioTauNonIsoSS",3,0.5,3.5,"","");     

  RatioTauIsoSS=HConfig.GetTH1D(Name+"_RatioTauIsoSS","RatioTauIsoSS",3,0.5,3.5,"","");     
  RatioTauIsoOS=HConfig.GetTH1D(Name+"_RatioTauIsoOS","RatioTauIsoOS",3,0.5,3.5,"","");     


  Met1=HConfig.GetTH1D(Name+"_Met1","Met1",50,0.,100,"Met1","");     
  Met2=HConfig.GetTH1D(Name+"_Met2","Met2",50,0.,100,"Met2","");     
  Met3=HConfig.GetTH1D(Name+"_Met3","Met3",50,0.,100,"Met3","");     
  Met4=HConfig.GetTH1D(Name+"_Met4","Met4",50,0.,100,"Met4","");     
  Met5=HConfig.GetTH1D(Name+"_Met5","Met5",50,0.,100,"Met5","");     
  Met6=HConfig.GetTH1D(Name+"_Met6","Met6",50,0.,100,"Met6","");     
  Met7=HConfig.GetTH1D(Name+"_Met7","Met7",50,0.,100,"Met7","");     
  Met8=HConfig.GetTH1D(Name+"_Met8","Met8",50,0.,100,"Met8","");     
  Met9=HConfig.GetTH1D(Name+"_Met9","Met9",50,0.,100,"Met9","");     



  RatioTauMet2IsoSS=HConfig.GetTH1D(Name+"_RatioTauMet2IsoSS","RatioTauMet2IsoSS",3,0.5,3.5,"","");     
  RatioTauMet2IsoOS=HConfig.GetTH1D(Name+"_RatioTauMet2IsoOS","RatioTauMet2IsoOS",3,0.5,3.5,"","");     

  RatioTauMet3IsoSS=HConfig.GetTH1D(Name+"_RatioTauMet3IsoSS","RatioTauMet3IsoSS",3,0.5,3.5,"","");     
  RatioTauMet3IsoOS=HConfig.GetTH1D(Name+"_RatioTauMet3IsoOS","RatioTauMet3IsoOS",3,0.5,3.5,"","");     

  RatioTauMet4IsoSS=HConfig.GetTH1D(Name+"_RatioTauMet4IsoSS","RatioTauMet4IsoSS",3,0.5,3.5,"","");     
  RatioTauMet4IsoOS=HConfig.GetTH1D(Name+"_RatioTauMet4IsoOS","RatioTauMet4IsoOS",3,0.5,3.5,"","");     

  RatioTauMet5IsoSS=HConfig.GetTH1D(Name+"_RatioTauMet5IsoSS","RatioTauMet5IsoSS",3,0.5,3.5,"","");     
  RatioTauMet5IsoOS=HConfig.GetTH1D(Name+"_RatioTauMet5IsoOS","RatioTauMet5IsoOS",3,0.5,3.5,"","");     


  RatioTauMet6IsoSS=HConfig.GetTH1D(Name+"_RatioTauMet6IsoSS","RatioTauMet6IsoSS",3,0.5,3.5,"","");     
  RatioTauMet6IsoOS=HConfig.GetTH1D(Name+"_RatioTauMet6IsoOS","RatioTauMet6IsoOS",3,0.5,3.5,"","");     

  RatioTauMet7IsoSS=HConfig.GetTH1D(Name+"_RatioTauMet7IsoSS","RatioTauMet7IsoSS",3,0.5,3.5,"","");     
  RatioTauMet7IsoOS=HConfig.GetTH1D(Name+"_RatioTauMet7IsoOS","RatioTauMet7IsoOS",3,0.5,3.5,"","");     

  RatioTauMet8IsoSS=HConfig.GetTH1D(Name+"_RatioTauMet8IsoSS","RatioTauMet8IsoSS",3,0.5,3.5,"","");     
  RatioTauMet8IsoOS=HConfig.GetTH1D(Name+"_RatioTauMet8IsoOS","RatioTauMet8IsoOS",3,0.5,3.5,"","");     

  RatioTauMet9IsoSS=HConfig.GetTH1D(Name+"_RatioTauMet9IsoSS","RatioTauMet9IsoSS",3,0.5,3.5,"","");     
  RatioTauMet9IsoOS=HConfig.GetTH1D(Name+"_RatioTauMet9IsoOS","RatioTauMet9IsoOS",3,0.5,3.5,"","");     


  FlightLengthMET=HConfig.GetTH2D(Name+"_FlightLengthMET","Significance of #phi angle rotation",60,0,2,50,0,100,"","");
  SignigMET=HConfig.GetTH2D(Name+"_SignigMET","Significance of #phi angle rotation",20,0,10,50,0,100,"","");

 
  Selection::ConfigureHistograms(); 
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}  


void  WJetCorrectionFactor::Store_ExtraDist(){

 Extradist1d.push_back(&MetSSTauIso);
 Extradist1d.push_back(&MetOSTauIso);

 Extradist1d.push_back(&MetSSTauNonIso);
 Extradist1d.push_back(&MetOSTauNonIso);

 Extradist1d.push_back(&MetSSTauIsoSig2);
 Extradist1d.push_back(&MetOSTauIsoSig2);
 Extradist1d.push_back(&MetSSTauNonIsoSig2);
 Extradist1d.push_back(&MetOSTauNonIsoSig2);
 
 Extradist1d.push_back(&MetSSTauIsoSig1);
 Extradist1d.push_back(&MetOSTauIsoSig1);
 Extradist1d.push_back(&MetSSTauNonIsoSig1);
 Extradist1d.push_back(&MetOSTauNonIsoSig1);
 
 Extradist1d.push_back(&MetSSTauIsoSig0);
 Extradist1d.push_back(&MetOSTauIsoSig0);
 Extradist1d.push_back(&MetSSTauNonIsoSig0); 
 Extradist1d.push_back(&MetOSTauNonIsoSig0);

 Extradist2d.push_back(&FlightLengthMET);
 Extradist2d.push_back(&SignigMET);

 Extradist1d.push_back(&RatioTauNonIsoSS);
 Extradist1d.push_back(&RatioTauNonIsoOS);

 Extradist1d.push_back(&RatioTauIsoSS);
 Extradist1d.push_back(&RatioTauIsoOS);

 Extradist1d.push_back(&Met1);
 Extradist1d.push_back(&Met2);
 Extradist1d.push_back(&Met3);
 Extradist1d.push_back(&Met4);
 Extradist1d.push_back(&Met5);
 Extradist1d.push_back(&Met6);
 Extradist1d.push_back(&Met7);
 Extradist1d.push_back(&Met8);
 Extradist1d.push_back(&Met9);


 Extradist1d.push_back(&RatioTauMet2IsoSS);
 Extradist1d.push_back(&RatioTauMet2IsoOS);

 Extradist1d.push_back(&RatioTauMet3IsoSS);
 Extradist1d.push_back(&RatioTauMet3IsoOS);

 Extradist1d.push_back(&RatioTauMet4IsoSS); 
 Extradist1d.push_back(&RatioTauMet4IsoOS); 

 Extradist1d.push_back(&RatioTauMet5IsoSS); 
 Extradist1d.push_back(&RatioTauMet5IsoOS); 
 
 
 Extradist1d.push_back(&RatioTauMet6IsoSS); 
 Extradist1d.push_back(&RatioTauMet6IsoOS); 

 Extradist1d.push_back(&RatioTauMet7IsoSS); 
 Extradist1d.push_back(&RatioTauMet7IsoOS); 

 Extradist1d.push_back(&RatioTauMet8IsoSS); 
 Extradist1d.push_back(&RatioTauMet8IsoOS); 

 Extradist1d.push_back(&RatioTauMet9IsoSS); 
 Extradist1d.push_back(&RatioTauMet9IsoOS); 



}

void  WJetCorrectionFactor::doEvent(){
 
   unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  
  
  // Apply Selection
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
  


  unsigned int TauCandidate;
  unsigned int MuonCandidate;
  
  
  ///////////////////////////////////
  // Fill vector of object candidate

  std::vector<unsigned int> tau_idx_good;
  unsigned int tau_idx(999);
  for(unsigned int i=0;i<Ntp->NPFTaus();i++){
    if(Ntp->PFTau_p4(i).Pt()> 20){
      if(fabs(Ntp->PFTau_p4(i).Eta())<2.3){
	if(Ntp->PFTau_hpsDecayMode(i)==10 && Ntp->PFTau_isHPSByDecayModeFinding(i) && Ntp->PFTau_isHPSAgainstMuonTight(i) && Ntp->PFTau_isHPSAgainstElectronsTight(i)  ){
		  //if( Ntp->PFTau_isHPSByDecayModeFinding(i) && Ntp->PFTau_isHPSAgainstMuonTight(i) && Ntp->PFTau_isHPSAgainstElectronsTight(i)  ){
	  if(  Ntp->PFTau_TIP_hassecondaryVertex(i))
	    {
	  tau_idx_good.push_back(i);
	  }
	}
      }
    }
  }
  double temptaupt(0);
  for(unsigned int i =0; i<tau_idx_good.size(); i++){
    if(Ntp->PFTau_p4(tau_idx_good.at(i)).Pt() >temptaupt ){
      temptaupt = Ntp->PFTau_p4(tau_idx_good.at(i)).Pt();
      TauCandidate= tau_idx_good.at(i);
    }
  }






  ////////////////////////////////////////////////////
  // Find the leading pt ones if more than one found 


  std::vector<unsigned int> mu_idx_good;
  unsigned mu_idx(999);
  for(unsigned int i=0;i<Ntp->NMuons();i++){
    if(Ntp->isGoodMuon(i) && 
       fabs(Ntp->Muon_p4(i).Eta())<2.1 && 
       Ntp->Muon_p4(i).Pt() > 20){
      mu_idx_good.push_back(i);
    }  
  }


  double tempmupt(0);
  for(unsigned int i =0; i<mu_idx_good.size(); i++){
    if(Ntp->Muon_p4(mu_idx_good.at(i)).Pt() > tempmupt){
      tempmupt = Ntp->Muon_p4(mu_idx_good.at(i)).Pt();
      MuonCandidate = mu_idx_good.at(i);
    }
  }
    


  value.at(PrimeVtx)=nGoodVtx;
  
  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("eta2p1_LooseIsoPFTau"))value.at(TriggerOk)+=1;

  value.at(hasTau) =  tau_idx_good.size();
  value.at(hasMuon) = mu_idx_good.size();
  bool addMuonID = false;
  if((value.at(hasTau)==cut.at(hasTau)) && (value.at(hasMuon) ==cut.at(hasMuon))){

    unsigned int PFJetIndex(0);
    double TauJetdR(999.);
    for(unsigned int ijet = 0; ijet < Ntp->NPFJets(); ijet++){
      if(Ntp->PFTau_p4(TauCandidate).DeltaR(Ntp->PFJet_p4(ijet))  <TauJetdR ){
	TauJetdR=Ntp->PFTau_p4(TauCandidate).DeltaR(Ntp->PFJet_p4(ijet));
	PFJetIndex = ijet;
      }
    }
    addMuonID = (dz(Ntp->Muon_p4(MuonCandidate),Ntp->Muon_Poca(MuonCandidate),Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate) ) < 0.5 &&
		 dxy(Ntp->Muon_p4(MuonCandidate),Ntp->Muon_Poca(MuonCandidate),Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate) ) < 0.2 );


    //    value.at(TauIsMediumIsolated) =   Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate);
    value.at(TauIsMediumIsolated) =   Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate);
    value.at(TauHPSID) = Ntp->PFTau_hpsDecayMode(TauCandidate);

    value.at(MuonIso) =  ((Ntp->Muon_sumChargedHadronPt04(MuonCandidate) + max(0.,Ntp->Muon_sumNeutralHadronEt04(MuonCandidate) + Ntp->Muon_sumPhotonEt04(MuonCandidate)  - 0.5*Ntp->Muon_sumPUPt04(MuonCandidate)))/ Ntp->Muon_p4(MuonCandidate).Pt());
    value.at(MuonID) =  Ntp->isMuonID(MuonCandidate); 

 
    value.at(FlightLenghtSignificance)=Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate),Ntp->PFTau_TIP_primaryVertex_cov(TauCandidate),Ntp->PFTau_a1_lvp(TauCandidate).Vertex(),Ntp->PFTau_a1_lvp(TauCandidate).VertexCov());
    value.at(MET) =   sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) );
    value.at(charge) =Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate);
    value.at(PUJetID) =  Ntp->PFJet_PUJetID_discr(PFJetIndex);
    value.at(Mass) =(Ntp->Muon_p4(MuonCandidate) + Ntp->PFTau_p4(TauCandidate)).M();


  }
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  pass.at(hasTau)=(value.at(hasTau)==cut.at(hasTau));
  pass.at(TauIsMediumIsolated)=(value.at(TauIsMediumIsolated) == cut.at(TauIsMediumIsolated));
  pass.at(TauHPSID)=(value.at(TauHPSID) ==cut.at(TauHPSID));

  pass.at(hasMuon)=(value.at(hasMuon) ==cut.at(hasMuon));
  pass.at(MuonIso)=true;//(value.at(MuonIso) <= cut.at(MuonIso));
  pass.at(MuonID) = (value.at(MuonID) == cut.at(MuonID) && addMuonID );

  pass.at(PUJetID) = (value.at(PUJetID) > cut.at(PUJetID));



  pass.at(FlightLenghtSignificance)=(value.at(FlightLenghtSignificance)>=cut.at(FlightLenghtSignificance));
  pass.at(MET)=true;//(value.at(MET)<=cut.at(MET));
  pass.at(charge)=true;//(Ntp->KFTau_Fit_charge(TauCandidate)*Ntp->Muon_Charge(FirstLeadingMuon) == -1);
  if(value.at(Mass) > 40 && value.at(Mass) < 90)  pass.at(Mass)=true;
  else  pass.at(Mass)=false;




    double wobs=1;
    double w;




    if(!Ntp->isData()){w = Ntp->PUWeightFineBins();}
    else{w=1;}
    bool status=AnalysisCuts(t,w,wobs); 
    if(status){

//       FlightLengthMET.at(t).Fill( (Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate) - Ntp->PFTau_a1_lvp(TauCandidate).Vertex()).Mag()  ,sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi()))) ,w);
     
//       SignigMET.at(t).Fill(Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate),Ntp->PFTau_TIP_primaryVertex_cov(TauCandidate),Ntp->PFTau_a1_lvp(TauCandidate).Vertex(),Ntp->PFTau_a1_lvp(TauCandidate).VertexCov())
// 		      ,sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi()))) ,w);


      double FlightLenghtSign = 1;//Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate),Ntp->PFTau_TIP_primaryVertex_cov(TauCandidate),Ntp->PFTau_a1_lvp(TauCandidate).Vertex(),Ntp->PFTau_a1_lvp(TauCandidate).VertexCov());
      double MissingTransverseMass = sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) );
      

      Met1.at(t).Fill(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) ),1);
      Met2.at(t).Fill(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pcT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pcT1Txy_phi())) ),1);
      Met3.at(t).Fill(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rt_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rt_phi())) ),1);
      Met4.at(t).Fill(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pcT1_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pcT1_phi())) ),1);
      Met5.at(t).Fill(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT1_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT1_phi())) ),1);
      Met6.at(t).Fill(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pc_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pc_phi())) ),1);
      Met7.at(t).Fill(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT1Txy_phi())) ),1);
      Met8.at(t).Fill(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrMVA_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrMVA_phi())) ),1);
      Met9.at(t).Fill(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrMVAMuTau_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrMVAMuTau_phi())) ),1);

      



      // PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2
      if(Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate)){
	
	if(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate)==1) {
	  MetSSTauIso.at(t).Fill(MissingTransverseMass,w);
	  

	  if(MissingTransverseMass <30)	  RatioTauIsoSS.at(t).Fill(1,1);
	  if(MissingTransverseMass >70)	  RatioTauIsoSS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) ) <30)	  RatioTauMet2IsoSS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) ) >70)	  RatioTauMet2IsoSS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rt_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rt_phi())) ) <30)	  RatioTauMet3IsoSS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rt_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rt_phi())) ) >70)	  RatioTauMet3IsoSS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pcT1_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pcT1_phi())) ) <30)	  RatioTauMet4IsoSS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pcT1_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pcT1_phi())) ) >70)	  RatioTauMet4IsoSS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT1_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT1_phi())) ) <30)	  RatioTauMet5IsoSS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT1_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT1_phi())) ) >70)	  RatioTauMet5IsoSS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pc_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pc_phi())) ) <30)	  RatioTauMet6IsoSS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pc_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pc_phi())) ) >70)	  RatioTauMet6IsoSS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT1Txy_phi())) ) <30)	  RatioTauMet7IsoSS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT1Txy_phi())) ) >70)	  RatioTauMet7IsoSS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrMVA_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrMVA_phi())) ) <30)	  RatioTauMet8IsoSS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrMVA_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrMVA_phi())) ) >70)	  RatioTauMet8IsoSS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrMVAMuTau_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrMVAMuTau_phi())) ) <30)	  RatioTauMet9IsoSS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrMVAMuTau_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrMVAMuTau_phi())) ) >70)	  RatioTauMet9IsoSS.at(t).Fill(2,1);

	  if(FlightLenghtSign>0)	MetSSTauIsoSig0.at(t).Fill(MissingTransverseMass,w);
	  if(FlightLenghtSign>1)	MetSSTauIsoSig1.at(t).Fill(MissingTransverseMass,w);
	  if(FlightLenghtSign>2.2)	MetSSTauIsoSig2.at(t).Fill(MissingTransverseMass,w);
	  
	}
	if(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate)==-1){

	  if(MissingTransverseMass <30)	  RatioTauIsoOS.at(t).Fill(1,1);
	  if(MissingTransverseMass >70)	  RatioTauIsoOS.at(t).Fill(2,1);



	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) ) <30)	  RatioTauMet2IsoOS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) ) >70)	  RatioTauMet2IsoOS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rt_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rt_phi())) ) <30)	  RatioTauMet3IsoOS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rt_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rt_phi())) ) >70)	  RatioTauMet3IsoOS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pcT1_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pcT1_phi())) ) <30)	  RatioTauMet4IsoOS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pcT1_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pcT1_phi())) ) >70)	  RatioTauMet4IsoOS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT1_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT1_phi())) ) <30)	  RatioTauMet5IsoOS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT1_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT1_phi())) ) >70)	  RatioTauMet5IsoOS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pc_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pc_phi())) ) <30)	  RatioTauMet6IsoOS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pc_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pc_phi())) ) >70)	  RatioTauMet6IsoOS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT1Txy_phi())) ) <30)	  RatioTauMet7IsoOS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT1Txy_phi())) ) >70)	  RatioTauMet7IsoOS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrMVA_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrMVA_phi())) ) <30)	  RatioTauMet8IsoOS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrMVA_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrMVA_phi())) ) >70)	  RatioTauMet8IsoOS.at(t).Fill(2,1);

	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrMVAMuTau_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrMVAMuTau_phi())) ) <30)	  RatioTauMet9IsoOS.at(t).Fill(1,1);
	  if(sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrMVAMuTau_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrMVAMuTau_phi())) ) >70)	  RatioTauMet9IsoOS.at(t).Fill(2,1);

	  MetOSTauIso.at(t).Fill(MissingTransverseMass,w);
	  
	  if(FlightLenghtSign>0)	MetOSTauIsoSig0.at(t).Fill(MissingTransverseMass,w);
	  if(FlightLenghtSign>1)	MetOSTauIsoSig1.at(t).Fill(MissingTransverseMass,w);
	  if(FlightLenghtSign>2.2)	MetOSTauIsoSig2.at(t).Fill(MissingTransverseMass,w);
	  
	  
	}
	

	
      }
      //     if(Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2(Ntp->KFTau_MatchedHPS_idx(FirstLeadingTau)))
      {
	
	if(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate)==1)  {
	  
	  if(MissingTransverseMass <30)	  RatioTauNonIsoSS.at(t).Fill(1,1);
	  if(MissingTransverseMass >70)	  RatioTauNonIsoSS.at(t).Fill(2,1);


	  MetSSTauNonIso.at(t).Fill(MissingTransverseMass,w);
	  
	  if(FlightLenghtSign>0)   MetSSTauNonIsoSig0.at(t).Fill(MissingTransverseMass,w);
	  if(FlightLenghtSign>1)   MetSSTauNonIsoSig1.at(t).Fill(MissingTransverseMass,w);
	  if(FlightLenghtSign>2.2)   MetSSTauNonIsoSig2.at(t).Fill(MissingTransverseMass,w);
	}
	if(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate)==-1) {
	  
	  if(MissingTransverseMass <30)	  RatioTauNonIsoOS.at(t).Fill(1,1);
	  if(MissingTransverseMass >70)	  RatioTauNonIsoOS.at(t).Fill(2,1);

	  MetOSTauNonIso.at(t).Fill(MissingTransverseMass,w);
	  
	  if(FlightLenghtSign>0)  MetOSTauNonIsoSig0.at(t).Fill(MissingTransverseMass,w);
	  if(FlightLenghtSign>1)  MetOSTauNonIsoSig1.at(t).Fill(MissingTransverseMass,w);
	  if(FlightLenghtSign>2.2)  MetOSTauNonIsoSig2.at(t).Fill(MissingTransverseMass,w);
	  
	}
      }
    }
}



double WJetCorrectionFactor::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs((-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt());
}

double WJetCorrectionFactor::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2));
}






 
