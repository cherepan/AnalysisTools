#include "BackgroundEstimator.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"
#include <algorithm>
BackgroundEstimator::BackgroundEstimator(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

BackgroundEstimator::~BackgroundEstimator(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "BackgroundEstimator::~BackgroundEstimator Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "BackgroundEstimator::~BackgroundEstimator()" << std::endl;
}

void  BackgroundEstimator::Configure(){
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
  
  aqcd  =HConfig.GetTH1D(Name+"_aqcd","aqcd",3,-1.5,1.5,"aqcd","");     
  bqcd  =HConfig.GetTH1D(Name+"_bqcd","bqcd",3,-1.5,1.5,"bqcd","");     
  c1qcd  =HConfig.GetTH1D(Name+"_c1qcd","c1qcd",3,-1.5,1.5,"c1qcd","");     
  d1qcd  =HConfig.GetTH1D(Name+"_d1qcd","d1qcd",3,-1.5,1.5,"d1qcd","");     
  c2qcd  =HConfig.GetTH1D(Name+"_c2qcd","c2qcd",3,-1.5,1.5,"c2qcd","");     
  d2qcd  =HConfig.GetTH1D(Name+"_d2qcd","d2qcd",3,-1.5,1.5,"d2qcd","");     

  
  wj65os=HConfig.GetTH1D(Name+"_wj65os","wj65os",3,-1.5,1.5,"wj65os","");     
  wj65ss=HConfig.GetTH1D(Name+"_wj65ss","wj65ss",3,-1.5,1.5,"wj65ss","");     

  wj25os=HConfig.GetTH1D(Name+"_wj25os","wj25os",3,-1.5,1.5,"wj25os","");     
  wj25ss=HConfig.GetTH1D(Name+"_wj25ss","wj25ss",3,-1.5,1.5,"wj25ss","");     
  
  HPSTauDiscriminators=HConfig.GetTH1D(Name+"_HPSTauDiscriminators","HPSTauDiscriminators",20,0.5,20.5,"HPSTauDiscriminators","");     
  
  ChargeAll=HConfig.GetTH1D(Name+"_ChargeAll","ChargeAll",3,-1.5,1.5,"pair charge","");     

  ChargeAllTauIso=HConfig.GetTH2D(Name+"_ChargeAllTauIso","",3,-1.5,1.5,2,-0.5,1.5,"pair charge","");     
  ChargeAllTauHPS=HConfig.GetTH2D(Name+"_ChargeAllTauHPS","",3,-1.5,1.5,11,0.5,10.5,"pair charge","");     
  ChargeAllMuonIso=HConfig.GetTH2D(Name+"_ChargeAllMuonIso","",3,-1.5,1.5,50,0,1,"pair charge","");     
  ChargeAllMuonID=HConfig.GetTH2D(Name+"_ChargeAllMuonID","",3,-1.5,1.5,5,-0.5,4.5,"pair charge","");     
  ChargeAllSign=HConfig.GetTH2D(Name+"_ChargeAllSign","",3,-1.5,1.5,50,0,10,"pair charge","");     
  ChargeAllMET=HConfig.GetTH2D(Name+"_ChargeAllMET","",3,-1.5,1.5,100,0,100,"pair charge","");     
  
  MissMass1=HConfig.GetTH1D(Name+"_MissMass1","MissMass1",100,0,100,"Mt, GeV","");     
  MissMass2=HConfig.GetTH1D(Name+"_MissMass2","MissMass2",100,0,100,"Mt, GeV","");     

  TauIsolTightEfficiency=HConfig.GetTH1D(Name+"_TauIsolTightEfficiency","TauIsolTightEfficiency",10,-0.5,10.5," ","");     

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
} 


void  BackgroundEstimator::Store_ExtraDist(){
  
  Extradist1d.push_back(&aqcd);  
  Extradist1d.push_back(&bqcd);  
  Extradist1d.push_back(&c1qcd);  
  Extradist1d.push_back(&d1qcd);  
  Extradist1d.push_back(&c2qcd);  
  Extradist1d.push_back(&d2qcd);  
			   
			   
  Extradist1d.push_back(&wj65os);
  Extradist1d.push_back(&wj65ss);
			   
  Extradist1d.push_back(&wj25os);
  Extradist1d.push_back(&wj25ss);
  Extradist1d.push_back(&HPSTauDiscriminators);
  Extradist1d.push_back(&ChargeAll);

  Extradist2d.push_back(&ChargeAllTauIso);
  Extradist2d.push_back(&ChargeAllTauHPS);
  Extradist2d.push_back(&ChargeAllMuonIso);
  Extradist2d.push_back(&ChargeAllMuonID);
  Extradist2d.push_back(&ChargeAllSign);
  Extradist2d.push_back(&ChargeAllMET);


  Extradist1d.push_back(&MissMass1);
  Extradist1d.push_back(&MissMass2);
  Extradist1d.push_back(&TauIsolTightEfficiency);
}

void  BackgroundEstimator::doEvent(){
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
	  //if(Ntp->PFTau_isHPSByDecayModeFinding(i) && Ntp->PFTau_isHPSAgainstMuonTight(i) && Ntp->PFTau_isHPSAgainstElectronsTight(i)  ){
	  if(  Ntp->PFTau_TIP_hassecondaryVertex(i) )
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
       Ntp->Muon_p4(i).Pt() > 20 ){
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
    value.at(TauIsMediumIsolated) =   Ntp->PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits(TauCandidate);
    value.at(TauHPSID) = Ntp->PFTau_hpsDecayMode(TauCandidate);

    //  if( Ntp->PFTau_hpsDecayMode(TauCandidate) == 1)    std::cout<<"------------------------> Tau HPS id    "<< Ntp->PFTau_hpsDecayMode(TauCandidate)  <<std::endl;

    value.at(MuonIso) =  ((Ntp->Muon_sumChargedHadronPt04(MuonCandidate) + max(0.,Ntp->Muon_sumNeutralHadronEt04(MuonCandidate) + Ntp->Muon_sumPhotonEt04(MuonCandidate)  - 0.5*Ntp->Muon_sumPUPt04(MuonCandidate)))/ Ntp->Muon_p4(MuonCandidate).Pt());
    value.at(MuonID) =  Ntp->isMuonID(MuonCandidate); 

    //    std::cout<<"  1  "<<std::endl; 
    value.at(FlightLenghtSignificance)=Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate),Ntp->PFTau_TIP_primaryVertex_cov(TauCandidate),Ntp->PFTau_a1_lvp(TauCandidate).Vertex(),Ntp->PFTau_a1_lvp(TauCandidate).VertexCov());
    value.at(MET) =   sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) );
    value.at(charge) =Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate);
    value.at(Mass) =(Ntp->Muon_p4(MuonCandidate) + Ntp->PFTau_p4(TauCandidate)).M();


    value.at(PUJetID) = Ntp->PFJet_PUJetID_discr(PFJetIndex);
    //    std::cout<<"  2  "<<std::endl; 
    ChargeAll.at(t).Fill(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate),1);
    //    std::cout<<"  3  "<<std::endl; 
    ChargeAllTauIso.at(t).Fill(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate),Ntp->PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits(TauCandidate),1);
    //    std::cout<<"  4  "<<std::endl; 
    ChargeAllTauHPS.at(t).Fill(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate),Ntp->PFTau_hpsDecayMode(TauCandidate),1);
    //    std::cout<<"  5  "<<std::endl; 
    ChargeAllMuonIso.at(t).Fill(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate),((Ntp->Muon_sumChargedHadronPt04(MuonCandidate) + max(0.,Ntp->Muon_sumNeutralHadronEt04(MuonCandidate) + Ntp->Muon_sumPhotonEt04(MuonCandidate)  - 0.5*Ntp->Muon_sumPUPt04(MuonCandidate)))/ Ntp->Muon_p4(MuonCandidate).Pt()),1);
    //   std::cout<<"  6  "<<std::endl; 
    ChargeAllMuonID.at(t).Fill(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate), Ntp->isMuonID(MuonCandidate),1);
    // ChargeAllSign.at(t).Fill(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate),Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate),Ntp->PFTau_TIP_primaryVertex_cov(TauCandidate),Ntp->PFTau_a1_lvp(TauCandidate).Vertex(),Ntp->PFTau_a1_lvp(TauCandidate).VertexCov()),1);
    //    std::cout<<"  7  "<<std::endl; 
    ChargeAllMET.at(t).Fill(Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate),sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) ),1);
  }

  //    std::cout<<"  3  "<<std::endl; 
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  pass.at(hasTau)=(value.at(hasTau)==cut.at(hasTau));
  pass.at(TauIsMediumIsolated)=(value.at(TauIsMediumIsolated) == cut.at(TauIsMediumIsolated));
  pass.at(TauHPSID)=(value.at(TauHPSID) ==cut.at(TauHPSID));

  pass.at(hasMuon)=(value.at(hasMuon) ==cut.at(hasMuon));
  pass.at(MuonIso)=true;//(value.at(MuonIso) <= cut.at(MuonIso));

  pass.at(MuonID) = (value.at(MuonID) == cut.at(MuonID) && addMuonID);

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
    //   std::cout<<"------------------------- charge tau      "<< Ntp->PFTau_Charge(TauCandidate)  << "------------------------- charge mmuon      "<< Ntp->Muon_Charge(MuonCandidate) <<std::endl;

    //-------------------------- ABCD with OS SS and Muon isolation ---------------------
    //-------------------------- gives an estimation of combined co ---------------------
    //-------------------------- contribution QCD and W+Jet backgrd ---------------------

    double MissingTransverseMass = sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0rtT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0rtT1Txy_phi())) );
    double MissingTransverseMass2 = sqrt(2*Ntp->Muon_p4(MuonCandidate).Pt()*Ntp->MET_CorrT0pcT1Txy_et()*(1  - cos(Ntp->Muon_p4(MuonCandidate).Phi() - Ntp->MET_CorrT0pcT1Txy_phi())) );

    MissMass1.at(t).Fill(MissingTransverseMass,w);
    MissMass2.at(t).Fill(MissingTransverseMass2,w);

    double PairCharge = Ntp->PFTau_Charge(TauCandidate)*Ntp->Muon_Charge(MuonCandidate);
    double MuonIsolation = ((Ntp->Muon_sumChargedHadronPt04(MuonCandidate) + max(0.,Ntp->Muon_sumNeutralHadronEt04(MuonCandidate) + Ntp->Muon_sumPhotonEt04(MuonCandidate)  - 0.5*Ntp->Muon_sumPUPt04(MuonCandidate)))/ Ntp->Muon_p4(MuonCandidate).Pt());
    //bool MID = Ntp->isMuonID(MuonCandidate); 
    
     //  region A
     if(MuonIsolation<0.12 && PairCharge == -1 && MissingTransverseMass < 30 &&  Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate) ){
       aqcd.at(t).Fill(1,w); 
       TauIsolTightEfficiency.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits(TauCandidate),1);
     }

     //  region B
     if(MuonIsolation<0.12 && PairCharge == 1 && MissingTransverseMass < 30 &&  Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate) ){
       bqcd.at(t).Fill(1,w); 
    
       //     std::cout<<"---- "<<std::endl;
     }


     //  region C2
     if(MuonIsolation>0.2 && PairCharge == -1 && MissingTransverseMass < 30 ){
       c2qcd.at(t).Fill(1,w); 
    
       if(Ntp->PFTau_isTightIsolationDBSumPtCorr(TauCandidate))HPSTauDiscriminators.at(t).Fill(1,w); 
       if(Ntp->PFTau_isMediumIsolationDBSumPtCorr(TauCandidate))HPSTauDiscriminators.at(t).Fill(2,w); 
       if(Ntp->PFTau_isLooseIsolationDBSumPtCorr(TauCandidate))HPSTauDiscriminators.at(t).Fill(3,w); 
       if(Ntp->PFTau_isVLooseIsolationDBSumPtCorr(TauCandidate))HPSTauDiscriminators.at(t).Fill(4,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection(TauCandidate))HPSTauDiscriminators.at(t).Fill(5,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection(TauCandidate))HPSTauDiscriminators.at(t).Fill(6,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection(TauCandidate))HPSTauDiscriminators.at(t).Fill(7,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection(TauCandidate))HPSTauDiscriminators.at(t).Fill(8,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits(TauCandidate))HPSTauDiscriminators.at(t).Fill(9,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate))HPSTauDiscriminators.at(t).Fill(10,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits(TauCandidate))HPSTauDiscriminators.at(t).Fill(11,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA(TauCandidate))HPSTauDiscriminators.at(t).Fill(12,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA(TauCandidate))HPSTauDiscriminators.at(t).Fill(13,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA(TauCandidate))HPSTauDiscriminators.at(t).Fill(14,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2(TauCandidate))HPSTauDiscriminators.at(t).Fill(15,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2(TauCandidate))HPSTauDiscriminators.at(t).Fill(16,w); 
       if(Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA2(TauCandidate))HPSTauDiscriminators.at(t).Fill(17,w); 
    
    
    
     }
     //  region C1
     if(MuonIsolation<0.2 && MuonIsolation>0.12 && PairCharge == -1 && MissingTransverseMass < 30 ){
       c1qcd.at(t).Fill(1,w); 

        }
     //  region D2
     if(MuonIsolation>0.2 && PairCharge == 1 && MissingTransverseMass < 30 ){
       d2qcd.at(t).Fill(1,w); 
     }

     //  region D1
     if(MuonIsolation<0.2 && MuonIsolation>0.12 && PairCharge == 1 && MissingTransverseMass < 30 ){
       d1qcd.at(t).Fill(1,w); 
     }
  
     //   OS High MET
    
     if(MuonIsolation<0.12 && PairCharge == -1 && MissingTransverseMass > 70){
       wj65os.at(t).Fill(1,w); 
     }


     //   SS High  MET
     if(MuonIsolation<0.12 && PairCharge == 1 && MissingTransverseMass > 70){
       wj65ss.at(t).Fill(1,w); 
     }

     //   OS Low  MET
  
     if(MuonIsolation<0.12 && PairCharge == -1 && MissingTransverseMass <30){
       wj25os.at(t).Fill(1,w); 
     }
  

     //   SS Low  MET
     if(MuonIsolation<0.12 && PairCharge == 1 && MissingTransverseMass < 30){
       wj25ss.at(t).Fill(1,w); 
     }
//---------------------------------------------------------
//     //  region A
//     if(MuonIsolation<0.12 && PairCharge == -1  ){
//       aqcd.at(t).Fill(1,1); 
      
//     }

//     //  region B
//     if(MuonIsolation<0.12 && PairCharge == 1 ){
//       bqcd.at(t).Fill(1,1); 
      
//       //     std::cout<<"---- "<<std::endl;
//     }

   
//     //  region C2
//     if(MuonIsolation>0.2 && PairCharge == -1  ){
//       c2qcd.at(t).Fill(1,1); 
      
//       if(Ntp->PFTau_isTightIsolationDBSumPtCorr(TauCandidate))HPSTauDiscriminators.at(t).Fill(1,1); 
//       if(Ntp->PFTau_isMediumIsolationDBSumPtCorr(TauCandidate))HPSTauDiscriminators.at(t).Fill(2,1); 
//       if(Ntp->PFTau_isLooseIsolationDBSumPtCorr(TauCandidate))HPSTauDiscriminators.at(t).Fill(3,1); 
//       if(Ntp->PFTau_isVLooseIsolationDBSumPtCorr(TauCandidate))HPSTauDiscriminators.at(t).Fill(4,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection(TauCandidate))HPSTauDiscriminators.at(t).Fill(5,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection(TauCandidate))HPSTauDiscriminators.at(t).Fill(6,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection(TauCandidate))HPSTauDiscriminators.at(t).Fill(7,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection(TauCandidate))HPSTauDiscriminators.at(t).Fill(8,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits(TauCandidate))HPSTauDiscriminators.at(t).Fill(9,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(TauCandidate))HPSTauDiscriminators.at(t).Fill(10,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits(TauCandidate))HPSTauDiscriminators.at(t).Fill(11,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA(TauCandidate))HPSTauDiscriminators.at(t).Fill(12,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA(TauCandidate))HPSTauDiscriminators.at(t).Fill(13,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA(TauCandidate))HPSTauDiscriminators.at(t).Fill(14,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2(TauCandidate))HPSTauDiscriminators.at(t).Fill(15,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2(TauCandidate))HPSTauDiscriminators.at(t).Fill(16,1); 
//       if(Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA2(TauCandidate))HPSTauDiscriminators.at(t).Fill(17,1); 
      
      
      
//     }

//     //  region C1
//     if(MuonIsolation<0.2 && MuonIsolation>0.12 && PairCharge == -1  ){
//       c1qcd.at(t).Fill(1,1); 

//        }
//     //  region D2
//     if(MuonIsolation>0.2 && PairCharge == 1  ){
//       d2qcd.at(t).Fill(1,1); 
//     }

//     //  region D1
//     if(MuonIsolation<0.2 && MuonIsolation>0.12 && PairCharge == 1  ){
//       d1qcd.at(t).Fill(1,1); 
//     }
    
//     //   OS High MET
      
//     if(MuonIsolation<0.12 && PairCharge == -1){
//       wj65os.at(t).Fill(1,1); 
//     }


//     //   SS High  MET
//     if(MuonIsolation<0.12 && PairCharge == 1){
//       wj65ss.at(t).Fill(1,1); 
//     }

//     //   OS Low  MET
    
//     if(MuonIsolation<0.12 && PairCharge == -1 ){
//       wj25os.at(t).Fill(1,1); 
//     }
    

//     //   SS Low  MET
//     if(MuonIsolation<0.12 && PairCharge == 1  ){
//       wj25ss.at(t).Fill(1,1); 
//     }


  }
}


double BackgroundEstimator::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs((-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt());
}

double BackgroundEstimator::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2));
}





 

