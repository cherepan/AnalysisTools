
#include "Skimmer.h"
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

Skimmer::Skimmer(TString Name_, TString id_):
  Selection(Name_,id_)
  ,muoneta(2.1)
  ,jeteta(2.3)
  ,TauTrackPtThreshold(5.0)
		//  ,ran()
{
  verbose=false;
}

Skimmer::~Skimmer(){
  for(int j=0; j<Npassed.size(); j++){
    //std::cout << "Skimmer::~Skimmer Selection Summary before: "
	// << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 //<< Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  //std::cout << "Skimmer::~Skimmer()" << std::endl;
}

void  Skimmer::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==hasMuon)	      cut.at(hasMuon)=1;
    if(i==hasTau)	      cut.at(hasTau)=1;

    //if(i==ZMassmax)           cut.at(ZMassmax)=80;
  //std::cout << "Setting cut no. i=" << i << std::endl;
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

  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");


  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");




  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    //std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
  //  ran = new TRandom();
}



void  Skimmer::Store_ExtraDist(){
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&NGoodVtx);



}

void  Skimmer::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
    std::cout << "Skimmer Ntp->GetMCID(): " << Ntp->GetMCID() << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){/* std::cout << "failed to find id" <<std::endl;*/ return;}

  if(verbose)std::cout << "void  Skimmer::doEvent() Mu" << std::endl;


//   // check out triggers
//   for(unsigned int itr = 0; itr <  Ntp->NHLTTriggers(); itr ++){
//     std::cout << "Trigger names  " <<Ntp->HTLTriggerName(itr) << std::endl;
//   }




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
    if(fabs(Ntp->Muon_p4(i).Eta())<2.1 && Ntp->Muon_p4(i).Pt() > 17 && Ntp->Muon_isGlobalMuon(i) && Ntp->Muon_isPFMuon(i)){
      mu_idx_good.push_back(i);
    }  
  }

  if(verbose)std::cout << "void  Skimmer::doEvent() MuA" << std::endl;
  //std::cout << "nmus = " << mu_idx_iso.size() << std::endl;
  value.at(hasMuon)=mu_idx_good.size();
  pass.at(hasMuon)=(value.at(hasMuon)==cut.at(hasMuon));



  if(verbose)std::cout << "void  Skimmer::doEvent() MuEnd - mu_idx = " << mu_idx << " NMuons = " << Ntp->NMuons() << std::endl;
  if(verbose)std::cout << "void  Skimmer::doEvent() Tau" << std::endl;
  
  std::vector<unsigned int> tau_idx_pt, tau_idx_eta, tau_idx_iso, tau_idx_good ;
  unsigned int tau_idx(999);
  for(unsigned int i=0;i<Ntp->NPFTaus();i++){
    if(Ntp->PFTau_p4(i).Pt()> 20){
      if(fabs(Ntp->PFTau_p4(i).Eta())<2.3){
	if( Ntp->PFTau_isHPSByDecayModeFinding(i)  && (Ntp->PFTau_hpsDecayMode(i) == 10 || Ntp->PFTau_hpsDecayMode(i) == 1 )){ 
	  tau_idx_good.push_back(i);
	}
      }
    }
  }
  if(verbose)std::cout << "void  Skimmer::doEvent() Tau -a" << std::endl;
  //event
  //should have only one element, shouldn't it? for(int i=0;i<tau_idx_iso.size();i++){
  

  if(verbose)std::cout << "void  Skimmer::doEvent() Tau -b" << std::endl;
//   std::cout << "tau_idx_eta.size    = " << tau_idx_eta.size() << std::endl;
//   std::cout << "tau_idx_pt.size     = " << tau_idx_pt.size() << std::endl;



  value.at(hasTau)=tau_idx_good.size();
  pass.at(hasTau)=(value.at(hasTau)==cut.at(hasTau));



//     std::cout << "tau_idx_iso.size: " << value.at(hasTau) << " mu_idx_iso.size: " << value.at(hasTau) << std::endl;

  //else{std::cout << "no Tau cuts made since no mu passed mu-cuts" << std::endl;}  
  if(verbose)std::cout << "void  Skimmer::doEvent() Tau -c " << std::endl;



  double wobs(1),w(1),w1(1);
   if(!Ntp->isData()){
     w1*=Ntp->PUWeight3D();
   }
   else{w1=1;}
 
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  
  ///////////////////////////////////////////////////////////
  // Add plots

  
  if(verbose)std::cout << "void  Skimmer::doEvent() -Tau done" << std::endl;
  if(status){
    if(verbose)std::cout << "MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w1);;
  }
}
