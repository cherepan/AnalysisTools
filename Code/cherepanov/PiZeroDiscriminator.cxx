
#include "PiZeroDiscriminator.h"
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

PiZeroDiscriminator::PiZeroDiscriminator(TString Name_, TString id_):
  Selection(Name_,id_)
  ,muoneta(2.1)
  ,jeteta(2.3)
  ,TauTrackPtThreshold(5.0)
		//  ,ran()
{
  verbose=true;
}

PiZeroDiscriminator::~PiZeroDiscriminator(){
  for(int j=0; j<Npassed.size(); j++){
    //std::cout << "PiZeroDiscriminator::~PiZeroDiscriminator Selection Summary before: "
	// << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 //<< Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  //std::cout << "PiZeroDiscriminator::~PiZeroDiscriminator()" << std::endl;
}

void  PiZeroDiscriminator::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==hasMuon)	      cut.at(hasMuon)=0;
    if(i==hasValidPhoton)     cut.at(hasValidPhoton)=0;
    if(i==hasTau)	      cut.at(hasTau)=0;

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

    else if(i==hasValidPhoton){
      title.at(i)="$Number of good photons$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of good photons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasValidPhoton_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasValidPhoton_",htitle,10,-0.5,9.5,hlabel,"Events"));
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

  A1Mass=HConfig.GetTH1D(Name+"_A1Mass","A1Mass",50,0.45 ,1.777,"A1 mass, GeV","");
  NPiZero=HConfig.GetTH1D(Name+"_NPiZero","NPiZero",3,-0.5,2.5,"Number of PiZeros","");

  NGamma=HConfig.GetTH1D(Name+"_NGamma","NGamma",20,-0.5,19.5,"Number of Gammas","");
  NPhotons=HConfig.GetTH1D(Name+"_NPhotons","NPhotons",20,-0.5,19.5,"Number of Photons","");

  PFTauNGamma=HConfig.GetTH1D(Name+"_PFTauNGamma","PFTauNGamma",20,-0.5,19.5,"Number of PF Gammas","");
  PiZeroMass=HConfig.GetTH1D(Name+"_Pi0Mass","Pi0Mass",50,0.05,0.25,"M_{#pi^{0}}, GeV","");
       
  Pi0EnergyResolution=HConfig.GetTH1D(Name+"_Pi0EnergyResolution","Pi0EnergyResolution",50,-20,20,"#Delta E(mc - reco), GeV","");
  Pi0DirectionResolution=HConfig.GetTH1D(Name+"_Pi0DirectionResolution","Pi0DirectionResolution",50,0,0.4,"#Delta R","");

  Pi0EnergyResolutionSingle=HConfig.GetTH1D(Name+"_Pi0EnergyResolutionSingle","Pi0EnergyResolutionSingle",50,-20,20,"#Delta E(mc - reco), GeV","");
  Pi0DirectionResolutionSingle=HConfig.GetTH1D(Name+"_Pi0DirectionResolutionSingle","Pi0DirectionResolutionSingle",50,0,0.4,"#Delta R","");


  Pi0EnergyResolutionPair=HConfig.GetTH1D(Name+"_Pi0EnergyResolutionPair","Pi0EnergyResolutionPair",50,-20,20,"#Delta E(mc - reco), GeV","");
  Pi0DirectionResolutionPair=HConfig.GetTH1D(Name+"_Pi0DirectionResolutionPair","Pi0DirectionResolutionPair",50,0,0.4,"#Delta R","");


  Pi0EnergyResolutionSingleID=HConfig.GetTH1D(Name+"_Pi0EnergyResolutionSingleID","Pi0EnergyResolutionSingleIDdddddddd",50,-20,20,"#Delta E(mc - reco), GeV","");
  Pi0DirectionResolutionSingleID=HConfig.GetTH1D(Name+"_Pi0DirectionResolutionSingleID","Pi0DirectionResolutionSingleID",50,0,0.4,"#Delta R","");


  Pi0EnergyResolutionPairID=HConfig.GetTH1D(Name+"_Pi0EnergyResolutionPairID","Pi0EnergyResolutionPairID",50,-20,20,"#Delta E(mc - reco), GeV","");
  Pi0DirectionResolutionPairID=HConfig.GetTH1D(Name+"_Pi0DirectionResolutionPairID","Pi0DirectionResolutionPairID",50,0,0.4,"#Delta R","");


  LeadingGammaEnergyResolution=HConfig.GetTH1D(Name+"_LeadingGammaEnergyResolution","LeadingGammaEnergyResolution",50,-20,20,"#Delta E(mc - reco), GeV","");
  LeadingGammaDirectionResolution=HConfig.GetTH1D(Name+"_LeadingGammaDirectionResolution","LeadingGammaDirectionResolution",50,-20,20,"#Delta R","");
  dRVsLeadingGammaEnergy=HConfig.GetTH2D(Name+"_dRVsLeadingGammaEnergy","dR vs Energy",100,0,0.6,100,0,25,"","");

    
  dRVsPi0Energy=HConfig.GetTH2D(Name+"_dRVsPi0Energy","dR vs Energy",100,0,0.6,100,0,25,"","");

  dRVsPi0EnergyPair=HConfig.GetTH2D(Name+"_dRVsPi0EnergyPair","dR vs Energy",100,0,0.6,100,0,25,"","");
  dRVsPi0EnergySingle=HConfig.GetTH2D(Name+"_dRVsPi0EnergySingle","dR vs Energy",100,0,0.6,100,0,25,"","");

  dRVsPi0EnergyPairCutLoose=HConfig.GetTH2D(Name+"_dRVsPi0EnergyPairCutLoose","dR vs Energy",100,0,0.6,100,0,25,"","");
  dRVsPi0EnergySingleCutLoose=HConfig.GetTH2D(Name+"_dRVsPi0EnergySingleCutLoose","dR vs Energy",100,0,0.6,100,0,25,"","");

  dRVsPi0EnergyPairCutTight=HConfig.GetTH2D(Name+"_dRVsPi0EnergyPairCutTight","dR vs Energy",100,0,0.6,100,0,25,"","");
  dRVsPi0EnergySingleCutTight=HConfig.GetTH2D(Name+"_dRVsPi0EnergySingleCutTight","dR vs Energy",100,0,0.6,100,0,25,"","");


  Pi0Energy=HConfig.GetTH1D(Name+"_Pi0Energy","Pi0 Energy",100,0,25,"E_{#pi^{0}}, GeV","");

  EnergySpectrSingle=HConfig.GetTH1D(Name+"_EnergySpectrSingle","EnergySpectrSingle",100,0,10,"E_{#pi^{0}}, GeV","");
  EnergySpectrPair=HConfig.GetTH1D(Name+"_EnergySpectrPair","EnergySpectrPair",100,0,10,"E_{#pi^{0}}, GeV","");


  LeadingGammaEnergy=HConfig.GetTH1D(Name+"_LeadingGammaEnergy","LeadingGammaEnergy",100,0,10,"E_{#gamma}, GeV","");
  PairSingle=HConfig.GetTH1D(Name+"_PairSingle","PairSingle",2,-0.5,1.5," 0 - single gamma, 1 - paired","");
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    //std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
  //  ran = new TRandom();
}

void  PiZeroDiscriminator::Store_ExtraDist(){
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&NGoodVtx);

  Extradist1d.push_back(&A1Mass);
  Extradist1d.push_back(&NPiZero);
  Extradist1d.push_back(&NGamma);
  Extradist1d.push_back(&NPhotons);

  Extradist1d.push_back(&PiZeroMass);
  Extradist1d.push_back(&Pi0EnergyResolution);
  Extradist1d.push_back(&Pi0DirectionResolution);

  Extradist1d.push_back(&Pi0EnergyResolutionPair);
  Extradist1d.push_back(&Pi0DirectionResolutionPair);

  Extradist1d.push_back(&Pi0EnergyResolutionSingle);
  Extradist1d.push_back(&Pi0DirectionResolutionSingle);

  Extradist1d.push_back(&Pi0EnergyResolutionPairID);
  Extradist1d.push_back(&Pi0DirectionResolutionPairID);

  Extradist1d.push_back(&Pi0EnergyResolutionSingleID);
  Extradist1d.push_back(&Pi0DirectionResolutionSingleID);



  Extradist2d.push_back(&dRVsPi0Energy);

  Extradist2d.push_back(&dRVsPi0EnergyPair);
  Extradist2d.push_back(&dRVsPi0EnergySingle);
  Extradist1d.push_back(&Pi0Energy);


  Extradist1d.push_back(&PairSingle);
  Extradist2d.push_back(&dRVsPi0EnergyPairCutLoose);
  Extradist2d.push_back(&dRVsPi0EnergySingleCutLoose);

  Extradist2d.push_back(&dRVsPi0EnergyPairCutTight);
  Extradist2d.push_back(&dRVsPi0EnergySingleCutTight);


  Extradist1d.push_back(&EnergySpectrSingle);
  Extradist1d.push_back(&EnergySpectrPair);
  Extradist1d.push_back(&LeadingGammaEnergy);

  Extradist1d.push_back(&LeadingGammaEnergyResolution);
  Extradist1d.push_back(&LeadingGammaDirectionResolution);
  Extradist2d.push_back(&dRVsLeadingGammaEnergy);



}

void  PiZeroDiscriminator::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "void  PiZeroDiscriminator::doEvent() 1 " << id<<std::endl;

  //  std::cout << "PiZeroDiscriminator Ntp->GetMCID(): " << Ntp->GetMCID() << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){/* std::cout << "failed to find id" <<std::endl;*/ return;}

  if(verbose)std::cout << "void  PiZeroDiscriminator::doEvent() Mu" << std::endl;


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
    if(Ntp->isGoodMuon(i) && fabs(Ntp->Muons_p4(i).Eta())<2.4 && Ntp->Muons_p4(i).Pt() > 17){
      mu_idx_good.push_back(i);
    }  
  }



  if(verbose)std::cout << "void  PiZeroDiscriminator::doEvent() MuA" << std::endl;
  //std::cout << "nmus = " << mu_idx_iso.size() << std::endl;
  value.at(hasMuon)=mu_idx_good.size();
  pass.at(hasMuon)=(value.at(hasMuon)>cut.at(hasMuon));



  if(verbose)std::cout << "void  PiZeroDiscriminator::doEvent() MuEnd - mu_idx = " << mu_idx << " NMuons = " << Ntp->NMuons() << std::endl;
  if(verbose)std::cout << "void  PiZeroDiscriminator::doEvent() Tau" << std::endl;
  
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
  if(verbose)std::cout << "void  PiZeroDiscriminator::doEvent() Tau -a" << std::endl;
  //event
  //should have only one element, shouldn't it? for(int i=0;i<tau_idx_iso.size();i++){

  std::vector<int> good_gammaidx;
  
  if(Ntp->PFTau_MatchedPFJetGammaSize(tau_idx_good.at(0))!=0) LeadingGammaEnergy.at(t).Fill(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E(),1);

  if(tau_idx_good.size()!=0){
    for(unsigned int ig = 0; ig <Ntp->PFTau_MatchedPFJetGammaSize(tau_idx_good.at(0));  ig++ ){
      //if(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),ig).E() > 1){
      if(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),ig).E() > 0.6){
	good_gammaidx.push_back(ig);
	
      }
    }
  }
  std::vector<int> gamma_sorted;
  while(good_gammaidx.size()!=0){
    double buff(-1.);
    int gindex(0);
    for(unsigned int j =0; j < good_gammaidx.size(); j++){
      if(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),good_gammaidx.at(j)).E() > buff){
	buff = Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),good_gammaidx.at(j)).E();
	gindex = j;
      }
    }
    gamma_sorted.push_back(good_gammaidx.at(gindex));
    good_gammaidx.erase(good_gammaidx.begin() + gindex);
  }
  
//   for(unsigned int ig =0; ig < gamma_sorted.size(); ig++){
//     std::cout<<" gamma energy   "<<  Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(ig)).E()<<std::endl;
//   }

  //std::cout<<" ngammas size  "<< good_gammaidx.size()<<std::endl;

  value.at(hasValidPhoton)=gamma_sorted.size();
  pass.at(hasValidPhoton)=(value.at(hasValidPhoton)>0);

  //std::cout<<" pass photon  "<<pass.at(hasValidPhoton) <<std::endl;

  if(verbose)std::cout << "void  PiZeroDiscriminator::doEvent() Tau -b" << std::endl;
//   std::cout << "tau_idx_eta.size    = " << tau_idx_eta.size() << std::endl;
//   std::cout << "tau_idx_pt.size     = " << tau_idx_pt.size() << std::endl;



  value.at(hasTau)=tau_idx_good.size();
  pass.at(hasTau)=(value.at(hasTau)>cut.at(hasTau));





  double wobs(1),w(1),w1(1);
   if(!Ntp->isData()){
     w1*=Ntp->EvtWeight3D();
   }
   else{w1=1;}
 
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  
  ///////////////////////////////////////////////////////////
  // Add plots

  
  if(verbose)std::cout << "void  PiZeroDiscriminator::doEvent() -Tau done" << std::endl;
  if(status){

    bool Pi0MadeOfPair = false;
    //std::cout<<" ngammas size  2  "<< gamma_sorted.size()<<std::endl;
    //  if(verbose)std::cout << "MC type: " << Ntp->GetMCID() <<std::endl;
    //std::cout << "MC type: " << Ntp->GetMCID() <<std::endl;

 std::cout<<" Npuons   "<<std::endl;
  std::cout<<" Npuons   "<< Ntp->NPions(tau_idx_good.at(0))<<std::endl;
 if(Ntp->NPions(tau_idx_good.at(0))==3){
   std::cout<<" SumPionsPx and charge   "<< (Ntp->PFTau_PionsP4(tau_idx_good.at(0),0) +Ntp->PFTau_PionsP4(tau_idx_good.at(0),2) + Ntp->PFTau_PionsP4(tau_idx_good.at(0),2) ).Px() <<"  " << Ntp->PFTau_PionsCharge(tau_idx_good.at(0),0) *Ntp->PFTau_PionsCharge(tau_idx_good.at(0),2) *Ntp->PFTau_PionsCharge(tau_idx_good.at(0),2)<< std::endl;
   std::cout<<" a1 charge   "<< Ntp->PFTau_Charge(tau_idx_good.at(0))<<std::endl;
   std::cout<<" a1 charge   "<< Ntp->PFTau_p4(tau_idx_good.at(0)).Px()<<std::endl;

 }
    //    std::cout<<" pions   "<<Ntp->NPions(tau_idx_good.at(0))<<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w1);;
   
    A1Mass.at(t).Fill(Ntp->PFTau_3PS_A1_LV(tau_idx_good.at(0)).M(),1);
    NPiZero.at(t).Fill(Ntp->PFTau_PiZeroSize(tau_idx_good.at(0)),1);
    NGamma.at(t).Fill(Ntp->PFTau_MatchedPFJetGammaSize(tau_idx_good.at(0)),1);
    PFTauNGamma.at(t).Fill(Ntp->PFTau_GammaSize(tau_idx_good.at(0)),1);

   
    TLorentzVector Pi0LV;
   
    //---------- make pair 
    double mpi0 = 0.137;
    int gammaToPair(0);
    double bufMass(999);

    TLorentzVector LeadingGamma = Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(0));

    for(unsigned int igamma =1; igamma < gamma_sorted.size(); igamma++){

      if(fabs((Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(0)) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)) ).M() - mpi0) < bufMass){
	bufMass = (Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)) ).M() - mpi0;
	gammaToPair = gamma_sorted.at(igamma);

	std::cout<<"gamma to pair index   "<< gammaToPair << "  lead  " << gamma_sorted.at(0)<<std::endl;
      }
    }

    if(fabs((Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(0)) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gammaToPair) ).M() - mpi0) < 0.2 ){
      Pi0LV = Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(0)) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gammaToPair);
      Pi0MadeOfPair=true;
    }
    else {
      Pi0LV = Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(0));
    }
    
    PairSingle.at(t).Fill(Pi0MadeOfPair,1);

    if(Pi0MadeOfPair)PiZeroMass.at(t).Fill(Pi0LV.M(),1);
    
    
    
    bool PassePairId = false;
    bool PasseSingleId = false;
    
    dRVsPi0Energy.at(t).Fill(Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Pi0LV.E(),1);
    dRVsLeadingGammaEnergy.at(t).Fill(LeadingGamma.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , LeadingGamma.E(),1);

    if(!Pi0MadeOfPair)   EnergySpectrSingle.at(t).Fill(Pi0LV.E(),1);
    if(Pi0MadeOfPair)     EnergySpectrPair.at(t).Fill(Pi0LV.E(),1);



    if(Pi0MadeOfPair)     dRVsPi0EnergyPair.at(t).Fill(Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Pi0LV.E(),1);
    if(!Pi0MadeOfPair)     dRVsPi0EnergySingle.at(t).Fill(Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Pi0LV.E(),1);

    if(Pi0MadeOfPair) {

      if((Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12) ){dRVsPi0EnergyPairCutTight.at(t).Fill(Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Pi0LV.E(),1);PassePairId = true;}

      if((Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12  && Pi0LV.E() >3 ) ||
	 (Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12  && Pi0LV.E() <3 ) ||
	 (Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) < 0.12  && Pi0LV.E() <3 ) )	dRVsPi0EnergyPairCutLoose.at(t).Fill(Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Pi0LV.E(),1); 

    }
    if(!Pi0MadeOfPair){



      if((Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12) ){dRVsPi0EnergySingleCutTight.at(t).Fill(Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Pi0LV.E(),1);PasseSingleId  = true;}

      if((Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12  && Pi0LV.E() >3 ) ||
	 (Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12  && Pi0LV.E() <3 ) ||
	 (Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) < 0.12  && Pi0LV.E() <3 ) )   dRVsPi0EnergySingleCutLoose.at(t).Fill(Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Pi0LV.E(),1); 

    }



    if(id == 10230833){
      TLorentzVector TruthPi0 = Ntp->GetTruthTauProductLV(8,111);
      
      std::cout<<"--------- NPhotons  "<< Ntp->GetTruthPhotonsFromPi0(8,22).size()<< " number of regular reco::Photon in tau cone" <<  Ntp->nPhotonInConeOf05(tau_idx_good.at(0)) <<std::endl;
      if(Ntp->nPhotonInConeOf05(tau_idx_good.at(0))!=0)std::cout<<"Energy of regular reco::Photon  "<< Ntp->PFTau_Photon_P4(tau_idx_good.at(0),0).E() <<std::endl;
      if(Ntp->GetTruthPhotonsFromPi0(8,22).size()==2){



	std::cout<<" Photon 0  E  "<< Ntp->GetTruthPhotonsFromPi0(8,22).at(0).E() <<std::endl;
	std::cout<<" Photon 1  E  "<< Ntp->GetTruthPhotonsFromPi0(8,22).at(1).E() <<std::endl;

	std::cout<<" Leading PhotEnergy  "<< LeadingGamma.E() <<std::endl;
	if(Pi0MadeOfPair)	std::cout<<" Paired  PhotEnergy  "<< Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gammaToPair).E() <<std::endl;



	std::cout<<"--------- Loop over pfjetgammas  "<< std::endl;
	double deltaRbuf1(999.);
	double deltaRbuf2(999.);
	unsigned int ClosestIndex1(0);
	unsigned int ClosestIndex2(0);
	for(unsigned int igamma =0; igamma < gamma_sorted.size(); igamma++){
	  //	  if(Ntp->GetTruthPhotonsFromPi0(8,22).at(0).DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma))) < deltaRbuf1){
	  if(fabs(Ntp->GetTruthPhotonsFromPi0(8,22).at(0).E()  - Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)).E())   < deltaRbuf1){
	    deltaRbuf1 = fabs(Ntp->GetTruthPhotonsFromPi0(8,22).at(0).E()  - Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)).E()) ; ClosestIndex1 = gamma_sorted.at(igamma);
	  }
	  
	  //	  if(Ntp->GetTruthPhotonsFromPi0(8,22).at(1).DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma))) < deltaRbuf2){
	  if(fabs(Ntp->GetTruthPhotonsFromPi0(8,22).at(1).E()  - Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)).E() )   < deltaRbuf2){
	    deltaRbuf2 =fabs(Ntp->GetTruthPhotonsFromPi0(8,22).at(1).E()  - Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)).E() ); ClosestIndex2 = gamma_sorted.at(igamma);
	  }
	  
	  
	  //	  	  std::cout<<" All   PhotEnergy  "<< Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)).E() << "  DeltaR with 1st  " <<Ntp->GetTruthPhotonsFromPi0(8,22).at(0).DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma))) << "   deltR with 2nd   " << Ntp->GetTruthPhotonsFromPi0(8,22).at(1).DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)))<<std::endl;
	  std::cout<<" All   PhotEnergy  "<< Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)).E() << "   Ivariant mass with leading " << (Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gamma_sorted.at(igamma)) + LeadingGamma).M() <<std::endl;
	}
	
	std::cout<<" Closest1  delta Energy "<< Ntp->GetTruthPhotonsFromPi0(8,22).at(0).DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),ClosestIndex1)) <<" energy   " <<Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),ClosestIndex1).E() <<std::endl;
	std::cout<<" Closest2  delta Energy "<< Ntp->GetTruthPhotonsFromPi0(8,22).at(1).DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),ClosestIndex2)) <<" energy   " <<Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),ClosestIndex2).E() <<std::endl;
	std::cout<<" InvariantMass of two above "<< (Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),ClosestIndex2) +  Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),ClosestIndex1)).M()<<std::endl;
	
      }
      
      Pi0EnergyResolution.at(t).Fill(TruthPi0.E() - Pi0LV.E(),1 );
      Pi0DirectionResolution.at(t).Fill(TruthPi0.DeltaR( Pi0LV),1 );
      LeadingGammaEnergyResolution.at(t).Fill(TruthPi0.E() - LeadingGamma.E(),1 );
      LeadingGammaDirectionResolution.at(t).Fill(TruthPi0.DeltaR( LeadingGamma),1 );


      if(Pi0MadeOfPair){
	Pi0EnergyResolutionPair.at(t).Fill(TruthPi0.E() - Pi0LV.E(),1 );
	Pi0DirectionResolutionPair.at(t).Fill(TruthPi0.DeltaR( Pi0LV),1 );

	if(!PassePairId){
	  Pi0EnergyResolutionPairID.at(t).Fill(TruthPi0.E() - Pi0LV.E(),1 );
	  Pi0DirectionResolutionPairID.at(t).Fill(TruthPi0.DeltaR( Pi0LV),1 );
	}

      }
      if(!Pi0MadeOfPair){
	Pi0EnergyResolutionSingle.at(t).Fill(TruthPi0.E() - Pi0LV.E(),1 );
	Pi0DirectionResolutionSingle.at(t).Fill(TruthPi0.DeltaR( Pi0LV),1 );
	if(!PasseSingleId){
	  Pi0EnergyResolutionSingleID.at(t).Fill(TruthPi0.E() - Pi0LV.E(),1 );
	  Pi0DirectionResolutionSingleID.at(t).Fill(TruthPi0.DeltaR( Pi0LV),1 );
	}
      }
    }



    Pi0Energy.at(t).Fill(Pi0LV.E(),1);
    
  }

}


// //     double SumGammaEnergy(0);
// //     double SumGammaEnergyInTauCone(0);
// //     for(unsigned int igamma = 0; igamma <Ntp->PFTau_MatchedPFJetGammaSize(tau_idx_good.at(0)); igamma++ ){
// //       SumGammaEnergy+=Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),igamma).E();

// //       if(Ntp->PFTau_isMatchedPFJetGammaSCValid(tau_idx_good.at(0),igamma))std::cout<<" energy  "<<  Ntp->PFTau_MatchedJetSCEnergy(tau_idx_good.at(0),igamma)<< "   id  "<< id << std::endl;
// //       if(Ntp->PFTau_3PS_A1_LV(tau_idx_good.at(0)).DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),igamma)) < 0.1){
// //       SumGammaEnergyInTauCone+=Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),igamma).E();

// //       }
// //     }
  
//     TotalGammaEnergy.at(t).Fill(SumGammaEnergy,1);
//     TotalGammaEnergyInTauCone.at(t).Fill(SumGammaEnergyInTauCone,1);
    
//     if(Ntp->PFTau_MatchedPFJetGammaSize(tau_idx_good.at(0)) > 0){EnergyOfFirst.at(t).Fill(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E(),1);
      
//       std::cout<<" energy of first Gamma  "<<  Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).Pt() << std::endl;
//       if(Ntp->PFTau_MatchedPFJetGammaSize(tau_idx_good.at(0)) > 1)std::cout<<" energy of second Gamma  "<<  Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),1).Pt() << std::endl;


//       double bufMass(999);
//       int gammaToPair(0);
//       double mpi0 = 0.137;
//       TLorentzVector Pi0LV;
//         for(unsigned int igamma = 1; igamma <Ntp->PFTau_MatchedPFJetGammaSize(tau_idx_good.at(0)); igamma++ ){


// 	  if(fabs((Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),igamma) ).M() - mpi0) < bufMass){

// 	    bufMass = fabs((Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),igamma) ).M() - mpi0);
// 	    gammaToPair = igamma;

// 	  }

// 	  if(fabs((Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gammaToPair) ).M() - mpi0) < 0.4 ){
// 	    Pi0LV = Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gammaToPair);}
// 	  else {
// 	    Pi0LV = Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0);
// 	  }


// 	  PiZeroMass.at(t).Fill((Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gammaToPair) ).M(),1);
// 	  std::cout<<"Invariant mass with all otehr photon "<< (Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),gammaToPair) ).M() <<std::endl;
// 	}
    
//       if(id == 10230833){
// 	TLorentzVector TruthPi0 = Ntp->GetTruthTauProductLV(8,111);
	
// 	TruthPi0GammaEnergy.at(t).Fill(TruthPi0.E() - Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E(),1 );
// 	TruthPi0GammaDealtR.at(t).Fill(TruthPi0.DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0)),1);
// 	std::cout<<" Truth pi0 momenta  "<< TruthPi0.Pt() << std::endl;
// 	std::cout<<" deltaR with gamma  "<< TruthPi0.DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0)) << std::endl;

// 	TruthPi0GammaEnergy2.at(t).Fill(TruthPi0.E() - Pi0LV.E(),1 );
// 	TruthPi0GammaDealtR2.at(t).Fill(TruthPi0.DeltaR( Pi0LV),1 );
//       }

//       EnergyOfPi0VsdR2.at(t).Fill(Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Pi0LV.E(),1);
//       Pi0Energy.at(t).Fill(Pi0LV.E(),1);

//       if( (Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12  && Pi0LV.E() >3 ) ||
// 	  (Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12  && Pi0LV.E() <3 ) ||
// 	  (Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) < 0.12  && Pi0LV.E() <3 ) ){

// 	  EnergyOfPi0VsdR3.at(t).Fill(Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Pi0LV.E(),1);

// 	  }

//       if((Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12 )){

// 	  EnergyOfPi0VsdR4.at(t).Fill(Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Pi0LV.E(),1);

// 	if(id == 10230833){
// 	  TLorentzVector TruthPi0 = Ntp->GetTruthTauProductLV(8,111);


// 	  TruthPi0GammaEnergyCut.at(t).Fill(TruthPi0.E() - Pi0LV.E(),1);
// 	  TruthPi0GammaDealtRCut.at(t).Fill(TruthPi0.DeltaR(Pi0LV),1);
// 	}
//       }

//       if((Pi0LV.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) < 0.12 )){
// 	if(id == 10230833){
// 	  TLorentzVector TruthPi0 = Ntp->GetTruthTauProductLV(8,111);


// 	  TruthPi0GammaEnergyCut.at(t).Fill(TruthPi0.E() - Pi0LV.E(),1);
// 	  TruthPi0GammaDealtRCut.at(t).Fill(TruthPi0.DeltaR(Pi0LV),1);
// 	}
 
//       }



//       EnergyOfFirstVsdR.at(t).Fill( Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E(),1);
//       if( (Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12 && Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E() > 3 ) || 
// 	  (Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) > 0.12 && Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E() < 3 ) || 
// 	  (Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) < 0.12 && Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E() < 3 )   ){


// 	EnergyOfFirstVsdRWithCut.at(t).Fill( Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))) , Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E(),1);


// 	if(id == 10230833){
// 	  TLorentzVector TruthPi0 = Ntp->GetTruthTauProductLV(8,111);
	  
// 	  TruthPi0GammaEnergyWithCut.at(t).Fill(TruthPi0.E() - Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E(),1 );
// 	  TruthPi0GammaDealtRWithCut.at(t).Fill(TruthPi0.DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0)),1);
// 	}
//       }
//       else {

// 	if(id == 10230833){
// 	  TLorentzVector TruthPi0 = Ntp->GetTruthTauProductLV(8,111);
	  
// 	  TruthPi0GammaEnergyWithInvertCut.at(t).Fill(TruthPi0.E() - Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E(),1 );
// 	  TruthPi0GammaDealtRWithInvertCut.at(t).Fill(TruthPi0.DeltaR(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0)),1);
// 	}

//       }



//       RatioFirstToAll.at(t).Fill(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E()/(SumGammaEnergy - Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E()) ,1 );



//     }
    
    
//     if(Ntp->PFTau_MatchedPFJetGammaSize(tau_idx_good.at(0)) > 1){EnergyOfFirstTwo.at(t).Fill(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E()+
// 											     Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),1).E(),1);

//       EnergyOfSecond.at(t).Fill(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E()+
// 				Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),1).E(),1);

//       TLorentzVector Pi0 = Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0) + Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),1).E();
//       EnergyOfPi0VsdR.at(t).Fill(Pi0.DeltaR(Ntp->PFTau_p4(tau_idx_good.at(0))),Pi0.E(),1);
//       Pi0Mass.at(t).Fill(Pi0.M(),1);
//       RatioFirstToSecond.at(t).Fill(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E()/Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),1).E() ,1);

//     }

    
//     if(Ntp->PFTau_MatchedPFJetGammaSize(tau_idx_good.at(0)) > 2)EnergyOfFirstThree.at(t).Fill(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E()+
// 											      Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),1).E()+
// 											      Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),2).E(),1);


    
//     if(Ntp->PFTau_MatchedPFJetGammaSize(tau_idx_good.at(0)) > 3)EnergyOfFirstFour.at(t).Fill(Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),0).E()+
// 											     Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),1).E()+
// 											     Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),2).E()+
// 											     Ntp->PFTau_MatchedPFJetGammaP4(tau_idx_good.at(0),3).E(),1);








//   }
// }
