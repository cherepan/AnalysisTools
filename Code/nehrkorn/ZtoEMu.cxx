#include "ZtoEMu.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

#include "Parameters.h"
#include "TMath.h"

#include <TFile.h>

ZtoEMu::ZtoEMu(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_pt(30)
  ,e_pt(30)
  ,mu_eta(2.1)
  ,e_eta(2.3)
  ,mu_reliso(0.12) //tight
  ,jet_pt(30)
  ,jet_eta(2.4)
  ,jet_sum(70)
  ,zmin(88)
  ,zmax(94)
{
    //verbose=true;
TFile* FRFile = new TFile("/net/scratch_cms/institut_3b/nehrkorn/workdirChargedHiggs_May_27_2013/Code/nehrkorn/FakeRates_2012_19ifb.root");
ElectronFakeRate = (TH2D*)(FRFile->Get("ElectronFakeRateHist"));
MuonFakeRate = (TH2D*)(FRFile->Get("MuonFakeRateHist"));
 
}

ZtoEMu::~ZtoEMu(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZtoEMu::~ZtoEMu Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZtoEMu::~ZtoEMu()" << std::endl;
}

void  ZtoEMu::Configure(){

  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==NMu)                cut.at(NMu)=1;
    if(i==NMuPt)              cut.at(NMuPt)=1;
    if(i==NMuEta)             cut.at(NMuEta)=1;
    if(i==NMuRelIso)          cut.at(NMuRelIso)=1;
    if(i==looseMuonVeto)      cut.at(looseMuonVeto)=0;
    if(i==NE)                 cut.at(NE)=1;
    if(i==NEPt)               cut.at(NEPt)=1;
    if(i==NEEta)              cut.at(NEEta)=1;
    if(i==SameVtx)            cut.at(SameVtx)=true;
    if(i==charge)             cut.at(charge)=0;
    if(i==ptBalance)          cut.at(ptBalance)=20;;
    if(i==ZMassmax)           cut.at(ZMassmax)=zmax;
    if(i==ZMassmin)           cut.at(ZMassmin)=zmin;
	if(i==jetVeto)            cut.at(jetVeto)=jet_sum;
	if(i==MtMu)               cut.at(MtMu)=40;
	if(i==drMuE)              cut.at(drMuE)=0.2;
	if(i==qualitycuts)        cut.at(qualitycuts)=true;
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
  
    if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>=$";
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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,17,-0.5,16.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,17,-0.5,16.5,hlabel,"Events"));
    }

    else if(i==NMu){
      title.at(i)="Number $\\mu =$";
      title.at(i)+=cut.at(NMu);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMu_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMu_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==NE){
      title.at(i)="Number $e =$";
      title.at(i)+=cut.at(NE);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of e";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NE_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NE_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==NMuPt){
      title.at(i)="Number of $\\mu$ [$P_{T}^{\\mu}>$";
      title.at(i)+=mu_pt;
      title.at(i)+="(GeV)] $>=$";
      title.at(i)+=cut.at(NMuPt);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NMuPt_","P_{T,#mu} (N-1 Distribution)",100,0,200,"P_{T,#mu} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_NMuPt_","P_{T,#mu} (Accumulative Distribution)",100,0,200,"P_{T,#mu} (GeV)","Events");
    }
    else if(i==NEPt){
      title.at(i)="Number of $e$ [$P_{T}^{e}>$";
      title.at(i)+=e_pt;
      title.at(i)+="(GeV)] $>=$";
      title.at(i)+=cut.at(NEPt);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of e";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NEPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NEPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NEPt_","P_{T,e} (N-1 Distribution)",100,0,200,"P_{T,e} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_NEPt_","P_{T,e} (Accumulative Distribution)",100,0,200,"P_{T,e} (GeV)","Events");
    }
    else if(i==NMuEta){
      title.at(i)="Number of $\\mu$ [$|\\eta^{\\mu}|<$";
      char buffer[50];
      sprintf(buffer,"%5.2f",mu_eta);
      title.at(i)+=buffer;
      title.at(i)+="] $>=$";
      title.at(i)+=cut.at(NMuEta);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NMuEta_","#eta{#mu} (N-1 Distribution)",100,-7,7,"#eta_{#mu} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accumdist_NMuEta_","#eta_{#mu} (Accumulative Distribution)",100,-7,7,"#eta_{#mu} (GeV)","Events");
    }
    else if(i==NEEta){
      title.at(i)="Number of $e$ [$|\\eta^{e}|<$";
      char buffer[50];
      sprintf(buffer,"%5.2f",e_eta);
      title.at(i)+=buffer;
      title.at(i)+="] $>=$";
      title.at(i)+=cut.at(NEEta);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of e";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NEEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NEEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NEEta_","#eta{e} (N-1 Distribution)",100,-7,7,"#eta_{e} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accumdist_NEEta_","#eta_{e} (Accumulative Distribution)",100,-7,7,"#eta_{e} (GeV)","Events");
    }
    else if(i==ptBalance){
      title.at(i)="$p_{t,e+\\mu} < $";
      title.at(i)+=cut.at(ptBalance);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="p_t balance / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ptBalance_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ptBalance_",htitle,40,0,200,hlabel,"Events"));
    }
    else if(i==charge){
      title.at(i)="$e-\\mu$ Charge = ";
      title.at(i)+=cut.at(charge);
      title.at(i)+="";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="e-#mu Charge";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,21,-5.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,21,-5.5,5.5,hlabel,"Events"));
    }
    else if(i==ZMassmax){
      title.at(i)="$M_{e,\\mu} < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(ZMassmax));
      title.at(i)+=buffer;
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{e,#mu} / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmax_",htitle,20,60,120,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmax_",htitle,20,60,120,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_ZMassmax_","M_{e,#mu} (N-1 Distribution)",60,0,130,"M_{e,#mu} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_ZMassmax_","M_{e,#mu} (Accumulative Distribution)",60,0,130,"M_{e,#mu} (GeV)","Events");
    }
    else if(i==ZMassmin){
      title.at(i)="$M_{e,\\mu} > $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(ZMassmin));
      title.at(i)+=buffer;
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{e,#mu} / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmin_",htitle,20,60,120,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmin_",htitle,20,60,120,hlabel,"Events"));
   } 
    else if(i==NMuRelIso){
    	title.at(i)="Relative isolation of $\\mu < $";
    	char buffer[50];
    	sprintf(buffer,"%5.2f",cut.at(NMuRelIso));
    	title.at(i)+=buffer;
    	title.at(i)+="";
    	htitle=title.at(i);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="relative isolation of #mu";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuRelIso_",htitle,20,0,1,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuRelIso_",htitle,20,0,1,hlabel,"Events"));
    }
    else if(i==looseMuonVeto){
        title.at(i)="Number of loose $\\mu = $";
        title.at(i)+=cut.at(looseMuonVeto);
        title.at(i)+="";
        htitle=title.at(i);
        htitle.ReplaceAll("$","");
        htitle.ReplaceAll("\\","#");
        hlabel="loose muon veto";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_looseMuonVeto_",htitle,40,0,20,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_looseMuonVeto_",htitle,40,0,20,hlabel,"Events"));
    }
	else if(i==jetVeto){
      title.at(i)="$P_{t}$ sum of 2 highest $P_{t}$-jets$ < $";
      title.at(i)+=cut.at(jetVeto);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="jet P_{T} / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_jetVeto_",htitle,39,2,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_jetVeto_",htitle,39,2,200,hlabel,"Events"));
    }
	else if(i==MtMu){
      title.at(i)="$m_{T}^{\\mu,Miss} < $";
      title.at(i)+=cut.at(MtMu);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="m_{T}^{#mu,Miss} / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MtMu_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MtMu_",htitle,40,0,200,hlabel,"Events"));
    }
	else if(i==drMuE){
	  title.at(i)="$dR(e,\\mu) > $";
	  char buffer[50];
	  sprintf(buffer,"%5.2f",cut.at(drMuE));
	  title.at(i)+=buffer;
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="dR(e,#mu)";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_drMuE_",htitle,400,0.,10.,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_drMuE_",htitle,400,0.,10.,hlabel,"Events"));
	}
	else if(i==SameVtx){
	  title.at(i)="$e$ and $\\mu$ from same vtx";
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SameVtx_",htitle,20,0,0.1,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SameVtx_",htitle,20,0,0.1,hlabel,"Events"));
	}
    else if(i==qualitycuts){
  	  title.at(i)="quality cuts";
  	  htitle=title.at(i);
  	  htitle.ReplaceAll("$","");
  	  htitle.ReplaceAll("\\","#");
  	  hlabel="";
  	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_qualitycuts_",htitle,20,0,0.1,hlabel,"Events"));
  	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_qualitycuts_",htitle,20,0,0.1,hlabel,"Events"));
  	}
    
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  RelIsoE=HConfig.GetTH1D(Name+"_RelIsoE","RelIsoE",20,0.,1.,"Relative Isolation of Electron");
  RelIsoMu=HConfig.GetTH1D(Name+"_RelIsoMu","RelIsoMu",20,0.,1.,"Relative Isolation of Muon");
  EPt=HConfig.GetTH1D(Name+"_PtE","PtE",40,0.,200.,"P_{t}^{e} / GeV");
  MuPt=HConfig.GetTH1D(Name+"_PtMu","PtMu",40,0.,200.,"P_{t}^{#mu} / GeV");
  mtMu=HConfig.GetTH1D(Name+"_mtMu","mtMu",40,0.,200.,"m_{t}^{#mu} / GeV");
  mtE=HConfig.GetTH1D(Name+"_mtE","mtE",40,0.,200.,"m_{t}^{e} / GeV");
  NJets=HConfig.GetTH1D(Name+"_NJets","NJets",20,0,20,"number of jets");
  pzeta=HConfig.GetTH1D(Name+"_pzeta","pzeta",40,-100.,100.,"pzeta");
  pzetaDQM=HConfig.GetTH1D(Name+"_pzetaDQM","pzetaDQM",40,-100.,100.,"pzetaDQM");
  invmass=HConfig.GetTH1D(Name+"_invmass","invmass",40,0.,200.,"m_{e,#mu} / GeV");
  etaMu=HConfig.GetTH1D(Name+"_etaMu","etaMu",40,-3.5,3.5,"#eta_{#mu}");
  etaE=HConfig.GetTH1D(Name+"_etaE","etaE",40,-3.5,3.5,"#eta_{e}");
  jetsum=HConfig.GetTH1D(Name+"_jetsum","jetsum",80,0.,400,"P_{t}^{jets} / GeV");
  chargesum=HConfig.GetTH1D(Name+"_chargesum","chargesum",21,-5.5,5.5,"charge sum");
  chargesumID=HConfig.GetTH1D(Name+"_chargesumID","chargesumID",21,-5.5,5.5,"charge sum of tight e & tight #mu");
  drmue=HConfig.GetTH1D(Name+"_drmue","drmue",20,0.,1.,"dR(e,#mu)");
  drmueID=HConfig.GetTH1D(Name+"_drmueID","drmueID",20,0.,1.,"dR(e,#mu) (tight)");
  deltaphi=HConfig.GetTH1D(Name+"_deltaphi","deltaphi",40,0.,2*TMath::Pi(),"#phi_{e,#mu}");
  deltaphiID=HConfig.GetTH1D(Name+"_deltaphiID","deltaphiID",40,0.,2*TMath::Pi(),"#phi_{e,#mu} tight");
  ptbal=HConfig.GetTH1D(Name+"_ptbal","ptbal",40,0.,200.,"p_{t}^{e+#mu} / GeV");
  chargesumsigned=HConfig.GetTH1D(Name+"_chargesumsigned","chargesumsigned",21,-5.5,5.5,"charge sum");
  chargesumIDsigned=HConfig.GetTH1D(Name+"_chargesumIDsigned","chargesumIDsigned",21,-5.5,5.5,"charge sum of tight e & tight #mu");
  FirstJetPt=HConfig.GetTH1D(Name+"_FirstJetPt","FirstJetPt",40,0.,200.,"P_{t}^{1st jet} / GeV");
  SecondJetPt=HConfig.GetTH1D(Name+"_SecondJetPt","SecondJetPt",40,0.,200.,"P_{t}^{2nd jet} / GeV");
  ThirdJetPt=HConfig.GetTH1D(Name+"_ThirdJetPt","ThirdJetPt",40,0.,200.,"P_{t}^{3rd jet} / GeV");
  FourthJetPt=HConfig.GetTH1D(Name+"_FourthJetPt","FourthJetPt",40,0.,200.,"P_{t}^{4th jet} / GeV");
  
  invmass_zmass=HConfig.GetTH1D(Name+"_invmass_zmass","invmass_zmass",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_ptbalance=HConfig.GetTH1D(Name+"_invmass_ptbalance","invmass_ptbalance",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_mtmu=HConfig.GetTH1D(Name+"_invmass_mtmu","invmass_mtmu",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_jetveto=HConfig.GetTH1D(Name+"_invmass_jetveto","invmass_jetveto",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_charge=HConfig.GetTH1D(Name+"_invmass_charge","invmass_charge",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_loosemuonveto=HConfig.GetTH1D(Name+"_invmass_loosemuonveto","invmass_loosemuonveto",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_dremu=HConfig.GetTH1D(Name+"_invmass_dremu","invmass_dremu",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_only_object_id=HConfig.GetTH1D(Name+"_invmass_only_object_id","invmass_only_object_id",20,60.,120.,"m_{e,#mu} / GeV");
  
  invmass_no_loosemuonveto=HConfig.GetTH1D(Name+"_invmass_no_loosemuonveto","invmass_no_loosemuonveto",20,60.,120.,"m_{e,#mu} / GeV");
  nm2_charge=HConfig.GetTH1D(Name+"_nm2_charge","nm2_charge",20,60.,120.,"m_{e,#mu} / GeV");
  nm2_jetveto=HConfig.GetTH1D(Name+"_nm2_jetveto","nm2_jetveto",20,60.,120.,"m_{e,#mu} / GeV");
  nm2_mtmu=HConfig.GetTH1D(Name+"_nm2_mtmu","nm2_mtmu",20,60.,120.,"m_{e,#mu} / GeV");
  nm2_ptbalance=HConfig.GetTH1D(Name+"_nm2_ptbalance","nm2_ptbalance",20,60.,120.,"m_{e,#mu} / GeV");
  nm2_drmue=HConfig.GetTH1D(Name+"_nm2_drmue","nm2_drmue",20,60.,120.,"m_{e,#mu} / GeV");
  
  mproj1=HConfig.GetTH1D(Name+"_mproj1","mproj1",9,50.,79.,"Upper left projected to mass");
  mproj2=HConfig.GetTH1D(Name+"_mproj2","mproj2",8,79.,104.,"Upper middle projected to mass");
  mproj3=HConfig.GetTH1D(Name+"_mproj3","mproj3",9,104.,130.,"Upper right projected to mass");
  mproj4=HConfig.GetTH1D(Name+"_mproj4","mproj4",9,50.,79.,"Middle left projected to mass");
  mproj5=HConfig.GetTH1D(Name+"_mproj5","mrpoj5",8,79.,104.,"Middle middle projected to mass");
  mproj6=HConfig.GetTH1D(Name+"_mproj6","mproj6",9,104.,130.,"Middle right projected to mass");
  mproj7=HConfig.GetTH1D(Name+"_mproj7","mproj7",9,50.,79.,"Lower left projected to mass");
  mproj8=HConfig.GetTH1D(Name+"_mproj8","mproj8",8,79.,104.,"Lower middle projected to mass");
  mproj9=HConfig.GetTH1D(Name+"_mproj9","mproj9",9,104.,130.,"Lower right projected to mass");
  phiproj1=HConfig.GetTH1D(Name+"_phiproj1","phiproj1",16,0.,2.5,"Upper left projected to phi");
  phiproj2=HConfig.GetTH1D(Name+"_phiproj2","phiproj2",16,0.,2.5,"Upper middle projected to phi");
  phiproj3=HConfig.GetTH1D(Name+"_phiproj3","phiproj3",16,0.,2.5,"Upper right projected to phi");
  phiproj4=HConfig.GetTH1D(Name+"_phiproj4","phiproj4",8,2.5,3.7,"Middle left projected to phi");
  phiproj5=HConfig.GetTH1D(Name+"_phiproj5","phiproj5",8,2.5,3.7,"Middle middle projected to phi");
  phiproj6=HConfig.GetTH1D(Name+"_phiproj6","phiproj6",8,2.5,3.7,"Middle right projected to phi");
  phiproj7=HConfig.GetTH1D(Name+"_phiproj7","phiproj7",17,3.7,2*TMath::Pi(),"Lower left projected to phi");
  phiproj8=HConfig.GetTH1D(Name+"_phiproj8","phiproj8",17,3.7,2*TMath::Pi(),"Lower middle projected to phi");
  phiproj9=HConfig.GetTH1D(Name+"_phiproj9","phiproj9",17,3.7,2*TMath::Pi(),"Lower right projected to phi");

  InvmassVsDeltaPhi=HConfig.GetTH2D(Name+"_InvmassVsDeltaPhi","m_{e,#mu} vs. #Delta#phi_{e,#mu}",20,60,120,40,0.,2*TMath::Pi(),"m_{e,#mu} / GeV","#Delta#phi_{e,#mu} / rad");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ZtoEMu::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist1d.push_back(&RelIsoE);
 Extradist1d.push_back(&RelIsoMu);
 Extradist1d.push_back(&EPt);
 Extradist1d.push_back(&MuPt);
 Extradist1d.push_back(&mtMu);
 Extradist1d.push_back(&mtE);
 Extradist1d.push_back(&NJets);
 Extradist1d.push_back(&pzeta);
 Extradist1d.push_back(&pzetaDQM);
 Extradist1d.push_back(&invmass);
 Extradist1d.push_back(&jetsum);
 Extradist1d.push_back(&chargesum);
 Extradist1d.push_back(&chargesumID);
 Extradist1d.push_back(&drmue);
 Extradist1d.push_back(&drmueID);
 Extradist1d.push_back(&deltaphi);
 Extradist1d.push_back(&deltaphiID);
 Extradist1d.push_back(&ptbal);
 Extradist1d.push_back(&chargesumsigned);
 Extradist1d.push_back(&chargesumIDsigned);
 Extradist1d.push_back(&FirstJetPt);
 Extradist1d.push_back(&SecondJetPt);
 Extradist1d.push_back(&ThirdJetPt);
 Extradist1d.push_back(&FourthJetPt);
 
 Extradist1d.push_back(&invmass_zmass);
 Extradist1d.push_back(&invmass_ptbalance);
 Extradist1d.push_back(&invmass_mtmu);
 Extradist1d.push_back(&invmass_jetveto);
 Extradist1d.push_back(&invmass_charge);
 Extradist1d.push_back(&invmass_loosemuonveto);
 Extradist1d.push_back(&invmass_dremu);
 Extradist1d.push_back(&invmass_only_object_id);
 
 Extradist1d.push_back(&invmass_no_loosemuonveto);
 Extradist1d.push_back(&nm2_charge);
 Extradist1d.push_back(&nm2_jetveto);
 Extradist1d.push_back(&nm2_mtmu);
 Extradist1d.push_back(&nm2_ptbalance);
 Extradist1d.push_back(&nm2_drmue);
 
 Extradist1d.push_back(&mproj1);
 Extradist1d.push_back(&mproj2);
 Extradist1d.push_back(&mproj3);
 Extradist1d.push_back(&mproj4);
 Extradist1d.push_back(&mproj5);
 Extradist1d.push_back(&mproj6);
 Extradist1d.push_back(&mproj7);
 Extradist1d.push_back(&mproj8);
 Extradist1d.push_back(&mproj9);
 Extradist1d.push_back(&phiproj1);
 Extradist1d.push_back(&phiproj2);
 Extradist1d.push_back(&phiproj3);
 Extradist1d.push_back(&phiproj4);
 Extradist1d.push_back(&phiproj5);
 Extradist1d.push_back(&phiproj6);
 Extradist1d.push_back(&phiproj7);
 Extradist1d.push_back(&phiproj8);
 Extradist1d.push_back(&phiproj9);
 
 Extradist2d.push_back(&InvmassVsDeltaPhi);

}

void  ZtoEMu::doEvent(){
  if(verbose)std::cout << "ZtoEMu::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  //id=1;
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}
  
  value.at(TriggerOk)=0;
  /*std::cout << "Number of HLT Triggers:" << Ntp->NHLTTriggers() << std::endl;
  for(unsigned int i=0; i<Ntp->NHLTTriggers(); i++){
	  std::cout << "######################################" << std::endl;
	  std::cout << "Trigger: " << Ntp->HTLTriggerName(i) << std::endl;
	  std::cout << "Acceptet (1) or not (0): " << Ntp->TriggerAccept(Ntp->HTLTriggerName(i)) << std::endl;
	  std::cout << "######################################" << std::endl;
  }*/
  if(Ntp->TriggerAccept("HLT_IsoMu24_eta2p1"))value.at(TriggerOk)=1;
  //if(Ntp->GetMCID()==40)value.at(TriggerOk)=1; // NEEDS TO BE FIXED! ONLY TEMPORARY
  pass.at(TriggerOk)=(value.at(TriggerOk)==cut.at(TriggerOk));

  // Apply Selection
  if(verbose) std::cout << "find primary vertex" << std::endl;
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));

  unsigned int pVtx(999);
  if(nGoodVtx>=cut.at(PrimeVtx))pVtx=nGoodVtx;
  
  ///////////////////////////////////////////////
  //
  // Quality cuts
  //
  if(verbose)std::cout << "Quality cuts" << std::endl;
  value.at(qualitycuts)=false;
  std::vector<int> qualitymuons;
  std::vector<int> qualityelectrons;
  qualitymuons.clear();
  qualityelectrons.clear();
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(Ntp->Muons_p4(i).Pt()>=mu_pt && Ntp->Muons_p4(i).Pt()<=70.){
		  qualitymuons.push_back(i);
	  }
  }
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(Ntp->Electron_p4(i).Pt()>=e_pt && Ntp->Electron_p4(i).Pt()<=70.){
		  qualityelectrons.push_back(i);
	  }
  }
  if(qualitymuons.size()==1 && qualityelectrons.size()==1)value.at(qualitycuts)=true;
  pass.at(qualitycuts)=(value.at(qualitycuts));
  
  ///////////////////////////////////////////////
  //
  // Blinding
  //
  bool blinddataonly = true;
  bool blindingactive = false;
  int blind = 1;
  if(pass.at(qualitycuts)){
	  double mass = (Ntp->Muons_p4(qualitymuons.at(0))+Ntp->Electron_p4(qualityelectrons.at(0))).M();
	  double delphi = Ntp->Muons_p4(qualitymuons.at(0)).DeltaPhi(Ntp->Electron_p4(qualityelectrons.at(0)));
	  if(delphi<0)delphi+=2*TMath::Pi();
	  if(blinddataonly){
		  if(Ntp->isData() && (mass>=80. && mass<=104. /*&& delphi>=2.5 && delphi<=3.7*/)){
			  blindingactive = true;
			  blind = 0;
			  //return;
		  }
	  }else{
		  if(mass>=80. && mass<=104. /*&& delphi>=2.5 && delphi<=3.7*/){
			  blindingactive = true;
			  blind = 0;
			  //return;
		  }
	  }
  }
  
  ///////////////////////////////////////////////
  //
  // Vertex constraint
  //
  if(verbose)std::cout << "Vertex constraint" << std::endl;
  int pos = 0;
  int posmu = 0;
  int pose = 0;
  int mutrack=0;
  int etrack=0;
  value.at(SameVtx)=false;
  if(qualitymuons.size()>0){
	  for(unsigned i=0;i<Ntp->NTracks();i++){
		  if(Ntp->Track_p4(i).DeltaR(Ntp->Muons_p4(qualitymuons.at(0)))<Ntp->Track_p4(mutrack).DeltaR(Ntp->Muons_p4(qualitymuons.at(0))))mutrack=i;
	  }
  }
  if(qualityelectrons.size()>0){
	  for(unsigned i=0;i<Ntp->NTracks();i++){
		  if(Ntp->Track_p4(i).DeltaR(Ntp->Electron_p4(qualityelectrons.at(0)))<Ntp->Track_p4(etrack).DeltaR(Ntp->Electron_p4(qualityelectrons.at(0)))){
			  if(Ntp->Track_p4(i).DeltaR(Ntp->Electron_p4(qualityelectrons.at(0)))<0.1 || Ntp->Track_p4(etrack).DeltaR(Ntp->Electron_p4(qualityelectrons.at(0)))<0.1){
				  etrack=i;
			  }
		  }
	  }
  }
  if(qualitymuons.size()>0 && qualityelectrons.size()>0 && mutrack==etrack){
	  if(verbose){
		  std::cout << "muon and electron from same track " << mutrack << std::endl;
		  std::cout << "dR(e,mu) = " << Ntp->Muons_p4(qualitymuons.at(0)).DeltaR(Ntp->Electron_p4(qualityelectrons.at(0))) << std::endl;
	  }
	  mutrack=-1;
	  etrack=-1;
  }
  if(verbose)std::cout << "looping over vertices" << std::endl;
  int muvtx=-1;
  int evtx=-1;
  for(unsigned i=0;i<Ntp->NVtx();i++){
	  for(unsigned j=0;j<Ntp->Vtx_Track_idx(i).size();j++){
		  if(mutrack>=0 && mutrack==Ntp->Vtx_Track_idx(i).at(j)){
			  muvtx=i;
		  }
		  if(etrack>=0 && etrack==Ntp->Vtx_Track_idx(i).at(j)){
			  evtx=i;
		  }
	  }
  }
  if(muvtx==evtx && muvtx!=-1 && evtx!=-1){
	  value.at(SameVtx)=true;
  }else if(muvtx==-1 && evtx!=-1){
	  if(mutrack!=-1 && qualitymuons.size()>0){
		  if(dxy(Ntp->Muons_p4(qualitymuons.at(0)),Ntp->Muon_Poca(qualitymuons.at(0)),Ntp->Vtx(evtx))<0.2){
			  if(dz(Ntp->Muons_p4(qualitymuons.at(0)),Ntp->Muon_Poca(qualitymuons.at(0)),Ntp->Vtx(evtx))<0.5){
				  muvtx=evtx;
				  value.at(SameVtx)=true;
			  }
		  }
	  }
  }
  else if(evtx==-1 && muvtx!=-1){
	  if(etrack!=-1 && qualityelectrons.size()>0){
		  if(dxy(Ntp->Electron_p4(qualityelectrons.at(0)),Ntp->Electron_Poca(qualityelectrons.at(0)),Ntp->Vtx(muvtx))<0.02){
			  if(dz(Ntp->Electron_p4(qualityelectrons.at(0)),Ntp->Electron_Poca(qualityelectrons.at(0)),Ntp->Vtx(muvtx))<0.1){
				  evtx=muvtx;
				  value.at(SameVtx)=true;
			  }
		  }
	  }
  }
  pass.at(SameVtx)=(value.at(SameVtx));
  if(value.at(SameVtx)){
	  pos=muvtx;
	  posmu=muvtx;
	  pose=evtx;
  }else{
	  pos=0;
	  if(muvtx!=-1){
		  posmu=muvtx;
	  }else{
		  posmu=0;
	  }
	  if(evtx!=-1){
		  pose=evtx;
	  }else{
		  pose=0;
	  }
  }

  ///////////////////////////////////////////////
  //
  // Mu Cuts
  //
  if(verbose) std::cout << "Muon cuts" << std::endl;
  std::vector<unsigned int> GoodMuons;
  bool MVA_ID = true;
  bool fakeMuon = false;

  for(unsigned i=0;i<qualitymuons.size();i++){
	  if(Ntp->isTightMuon(qualitymuons.at(i),posmu) &&
			  dxy(Ntp->Muons_p4(qualitymuons.at(i)),Ntp->Muon_Poca(qualitymuons.at(i)),Ntp->Vtx(posmu))<0.02 &&
			  dz(Ntp->Muons_p4(qualitymuons.at(i)),Ntp->Muon_Poca(qualitymuons.at(i)),Ntp->Vtx(posmu))<0.2
			  ){
		  if(fabs(Ntp->Muons_p4(qualitymuons.at(i)).Eta())<1.479 && Ntp->Muon_RelIso(qualitymuons.at(i))<0.15){
			  GoodMuons.push_back(qualitymuons.at(i));
		  }else if(fabs(Ntp->Muons_p4(qualitymuons.at(i)).Eta())>=1.479 && Ntp->Muon_RelIso(qualitymuons.at(i))<0.1){
			  GoodMuons.push_back(qualitymuons.at(i));
		  }
		continue;
	  }else if(Ntp->isFakeMuon(qualitymuons.at(i),posmu) && Ntp->isData() && MVA_ID){
		  fakeMuon = true;
		  GoodMuons.push_back(qualitymuons.at(i));
		  fakeRateMu = Fakerate(Ntp->Muons_p4(qualitymuons.at(i)),MuonFakeRate,"muon");
	  }
  }

  for(unsigned i=0;i<GoodMuons.size();i++){
    dist.at(NMuPt).push_back(Ntp->Muons_p4(GoodMuons.at(i)).Pt());
    if(Ntp->Muons_p4(GoodMuons.at(i)).Pt()<mu_pt){
      GoodMuons.erase(GoodMuons.begin()+i);
      i--;
    }
  }
  
  value.at(NMuPt)=GoodMuons.size();
  pass.at(NMuPt)=(value.at(NMuPt)>=cut.at(NMuPt));
  
  for(unsigned i=0;i<GoodMuons.size();i++){
    dist.at(NMuEta).push_back(Ntp->Muons_p4(GoodMuons.at(i)).Eta());
    if(fabs(Ntp->Muons_p4(GoodMuons.at(i)).Eta())>mu_eta){
      GoodMuons.erase(GoodMuons.begin()+i);
      i--;
    }
  }
  
  value.at(NMuEta)=GoodMuons.size();
  pass.at(NMuEta)=(value.at(NMuEta)>=cut.at(NMuEta));
  
  for(unsigned i=0;i<GoodMuons.size();i++){
	  dist.at(NMuRelIso).push_back(Ntp->Muon_RelIso(GoodMuons.at(i)));
	  if(Ntp->Muon_RelIso(GoodMuons.at(i))>mu_reliso){
		  GoodMuons.erase(GoodMuons.begin()+i);
		  i--;
	  }
  }
  
  value.at(NMuRelIso)=GoodMuons.size();
  pass.at(NMuRelIso)=(value.at(NMuRelIso)>=cut.at(NMuRelIso));
  
  value.at(NMu)=GoodMuons.size();
  pass.at(NMu)=(value.at(NMu)==cut.at(NMu));
  
  unsigned int muidx1(999),muInfoidx1(999);
  if(GoodMuons.size()>=1){muidx1=GoodMuons.at(0);muInfoidx1=muidx1;}
  if(verbose)std::cout << "void  ZtoEMu::doEvent() E " << muidx1 <<" "<< muInfoidx1 << std::endl;

  ///////////////////////////////////////////////
  //
  // E Cuts
  //
  if(verbose) std::cout << "electrons cuts" << std::endl;
  std::vector<unsigned int> GoodElectrons;
  bool fakeElectron = false;
  
  for(unsigned i=0;i<qualityelectrons.size();i++){
	  if(!Ntp->Electron_HasMatchedConversions(qualityelectrons.at(i)) &&
			  Ntp->Electron_numberOfMissedHits(qualityelectrons.at(i))==0 &&
			  dz(Ntp->Electron_p4(qualityelectrons.at(i)),Ntp->Electron_Poca(qualityelectrons.at(i)),Ntp->Vtx(pose))<0.2 &&
			  dxy(Ntp->Electron_p4(qualityelectrons.at(i)),Ntp->Electron_Poca(qualityelectrons.at(i)),Ntp->Vtx(pose))<0.02
			  ){
		  if(MVA_ID){
			  if(fabs(Ntp->Electron_supercluster_eta(qualityelectrons.at(i)))<0.8 && Ntp->Electron_RelIso(qualityelectrons.at(i))<0.15){
				  if(Ntp->Electron_p4(qualityelectrons.at(i)).Pt()<20 && Ntp->Electron_MVA_discriminator(qualityelectrons.at(i))>0.925){
					  GoodElectrons.push_back(qualityelectrons.at(i));
				  }else if(Ntp->Electron_p4(qualityelectrons.at(i)).Pt()>=20 && Ntp->Electron_MVA_discriminator(qualityelectrons.at(i))>0.905){
					  GoodElectrons.push_back(qualityelectrons.at(i));
				  }else if(Ntp->isFakeElectron(qualityelectrons.at(i),pose) && Ntp->isData()){
					  fakeElectron = true;
					  GoodElectrons.push_back(qualityelectrons.at(i));
					  fakeRateE = Fakerate(Ntp->Electron_p4(qualityelectrons.at(i)),ElectronFakeRate,"electron");
				  }
			  }else if(fabs(Ntp->Electron_supercluster_eta(qualityelectrons.at(i)))>=0.8 && fabs(Ntp->Electron_supercluster_eta(qualityelectrons.at(i)))<1.479 && Ntp->Electron_RelIso(qualityelectrons.at(i))<0.15){
				  if(Ntp->Electron_p4(qualityelectrons.at(i)).Pt()<20 && Ntp->Electron_MVA_discriminator(qualityelectrons.at(i))>0.915){
					  GoodElectrons.push_back(qualityelectrons.at(i));
				  }else if(Ntp->Electron_p4(qualityelectrons.at(i)).Pt()>=20 && Ntp->Electron_MVA_discriminator(qualityelectrons.at(i))>0.955){
					  GoodElectrons.push_back(qualityelectrons.at(i));
				  }else if(Ntp->isFakeElectron(qualityelectrons.at(i),pose) && Ntp->isData()){
					  fakeElectron = true;
					  GoodElectrons.push_back(qualityelectrons.at(i));
					  fakeRateE = Fakerate(Ntp->Electron_p4(qualityelectrons.at(i)),ElectronFakeRate,"electron");
				  }
			  }else if(fabs(Ntp->Electron_supercluster_eta(qualityelectrons.at(i)))>=1.479 && Ntp->Electron_RelIso(qualityelectrons.at(i))<0.1){
				  if(Ntp->Electron_p4(qualityelectrons.at(i)).Pt()<20 && Ntp->Electron_MVA_discriminator(qualityelectrons.at(i))>0.965){
					  GoodElectrons.push_back(qualityelectrons.at(i));
				  }else if(Ntp->Electron_p4(qualityelectrons.at(i)).Pt()>=20 && Ntp->Electron_MVA_discriminator(qualityelectrons.at(i))>0.975){
					  GoodElectrons.push_back(qualityelectrons.at(i));
				  }else if(Ntp->isFakeElectron(qualityelectrons.at(i),pose) && Ntp->isData()){
					  fakeElectron = true;
					  GoodElectrons.push_back(qualityelectrons.at(i));
					  fakeRateE = Fakerate(Ntp->Electron_p4(qualityelectrons.at(i)),ElectronFakeRate,"electron");
				  }
			  }
		  }else{
			  if(Ntp->isTightElectron(qualityelectrons.at(i),pose))GoodElectrons.push_back(qualityelectrons.at(i));
		  }
	  }
  }

  if(Ntp->Electron_NMVA()!=Ntp->NElectrons()){
	  std::cout << "#########" << std::endl;
	  std::cout << "nmva = " << Ntp->Electron_NMVA() << ", nele = " << Ntp->NElectrons() << std::endl; 
  }
  
  for(unsigned i=0;i<GoodElectrons.size();i++){
    dist.at(NEPt).push_back(Ntp->Electron_p4(GoodElectrons.at(i)).Pt());
    if(Ntp->Electron_p4(GoodElectrons.at(i)).Pt()<e_pt){
      GoodElectrons.erase(GoodElectrons.begin()+i);
      i--;
    }
  }
  value.at(NEPt)=GoodElectrons.size();
  pass.at(NEPt)=(value.at(NEPt)>=cut.at(NEPt));
  
  for(unsigned i=0;i<GoodElectrons.size();i++){
    dist.at(NEEta).push_back(Ntp->Electron_p4(GoodElectrons.at(i)).Eta());
    if(fabs(Ntp->Electron_supercluster_eta(GoodElectrons.at(i)))>e_eta){
      GoodElectrons.erase(GoodElectrons.begin()+i);
      i--;
    }
  }
  value.at(NEEta)=GoodElectrons.size();
  pass.at(NEEta)=(value.at(NEEta)>=cut.at(NEEta));
  
  value.at(NE)=GoodElectrons.size();
  pass.at(NE)=(value.at(NE)==cut.at(NE));
  
  unsigned int eidx1(999),eInfoidx1(999);
  if(GoodElectrons.size()>=1){eidx1=GoodElectrons.at(0);eInfoidx1=eidx1;}
  if(verbose)std::cout << "void  ZtoEMu::doEvent() E " << eidx1 <<" "<< eInfoidx1 << std::endl;
  
  ///////////////////////////////////////////////
  //
  // dR cleaning of e & mu
  //
  if(verbose)std::cout << "dR cleaning" << std::endl;
  
  value.at(drMuE)=10.;
  if(muidx1!=999 && eidx1!=999){
	  if(Ntp->Muon_Track_idx(GoodMuons.at(0))==Ntp->Electron_Track_idx(GoodElectrons.at(0))){
		  value.at(drMuE)=10.;
	  }
	  value.at(drMuE)=Ntp->Muons_p4(GoodMuons.at(0)).DeltaR(Ntp->Electron_p4(GoodElectrons.at(0)));
  }
  pass.at(drMuE)=(value.at(drMuE)>cut.at(drMuE));
  
  ///////////////////////////////////////////////
  //
  // Loose Muon Veto
  //
  if(verbose)std::cout << "Loose muon veto" << std::endl;
  int LooseMuons = 0;
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(Ntp->isLooseMuon(i)
			  && !Ntp->isTightMuon(i,posmu)
			  && dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(posmu))<0.02
			  && dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(posmu))<0.2
			  && Ntp->Muons_p4(i).Pt()>26.
			  ){
		  LooseMuons++;
	  }
  }
  value.at(looseMuonVeto)=LooseMuons;
  pass.at(looseMuonVeto)=(value.at(looseMuonVeto)==cut.at(looseMuonVeto));

  ///////////////////////////////////////////////
  //
  // charge cut
  //
  if(verbose) std::cout << "Charge cut" << std::endl;
  value.at(charge)=-5;
  if(eidx1!=999 && muidx1!=999){
	  value.at(charge)=Ntp->Electron_Charge(GoodElectrons.at(0))+Ntp->Muon_Charge(GoodMuons.at(0));
  }
  pass.at(charge)=(value.at(charge)==0);
  
  ///////////////////////////////////////////////
  //
  // jet veto
  //
  if(verbose)std::cout << "jet veto" << std::endl;
  if(verbose)std::cout << "Cleaning jets" << std::endl;
  std::vector<int> jetidx;
  bool etrackjet;
  bool mutrackjet;
  jetidx.clear();
  if(eidx1!=999 && muidx1!=999){
	  for(unsigned i=0;i<Ntp->NPFJets();i++){
		  etrackjet=false;
		  mutrackjet=false;
		  for(unsigned j=0;j<Ntp->PFJet_nTrk(i);j++){
			  if(Ntp->PFJet_TracksP4(i,j).DeltaR(Ntp->Muons_p4(GoodMuons.at(0)))<0.001){
				  mutrackjet=true;
			  }
			  if(Ntp->PFJet_TracksP4(i,j).DeltaR(Ntp->Electron_p4(GoodElectrons.at(0)))<0.1){
				  etrackjet=true;
			  }
		  }
		  if(!mutrackjet && !etrackjet){
			  jetidx.push_back(i);
		  }
	  }
  }
  
  if(verbose)std::cout << "Looking for jets from Vtx" << std::endl;
  int leadingtrack;
  bool jetfromvtx;
  std::vector<int> jetsfromvtx;
  std::vector<int> leadingtracks;
  jetsfromvtx.clear();
  leadingtracks.clear();
  for(unsigned i=0;i<jetidx.size();i++){
	  leadingtrack = 0;
	  jetfromvtx = false;
	  int counter = 0;
	  for(unsigned j=0;j<Ntp->PFJet_nTrk(jetidx.at(i));j++){
		  if(Ntp->PFJet_TracksP4(jetidx.at(i),j).Pt()>Ntp->PFJet_TracksP4(jetidx.at(i),leadingtrack).Pt()){
			  leadingtrack=j;
		  }
	  }
	  if(Ntp->PFJet_nTrk(jetidx.at(i))==0)leadingtrack = -1;
	  for(unsigned j=0;j<Ntp->Vtx_nTrk(pos);j++){
		  if(leadingtrack>=0 && Ntp->PFJet_TracksP4(jetidx.at(i),leadingtrack).DeltaR(Ntp->Vtx_TracksP4(pos,j))<0.00001){
			  jetfromvtx=true;
			  counter++;
		  }
	  }
	  if(counter>1)std::cout << "More than one vtx track associated to leading track from jet" << std::endl;
	  if(jetfromvtx){
		  jetsfromvtx.push_back(jetidx.at(i));
		  leadingtracks.push_back(leadingtrack);
	  }
  }
  
  if(verbose)std::cout<< "Find two highest pt jets" << std::endl;
  int firstjet_idx=-1;
  int secondjet_idx=-1;
  double initialpt=0.;
  for(unsigned i=0;i<jetsfromvtx.size();i++){
	  if(Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()>initialpt){
		  initialpt=Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt();
		  firstjet_idx=jetsfromvtx.at(i);
	  }
  }
  initialpt=0.;
  for(unsigned i=0;i<jetsfromvtx.size();i++){
	  if(jetsfromvtx.size()>1 && firstjet_idx!=-1 && Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()>initialpt && Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()<Ntp->PFJet_p4(firstjet_idx).Pt()){
		  initialpt=Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt();
		  secondjet_idx=jetsfromvtx.at(i);
	  }
  }
  
  if(verbose)std::cout << "applying veto" << std::endl;
  value.at(jetVeto)=0;
  if(jetsfromvtx.size()>1){
	  value.at(jetVeto)=Ntp->PFJet_p4(firstjet_idx).Pt()+Ntp->PFJet_p4(secondjet_idx).Pt();
  }
  pass.at(jetVeto)=(value.at(jetVeto)<cut.at(jetVeto));
  
  ///////////////////////////////////////////////
  //
  // Mt Mu cut
  //
  if(verbose) std::cout << "Mt Mu cut" << std::endl;
  value.at(MtMu)=999;
  if(muidx1!=999){
	  value.at(MtMu)=sqrt(2*Ntp->Muons_p4(GoodMuons.at(0)).Pt()*Ntp->MET_Corr_et()*(1-cosphi2d(Ntp->Muons_p4(GoodMuons.at(0)).Px(),Ntp->Muons_p4(GoodMuons.at(0)).Py(),Ntp->MET_Corr_ex(),Ntp->MET_Corr_ey())));
  }
  pass.at(MtMu)=(value.at(MtMu)<cut.at(MtMu));
  
  ///////////////////////////////////////////////
  //
  // Pt balance cut
  //
  if(verbose) std::cout << "pt balance cut" << std::endl;
  value.at(ptBalance)=999.;
  for(unsigned i=0;i<GoodMuons.size();i++){
	  for(unsigned j=0;j<GoodElectrons.size();j++){
		  value.at(ptBalance) = (Ntp->Muons_p4(GoodMuons.at(i))+Ntp->Electron_p4(GoodElectrons.at(j))).Pt();
	  }
  }
  pass.at(ptBalance)=(value.at(ptBalance)<cut.at(ptBalance));
  
  ///////////////////////////////////////////////
  //
  // Invariant mass cut
  //
  if(verbose) std::cout << "Invariant mass cut" << std::endl;
  value.at(ZMassmax)=zmax+1.;
  value.at(ZMassmin)=zmin-1.;
  if(eidx1!=999 && muidx1!=999){
	  value.at(ZMassmax)=(Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M();
	  value.at(ZMassmin)=(Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M();
  }
  pass.at(ZMassmax)=(value.at(ZMassmax)<cut.at(ZMassmax));
  pass.at(ZMassmin)=(value.at(ZMassmin)>cut.at(ZMassmin));
  
  ////////////////////////////////////////////////
  //
  // QCD
  //
  if(MVA_ID){
	  fakeRate = 0.;
	  if(!pass.at(charge)
		  && Ntp->isData()
		  ){
		  if(fakeMuon && !fakeElectron){
			  fakeRate = fakeRateMu;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge) = true;
			  blind = 1;
		  }else if(fakeElectron && !fakeMuon){
			  fakeRate = fakeRateE;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge) = true;
			  blind = 1;
		  }else{
			  fakeRate = 1.;
		  }
	  }else if(pass.at(charge)
		  && Ntp->isData()
		  && !fakeMuon
		  && !fakeElectron
		  ){
		  fakeRate = 1.;
	  }
  }else{
	fakeRate = 1.;
  }
  /*if(!pass.at(charge)){
    if(Ntp->isData()){
      if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
      pass.at(charge)=true;
    }
  }*/
  
  //////////////////////////////////////////////////////////
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1);
  if(!Ntp->isData()){
    w*=Ntp->EvtWeight3D()*blind;
    if(verbose)std::cout << "void  ZtoEMu::doEvent() k" << w << " " << wobs << std::endl;
  }
  else{w=1*fakeRate*blind;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose) std::cout << "add plots" << std::endl;
  blindingactive = false;
  if(!blindingactive){
  if(status){
    if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);
  }
  if(verbose)std::cout << "filling NJets" << std::endl;
  
  NJets.at(t).Fill(Ntp->NPFJets(),w);
  
  if(jetsfromvtx.size()>1){
	  jetsum.at(t).Fill(Ntp->PFJet_p4(firstjet_idx).Pt()+Ntp->PFJet_p4(secondjet_idx).Pt(),w);
  }
  
  if(jetsfromvtx.size()>0)FirstJetPt.at(t).Fill(Ntp->PFJet_p4(firstjet_idx).Pt(),w);
  if(jetsfromvtx.size()>1)SecondJetPt.at(t).Fill(Ntp->PFJet_p4(secondjet_idx).Pt(),w);
  if(jetsfromvtx.size()>2)ThirdJetPt.at(t).Fill(Ntp->PFJet_p4(jetsfromvtx.at(2)).Pt(),w);
  if(jetsfromvtx.size()>3)FourthJetPt.at(t).Fill(Ntp->PFJet_p4(jetsfromvtx.at(3)).Pt(),w);
  
  
  if(verbose)std::cout << "looping over good electrons" << std::endl;
  for(unsigned int i=0;i<GoodElectrons.size();i++){
		EPt.at(t).Fill(Ntp->Electron_p4(GoodElectrons.at(i)).Pt(),w);
		etaE.at(t).Fill(Ntp->Electron_p4(GoodElectrons.at(i)).Eta(),w);
		mtE.at(t).Fill(sqrt(2*Ntp->Electron_p4(GoodElectrons.at(i)).Pt()*Ntp->MET_Corr_et()*(1-cosphi2d(Ntp->Electron_p4(GoodElectrons.at(i)).Px(),Ntp->Electron_p4(GoodElectrons.at(i)).Py(),Ntp->MET_Corr_ex(),Ntp->MET_Corr_ey()))),w);
		RelIsoE.at(t).Fill(Ntp->Electron_RelIso(GoodElectrons.at(i)),w);
  }
  if(verbose)std::cout << "looping over good muons" << std::endl;	
  for(unsigned int i=0;i<GoodMuons.size();i++){
	  MuPt.at(t).Fill(Ntp->Muons_p4(GoodMuons.at(i)).Pt(),w);
	  etaMu.at(t).Fill(Ntp->Muons_p4(GoodMuons.at(i)).Eta(),w);
	  mtMu.at(t).Fill(sqrt(2*Ntp->Muons_p4(GoodMuons.at(i)).Pt()*Ntp->MET_Corr_et()*(1-cosphi2d(Ntp->Muons_p4(GoodMuons.at(i)).Px(),Ntp->Muons_p4(GoodMuons.at(i)).Py(),Ntp->MET_Corr_ex(),Ntp->MET_Corr_ey()))),w);
	  RelIsoMu.at(t).Fill(Ntp->Muon_RelIso(GoodMuons.at(i)),w);
  }
  
  if(verbose)std::cout << "Calculating pzeta & poca difference" << std::endl;
  if(verbose)std::cout << "Calculating charges and dr(e,mu)" << std::endl;
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  for(unsigned j=0;j<Ntp->NElectrons();j++){
		  double dp = Ntp->Muons_p4(i).DeltaPhi(Ntp->Electron_p4(j));
		  if(dp<0)dp+=2*TMath::Pi();
		  deltaphi.at(t).Fill(dp,w);
		  drmue.at(t).Fill(Ntp->Muons_p4(i).DeltaR(Ntp->Electron_p4(j)),w);
		  chargesum.at(t).Fill(fabs(Ntp->Muon_Charge(i)+Ntp->Electron_Charge(j)),w);
		  chargesumsigned.at(t).Fill(Ntp->Muon_Charge(i)+Ntp->Electron_Charge(j),w);
	  }
  }

  if(pass.at(qualitycuts)){
	  // projections from 2d plot invmass vs deltaphi
	  double mass = (Ntp->Muons_p4(qualitymuons.at(0))+Ntp->Electron_p4(qualityelectrons.at(0))).M();
	  double qualityphi = Ntp->Muons_p4(qualitymuons.at(0)).DeltaPhi(Ntp->Electron_p4(qualityelectrons.at(0)));
	  if(qualityphi<0)qualityphi+=2*TMath::Pi();
	  if(status){
		  if(mass<79 && qualityphi<2.5){
			  mproj1.at(t).Fill(mass,w);
			  phiproj1.at(t).Fill(qualityphi,w);
		  }else if(mass>=79 && mass<=104 && qualityphi<2.5){
			  mproj2.at(t).Fill(mass,w);
			  phiproj2.at(t).Fill(qualityphi,w);
		  }else if(mass>104 && qualityphi<2.5){
			  mproj3.at(t).Fill(mass,w);
			  phiproj3.at(t).Fill(qualityphi,w);
		  }else if(mass<79 && qualityphi>=2.5 && qualityphi<=3.7){
			  mproj4.at(t).Fill(mass,w);
			  phiproj4.at(t).Fill(qualityphi,w);
		  }else if(mass>=79 && mass<=104 && qualityphi>=2.5 && qualityphi<=3.7){
			  mproj5.at(t).Fill(mass,w);
			  phiproj5.at(t).Fill(qualityphi,w);
		  }else if(mass>104 && qualityphi>=2.5 && qualityphi<=3.7){
			  mproj6.at(t).Fill(mass,w);
			  phiproj6.at(t).Fill(qualityphi,w);
		  }else if(mass<79 && qualityphi>3.7){
			  mproj7.at(t).Fill(mass,w);
			  phiproj7.at(t).Fill(qualityphi,w);
		  }else if(mass>=79 && mass<=104 && qualityphi>3.7){
			  mproj8.at(t).Fill(mass,w);
			  phiproj8.at(t).Fill(qualityphi,w);
		  }else if(mass>104 && qualityphi>3.7){
			  mproj9.at(t).Fill(mass,w);
			  phiproj9.at(t).Fill(qualityphi,w);
		  }
	  }
  }

  for(unsigned i=0;i<GoodMuons.size();i++){
	  for(unsigned j=0;j<GoodElectrons.size();j++){
		  pzeta.at(t).Fill(calculatePzeta(i,j,GoodElectrons,GoodMuons),w);
		  pzetaDQM.at(t).Fill(calculatePzetaDQM(i,j,GoodElectrons,GoodMuons),w);
		  ptbal.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(i))+Ntp->Electron_p4(GoodElectrons.at(j))).Pt(),w);
		  double dpID = Ntp->Muons_p4(GoodMuons.at(i)).DeltaPhi(Ntp->Electron_p4(GoodElectrons.at(j)));
		  if(dpID<0)dpID+=2*TMath::Pi();
		  deltaphiID.at(t).Fill(dpID,w);
		  if(Ntp->Muons_p4(GoodMuons.at(i)).DeltaR(Ntp->Electron_p4(GoodElectrons.at(j)))>0.2){
			if(pass.at(qualitycuts))InvmassVsDeltaPhi.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(i))+Ntp->Electron_p4(GoodElectrons.at(j))).M(),dpID,w);
			if(pass.at(qualitycuts))invmass.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(i))+Ntp->Electron_p4(GoodElectrons.at(j))).M(),w);
			chargesumID.at(t).Fill(fabs(Ntp->Muon_Charge(GoodMuons.at(i))+Ntp->Electron_Charge(GoodElectrons.at(j))),w);
			chargesumIDsigned.at(t).Fill(Ntp->Muon_Charge(GoodMuons.at(i))+Ntp->Electron_Charge(GoodElectrons.at(j)),w);
		  }
		  drmueID.at(t).Fill(Ntp->Muons_p4(GoodMuons.at(i)).DeltaR(Ntp->Electron_p4(GoodElectrons.at(j))),w);
	  }
  }
  
  if(verbose)std::cout << "invariant mass with different cuts applied" << std::endl;
  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(qualitycuts)
		  && pass.at(SameVtx)
		  && pass.at(NMuPt)
		  && pass.at(NMuEta)
		  && pass.at(NEPt)
		  && pass.at(NEEta)
		  && pass.at(NMu)
		  && pass.at(NE)
		  ){
	  invmass_only_object_id.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
	  if(pass.at(drMuE)){
		  invmass_dremu.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  //if(pass.at(looseMuonVeto)){
			  //invmass_loosemuonveto.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
			  if(pass.at(charge)){
				  invmass_charge.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
				  if(pass.at(jetVeto)){
					  invmass_jetveto.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
					  if(pass.at(MtMu)){
						  invmass_mtmu.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
						  if(pass.at(ptBalance)){
							  invmass_ptbalance.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
							  invmass_no_loosemuonveto.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
							  if(pass.at(ZMassmax)
									  && pass.at(ZMassmin)){
								  invmass_zmass.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
							  }
						  }
					  }
				  }
			  }
		  //}
		  if(pass.at(looseMuonVeto))invmass_loosemuonveto.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  if(pass.at(charge)
				  && pass.at(jetVeto)
				  && pass.at(MtMu)
				  ){
			  nm2_ptbalance.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(charge)
				  && pass.at(jetVeto)
				  && pass.at(ptBalance)
				  ){
			  nm2_mtmu.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(charge)
				  && pass.at(MtMu)
				  && pass.at(ptBalance)
				  ){
			  nm2_jetveto.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(jetVeto)
				  && pass.at(MtMu)
				  && pass.at(ptBalance)
				  ){
			  nm2_charge.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
	  }
	  if(pass.at(charge)
			  && pass.at(jetVeto)
			  && pass.at(MtMu)
			  && pass.at(ptBalance)
			  ){
		  nm2_drmue.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
	  }
  }

  }

  if(verbose)std::cout << "ZtoEMu::doEvent() doEvent END" << std::endl;
}

double ZtoEMu::calculatePzeta(int muiterator, int eiterator,std::vector<unsigned int> vec1,std::vector<unsigned int> vec2){
  pex=Ntp->Electron_p4(vec1.at(eiterator)).Px();
  pey=Ntp->Electron_p4(vec1.at(eiterator)).Py();
  pmux=Ntp->Muons_p4(vec2.at(muiterator)).Px();
  pmuy=Ntp->Muons_p4(vec2.at(muiterator)).Py();
  phie=Ntp->Electron_p4(vec1.at(eiterator)).Phi();
  phimu=Ntp->Muons_p4(vec2.at(muiterator)).Phi();
  combpt=TMath::Sqrt(pow(pex+pmux,2)+pow(pey+pmuy,2));
  aemu=TMath::ACos(pmux*pex+pmuy*pey/(Ntp->Muons_p4(vec2.at(muiterator)).Pt()*Ntp->Electron_p4(vec1.at(eiterator)).Pt()));
  if(phie<phimu && fabs(phie-phimu)<TMath::Pi())phi1=phie;
  else if(phimu<phie && fabs(phie-phimu)>TMath::Pi())phi1=phie;
  else if(phie<phimu && fabs(phie-phimu)>TMath::Pi())phi1=phimu;
  else if(phimu<phie && fabs(phie-phimu)<TMath::Pi())phi1=phimu;
  beta=TMath::ACos(((pex+pmux)*TMath::Cos(phi1+0.5*aemu)+(pey+pmuy)*TMath::Sin(phi1+0.5*aemu))/combpt);
  gamma=TMath::ACos((Ntp->MET_Corr_ex()*TMath::Cos(phi1+0.5*aemu)+Ntp->MET_Corr_ey()*TMath::Sin(phi1+0.5*aemu))/Ntp->MET_Corr_et());
  if(Ntp->MET_Corr_phi()>(phi1+0.5*aemu+0.5*TMath::Pi()) && Ntp->MET_Corr_phi()<(phi1+0.5*aemu+1.5*TMath::Pi()))gamma*=-1;
  pvis=TMath::Sin(beta)*combpt;
  pmiss=TMath::Sin(gamma)*Ntp->MET_Corr_et();
  return pmiss-pvis;
}

double ZtoEMu::calculatePzetaDQM(int muiterator, int eiterator,std::vector<unsigned int> vec1,std::vector<unsigned int> vec2){
	double cosPhi1 = TMath::Cos(Ntp->Electron_p4(vec1.at(eiterator)).Phi());
	double sinPhi1 = TMath::Sin(Ntp->Electron_p4(vec1.at(eiterator)).Phi());
	double cosPhi2 = TMath::Cos(Ntp->Muons_p4(vec2.at(muiterator)).Phi());
	double sinPhi2 = TMath::Sin(Ntp->Muons_p4(vec2.at(muiterator)).Phi());
	double zetaX = cosPhi1 + cosPhi2;
	double zetaY = sinPhi1 + sinPhi2;
	double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
	if(zetaR>0.){
		zetaX/=zetaR;
		zetaY/=zetaR;
	}
	double pxVis=Ntp->Electron_p4(vec1.at(eiterator)).Px()+Ntp->Muons_p4(vec2.at(muiterator)).Px();
	double pyVis=Ntp->Electron_p4(vec1.at(eiterator)).Py()+Ntp->Muons_p4(vec2.at(muiterator)).Py();
	double pZetaVis=pxVis*zetaX+pyVis*zetaY;
	double px=pxVis+Ntp->MET_Corr_ex();
	double py=pyVis+Ntp->MET_Corr_ey();
	double pZeta=px*zetaX+py*zetaY;
	return pZeta-1.5*pZetaVis;
}

double ZtoEMu::cosphi2d(double px1, double py1, double px2, double py2){
	return (px1*px2+py1*py2)/(sqrt(pow(px1,2)+pow(py1,2))*sqrt(pow(px2,2)+pow(py2,2)));
}

double ZtoEMu::cosphi3d(TVector3 vec1, TVector3 vec2){
	return (vec1.Dot(vec2))/vec1.Mag()/vec2.Mag();
}

double ZtoEMu::jetdxy(int vtx_idx, int leadingtrack_idx, int jet_idx){
	return fabs((-(Ntp->Track_Poca(leadingtrack_idx).X()-Ntp->Vtx(vtx_idx).X())*Ntp->PFJet_p4(jet_idx).Py()+(Ntp->Track_Poca(leadingtrack_idx).Y()-Ntp->Vtx(vtx_idx).Y())*Ntp->PFJet_p4(jet_idx).Px())/Ntp->PFJet_p4(jet_idx).Pt());
}

double ZtoEMu::jetdz(int vtx_idx, int leadingtrack_idx, int jet_idx){
	return fabs(Ntp->Track_Poca(leadingtrack_idx).Z()-Ntp->Vtx(vtx_idx).Z()-((Ntp->Track_Poca(leadingtrack_idx).X()-Ntp->Vtx(vtx_idx).X())*Ntp->PFJet_p4(jet_idx).Px()+(Ntp->Track_Poca(leadingtrack_idx).Y()-Ntp->Vtx(vtx_idx).Y())*Ntp->PFJet_p4(jet_idx).Py())*Ntp->PFJet_p4(jet_idx).Pz()/pow(Ntp->PFJet_p4(jet_idx).Pt(),2));
}

bool ZtoEMu::jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx){
	for(unsigned i=0;i<vtx_track_idx.size();i++){
		if(vtx_track_idx.at(i)==leadingtrack_idx)return true;
	}
	return false;
}

double ZtoEMu::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs((-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt());
}

double ZtoEMu::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2));
}

double ZtoEMu::Fakerate(TLorentzVector vec, TH2D *fakeRateHist, std::string type){
	
	double eta1, eta2;
	int ptbin = 0;
	int etabin = 0;
	double FakePt = vec.Pt();
	double FakeEta = vec.Eta();
	double fakerate = 0.;
	
	if(type=="muon"){
		eta1 = 2.1;
		eta2 = 1.2;
	}else if(type=="electron"){
		eta1 = 2.;
		eta2 = 1.479;
	}
	
	if(FakePt<15.)ptbin = 1;
	else if(FakePt>=15. && FakePt<20.)ptbin = 2;
	else if(FakePt>=20. && FakePt<25.)ptbin = 3;
	else if(FakePt>=25. && FakePt<30.)ptbin = 4;
	else if(FakePt>=30.)ptbin = 5;
	
	if(FakeEta<-eta1)etabin = 1;
	else if(FakeEta>=-eta1 && FakeEta<-eta2)etabin = 2;
	else if(FakeEta>=-eta2 && FakeEta<-0.8)etabin = 3;
	else if(FakeEta>=-0.8 && FakeEta<0.)etabin = 4;
	else if(FakeEta>=0. && FakeEta<0.8)etabin = 5;
	else if(FakeEta>=0.8 && FakeEta<eta2)etabin = 6;
	else if(FakeEta>=eta2 && FakeEta<eta1)etabin = 7;
	else if(FakeEta>=eta1)etabin = 8;
	
	if(ptbin==0 || etabin==0){
		fakerate = 0;
	}else{
		fakerate = fakeRateHist->GetBinContent(ptbin,etabin);
	}
	
	return fakerate;
}

// the number 1129 is arbitrary. NEEDS TO BE FIXED!!!
void ZtoEMu::Finish(){
	unsigned int k;
	if(Nminus0.at(0).at(k).Integral()!=0)if(HConfig.GetHisto(false,DataMCType::QCD,k))ScaleAllHistOfType(k,1./Nminus0.at(0).at(k).Integral());
	Selection::Finish();
}