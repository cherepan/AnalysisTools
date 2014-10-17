
#include "FitTest.h"
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

FitTest::FitTest(TString Name_, TString id_):
  Selection(Name_,id_)
  ,muoneta(2.1)
  ,jeteta(2.3)
  ,TauTrackPtThreshold(5.0)
		//  ,ran()
{
  verbose=false;
}

FitTest::~FitTest(){
  for(int j=0; j<Npassed.size(); j++){
    //std::cout << "FitTest::~FitTest Selection Summary before: "
	// << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 //<< Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  //std::cout << "FitTest::~FitTest()" << std::endl;
}

void  FitTest::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==hasTag)	      cut.at(hasTag)=1;
    if(i==TagPtMin)           cut.at(TagPtMin)=20;
    if(i==TagIso)             cut.at(TagIso)=0.2;
//    if(i==numIsoTags)         cut.at(numIsoTags)=1;
    if(i==TauPt)              cut.at(TauPt)=20;
    if(i==TauEta)             cut.at(TauEta)=2.0;
    if(i==TauIsIsolated)      cut.at(TauIsIsolated)=true;
    if(i==TauFit)             cut.at(TauFit)=1;
    if(i==deltaPhi)           cut.at(deltaPhi)=TMath::Pi()*3.0/4.0;
    if(i==EventId)            cut.at(EventId)=1;
    if(i==Charge)             cut.at(Charge)=0.0;
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
    if(i==EventId){
      title.at(i)="EventId ";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#mu} passing the Trigger";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_EventId_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_EventId_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }        




    else if(i==hasTag){
      title.at(i)="$Number of good muons$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of good muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TagPtMin){
      title.at(i)="$N_{\\mu} with P_{T}>$";
      title.at(i)+=cut.at(TagPtMin);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#mu}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagPtMin_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagPtMin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TagIso){
      title.at(i)="$N_{\\mu} with rel. isolation <=$";
      title.at(i)+=cut.at(TagIso);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#mu}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagIso_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagIso_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TauPt){
      title.at(i)="$N_{\\tau} with P_{T} >$";
      title.at(i)+=cut.at(TauPt);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TauEta){
      title.at(i)="$N_{\\tau} | \\eta <$";
      title.at(i)+=cut.at(TauEta);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TauIsIsolated){
      title.at(i)="$N_{\\tau} | isolated$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsolated_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsolated_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TauFit){
      title.at(i)="$N_{\\tau} | fitted good$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauFit_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauFit_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==deltaPhi){
      title.at(i)="$\\Delta \\phi distribution >$";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(deltaPhi));
      title.at(i)+=buffer;
      title.at(i)+="(rad)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#Delta #phi (rad)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaPhi_",htitle,32,0,TMath::Pi(),hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaPhi_",htitle,32,0,TMath::Pi(),hlabel,"Events"));
    }
    else if(i==Charge){
      title.at(i)="$|Q_{\\tau}+Q_{\\mu}|-0.5<$";
      title.at(i)+=cut.at(Charge);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="|Q_{#tau}+Q_{#mu}|";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Charge_",htitle,44,0,2.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Charge_",htitle,44,0,2.2,hlabel,"Events"));
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
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  TauFlightLength=HConfig.GetTH1D(Name+"_TauFlightLength","TauFlightLength",62,-0.55,2.55,"L (cm)","Events");
  TauFlightLengthTransverse=HConfig.GetTH1D(Name+"_TauFlightLengthTransverse","TauFlightLengthTransverse",62,-0.55,2.55,"L_{T} (cm)","Events");
  TauMomentum=HConfig.GetTH1D(Name+"_TauMomentum","TauMomentum",50,0,400,"|P| (GeV)","Events");
  TauMomentumTransverse=HConfig.GetTH1D(Name+"_TauMomentumTransverse","TauMomentum",50,0,400,"|P_{T}| (GeV)","Events");
  TauLife=HConfig.GetTH1D(Name+"_TauLife","TauLife",1000,-0.5,2,"#tau_{#tau} (ps)","Events");
  TauLifeTransverse=HConfig.GetTH1D(Name+"_TauLifeTransverse","TauLifeTransverse",1000,-0.5,2,"transverse #tau_{#tau} (ps)","Events");
  ResTauFlightLength=HConfig.GetTH1D(Name+"_ResTauFlightLength","ResTauFlightLength",100,-2,2,"L (cm)","Events");
  ResTauFlightLengthTransverse=HConfig.GetTH1D(Name+"_ResTauFlightLengthTransverse","ResTauFlightLengthTransverse",100,-2,2,"L_{T} (cm)","Events");
  ResTauMomentum=HConfig.GetTH1D(Name+"_ResTauMomentum","ResTauMomentum",80,-200,200,"|P| (GeV)","Events");
  ResTauMomentumTransverse=HConfig.GetTH1D(Name+"_ResTauMomentumTransverse","ResTauMomentum",80,-200,200,"|P_{T}| (GeV)","Events");
  ResTauLife=HConfig.GetTH1D(Name+"_ResTauLife","ResTauLife",1000,-1,1,"#tau_{#tau} (ps)","Events");
  ResTauLifeTransverse=HConfig.GetTH1D(Name+"_ResTauLifeTransverse","ResTauLifeTransverse",1000,-1,1,"transverse #tau_{#tau} (ps)","Events");

  //
  TauMomentumAvg=HConfig.GetTH1D(Name+"_TauMomentumAvg","TauMomentumAvg",50,0,400,"|P| (GeV)","Events");
  TauMomentumTransverseAvg=HConfig.GetTH1D(Name+"_TauMomentumTransverseAvg","TauMomentumAvg",50,0,400,"|P_{T}| (GeV)","Events");
  TauLifeAvg=HConfig.GetTH1D(Name+"_TauLifeAvg","TauLifeAvg",1000,-0.5,2,"#tau_{#tau} (ps)","Events");
  TauLifeTransverseAvg=HConfig.GetTH1D(Name+"_TauLifeTransverseAvg","TauLifeTransverseAvg",1000,-0.5,2,"transverse #tau_{#tau} (ps)","Events");
  ResTauMomentumAvg=HConfig.GetTH1D(Name+"_ResTauMomentumAvg","ResTauMomentumAvg",80,-200,200,"|P| (GeV)","Events");
  ResTauMomentumTransverseAvg=HConfig.GetTH1D(Name+"_ResTauMomentumTransverseAvg","ResTauMomentumAvg",80,-200,200,"|P_{T}| (GeV)","Events");
  ResTauLifeAvg=HConfig.GetTH1D(Name+"_ResTauLifeAvg","ResTauLifeAvg",1000,-1,1,"#tau_{#tau} (ps)","Events");
  ResTauLifeTransverseAvg=HConfig.GetTH1D(Name+"_ResTauLifeTransverseAvg","ResTauLifeTransverseAvg",1000,-1,1,"transverse #tau_{#tau} (ps)","Events");

  ZMass0=HConfig.GetTH1D(Name+"_ZMass0","ZMass0",80,85,95,"M_{Z}, (GeV)","Events");
  ZMass1=HConfig.GetTH1D(Name+"_ZMass1","ZMass1",80,85,95,"M_{Z}, (GeV)","Events");
  ZMass2=HConfig.GetTH1D(Name+"_ZMass2","ZMass2",80,85,95,"M_{Z}, (GeV)","Events");

  TauA1PtResol0=HConfig.GetTH1D(Name+"_TauA1PtResol0","TauA1PtResol0",80,-200,200,"#Delta Pt (GeV)","Events");
  TauA1PtResol1=HConfig.GetTH1D(Name+"_TauA1PtResol1","TauA1PtResol1",80,-200,200,"#Delta Pt (GeV)","Events");
  TauA1PtResol2=HConfig.GetTH1D(Name+"_TauA1PtResol2","TauA1PtResol2",80,-200,200,"#Delta Pt (GeV)","Events");

  TauMuPtResol0=HConfig.GetTH1D(Name+"_TauMuPtResol0","TauMuPtResol0",80,-200,200,"#Delta Pt (GeV)","Events");
  TauMuPtResol1=HConfig.GetTH1D(Name+"_TauMuPtResol1","TauMuPtResol1",80,-200,200,"#Delta Pt (GeV)","Events");
  TauMuPtResol2=HConfig.GetTH1D(Name+"_TauMuPtResol2","TauMuPtResol2",80,-200,200,"#Delta Pt (GeV)","Events");


  TauA1PhiResol0=HConfig.GetTH1D(Name+"_TauA1PhiResol0","TauA1PhiResol0",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauA1PhiResol1=HConfig.GetTH1D(Name+"_TauA1PhiResol1","TauA1PhiResol1",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauA1PhiResol2=HConfig.GetTH1D(Name+"_TauA1PhiResol2","TauA1PhiResol2",40,-1.5,1.5,"#Delta #phi (rad)","Events");

  TauMuPhiResol0=HConfig.GetTH1D(Name+"_TauMuPhiResol0","TauMuPhiResol0",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauMuPhiResol1=HConfig.GetTH1D(Name+"_TauMuPhiResol1","TauMuPhiResol1",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauMuPhiResol2=HConfig.GetTH1D(Name+"_TauMuPhiResol2","TauMuPhiResol2",40,-1.5,1.5,"#Delta #phi (rad)","Events");


  TauA1EtaResol0=HConfig.GetTH1D(Name+"_TauA1EtaResol0","TauA1EtaResol0",40,-2.5,2.5,"#Delta #eta","Events");
  TauA1EtaResol1=HConfig.GetTH1D(Name+"_TauA1EtaResol1","TauA1EtaResol1",40,-2.5,2.5,"#Delta #eta","Events");
  TauA1EtaResol2=HConfig.GetTH1D(Name+"_TauA1EtaResol2","TauA1EtaResol2",40,-2.5,2.5,"#Delta #eta","Events");

  TauMuEtaResol0=HConfig.GetTH1D(Name+"_TauMuEtaResol0","TauMuEtaResol0",40,-2.5,2.5,"#Delta #eta","Events");
  TauMuEtaResol1=HConfig.GetTH1D(Name+"_TauMuEtaResol1","TauMuEtaResol1",40,-2.5,2.5,"#Delta #eta","Events");
  TauMuEtaResol2=HConfig.GetTH1D(Name+"_TauMuEtaResol2","TauMuEtaResol2",40,-2.5,2.5,"#Delta #eta","Events");


  Niter0=HConfig.GetTH1D(Name+"_Niter0","Niter0",100,0,100,"Niter0","Events");
  Niter1=HConfig.GetTH1D(Name+"_Niter1","Niter1",100,0,100,"Niter1","Events");
  Niter2=HConfig.GetTH1D(Name+"_Niter2","Niter2",100,0,100,"Niter2","Events");

  TauMuEnergyPull0=HConfig.GetTH1D(Name+"_TauMuEnergyPull0","TauMuEnergyPull0",80,-40,40,"","Events");
  TauA1EnergyPull0=HConfig.GetTH1D(Name+"_TauA1EnergyPull0","TauA1EnergyPull0",80,-40,40,"","Events");

  TauMuEnergyPull1=HConfig.GetTH1D(Name+"_TauMuEnergyPull1","TauMuEnergyPull1",80,-40,40,"","Events");
  TauA1EnergyPull1=HConfig.GetTH1D(Name+"_TauA1EnergyPull1","TauA1EnergyPull1",80,-40,40,"","Events");
 
  TauMuEnergyPull2=HConfig.GetTH1D(Name+"_TauMuEnergyPull2","TauMuEnergyPull2",80,-40,40,"","Events");
  TauA1EnergyPull2=HConfig.GetTH1D(Name+"_TauA1EnergyPull2","TauA1EnergyPull2",80,-40,40,"","Events");

  Probability0=HConfig.GetTH1D(Name+"_Probability0","Probability0",20,0,1,"Probability","Events");
  Probability1=HConfig.GetTH1D(Name+"_Probability1","Probability1",20,0,1,"Probability","Events");
  Probability2=HConfig.GetTH1D(Name+"_Probability2","Probability2",20,0,1,"Probability","Events");




  TruthA1Pt=HConfig.GetTH1D(Name+"_TruthA1Pt","TruthA1Pt",50,10,60,"","Events");
  TruthA1Phi=HConfig.GetTH1D(Name+"_TruthA1Phi","TruthA1Phi",30,-3.14,3.14,"","Events");
  TruthA1Eta=HConfig.GetTH1D(Name+"_TruthA1Eta","TruthA1Eta",20,-2,2,"","Events");

  TruthMuPt=HConfig.GetTH1D(Name+"_TruthMuPt","TruthMuPt",50,10,60,"","Events");
  TruthMuPhi=HConfig.GetTH1D(Name+"_TruthMuPhi","TruthMuPhi",30,-3.14,3.14,"","Events");
  TruthMuEta=HConfig.GetTH1D(Name+"_TruthMuEta","TruthMuEta",20,-2,2,"","Events");


  TruthA1PtAfter=HConfig.GetTH1D(Name+"_TruthA1PtAfter","TruthA1PtAfter",50,10,60,"","Events");
  TruthA1PhiAfter=HConfig.GetTH1D(Name+"_TruthA1PhiAfter","TruthA1PhiAfter",30,-3.14,3.14,"","Events");
  TruthA1EtaAfter=HConfig.GetTH1D(Name+"_TruthA1EtaAfter","TruthA1EtaAfter",20,-2,2,"","Events");

  TruthMuPtAfter=HConfig.GetTH1D(Name+"_TruthMuPtAfter","TruthMuPtAfter",50,10,60,"","Events");
  TruthMuPhiAfter=HConfig.GetTH1D(Name+"_TruthMuPhiAfter","TruthMuPhiAfter",30,-3.14,3.14,"","Events");
  TruthMuEtaAfter=HConfig.GetTH1D(Name+"_TruthMuEtaAfter","TruthMuEtaAfter",20,-2,2,"","Events");



  TruthA1PtAfterQC=HConfig.GetTH1D(Name+"_TruthA1PtAfterQC","TruthA1PtAfterQC",50,10,60,"","Events");
  TruthA1PhiAfterQC=HConfig.GetTH1D(Name+"_TruthA1PhiAfterQC","TruthA1PhiAfterQC",30,-3.14,3.14,"","Events");
  TruthA1EtaAfterQC=HConfig.GetTH1D(Name+"_TruthA1EtaAfterQC","TruthA1EtaAfterQC",20,-2,2,"","Events");

  TruthMuPtAfterQC=HConfig.GetTH1D(Name+"_TruthMuPtAfterQC","TruthMuPtAfterQC",50,10,60,"","Events");
  TruthMuPhiAfterQC=HConfig.GetTH1D(Name+"_TruthMuPhiAfterQC","TruthMuPhiAfterQC",30,-3.14,3.14,"","Events");
  TruthMuEtaAfterQC=HConfig.GetTH1D(Name+"_TruthMuEtaAfterQC","TruthMuEtaAfterQC",20,-2,2,"","Events");



  Niter=HConfig.GetTH1D(Name+"_Niter","Niter",20,0,20,"Niter","Events");
  Probability=HConfig.GetTH1D(Name+"_Probability","Probability",20,0,1,"Probability","Events");
  ProbabilityOfCorrectZPt=HConfig.GetTH2D(Name+"_ProbabilityOfCorrectZPt","ProbabilityOfCorrectZPt",20,0,1,40,0,20,"Probability","Events");
  ProbabilityOfCorrectPtBalance=HConfig.GetTH2D(Name+"_ProbabilityOfCorrectPtBalance","ProbabilityOfCorrectPtBalance",20,0,1,60,-45,45,"Probability","Events");

  csum2Dim=HConfig.GetTH2D(Name+"_csum2Dim","csum2Dim",200,0,1,200,0,1,"csum, GeV","Events");
  csum2Prob=HConfig.GetTH2D(Name+"_csum2Prob","csum2Prob",20,0,1,200,0,1,"csum2Prob, GeV","Events");
  csumZero=HConfig.GetTH1D(Name+"_csumZero","csumZero",200,0,1,"csumZero, GeV","Events");
  csumGF=HConfig.GetTH1D(Name+"_csumGF","csumGF",200,0,1,"csum, GeV","Events");

  TauA1PhiResol=HConfig.GetTH1D(Name+"_TauA1PhiResol","TauA1PhiResol",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauMuPhiResol=HConfig.GetTH1D(Name+"_TauMuPhiResol","TauMuPhiResol",40,-1.5,1.5,"#Delta #phi (rad)","Events");

  TauA1EtaResol=HConfig.GetTH1D(Name+"_TauA1EtaResol","TauA1EtaResol",40,-2.5,2.5,"#Delta #eta","Events");
  TauMuEtaResol=HConfig.GetTH1D(Name+"_TauMuEtaResol","TauMuEtaResol",40,-2.5,2.5,"#Delta #eta","Events");

  TauA1PtResol=HConfig.GetTH1D(Name+"_TauA1PtResol","TauA1PtResol",80,-200,200,"#Delta Pt (GeV)","Events");
  TauMuPtResol=HConfig.GetTH1D(Name+"_TauMuPtResol","TauMuPtResol",80,-200,200,"#Delta Pt (GeV)","Events");

  TauA1EResol=HConfig.GetTH1D(Name+"_TauA1EResol","TauA1EResol",80,-200,200,"|E| (GeV)","Events");
  TauMuEResol=HConfig.GetTH1D(Name+"_TauMuEResol","TauMuEResol",80,-200,200,"|E| (GeV)","Events");


  TauA1PhiResolZPtCut=HConfig.GetTH1D(Name+"_TauA1PhiResolZPtCut","TauA1PhiResolZPtCut",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauMuPhiResolZPtCut=HConfig.GetTH1D(Name+"_TauMuPhiResolZPtCut","TauMuPhiResolZPtCut",40,-1.5,1.5,"#Delta #phi (rad)","Events");

  TauA1EtaResolZPtCut=HConfig.GetTH1D(Name+"_TauA1EtaResolZPtCut","TauA1EtaResolZPtCut",40,-2.5,2.5,"#Delta #eta","Events");
  TauMuEtaResolZPtCut=HConfig.GetTH1D(Name+"_TauMuEtaResolZPtCut","TauMuEtaResolZPtCut",40,-2.5,2.5,"#Delta #eta","Events");

  TauA1PtResolZPtCut=HConfig.GetTH1D(Name+"_TauA1PtResolZPtCut","TauA1PtResolZPtCut",80,-200,200,"#Delta Pt (GeV)","Events");
  TauMuPtResolZPtCut=HConfig.GetTH1D(Name+"_TauMuPtResolZPtCut","TauMuPtResolZPtCut",80,-200,200,"#Delta Pt (GeV)","Events");

  TauA1EResolZPtCut=HConfig.GetTH1D(Name+"_TauA1EResolZPtCut","TauA1EResolZPtCut",80,-200,200,"|E| (GeV)","Events");
  TauMuEResolZPtCut=HConfig.GetTH1D(Name+"_TauMuEResolZPtCut","TauMuEResolZPtCut",80,-200,200,"|E| (GeV)","Events");



  TauA1PhiResolCorrectAmbiga=HConfig.GetTH1D(Name+"_TauA1PhiResolCorrectAmbiga","TauA1PhiResolCorrectAmbiga",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauMuPhiResolCorrectAmbiga=HConfig.GetTH1D(Name+"_TauMuPhiResolCorrectAmbiga","TauMuPhiResolCorrectAmbiga",40,-1.5,1.5,"#Delta #phi (rad)","Events");

  TauA1EtaResolCorrectAmbiga=HConfig.GetTH1D(Name+"_TauA1EtaResolCorrectAmbiga","TauA1EtaResolCorrectAmbiga",40,-2.5,2.5,"#Delta #eta","Events");
  TauMuEtaResolCorrectAmbiga=HConfig.GetTH1D(Name+"_TauMuEtaResolCorrectAmbiga","TauMuEtaResolCorrectAmbiga",40,-2.5,2.5,"#Delta #eta","Events");

  TauA1PtResolCorrectAmbiga=HConfig.GetTH1D(Name+"_TauA1PtResolCorrectAmbiga","TauA1PtResolCorrectAmbiga",80,-200,200,"#Delta Pt (GeV)","Events");
  TauMuPtResolCorrectAmbiga=HConfig.GetTH1D(Name+"_TauMuPtResolCorrectAmbiga","TauMuPtResolCorrectAmbiga",80,-200,200,"#Delta Pt (GeV)","Events");

  TauA1EResolCorrectAmbiga=HConfig.GetTH1D(Name+"_TauA1EResolCorrectAmbiga","TauA1EResolCorrectAmbiga",80,-200,200,"|E| (GeV)","Events");
  TauMuEResolCorrectAmbiga=HConfig.GetTH1D(Name+"_TauMuEResolCorrectAmbiga","TauMuEResolCorrectAmbiga",80,-200,200,"|E| (GeV)","Events");


  TauA1PhiResolZero=HConfig.GetTH1D(Name+"_TauA1PhiResolZero","TauA1PhiResolZero",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauMuPhiResolZero=HConfig.GetTH1D(Name+"_TauMuPhiResolZero","TauMuPhiResolZero",40,-1.5,1.5,"#Delta #phi (rad)","Events");

  TauA1EtaResolZero=HConfig.GetTH1D(Name+"_TauA1EtaResolZero","TauA1EtaResolZero",40,-2.5,2.5,"#Delta #eta","Events");
  TauMuEtaResolZero=HConfig.GetTH1D(Name+"_TauMuEtaResolZero","TauMuEtaResolZero",40,-2.5,2.5,"#Delta #eta","Events");

  TauA1PtResolZero=HConfig.GetTH1D(Name+"_TauA1PtResolZero","TauA1PtResolZero",80,-200,200,"#Delta Pt (GeV)","Events");
  TauMuPtResolZero=HConfig.GetTH1D(Name+"_TauMuPtResolZero","TauMuPtResolZero",80,-200,200,"#Delta Pt (GeV)","Events");

  TauA1EResolZero=HConfig.GetTH1D(Name+"_TauA1EResolZero","TauA1EResolZero",80,-200,200,"|E| (GeV)","Events");
  TauMuEResolZero=HConfig.GetTH1D(Name+"_TauMuEResolZero","TauMuEResolZero",80,-200,200,"|E| (GeV)","Events");


  TauA1PtZero=HConfig.GetTH1D(Name+"_TauA1PtZero","TauA1PtZero",50,10,60,"","Events");
  TauA1PhiZero=HConfig.GetTH1D(Name+"_TauA1PhiZero","TauA1PhiZero",30,-3.14,3.14,"","Events");
  TauA1EtaZero=HConfig.GetTH1D(Name+"_TauA1EtaZero","TauA1EtaZero",20,-2,2,"","Events");

  TauMuPtZero=HConfig.GetTH1D(Name+"_TauMuPtZero","TauMuPtZero",50,10,60,"","Events");
  TauMuPhiZero=HConfig.GetTH1D(Name+"_TauMuPhiZero","TauMuPhiZero",30,-3.14,3.14,"","Events");
  TauMuEtaZero=HConfig.GetTH1D(Name+"_TauMuEtaZero","TauMuEtaZero",20,-2,2,"","Events");



  NiterQC=HConfig.GetTH1D(Name+"_NiterQC","NiterQC",20,0,20,"NiterQC","Events");
  ProbabilityQC=HConfig.GetTH1D(Name+"_ProbabilityQC","ProbabilityQC",20,0,1,"ProbabilityQC","Events");


  TauA1PhiResolQC=HConfig.GetTH1D(Name+"_TauA1PhiResolQC","TauA1PhiResolQC",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauMuPhiResolQC=HConfig.GetTH1D(Name+"_TauMuPhiResolQC","TauMuPhiResolQC",40,-1.5,1.5,"#Delta #phi (rad)","Events");

  TauA1EtaResolQC=HConfig.GetTH1D(Name+"_TauA1EtaResolQC","TauA1EtaResolQC",40,-2.5,2.5,"#Delta #eta","Events");
  TauMuEtaResolQC=HConfig.GetTH1D(Name+"_TauMuEtaResolQC","TauMuEtaResolQC",40,-2.5,2.5,"#Delta #eta","Events");

  TauA1PtResolQC=HConfig.GetTH1D(Name+"_TauA1PtResolQC","TauA1PtResolQC",80,-200,200,"#Delta Pt (GeV)","Events");
  TauMuPtResolQC=HConfig.GetTH1D(Name+"_TauMuPtResolQC","TauMuPtResolQC",80,-200,200,"#Delta Pt (GeV)","Events");

  TauA1EResolQC=HConfig.GetTH1D(Name+"_TauA1EResolQC","TauA1EResolQC",80,-200,200,"#Delta E (GeV)","Events");
  TauMuEResolQC=HConfig.GetTH1D(Name+"_TauMuEResolQC","TauMuEResolQC",80,-200,200,"#Delta E (GeV)","Events");

  TauA1PzResolQC=HConfig.GetTH1D(Name+"_TauA1PzResolQC","TauA1PzResolQC",80,-200,200,"#Delta E (GeV)","Events");
  TauMuPzResolQC=HConfig.GetTH1D(Name+"_TauMuPzResolQC","TauMuPzResolQC",80,-200,200,"#Delta E (GeV)","Events");



  TauA1PtResolKFP=HConfig.GetTH1D(Name+"_TauA1PtResolKFP","TauA1PtResolKFP",80,-200,200,"#Delta Pt (GeV)","Events");
  TauA1PhiResolKFP=HConfig.GetTH1D(Name+"_TauA1PhiResolKFP","TauA1PhiResolKFP",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauA1EtaResolKFP=HConfig.GetTH1D(Name+"_TauA1EtaResolKFP","TauA1EtaResolKFP",40,-2.5,2.5,"#Delta #eta","Events");
  
  TauA1PtResolKFM=HConfig.GetTH1D(Name+"_TauA1PtResolKFM","TauA1PtResolKFM",80,-200,200,"#Delta Pt (GeV)","Events");
  TauA1PhiResolKFM=HConfig.GetTH1D(Name+"_TauA1PhiResolKFM","TauA1PhiResolKFM",40,-1.5,1.5,"#Delta #phi (rad)","Events");
  TauA1EtaResolKFM=HConfig.GetTH1D(Name+"_TauA1EtaResolKFM","TauA1EtaResolKFM",40,-2.5,2.5,"#Delta #eta","Events");
  
  ZPtResol=HConfig.GetTH1D(Name+"_ZPtResol","ZPtResol",80,-200,200,"#Delta Pt (GeV)","Events");


  ZPzResolRelative=HConfig.GetTH1D(Name+"_ZPzResolRelative","ZPzResolRelative",100,-10,10,"#Delta Pz/Pz ","Events");
  ZEResolRelative=HConfig.GetTH1D(Name+"_ZEResolRelative","ZEResolRelative",100,-10,10,"#Delta E ","Events");


  ZPtResolZPtCut=HConfig.GetTH1D(Name+"_ZPtResolZPtCut","ZPtResolZPtCut",80,-200,200,"#Delta Pt (GeV)","Events");


  ZPzResolRelativeZPtCut=HConfig.GetTH1D(Name+"_ZPzResolRelativeZPtCut","ZPzResolRelativeZPtCut",100,-10,10,"#Delta Pz/Pz ","Events");
  ZEResolRelativeZPtCut=HConfig.GetTH1D(Name+"_ZEResolRelativeZPtCut","ZEResolRelativeZPtCut",100,-10,10,"#Delta E ","Events");



  ZPtTruth=HConfig.GetTH1D(Name+"_ZPtTruth","ZPtTruth",100,0,50,"Pt^{Z}, GeV ","Events");
  ZPtTruthTauMuPtResolution=HConfig.GetTH2D(Name+"_ZPtTruthTauMuPtResolution","ZPtTruthTauMuPtResolution",100,0,50,100,-100,100,"Pt^{Z}, GeV ","TauMuPt resolution");


  ZMass=HConfig.GetTH1D(Name+"_ZMass","ZMass",80,85,96,"|M| (GeV)","Events");

  pvsvsignificance=HConfig.GetTH1D(Name+"_pvsvsignificance","pvsvsignificance",60,0,20,"PV-SV significance","Events");
  pvMusvsignificance=HConfig.GetTH1D(Name+"_pvMusvsignificance","pvMusvsignificance",60,0,20,"PV-MUSV significance","Events");


  Chi2Dim=HConfig.GetTH2D(Name+"_Chi2Dim","Chi2Dim",20,0,1,20,0,1,"Probability of correct ambiguity point"," probability of opposite ambiguity");
  Chi2DimPlusMinus=HConfig.GetTH2D(Name+"_Chi2DimPlusMinus","Chi2DimPlusMinus",20,0,1,20,0,1,"Probability of Minus"," Probability of Plus");

  TauMuEResolVsSignificance=HConfig.GetTH2D(Name+"_TauMuEResolVsSignificance","TauMuEResolVsSignificance",100,-100,100,10,0,10,"Probability of correct ambiguity point"," significance");
  TauMuEResolVsSignificance2=HConfig.GetTH2D(Name+"_TauMuEResolVsSignificance2","TauMuEResolVsSignificance2",100,-100,100,10,0,10,"Probability of correct ambiguity point"," significance");
  TauMuEResolVsSignificanceSlices=HConfig.GetTH1D(Name+"_TauMuEResolVsSignificanceSlices","TauMuEResolVsSignificanceSlices",10,0,10,"Slices"," significance");



  TauMuEResolVsChi2Probability=HConfig.GetTH2D(Name+"_TauMuEResolVsChi2Probability","TauMuEResolVsChi2Probability",100,-100,100,20,0,1,"Probability of correct ambiguity point"," probability");


  MissingEnergyCheckCorrectAmbiga=HConfig.GetTH2D(Name+"_MissingEnergyCheckCorrectAmbiga","MissingEnergyCheckCorrectAmbiga",100,0,100,100,0,100," Correct Ambiguity  Et_{(#Sigma(3#nu))},  GeV"," MET");
  MissingEnergyCheckInCorrectAmbiga=HConfig.GetTH2D(Name+"_MissingEnergyCheckInCorrectAmbiga","MissingEnergyCheckInCorrectAmbiga",100,0,100,100,0,100,"Incorrect Ambiguity Et_{(#Sigma(3#nu))},  GeV"," MET");


  Ambiguity=HConfig.GetTH1D(Name+"_Ambiguity","Ambiguity",2,0.5,2.5,"Ambiguity ","Events");

  idEff=HConfig.GetTH1D(Name+"_idEff","idEff",4,0.5,4.5,"idEff ","Events");
  idPassHPS=HConfig.GetTH1D(Name+"_idPassHPS","idPassHPS",4,0.5,4.5," ","Events");
  idPassKFMinus=HConfig.GetTH1D(Name+"_idPassKFMinus","idPassKFMinus",4,0.5,4.5," ","Events");
  idPassKFPlus=HConfig.GetTH1D(Name+"_idPassKFPlus","idPassKFPlus",4,0.5,4.5," ","Events");
  idPassGF=HConfig.GetTH1D(Name+"_idPassGF","idPassGF",4,0.5,4.5," ","Events");

  idPassGFZero=HConfig.GetTH1D(Name+"_idPassGFZero","idPassGFZero",4,0.5,4.5," ","Events");

  AmbiguityIsCorrect=HConfig.GetTH1D(Name+"_AmbiguityIsCorrect","AmbiguityIsCorrect",2,-0.5,1.5,"AmbiguityIsCorrect ","Events");

  RightlyTakenAmbiguity=HConfig.GetTH1D(Name+"_RightlyTakenAmbiguity","RightlyTakenAmbiguity",2,0.5,2.5,"Correct Ambiguity ","Events");
  WronglyTakenAmbiguity=HConfig.GetTH1D(Name+"_WronglyTakenAmbiguity","WronglyTakenAmbiguity",2,0.5,2.5,"Incorrect Ambiguity ","Events");
  CorrectAmbiguityTruth=HConfig.GetTH1D(Name+"_CorrectAmbiguityTruth","CorrectAmbiguityTruth",2,0.5,2.5,"Correct Ambiguity Truth ","Events");

  AmbiguityEfficiency=HConfig.GetTH1D(Name+"_AmbiguityEfficiency","AmbiguityEfficiency",2,0.5,2.5,"Ambiguity Efficiency","Events");

  deltaPhiSignificance=HConfig.GetTH1D(Name+"_deltaPhiSignificance","deltaPhiSignificance",20,0,20,"delta Phi significance","Events");

  ProbabilityOfCorrect=HConfig.GetTH1D(Name+"_ProbabilityOfCorrect","ProbabilityOfCorrect",20,0,1,"#chi^{2} propability of correct ambiguity","Events");



  EtaTauMuEtaTauA1=HConfig.GetTH2D(Name+"_EtaTauMuEtaTauA1","EtaTauMuEtaTauA1",40,-2,2,40,-2,2,"#eta_{#tau_{#mu}}"," #eta_{#tau_{a1}}");


  ResTauMuPtEtaTauMu=HConfig.GetTH2D(Name+"_ResTauMuPtEtaTauMu","ResTauMuPtEtaTauMu",80,-80,80,40,-2,2,"#tau_{#mu}  pt resolution"," #eta_{#tau_{#mu}}");
  ResTauMuEEtaTauMu=HConfig.GetTH2D(Name+"_ResTauMuEEtaTauMu","ResTauMuEEtaTauMu",80,-80,80,40,-2,2,"#tau_{#mu}  e resolution, GeV"," #eta_{#tau_{#mu}}");
  ResTauMuPtPtTauMu=HConfig.GetTH2D(Name+"_ResTauMuPtPtTauMu","ResTauMuPtPtTauMu",80,-80,80,80,10,70,"#tau_{#mu}  pt resolution, GeV"," pt_{#tau_{#mu}}, GeV");


  ResTauA1PtEtaTauA1=HConfig.GetTH2D(Name+"_ResTauA1PtEtaTauA1","ResTauA1PtEtaTauA1",80,-80,80,40,-2,2,"#tau_{a1}  pt resolution"," #eta_{#tau_{a1}}");
  ResTauA1EEtaTauA1=HConfig.GetTH2D(Name+"_ResTauA1EEtaTauA1","ResTauA1EEtaTauA1",80,-80,80,40,-2,2,"#tau_{a1}  e resolution, GeV"," #eta_{#tau_{a1}}");
  ResTauA1PtPtTauA1=HConfig.GetTH2D(Name+"_ResTauA1PtPtTauA1","ResTauA1PtPtTauA1",80,-80,80,80,10,70,"#tau_{a1}  pt resolution, GeV"," pt_{#tau_{a1}}, GeV");

  PhiRotationSignificancePhysicalTaus=HConfig.GetTH1D(Name+"_PhiRotationSignificancePhysicalTaus","PhiRotationSignificancePhysicalTaus",20,0,1,"Significance of phi rotation for physical taus","Events");
  PhiRotationSignificanceUnPhysicalTaus=HConfig.GetTH1D(Name+"_PhiRotationSignificanceUnPhysicalTaus","PhiRotationSignificanceUnPhysicalTaus",20,0,10,"Significance of phi rotation for unpysical taus","Events");



  ResTauMuPtEtaTauMuProfile=HConfig.GetTH1D(Name+"_ResTauMuPtEtaTauMuProfile","ResTauMuPtEtaTauMuProfile",40,-2,2,"#eta_{#tau_{#mu}}","Pt resolution");
  ResTauMuEEtaTauMuProfile=HConfig.GetTH1D(Name+"_ResTauMuEEtaTauMuProfile","ResTauMuEEtaTauMuProfile",40,-2,2,"#eta_{#tau_{#mu}}","E resolution");
  ResTauMuPtPtTauMuProfile=HConfig.GetTH1D(Name+"_ResTauMuPtPtTauMuProfile","ResTauMuPtPtTauMuProfile",80,10,70,"Pt_{#mu}, GeV","Pt resolution");

  ResTauA1PtEtaTauA1Profile=HConfig.GetTH1D(Name+"_ResTauA1PtEtaTauA1Profile","ResTauA1PtEtaTauA1Profile",40,-2,2,"#eta_{#tau_{a1}}","pt resolution");
  ResTauA1EEtaTauA1Profile=HConfig.GetTH1D(Name+"_ResTauA1EEtaTauA1Profile","ResTauA1EEtaTauA1Profile",40,-2,2,"#eta_{#tau_{a1}}","E resolution");
  ResTauA1PtPtTauA1Profile=HConfig.GetTH1D(Name+"_ResTauA1PtPtTauA1Profile","ResTauA1PtPtTauA1Profile",80,10,70,"Pt_{a1}, GeV","pt resolution");


  ResTauMuPtProb=HConfig.GetTH2D(Name+"_ResTauMuPtProb","ResTauMuPtProb",80,-80,80,20,0,1,"#chi^{2} Probability","pt resolution ");
  ResTauMuEProb=HConfig.GetTH2D(Name+"_ResTauMuEProb","ResTauMuEProb",80,-80,80,20,0,1,"#chi^{2} Probability","pt resolution ");
  ResTauA1PtProb=HConfig.GetTH2D(Name+"_ResTauA1PtProb","ResTauA1PtProb",80,-80,80,20,0,1,"#chi^{2} Probability","pt resolution ");
  ResTauA1EProb=HConfig.GetTH2D(Name+"_ResTauA1EProb","ResTauA1EProb",80,-80,80,20,0,1,"#chi^{2} Probability","pt resolution ");
  
  ResTauMuPtProbProfile=HConfig.GetTH1D(Name+"_ResTauMuPtProbProfile","ResTauMuPtProbProfile",20,0,1,"#chi^{2} Probability","pt resolution");
  ResTauMuEProbProfile=HConfig.GetTH1D(Name+"_ResTauMuEProbProfile","ResTauMuEProbProfile",20,0,1,"#chi^{2} Probability","E resolution");

  ResTauA1PtProbProfile=HConfig.GetTH1D(Name+"_ResTauA1PtProbProfile","ResTauA1PtProbProfile",20,0,1,"#chi^{2} Probability","pt resolution");
  ResTauA1EProbProfile=HConfig.GetTH1D(Name+"_ResTauA1EProbProfile","ResTauA1EProbProfile",20,0,1,"#chi^{2} Probability","E resolution");
  MeasuredPhi=HConfig.GetTH1D(Name+"_MeasuredPhi","MeasuredPhi",40,-3.14,3.14,"  measured phi ","E resolution");

  KFTauPhi_Minus=HConfig.GetTH1D(Name+"_KFTauPhi_Minus","KFTauPhi_Minus",40,-3.14,3.14,"KFTauPhi_Minus ","E resolution");
  KFTauPhi_Plus=HConfig.GetTH1D(Name+"_KFTauPhi_Plus","KFTauPhi_Plus",40,-3.14,3.14,"KFTauPhi_Plus ","E resolution");
  KFTauPhi_Ambiga=HConfig.GetTH1D(Name+"_KFTauPhi_Ambiga","KFTauPhi_Ambiga",40,-3.14,3.14," KFTauPhi_Ambiga ","E resolution");
  KFTauPhi_Ambiga_DeltaPhiSignificanceCut=HConfig.GetTH1D(Name+"_KFTauPhi_Ambiga_DeltaPhiSignificanceCut","KFTauPhi_Ambiga_DeltaPhiSignificanceCut",40,-3.14,3.14," KFTauPhi_Ambiga deltaPhiSignificance Cut ","E resolution");

  KFTauPhi_Ambiga_vsDeltaPhiSignificance=HConfig.GetTH2D(Name+"_KFTauPhi_Ambiga_vsDeltaPhiSignificance","KFTauPhi_Ambiga_vsDeltaPhiSignificance",40,-3.14,3.14,50,0,10,"#chi^{2} Probability","pt resolution ");
  KFTauPhi_Ambiga_vsDeltaPhi=HConfig.GetTH2D(Name+"_KFTauPhi_Ambiga_vsDeltaPhi","KFTauPhi_Ambiga_vsDeltaPhi",40,-3.14,3.14,100,-0.3,0.3,"#chi^{2} Probability","pt resolution ");

  PhotonEnergyFraction1=HConfig.GetTH1D(Name+"_PhotonEnergyFraction1","PhotonEnergyFraction1",60,0,1,"PhotonEnergyFraction1","");
  PhotonEnergyFraction2=HConfig.GetTH1D(Name+"_PhotonEnergyFraction2","PhotonEnergyFraction2",60,0,1,"PhotonEnergyFraction2","");

  A1Mass=HConfig.GetTH1D(Name+"_A1Mass","A1Mass",60,0.8,1.5,"A1Mass","");
  deltaAmbigaEnergyForWrongAmb=HConfig.GetTH1D(Name+"_deltaAmbigaEnergyForWrongAmb","deltaAmbigaEnergyForWrongAmb",80,0,80,"deltaAmbigaEnergyForWrongAmb","");
  deltaAmbigaEnergyForCorreAmb=HConfig.GetTH1D(Name+"_deltaAmbigaEnergyForCorreAmb","deltaAmbigaEnergyForCorreAmb",80,0,80,"deltaAmbigaEnergyForCorreAmb","");




  TauA1PhiResolVsProb=HConfig.GetTH2D(Name+"_TauA1PhiResolVsProb","TauA1PhiResolVsProb",40,0,1,40,-1.5,1.5,"Prob","#Delta #phi (rad)");
  TauMuPhiResolVsProb=HConfig.GetTH2D(Name+"_TauMuPhiResolVsProb","TauMuPhiResolVsProb",40,0,1,40,-1.5,1.5,"Prob","#Delta #phi (rad)");

  TauA1EtaResolVsProb=HConfig.GetTH2D(Name+"_TauA1EtaResolVsProb","TauA1EtaResolVsProb",40,0,1,40,-2.5,2.5,"Prob","#Delta #eta");
  TauMuEtaResolVsProb=HConfig.GetTH2D(Name+"_TauMuEtaResolVsProb","TauMuEtaResolVsProb",40,0,1,40,-2.5,2.5,"Prob","#Delta #eta");

  TauA1PtResolVsProb=HConfig.GetTH2D(Name+"_TauA1PtResolVsProb","TauA1PtResolVsProb",40,0,1,80,-200,200,"Prob","#Delta Pt (GeV)");
  TauMuPtResolVsProb=HConfig.GetTH2D(Name+"_TauMuPtResolVsProb","TauMuPtResolVsProb",40,0,1,80,-200,200,"Prob","#Delta Pt (GeV)");

  TauA1EResolVsProb=HConfig.GetTH2D(Name+"_TauA1EResolVsProb","TauA1EResolVsProb",40,0,1,80,-200,200,"Prob","|E| (GeV)");
  TauMuEResolVsProb=HConfig.GetTH2D(Name+"_TauMuEResolVsProb","TauMuEResolVsProb",40,0,1,80,-200,200,"Prob","|E| (GeV)");


  AmbiguitySolverH=HConfig.GetTH1D(Name+"_AmbiguitySolverH","AmbiguitySolverH",2,-0.5,1.5,"Physical region = 1, Unphysical = 0","");
  AmbiguitySolverHvsSign=HConfig.GetTH2D(Name+"_AmbiguitySolverHvsSign","AmbiguitySolverHvsSign",2,-0.5,1.5,100,0,10,"Physical region = 1, Unphysical = 0","");


  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    //std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
  //  ran = new TRandom();
}



void  FitTest::Store_ExtraDist(){
//  Extradist1d.push_back(&NVtx);
//  Extradist1d.push_back(&NGoodVtx);
//  Extradist1d.push_back(&TauFlightLength);
//  Extradist1d.push_back(&TauFlightLengthTransverse);
//  Extradist1d.push_back(&TauMomentum);
//  Extradist1d.push_back(&TauMomentumTransverse);
//  Extradist1d.push_back(&TauLife);
//  Extradist1d.push_back(&TauLifeTransverse);
//  Extradist1d.push_back(&ResTauFlightLength);
//  Extradist1d.push_back(&ResTauFlightLengthTransverse);
//  Extradist1d.push_back(&ResTauMomentum);
//  Extradist1d.push_back(&ResTauMomentumTransverse);
//  Extradist1d.push_back(&ResTauLife);
//  Extradist1d.push_back(&ResTauLifeTransverse);

//  Extradist1d.push_back(&TauMomentumAvg);
//  Extradist1d.push_back(&TauMomentumTransverseAvg);
//  Extradist1d.push_back(&ResTauMomentumAvg);
//  Extradist1d.push_back(&ResTauMomentumTransverseAvg);
//  Extradist1d.push_back(&TauLifeAvg);
//  Extradist1d.push_back(&TauLifeTransverseAvg);
//  Extradist1d.push_back(&ResTauLifeAvg);
//  Extradist1d.push_back(&ResTauLifeTransverseAvg);


 Extradist1d.push_back(&ZMass0);
 Extradist1d.push_back(&ZMass1);
 Extradist1d.push_back(&ZMass2);

 Extradist1d.push_back(&TauA1PtResol0);
 Extradist1d.push_back(&TauA1PtResol1);
 Extradist1d.push_back(&TauA1PtResol2);

 Extradist1d.push_back(&TauMuPtResol0);
 Extradist1d.push_back(&TauMuPtResol1);
 Extradist1d.push_back(&TauMuPtResol2);

 Extradist1d.push_back(&Probability0);
 Extradist1d.push_back(&Probability1);
 Extradist1d.push_back(&Probability2);

 Extradist1d.push_back(&TauA1PhiResol0);
 Extradist1d.push_back(&TauA1PhiResol1);
 Extradist1d.push_back(&TauA1PhiResol2);

 Extradist1d.push_back(&TauMuPhiResol0);
 Extradist1d.push_back(&TauMuPhiResol1);
 Extradist1d.push_back(&TauMuPhiResol2);

 Extradist1d.push_back(&TauA1EtaResol0);
 Extradist1d.push_back(&TauA1EtaResol1);
 Extradist1d.push_back(&TauA1EtaResol2);

 Extradist1d.push_back(&TauMuEtaResol0);
 Extradist1d.push_back(&TauMuEtaResol1);
 Extradist1d.push_back(&TauMuEtaResol2);

 Extradist1d.push_back(&TauMuEnergyPull0);
 Extradist1d.push_back(&TauA1EnergyPull0);

 Extradist1d.push_back(&TauMuEnergyPull1);
 Extradist1d.push_back(&TauA1EnergyPull1);

 Extradist1d.push_back(&TauMuEnergyPull2);
 Extradist1d.push_back(&TauA1EnergyPull2);

 Extradist1d.push_back(&Niter0);
 Extradist1d.push_back(&Niter1);
 Extradist1d.push_back(&Niter2);



 Extradist1d.push_back(&TruthA1Pt);
 Extradist1d.push_back(&TruthA1Phi);
 Extradist1d.push_back(&TruthA1Eta);

 Extradist1d.push_back(&TruthMuPt);
 Extradist1d.push_back(&TruthMuPhi);
 Extradist1d.push_back(&TruthMuEta);

 Extradist1d.push_back(&TruthA1PtAfter);
 Extradist1d.push_back(&TruthA1PhiAfter);
 Extradist1d.push_back(&TruthA1EtaAfter);

 Extradist1d.push_back(&TruthMuPtAfter);
 Extradist1d.push_back(&TruthMuPhiAfter);
 Extradist1d.push_back(&TruthMuEtaAfter);



 Extradist1d.push_back(&TruthA1PtAfterQC);
 Extradist1d.push_back(&TruthA1PhiAfterQC);
 Extradist1d.push_back(&TruthA1EtaAfterQC);

 Extradist1d.push_back(&TruthMuPtAfterQC);
 Extradist1d.push_back(&TruthMuPhiAfterQC);
 Extradist1d.push_back(&TruthMuEtaAfterQC);


 Extradist2d.push_back(&Chi2Dim);

 Extradist1d.push_back(&pvsvsignificance);
 Extradist1d.push_back(&pvMusvsignificance);


 Extradist1d.push_back(&Niter);
 Extradist1d.push_back(&Probability);

 Extradist1d.push_back(&TauA1PhiResol);
 Extradist1d.push_back(&TauMuPhiResol);

 Extradist1d.push_back(&TauA1EtaResol);
 Extradist1d.push_back(&TauMuEtaResol);

 Extradist1d.push_back(&TauA1PtResol);
 Extradist1d.push_back(&TauMuPtResol);

 Extradist1d.push_back(&TauA1EResol);
 Extradist1d.push_back(&TauMuEResol);


 Extradist1d.push_back(&TauA1PhiResolZPtCut);
 Extradist1d.push_back(&TauMuPhiResolZPtCut);

 Extradist1d.push_back(&TauA1EtaResolZPtCut);
 Extradist1d.push_back(&TauMuEtaResolZPtCut);

 Extradist1d.push_back(&TauA1PtResolZPtCut);
 Extradist1d.push_back(&TauMuPtResolZPtCut);

 Extradist1d.push_back(&TauA1EResolZPtCut);
 Extradist1d.push_back(&TauMuEResolZPtCut);


 Extradist1d.push_back(&ZMass);
 Extradist1d.push_back(&ZPtResol);
 Extradist1d.push_back(&ZPzResolRelative);
 Extradist1d.push_back(&ZEResolRelative);

 Extradist1d.push_back(&ZPtResolZPtCut);
 Extradist1d.push_back(&ZPzResolRelativeZPtCut);
 Extradist1d.push_back(&ZEResolRelativeZPtCut);


 Extradist1d.push_back(&NiterQC);
 Extradist1d.push_back(&ProbabilityQC);

 Extradist1d.push_back(&TauA1PhiResolQC);
 Extradist1d.push_back(&TauMuPhiResolQC);

 Extradist1d.push_back(&TauA1EtaResolQC);
 Extradist1d.push_back(&TauMuEtaResolQC);

 Extradist1d.push_back(&TauA1PtResolQC);
 Extradist1d.push_back(&TauMuPtResolQC);

 Extradist1d.push_back(&TauA1EResolQC);
 Extradist1d.push_back(&TauMuEResolQC);


 Extradist1d.push_back(&TauA1PzResolQC);
 Extradist1d.push_back(&TauMuPzResolQC);

 Extradist1d.push_back(&TauA1PtResolKFP);
 Extradist1d.push_back(&TauA1PhiResolKFP);
 Extradist1d.push_back(&TauA1EtaResolKFP);
 Extradist1d.push_back(&ZPtTruth);



 Extradist1d.push_back(&TauA1PtResolKFM);
 Extradist1d.push_back(&TauA1PhiResolKFM);
 Extradist1d.push_back(&TauA1EtaResolKFM);

 Extradist2d.push_back(&TauMuEResolVsSignificance);
 Extradist2d.push_back(&TauMuEResolVsChi2Probability);
 Extradist2d.push_back(&TauMuEResolVsSignificance2);
 Extradist2d.push_back(&ZPtTruthTauMuPtResolution);
 Extradist2d.push_back(&Chi2DimPlusMinus);

 Extradist1d.push_back(&TauA1PhiResolCorrectAmbiga);
 Extradist1d.push_back(&TauMuPhiResolCorrectAmbiga);
	    
 Extradist1d.push_back(&TauA1EtaResolCorrectAmbiga);
 Extradist1d.push_back(&TauMuEtaResolCorrectAmbiga);
	    
 Extradist1d.push_back(&TauA1PtResolCorrectAmbiga);
 Extradist1d.push_back(&TauMuPtResolCorrectAmbiga);
	    
 Extradist1d.push_back(&TauA1EResolCorrectAmbiga);
 Extradist1d.push_back(&TauMuEResolCorrectAmbiga);

 Extradist1d.push_back(&Ambiguity);
 Extradist1d.push_back(&AmbiguityIsCorrect);
 Extradist1d.push_back(&AmbiguityEfficiency);

 Extradist1d.push_back(&RightlyTakenAmbiguity);
 Extradist1d.push_back(&WronglyTakenAmbiguity);
 Extradist1d.push_back(&CorrectAmbiguityTruth);


 Extradist2d.push_back(&MissingEnergyCheckCorrectAmbiga);
 Extradist2d.push_back(&MissingEnergyCheckInCorrectAmbiga);

 Extradist1d.push_back(&deltaPhiSignificance);

 Extradist1d.push_back(&idEff);
 Extradist1d.push_back(&idPassHPS);
 Extradist1d.push_back(&idPassKFMinus);
 Extradist1d.push_back(&idPassKFPlus);
 Extradist1d.push_back(&idPassGF);


 Extradist1d.push_back(&ProbabilityOfCorrect);
 Extradist1d.push_back(&TauMuEResolVsSignificanceSlices);



 Extradist2d.push_back(&EtaTauMuEtaTauA1);


 Extradist2d.push_back(&ResTauMuPtEtaTauMu);
 Extradist2d.push_back(&ResTauMuEEtaTauMu);
 Extradist2d.push_back(&ResTauMuPtPtTauMu);


 Extradist2d.push_back(&ResTauA1PtEtaTauA1);
 Extradist2d.push_back(&ResTauA1PtPtTauA1);
 Extradist2d.push_back(&ResTauA1EEtaTauA1);

 Extradist2d.push_back(&ResTauMuPtProb);
 Extradist2d.push_back(&ResTauMuEProb);
 Extradist2d.push_back(&ResTauA1PtProb);
 Extradist2d.push_back(&ResTauA1EProb);



 Extradist1d.push_back(&ResTauMuPtEtaTauMuProfile);
 Extradist1d.push_back(&ResTauMuEEtaTauMuProfile);
 Extradist1d.push_back(&ResTauMuPtPtTauMuProfile);

 Extradist1d.push_back(&ResTauA1PtEtaTauA1Profile);
 Extradist1d.push_back(&ResTauA1PtPtTauA1Profile);
 Extradist1d.push_back(&ResTauA1EEtaTauA1Profile);

 Extradist1d.push_back(&ResTauMuPtProbProfile);
 Extradist1d.push_back(&ResTauMuEProbProfile);
 Extradist1d.push_back(&ResTauA1PtProbProfile);
 Extradist1d.push_back(&ResTauA1EProbProfile);

 Extradist1d.push_back(&MeasuredPhi);

 Extradist1d.push_back(&KFTauPhi_Minus);
 Extradist1d.push_back(&KFTauPhi_Plus);
 Extradist1d.push_back(&KFTauPhi_Ambiga);

 Extradist1d.push_back(&KFTauPhi_Ambiga_DeltaPhiSignificanceCut);

 Extradist2d.push_back(&KFTauPhi_Ambiga_vsDeltaPhiSignificance);
 Extradist2d.push_back(&KFTauPhi_Ambiga_vsDeltaPhi);



 Extradist1d.push_back(&PhotonEnergyFraction1);
 Extradist1d.push_back(&PhotonEnergyFraction2);
 Extradist1d.push_back(&A1Mass);

 Extradist1d.push_back(&deltaAmbigaEnergyForWrongAmb);
 Extradist1d.push_back(&deltaAmbigaEnergyForCorreAmb);



 Extradist1d.push_back(&TauA1PhiResolZero);
 Extradist1d.push_back(&TauMuPhiResolZero);
 
 Extradist1d.push_back(&TauA1EtaResolZero);
 Extradist1d.push_back(&TauMuEtaResolZero);
 
 Extradist1d.push_back(&TauA1PtResolZero);
 Extradist1d.push_back(&TauMuPtResolZero);
	    
 Extradist1d.push_back(&TauA1EResolZero);
 Extradist1d.push_back(&TauMuEResolZero);

 Extradist1d.push_back(&TauA1PhiZero);
 Extradist1d.push_back(&TauA1EtaZero);
 Extradist1d.push_back(&TauA1PtZero);
	    
 Extradist1d.push_back(&TauMuPhiZero);
 Extradist1d.push_back(&TauMuEtaZero);
 Extradist1d.push_back(&TauMuPtZero);

 Extradist1d.push_back(&idPassGFZero);

 Extradist2d.push_back(&csum2Dim);
 Extradist1d.push_back(&csumGF);
 Extradist1d.push_back(&csumZero);

 Extradist2d.push_back(&ProbabilityOfCorrectPtBalance);
 Extradist2d.push_back(&ProbabilityOfCorrectZPt);



 Extradist2d.push_back(&TauA1PhiResolVsProb);
 Extradist2d.push_back(&TauMuPhiResolVsProb);

 Extradist2d.push_back(&TauA1EtaResolVsProb);
 Extradist2d.push_back(&TauMuEtaResolVsProb);

 Extradist2d.push_back(&TauA1PtResolVsProb);
 Extradist2d.push_back(&TauMuPtResolVsProb);

 Extradist2d.push_back(&TauA1EResolVsProb);
 Extradist2d.push_back(&TauMuEResolVsProb);
 Extradist1d.push_back(&AmbiguitySolverH);
 Extradist2d.push_back(&AmbiguitySolverHvsSign);

 Extradist2d.push_back(&csum2Prob);


 Extradist1d.push_back(&PhiRotationSignificancePhysicalTaus);
 Extradist1d.push_back(&PhiRotationSignificanceUnPhysicalTaus);


}

void  FitTest::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
    std::cout << "FitTest Ntp->GetMCID(): " << Ntp->GetMCID() << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){/* std::cout << "failed to find id" <<std::endl;*/ return;}

  if(verbose)std::cout << "void  FitTest::doEvent() Mu" << std::endl;


//   // check out triggers
//   for(unsigned int itr = 0; itr <  Ntp->NHLTTriggers(); itr ++){
//     std::cout << "Trigger names  " <<Ntp->HTLTriggerName(itr) << std::endl;
//   }

  int idToFill(0);
  if(id ==998)idToFill = 1;
  if(id ==10230833)idToFill = 2;
  if(id ==10231433)idToFill = 3;
  if(id ==10231833)idToFill = 4;
  idEff.at(t).Fill(idToFill,1);




  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("eta2p1_LooseIsoPFTau"))value.at(TriggerOk)+=1;
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);

  // Apply Selection
  std::vector<unsigned int> mu_idx_good, mu_idx_pt, mu_idx_iso;
  unsigned mu_idx(999);
  //double mu_pt(0);
  for(unsigned int i=0;i<Ntp->NMuons();i++){
    //std::cout << "Checking if muon number = " << i << " is good..." << std::endl;
    if(Ntp->isGoodMuon(i) && fabs(Ntp->Muons_p4(i).Eta())<muoneta){
      mu_idx_good.push_back(i);
    }  
  }
  for(unsigned int i=0;i<mu_idx_good.size();i++){	  
    if(Ntp->Muons_p4(mu_idx_good.at(i)).Pt() > cut.at(TagPtMin)){
      mu_idx_pt.push_back(mu_idx_good.at(i));
    }
  }	
  for(unsigned int i=0;i<mu_idx_pt.size();i++){
    if((Ntp->Muon_emEt05(mu_idx_pt.at(i)) + Ntp->Muon_hadEt05(mu_idx_pt.at(i)) + Ntp->Muon_sumPt05(mu_idx_pt.at(i)))/Ntp->Muons_p4(mu_idx_pt.at(i)).Pt()<=cut.at(TagIso)){
      mu_idx=mu_idx_pt.at(i);//mu_idx is the idx of the effectively selected muon.
      mu_idx_iso.push_back(mu_idx);
    }
  }
  if(verbose)std::cout << "void  FitTest::doEvent() MuA" << std::endl;
  //std::cout << "nmus = " << mu_idx_iso.size() << std::endl;
  value.at(hasTag)=mu_idx_good.size();
  pass.at(hasTag)=(value.at(hasTag)>=cut.at(hasTag));
  value.at(TagPtMin)=mu_idx_pt.size();
  pass.at(TagPtMin)=(value.at(TagPtMin)>0);
  value.at(TagIso)=mu_idx_iso.size();
  pass.at(TagIso)=(value.at(TagIso)==1);


  if(verbose)std::cout << "void  FitTest::doEvent() MuEnd - mu_idx = " << mu_idx << " NMuons = " << Ntp->NMuons() << std::endl;
  if(verbose)std::cout << "void  FitTest::doEvent() Tau" << std::endl;
  
  std::vector<unsigned int> tau_idx_pt, tau_idx_eta, tau_idx_iso, tau_idx_fit, tau_idx_phi, tau_idx_charge;
  unsigned int tau_idx(999);
  if(pass.at(TagIso)==true){
    for(unsigned int i=0;i<Ntp->NPFTaus();i++){
      if(Ntp->PFTau_p4(i).Pt()>cut.at(TauPt)){
        tau_idx_pt.push_back(i);
      }
    }
    for(unsigned int i=0;i<tau_idx_pt.size();i++){
      if(fabs(Ntp->PFTau_p4(tau_idx_pt.at(i)).Eta())<cut.at(TauEta)){
        tau_idx_eta.push_back(tau_idx_pt.at(i));
      }
    }
    for(unsigned int i=0;i<tau_idx_eta.size();i++){
      if(Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA(tau_idx_eta.at(i))){
        tau_idx_iso.push_back(tau_idx_eta.at(i));
      }
    }    
    for(unsigned int i=0;i<tau_idx_iso.size();i++){
      if(Ntp->PFTau_hpsDecayMode(tau_idx_iso.at(i))==10 && 
	 Ntp->PFTau_isHPSByDecayModeFinding(tau_idx_iso.at(i)) && 
	 Ntp->PFTau_TIP_hassecondaryVertex(tau_idx_iso.at(i)) &&
	 Ntp->PFTau_TIP_hasA1Momentum(tau_idx_iso.at(i))){
        tau_idx=tau_idx_iso.at(i);
        tau_idx_fit.push_back(tau_idx);
      }
    }
    double delta_phi(0), charge(999);
    if(verbose)std::cout << "void  FitTest::doEvent() Tau -a" << std::endl;
    //event
    //should have only one element, shouldn't it? for(int i=0;i<tau_idx_iso.size();i++){
    if(tau_idx_fit.size()>=1 && mu_idx_iso.size()>=1){
      //std::cout << "The index of the tau is " << tau_idx_fit.at(0) << " = " << tau_idx << std::endl;
      if(fabs(Tools::DeltaPhi(Ntp->Muons_p4(mu_idx),Ntp->PFTau_p4(tau_idx)))>=cut.at(deltaPhi)){
        delta_phi=fabs(Tools::DeltaPhi(Ntp->Muons_p4(mu_idx),Ntp->PFTau_p4(tau_idx)));
        tau_idx_phi.push_back(tau_idx_fit.at(0));
      }
      if(fabs(Ntp->PFTau_Charge(tau_idx) + Ntp->Muon_Charge(mu_idx))<0.5){ 
        charge=Ntp->PFTau_Charge(tau_idx) + Ntp->Muon_Charge(mu_idx);
        tau_idx_charge.push_back(tau_idx_fit.at(0));      
      }
    }

    if(verbose)std::cout << "void  FitTest::doEvent() Tau -b" << std::endl;
    /*std::cout << "tau_idx_charge.size = " << tau_idx_charge.size() << std::endl;
      std::cout << "tau_idx_phi.size    = " << tau_idx_phi.size() << std::endl;
      std::cout << "tau_idx_iso.size    = " << tau_idx_iso.size() << std::endl;
      std::cout << "tau_idx_eta.size    = " << tau_idx_eta.size() << std::endl;
      std::cout << "tau_idx_pt.size     = " << tau_idx_pt.size() << std::endl;*/



    bool passId = false;
    if(id == 998) passId = true;
    value.at(TauPt)=tau_idx_pt.size();
    value.at(TauEta)=tau_idx_eta.size();
    value.at(TauIsIsolated)=tau_idx_iso.size();
    value.at(TauFit)=tau_idx_fit.size();
    value.at(deltaPhi)=delta_phi;
    value.at(Charge)=charge;

    pass.at(EventId)=true;//passId;
    pass.at(TauPt)=true;//(value.at(TauPt)>0);
    pass.at(TauEta)=true;//(value.at(TauEta)>0);
    pass.at(TauIsIsolated)=true;//(value.at(TauIsIsolated)>0);
    pass.at(TauFit)=true;//(value.at(TauFit)==cut.at(TauFit));
    pass.at(deltaPhi)=true;//(value.at(deltaPhi)>=cut.at(deltaPhi));
    pass.at(Charge)=true;//(fabs(charge)<0.5);
  }
  //else{std::cout << "no Tau cuts made since no mu passed mu-cuts" << std::endl;}  
  if(verbose)std::cout << "void  FitTest::doEvent() Tau -c " << std::endl;
  if(mu_idx==999 || tau_idx==999) { 
    //std::cout << "Either mu or tau did not pass" << std::endl;
    value.at(deltaPhi)=-1;
    pass.at(deltaPhi)=false;
    value.at(Charge)=-10.0;
    pass.at(Charge)=false;
  }
  //////////////////////////////////////////////////////////////////////////////////
  // QCD Control sample
  /*  if(!pass.at(Charge) && fabs(value.at(Charge))==1){
    if(Ntp->isData()){
      if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
      pass.at(Charge)=true;
    }
  }
  */


  //>>>>>>>>>>>>>>>>>>>>>>>Fill histos for efficiency 
  if(pass.at(EventId)){
    if(Ntp->CheckDecayID(2,5)){
      TLorentzVector TruthA1 = Ntp->GetTruthTauProductLV(5, 20213);
      TLorentzVector TruthMu = Ntp->GetTruthTauProductLV(2, 13);

      TruthA1Pt.at(t).Fill(TruthA1.Pt(),1);
      TruthA1Phi.at(t).Fill(TruthA1.Phi(),1);
      TruthA1Eta.at(t).Fill(TruthA1.Eta(),1);
      TruthMuPt.at(t).Fill(TruthMu.Pt(),1);
      TruthMuPhi.at(t).Fill(TruthMu.Phi(),1);
      TruthMuEta.at(t).Fill(TruthMu.Eta(),1);
    }
  }

  //<<<<<<<<<<<<<<<<<<<<<<<Fill histos for efficiency 

  double wobs(1),w(1);
  if(!Ntp->isData()){
    w*=Ntp->EvtWeight3D();
  }
  else{w=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  
  ///////////////////////////////////////////////////////////
  // Add plots


  if(verbose)std::cout << "void  FitTest::doEvent() -Tau done" << std::endl;
  if(status){
    if(verbose)std::cout << "MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);;
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


      std::cout<<" ID-->  "<< id  <<std::endl;
      std::cout<<" Pi0Size  "<< Ntp->PFTau_PiZeroSize(tau_idx)  <<std::endl;
      if(Ntp->PFTau_PiZeroSize(tau_idx)!=0){std::cout<<" Pi0Size  "<< Ntp->PFTau_PiZeroP4(tau_idx,0).E()  <<std::endl;}


      PhotonEnergyFraction1.at(t).Fill(Ntp->PFJet_photonEnergy(matchedJet)/Ntp->PFJet_p4(matchedJet).E());
      PhotonEnergyFraction2.at(t).Fill(Ntp->PFJet_neutralEmEnergyFraction(matchedJet));

      //<<<<<<<<< search for a matched jet
      
      int idToFill(0);
      if(id ==998)idToFill = 1;
      if(id ==10230833)idToFill = 2;
      if(id ==10231333)idToFill = 3;
      if(id ==10231833)idToFill = 4;
      if(Ntp->PFTau_hpsDecayMode(tau_idx) == 10)idPassHPS.at(t).Fill(idToFill,1);
      std::cout << "PFTau_hpsDecayMode(KFTau_MatchedHPS_idx(tau_idx))  && id " << Ntp->PFTau_hpsDecayMode(tau_idx)<<"  " <<id<<std::endl;
      if(verbose)std::cout << "void  FitTest::doEvent() A" << id<<std::endl;
      std::cout << " id  " << id <<std::endl;
      std::vector<bool> Tau_FitOk; 
      std::vector<TLorentzVector> Tau_sol;
      std::vector<double> Tau_FitChi2;

      std::vector<bool> EventFit_Ok;
      std::vector<TLorentzVector> TauA1_TPF;
      std::vector<double> DeltaPhiSign;

      std::vector<double> deltaphiRotation;
      std::vector<double> SignificanceOfTauRotation;
      std::vector<TLorentzVector> TauA1Neutrino_TPF;
      std::vector<TLorentzVector> TauA1_EF;
      std::vector<TLorentzVector> TauMu_EF;
      std::vector<TLorentzVector> TauMuNetrino_EF;
      std::vector<TLorentzVector> Z_EF;
      std::vector<double> Chi2EventFit;
      std::vector<double> Chi2ProbabilityEventFit;
      std::vector<double> NiterationsEventFit;
      std::vector<double> csum_GF;
      std::vector<double> SigmaEA1;
      std::vector<double> SigmaEMu;
      double flightsignificance;
      double flightsignificanceMVPV;
      TVector3 pv;
      TMatrixTSym<double> PVcov;
      TVector3 sv;
      TMatrixTSym<double> SVcov;
      //  Generate Random value around z mass
      //      double ZM = ran.Gaus(91.18,2.5);


    
      TLorentzVector A1Reco = Ntp->PFTau_a1_lvp(tau_idx).LV();
      TLorentzVector MuReco = Ntp->Muons_p4(mu_idx);
      A1Mass.at(t).Fill(A1Reco.M(),1);
	    //------------------------------ check MC ambiguity 
	     if(Ntp->CheckDecayID(2,5)){

	    TLorentzVector TruthA1 = Ntp->GetTruthTauProductLV(5, 20213);
	    TLorentzVector TruthMu = Ntp->GetTruthTauProductLV(2, 13);
	    
	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);


	    double angle = TruthTauA1.Angle(TruthA1.Vect());
	    double rootAmb=sqrt((TruthA1.M()*TruthA1.M()+TruthA1.P()*TruthA1.P())*(pow(TruthA1.M()*TruthA1.M()-TruthTauA1.M()*TruthTauA1.M(),2)-4*TruthTauA1.M()*TruthTauA1.M()*TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle)));
	    double PtruthMinus  = ((TruthA1.M()*TruthA1.M()+TruthTauA1.M()*TruthTauA1.M() )*TruthA1.P()*cos(angle)  -rootAmb )/(TruthA1.M()*TruthA1.M() + TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle))/2;
	    double PtruthPlus  = ((TruthA1.M()*TruthA1.M()+TruthTauA1.M()*TruthTauA1.M() )*TruthA1.P()*cos(angle)  +rootAmb )/(TruthA1.M()*TruthA1.M() + TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle))/2;

	    int TruthAmibiga =0;


	     

	      if( fabs(TruthTauA1.P()  -PtruthMinus ) < fabs(TruthTauA1.P()  - PtruthPlus))TruthAmibiga =1;
	      if( fabs(TruthTauA1.P()  -PtruthMinus ) > fabs(TruthTauA1.P()  - PtruthPlus))TruthAmibiga =2;
	      CorrectAmbiguityTruth.at(t).Fill(TruthAmibiga,1.);


	     }
	    //------------------------------ check MC ambiguity 


	     std::cout<<"deba 1"<<std::endl;

      for(unsigned int j=0;j< MultiProngTauSolver::NAmbiguity;j++){
	LorentzVectorParticle theTau;
	LorentzVectorParticle theZ;
	std::vector<LorentzVectorParticle> daughter;
	std::vector<LorentzVectorParticle> theZdaughter;

	TLorentzVector NeutrinoA1(0,0,0,0);
	TLorentzVector NeutrinoMu(0,0,0,0);

	TLorentzVector TauA1EventFit(0,0,0,0);
	TLorentzVector TauA1ThreeProngFit(0,0,0,0);
	TLorentzVector TauMuEventFit(0,0,0,0);
	TLorentzVector Z_sol(0,0,0,0);
	double deltaphi(0);

	double LC_Eventchi2(0);
	double LC_Eventchi2Probability(0);
	double LC_chi2(0);
	double phisign(0);
	int NiterationsEF(0);
	double csum(0);
	double SEA1(0);
	double SEMU(0);
	//if(verbose)std::cout << "void  FitTest::doEvent() A1" << std::endl;
	bool  ThreeProngFitSuccess =false;
	bool  EventFitSuccess =false;
	ThreeProngFitSuccess=Ntp->ThreeProngTauFit(tau_idx,j,theTau,daughter,LC_chi2,phisign);
	MeasuredPhi.at(t).Fill(Ntp->MeasuredTauDirection(tau_idx).Phi(),1);
 	if(ThreeProngFitSuccess){


	  std::cout<<" checkTauCovariance "<< theTau.Covariance(3,3)<<std::endl;
	  std::cout<<" checkTauCovariance "<< theTau.Covariance(4,4)<<std::endl;
	  std::cout<<" checkTauCovariance "<< theTau.Covariance(5,5)<<std::endl;


	  int idToFill(0);
	  if(id ==998)idToFill = 1;
	  if(id ==10230833)idToFill = 2;
	  if(id ==10231433)idToFill = 3;
	  if(id ==10231833)idToFill = 4;

	  if(j==1)idPassKFMinus.at(t).Fill(idToFill,1);
	  if(j==2)idPassKFPlus.at(t).Fill(idToFill,1);
     
	  if(j==1) KFTauPhi_Minus.at(t).Fill(theTau.LV().Phi(),1);
	  if(j==2) KFTauPhi_Plus.at(t).Fill(theTau.LV().Phi(),1);

	  //>>>>>>>>>>> Fill Resol plots dor TPF
	  if(Ntp->CheckDecayID(2,5)){
	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	    if(j==1){
	      TauA1PtResolKFM.at(t).Fill(TruthTauA1.Pt() - theTau.LV().Pt(),1);
	      TauA1PhiResolKFM.at(t).Fill(TruthTauA1.Phi() - theTau.LV().Phi(),1);
	      TauA1EtaResolKFM.at(t).Fill(TruthTauA1.Eta() - theTau.LV().Eta(),1);
	    }
	    if(j==2){
	      TauA1PtResolKFP.at(t).Fill(TruthTauA1.Pt() -  theTau.LV().Pt(),1);
	      TauA1PhiResolKFP.at(t).Fill(TruthTauA1.Phi() - theTau.LV().Phi(),1);
	      TauA1EtaResolKFP.at(t).Fill(TruthTauA1.Eta() - theTau.LV().Eta(),1);
	    }
	  }

	  //<<<<<<<<<<< Fill Resol plots dor TPF
	  //std::cout<<"significance  "<<Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(tau_idx),Ntp->PFTau_TIP_primaryVertex_cov(tau_idx),theTau.Vertex(),theTau.VertexCov()) <<std::endl;
	  flightsignificance=Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(tau_idx),Ntp->PFTau_TIP_primaryVertex_cov(tau_idx),theTau.Vertex(),theTau.VertexCov());
	  pvsvsignificance.at(t).Fill(Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(tau_idx),Ntp->PFTau_TIP_primaryVertex_cov(tau_idx),theTau.Vertex(),theTau.VertexCov()),1);
	  significan=Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(tau_idx),Ntp->PFTau_TIP_primaryVertex_cov(tau_idx),theTau.Vertex(),theTau.VertexCov());
	  NeutrinoA1=daughter.at(1).LV();
	  
	  std::cout<<"debosa 1"<<std::endl;
	  TauA1ThreeProngFit=theTau.LV();
	  EventFitSuccess = Ntp->EventFit(tau_idx,mu_idx,theTau,theZ,theZdaughter,LC_Eventchi2,NiterationsEF,csum);

	  if(EventFitSuccess){
	    std::cout<<"debosa 3"<<std::endl;
	  LC_Eventchi2Probability = TMath::Prob(LC_Eventchi2,1);

	  TauA1EventFit =theZdaughter.at(0).LV();
	  TauMuEventFit =theZdaughter.at(1).LV();

	  std::cout<<"debosa 4"<<std::endl;




	  TVector3 pv=Ntp->PFTau_TIP_primaryVertex_pos(tau_idx);
	  TMatrixTSym<double> pvcov=Ntp->PFTau_TIP_primaryVertex_cov(tau_idx);

	  phisign=Ntp->DeltaPhiSignificance(pv.X(),pv.Y(),pvcov(0,0),pvcov(1,1),theTau.Vertex().X(),theTau.Vertex().Y(), theTau.VertexCov()(0,0),theTau.VertexCov()(1,1));

	  deltaphi = atan2(theTau.Vertex().Y() - pv.Y(),theTau.Vertex().X() - pv.X() ) - theTau.LV().Phi();


	  NeutrinoMu=TauMuEventFit-Ntp->Muons_p4(mu_idx);



	  flightsignificanceMVPV=Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(tau_idx),Ntp->PFTau_TIP_primaryVertex_cov(tau_idx),theZdaughter.at(1).Vertex(),theZdaughter.at(1).VertexCov());
	  pvMusvsignificance.at(t).Fill(flightsignificanceMVPV,1);
	  std::cout<<"theZdaughter.at(1).VertexCov()-----"<<std::endl;
	  for(int str =0; str <theZdaughter.at(1).VertexCov().GetNrows(); str++){
	    for(int kol =0; kol < theZdaughter.at(1).VertexCov().GetNcols(); kol++){
	      std::cout<<"  "<< theZdaughter.at(1).VertexCov()(str,kol)<<"  ";
	    }    
	    std::cout<<std::endl;
	  }



	  Z_sol=theZ.LV();
	  SEA1 = sqrt((pow(TauA1EventFit.Px(),2)*theZdaughter.at(0).Covariance(3,3)+pow(TauA1EventFit.Py(),2)*theZdaughter.at(0).Covariance(4,4)+ pow(TauA1EventFit.Pz(),2)*theZdaughter.at(0).Covariance(5,5))/TauA1EventFit.P()/TauA1EventFit.P());
	  SEMU= sqrt((pow(TauMuEventFit.Px(),2)*theZdaughter.at(1).Covariance(3,3)+pow(TauMuEventFit.Py(),2)*theZdaughter.at(1).Covariance(4,4)+ pow(TauMuEventFit.Pz(),2)*theZdaughter.at(1).Covariance(5,5))/TauMuEventFit.P()/TauMuEventFit.P());
	  }

	}
	SignificanceOfTauRotation.push_back(phisign);
	TauA1Neutrino_TPF.push_back(NeutrinoA1);
	TauMuNetrino_EF.push_back(NeutrinoMu);

	Tau_FitOk.push_back(ThreeProngFitSuccess);
	Tau_FitChi2.push_back(LC_chi2);
	Chi2ProbabilityEventFit.push_back(LC_Eventchi2Probability);
	EventFit_Ok.push_back(EventFitSuccess);
  	Z_EF.push_back(Z_sol);
 	Chi2EventFit.push_back(LC_Eventchi2);
 	TauA1_EF.push_back(TauA1EventFit);
 	TauMu_EF.push_back(TauMuEventFit);
	NiterationsEventFit.push_back(NiterationsEF);
	csum_GF.push_back(csum);
	DeltaPhiSign.push_back(phisign);   
	SigmaEA1.push_back(SEA1);
	SigmaEMu.push_back(SEMU);
	deltaphiRotation.push_back(deltaphi);
      }
      std::cout<<"sizes "<< Tau_FitOk.size() <<"  "<<EventFit_Ok.size()<<std::endl;
      if(Tau_FitOk.size()==3)	std::cout<<" Tau_FitOk "<<Tau_FitOk.at(0)<<" "<< Tau_FitOk.at(1)<<"  "<<Tau_FitOk.at(2)<<std::endl;
      if(EventFit_Ok.size()==3){	std::cout<<" EventFit_Ok "<<EventFit_Ok.at(0)<<" "<< EventFit_Ok.at(1)<<"  "<<EventFit_Ok.at(2)<<std::endl;
	std::cout<<" SEMU "<<SigmaEMu.at(0)<<" "<<SigmaEMu.at(1) <<"  "<<SigmaEMu.at(2)<<std::endl;}
      std::cout<<" SignificanceOfTauRotation size "<<SignificanceOfTauRotation.size()<<std::endl;
      if(SignificanceOfTauRotation.size()==3)	std::cout<<" SignificanceOfTauRotation "<<SignificanceOfTauRotation.at(0)<<" "<<SignificanceOfTauRotation.at(1) <<"  "<<SignificanceOfTauRotation.at(2)<<std::endl;
      PhiRotationSignificancePhysicalTaus.at(t).Fill(SignificanceOfTauRotation.at(1));
      PhiRotationSignificanceUnPhysicalTaus.at(t).Fill(SignificanceOfTauRotation.at(0));
    
      if(Tau_FitOk.at(0) && Tau_FitOk.at(1) && Tau_FitOk.at(2)){
	if(EventFit_Ok.at(0) && EventFit_Ok.at(1) && EventFit_Ok.at(2)){
	  if(Ntp->CheckDecayID(2,5)){
	    
	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);

	    //   	      std::cout<<"check truth taus are available  "<<TruthTauA1.Pt()<<"  "<<TruthTauMu.Pt()<<std::endl;
	    //   	      std::cout<<"TauMuKinematics 0 " <<TauMu_EF.at(0).Px()<<"  "<<TauMu_EF.at(0).Py()<<"  "<<TauMu_EF.at(0).Pz()<<"  "<<TauMu_EF.at(0).M() <<"   " <<TauMu_EF.at(0).Pt() <<std::endl;
	    //   	      std::cout<<"TauMuKinematics 1 " <<TauMu_EF.at(1).Px()<<"  "<<TauMu_EF.at(1).Py()<<"  "<<TauMu_EF.at(1).Pz()<<"  "<<TauMu_EF.at(1).M() <<"   " << TauMu_EF.at(1).Pt() <<std::endl;
	    //   	      std::cout<<"TauMuKinematics 2 " <<TauMu_EF.at(2).Px()<<"  "<<TauMu_EF.at(2).Py()<<"  "<<TauMu_EF.at(2).Pz()<<"  "<<TauMu_EF.at(2).M() <<"   " <<TauMu_EF.at(2).Pt()  <<std::endl;

	    //   	      std::cout<<"TauA1Kinematics 0 " <<TauA1_EF.at(0).Pt()<<"  "<<TauMu_EF.at(0).P()<<std::endl;
	    //   	      std::cout<<"TauA1Kinematics 1 " <<TauA1_EF.at(1).Pt()<<"  "<<TauMu_EF.at(1).P()<<std::endl;
	    //   	      std::cout<<"TauA1Kinematics 2 " <<TauA1_EF.at(2).Pt()<<"  "<<TauMu_EF.at(2).P()<<std::endl;

	    

	    //   	      std::cout<<"Z massa"<<Z_EF.at(0).M()<<"  "<<Z_EF.at(1).M()<<"  "<<Z_EF.at(2).M()<<std::endl;
	    
	    ZMass0.at(t).Fill(Z_EF.at(0).M(),1);
	    ZMass1.at(t).Fill(Z_EF.at(1).M(),1);
	    ZMass2.at(t).Fill(Z_EF.at(2).M(),1);

	    TauA1PtResol0.at(t).Fill(TauA1_EF.at(0).Pt() - TruthTauA1.Pt(),1);
	    TauA1PtResol1.at(t).Fill(TauA1_EF.at(1).Pt() - TruthTauA1.Pt(),1);
	    TauA1PtResol2.at(t).Fill(TauA1_EF.at(2).Pt() - TruthTauA1.Pt(),1);

	    TauMuPtResol0.at(t).Fill(TauMu_EF.at(0).Pt() - TruthTauMu.Pt(),1); 
	    TauMuPtResol1.at(t).Fill(TauMu_EF.at(1).Pt() - TruthTauMu.Pt(),1);
	    TauMuPtResol2.at(t).Fill(TauMu_EF.at(2).Pt() - TruthTauMu.Pt(),1);
	    
	    Niter0.at(t).Fill(NiterationsEventFit.at(0),1);
	    Niter1.at(t).Fill(NiterationsEventFit.at(1),1);
	    Niter2.at(t).Fill(NiterationsEventFit.at(2),1);

	    TauA1PhiResol0.at(t).Fill(TauA1_EF.at(0).Phi() - TruthTauA1.Phi(),1);
	    TauA1PhiResol1.at(t).Fill(TauA1_EF.at(1).Phi() - TruthTauA1.Phi(),1);
	    TauA1PhiResol2.at(t).Fill(TauA1_EF.at(2).Phi() - TruthTauA1.Phi(),1);

	    TauMuPhiResol0.at(t).Fill(TauMu_EF.at(0).Phi() - TruthTauMu.Phi(),1);
	    TauMuPhiResol1.at(t).Fill(TauMu_EF.at(1).Phi() - TruthTauMu.Phi(),1);
	    TauMuPhiResol2.at(t).Fill(TauMu_EF.at(2).Phi() - TruthTauMu.Phi(),1);
	      
	    TauA1EtaResol0.at(t).Fill(TauA1_EF.at(0).Eta() - TruthTauA1.Eta(),1);
	    TauA1EtaResol1.at(t).Fill(TauA1_EF.at(1).Eta() - TruthTauA1.Eta(),1);
	    TauA1EtaResol2.at(t).Fill(TauA1_EF.at(2).Eta() - TruthTauA1.Eta(),1);

	    TauMuEtaResol0.at(t).Fill(TauMu_EF.at(0).Eta() - TruthTauMu.Eta(),1);
	    TauMuEtaResol1.at(t).Fill(TauMu_EF.at(1).Eta() - TruthTauMu.Eta(),1);
	    TauMuEtaResol2.at(t).Fill(TauMu_EF.at(2).Eta() - TruthTauMu.Eta(),1);

	    
	    std::cout<<"TauMu  " <<TauMu_EF.at(0).Pt() - TruthTauMu.Pt()<<std::endl;
	    std::cout<<"TauMu  " <<TauMu_EF.at(1).Pt() - TruthTauMu.Pt()<<std::endl;
	    std::cout<<"TauMu  " <<TauMu_EF.at(2).Pt() - TruthTauMu.Pt()<<std::endl;

	    std::cout<<"TauA1  " <<TauA1_EF.at(0).Pt() - TruthTauA1.Pt()<<std::endl;
	    std::cout<<"TauA1  " <<TauA1_EF.at(1).Pt() - TruthTauA1.Pt()<<std::endl;
	    std::cout<<"TauA1  " <<TauA1_EF.at(2).Pt() - TruthTauA1.Pt()<<std::endl;
	    std::cout<<"---------  deb 1   "<<std::endl;
	    Probability0.at(t).Fill(TMath::Prob(Chi2EventFit.at(0),1),1);
	    Probability1.at(t).Fill(TMath::Prob(Chi2EventFit.at(1),1),1);
	      Probability2.at(t).Fill(TMath::Prob(Chi2EventFit.at(2),1),1);
	      

	      TauMuEnergyPull0.at(t).Fill((TauMu_EF.at(0).P() - TruthTauMu.P())/SigmaEMu.at(0),1);
	      TauA1EnergyPull0.at(t).Fill((TauA1_EF.at(0).P() - TruthTauMu.P())/SigmaEA1.at(0),1);

	      TauMuEnergyPull1.at(t).Fill((TauMu_EF.at(1).P() - TruthTauMu.P())/SigmaEMu.at(0),1);
	      TauA1EnergyPull1.at(t).Fill((TauA1_EF.at(1).P() - TruthTauMu.P())/SigmaEA1.at(0),1);

	      TauMuEnergyPull2.at(t).Fill((TauMu_EF.at(2).P() - TruthTauMu.P())/SigmaEMu.at(0),1);
	      TauA1EnergyPull2.at(t).Fill((TauA1_EF.at(2).P() - TruthTauMu.P())/SigmaEA1.at(0),1);
	      std::cout<<"---------  deb 2   "<<std::endl;

	      //----------------
	      //solve an ambiguity
	      int CorrectAmbiguity; int IncorrectAmbiguity;
	      if(fabs(TauA1_EF.at(1).E() - TruthTauA1.E()) < fabs(TauA1_EF.at(2).E() - TruthTauA1.E())){  CorrectAmbiguity=1; IncorrectAmbiguity =2;}
	      if(fabs(TauA1_EF.at(2).E() - TruthTauA1.E()) < fabs(TauA1_EF.at(1).E() - TruthTauA1.E())){  CorrectAmbiguity=2;IncorrectAmbiguity =1; }
	      //----------------
	      std::cout<<"---------  deb 3   "<<Chi2EventFit.size() <<std::endl;
	      Chi2Dim.at(t).Fill(TMath::Prob(Chi2EventFit.at(CorrectAmbiguity),1),TMath::Prob(Chi2EventFit.at(IncorrectAmbiguity),1),1);

	      csum2Dim.at(t).Fill(csum_GF.at(CorrectAmbiguity),csum_GF.at(IncorrectAmbiguity),1);

	      std::cout<<"---------  deb 4   "<<std::endl;
	      //  std::cout<<"---------  prob 1   "<<TMath::Prob(Chi2EventFit.at(1),1)<< " prob 2  " <<TMath::Prob(Chi2EventFit.at(2),1) << std::endl;

	      Chi2DimPlusMinus.at(t).Fill(TMath::Prob(Chi2EventFit.at(1),1),TMath::Prob(Chi2EventFit.at(2),1),1);
	      std::cout<<"---------  deb 5   "<<std::endl;

	     }
 
	  }
	}

    

	int AmbiguitySolution =0;
	bool AmbigPoint;
	if(Ntp->AmbiguitySolver(Tau_FitOk,EventFit_Ok,Chi2ProbabilityEventFit,AmbiguitySolution, AmbigPoint)){
	  isInPhysicalRegion=1;
	  std::cout<<"------------ picked solution <<  "<<AmbiguitySolution<<std::endl;


	  int idToFill(0);
	  if(id ==998)idToFill = 1;
	  if(id ==10230833)idToFill = 2;
	  if(id ==10231433)idToFill = 3;
	  if(id ==10231833)idToFill = 4;

	  idPassGF.at(t).Fill(idToFill,1);


	  if(Ntp->CheckDecayID(2,5)){

	    deltaPhiSignificance.at(t).Fill(DeltaPhiSign.at(AmbiguitySolution),1);
	    TLorentzVector TruthA1 = Ntp->GetTruthTauProductLV(5, 20213);
	    TLorentzVector TruthMu = Ntp->GetTruthTauProductLV(2, 13);
	    
	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	    std::cout<<"deltaR Mu  "<< TruthMu.DeltaR(MuReco) <<std::endl;
	    std::cout<<"deltaR A1 "<< TruthA1.DeltaR(A1Reco) <<std::endl;

	    KFTauPhi_Ambiga.at(t).Fill(TauA1_EF.at(AmbiguitySolution).Phi(),1);
	    if(DeltaPhiSign.at(AmbiguitySolution) < 0.1)KFTauPhi_Ambiga_DeltaPhiSignificanceCut.at(t).Fill(TauA1_EF.at(AmbiguitySolution).Phi(),1);


	    KFTauPhi_Ambiga_vsDeltaPhiSignificance.at(t).Fill(TauA1_EF.at(AmbiguitySolution).Phi(),DeltaPhiSign.at(AmbiguitySolution),1);
	    KFTauPhi_Ambiga_vsDeltaPhi.at(t).Fill(TauA1_EF.at(AmbiguitySolution).Phi(),deltaphiRotation.at(AmbiguitySolution),1);

	    TruthA1PtAfter.at(t).Fill(TruthA1.Pt(),1);
	    TruthA1PhiAfter.at(t).Fill(TruthA1.Phi(),1);
	    TruthA1EtaAfter.at(t).Fill(TruthA1.Eta(),1);
	    TruthMuPtAfter.at(t).Fill(TruthMu.Pt(),1);
	    TruthMuPhiAfter.at(t).Fill(TruthMu.Phi(),1);
	    TruthMuEtaAfter.at(t).Fill(TruthMu.Eta(),1);

	    TLorentzVector TZ1 = Ntp->MCSignalParticle_p4(0);
	    TLorentzVector TZ2 = Ntp->MCSignalParticle_p4(1);

	    ZPtResol.at(t).Fill((TauA1_EF.at(AmbiguitySolution) + TauMu_EF.at(AmbiguitySolution)).Pt() - TZ1.Pt(),1);
	    ZPzResolRelative.at(t).Fill(((TauA1_EF.at(AmbiguitySolution) + TauMu_EF.at(AmbiguitySolution)).Pz() - TZ1.Pz())/TZ1.Pz(),1);
	    ZEResolRelative.at(t).Fill(((TauA1_EF.at(AmbiguitySolution) + TauMu_EF.at(AmbiguitySolution)).E() - TZ1.E())/TZ1.E(),1);
	    ZPtTruth.at(t).Fill(TZ1.Pt(),1);

	    std::cout<<" <><><>  TruthTauMu Theta   "<<TruthTauMu.Theta()<< " Truth TauMu Pz  " <<TruthTauMu.Pz() <<std::endl;

	    std::cout<<"SP poca 1  (x,y,z)  "<<Ntp->MCSignalParticle_Poca(0).X()  << "    " <<Ntp->MCSignalParticle_Poca(0).Y()  << "  " << Ntp->MCSignalParticle_Poca(0).Z() <<std::endl;
	    std::cout<<"SP poca 2  (x,y,z)  "<<Ntp->MCSignalParticle_Poca(1).X()  << "    " <<Ntp->MCSignalParticle_Poca(1).Y()  << "  " << Ntp->MCSignalParticle_Poca(1).Z() <<std::endl;



	    TauMuEResolVsSignificance.at(t).Fill(TruthTauMu.E() - TauMu_EF.at(AmbiguitySolution).E(),flightsignificance,1);
	    TauMuEResolVsSignificance2.at(t).Fill(TruthTauMu.E() - TauMu_EF.at(AmbiguitySolution).E(),flightsignificanceMVPV,1);
	    TauMuEResolVsChi2Probability.at(t).Fill(TruthTauMu.E() - TauMu_EF.at(AmbiguitySolution).E(),Chi2ProbabilityEventFit.at(AmbiguitySolution),1);
	    ZPtTruthTauMuPtResolution.at(t).Fill(TZ1.Pt(),TruthTauMu.Pt() -TauMu_EF.at(AmbiguitySolution).Pt(),1);
	    
	    if(flightsignificance > 1  && Chi2ProbabilityEventFit.at(AmbiguitySolution) > 0.04){
	      TruthA1PtAfterQC.at(t).Fill(TruthA1.Pt(),1);
	      TruthA1PhiAfterQC.at(t).Fill(TruthA1.Phi(),1);
	      TruthA1EtaAfterQC.at(t).Fill(TruthA1.Eta(),1);
	      TruthMuPtAfterQC.at(t).Fill(TruthMu.Pt(),1);
	      TruthMuPhiAfterQC.at(t).Fill(TruthMu.Phi(),1);
	      TruthMuEtaAfterQC.at(t).Fill(TruthMu.Eta(),1);
	    }

	    
	    Niter.at(t).Fill(NiterationsEventFit.at(AmbiguitySolution),1);
	    Probability.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),1);
	    csumGF.at(t).Fill(csum_GF.at(AmbiguitySolution),1);



	    TauA1PhiResol.at(t).Fill(TruthTauA1.Phi() - TauA1_EF.at(AmbiguitySolution).Phi(),1);
	    TauMuPhiResol.at(t).Fill(TruthTauMu.Phi() -TauMu_EF.at(AmbiguitySolution).Phi(),1);
	    
	    TauA1EtaResol.at(t).Fill(TruthTauA1.Eta() - TauA1_EF.at(AmbiguitySolution).Eta(),1);
	    TauMuEtaResol.at(t).Fill(TruthTauMu.Eta() -TauMu_EF.at(AmbiguitySolution).Eta(),1);
	    
	    TauA1PtResol.at(t).Fill(TruthTauA1.Pt() - TauA1_EF.at(AmbiguitySolution).Pt(),1);
	    TauMuPtResol.at(t).Fill(TruthTauMu.Pt() -TauMu_EF.at(AmbiguitySolution).Pt(),1);
	    
	    TauA1EResol.at(t).Fill(TruthTauA1.E() - TauA1_EF.at(AmbiguitySolution).E(),1);
	    TauMuEResol.at(t).Fill(TruthTauMu.E() - TauMu_EF.at(AmbiguitySolution).E(),1);
		




	    if(TZ1.Pt() < 5){
	      TauA1PhiResolZPtCut.at(t).Fill(TruthTauA1.Phi() - TauA1_EF.at(AmbiguitySolution).Phi(),1);
	      TauMuPhiResolZPtCut.at(t).Fill(TruthTauMu.Phi() -TauMu_EF.at(AmbiguitySolution).Phi(),1);
	      
	      TauA1EtaResolZPtCut.at(t).Fill(TruthTauA1.Eta() - TauA1_EF.at(AmbiguitySolution).Eta(),1);
	      TauMuEtaResolZPtCut.at(t).Fill(TruthTauMu.Eta() -TauMu_EF.at(AmbiguitySolution).Eta(),1);
	    
	      TauA1PtResolZPtCut.at(t).Fill(TruthTauA1.Pt() - TauA1_EF.at(AmbiguitySolution).Pt(),1);
	      TauMuPtResolZPtCut.at(t).Fill(TruthTauMu.Pt() -TauMu_EF.at(AmbiguitySolution).Pt(),1);
	      
	      TauA1EResolZPtCut.at(t).Fill(TruthTauA1.E() - TauA1_EF.at(AmbiguitySolution).E(),1);
	      TauMuEResolZPtCut.at(t).Fill(TruthTauMu.E() - TauMu_EF.at(AmbiguitySolution).E(),1);

	      ZPtResolZPtCut.at(t).Fill((TauA1_EF.at(AmbiguitySolution) + TauMu_EF.at(AmbiguitySolution)).Pt() - TZ1.Pt(),1);
	      ZPzResolRelativeZPtCut.at(t).Fill(((TauA1_EF.at(AmbiguitySolution) + TauMu_EF.at(AmbiguitySolution)).Pz() - TZ1.Pz())/TZ1.Pz(),1);
	      ZEResolRelativeZPtCut.at(t).Fill(((TauA1_EF.at(AmbiguitySolution) + TauMu_EF.at(AmbiguitySolution)).E() - TZ1.E())/TZ1.E(),1);


	    }



	    if(flightsignificance > 2){
	      NiterQC.at(t).Fill(NiterationsEventFit.at(AmbiguitySolution),1);
	      ProbabilityQC.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),1);
	      
	      TauA1PhiResolQC.at(t).Fill(TruthTauA1.Phi() - TauA1_EF.at(AmbiguitySolution).Phi(),1);
	      TauMuPhiResolQC.at(t).Fill(TruthTauMu.Phi() -TauMu_EF.at(AmbiguitySolution).Phi(),1);
	    
	      TauA1EtaResolQC.at(t).Fill(TruthTauA1.Eta() - TauA1_EF.at(AmbiguitySolution).Eta(),1);
	      TauMuEtaResolQC.at(t).Fill(TruthTauMu.Eta() -TauMu_EF.at(AmbiguitySolution).Eta(),1);
	    
	      TauA1PtResolQC.at(t).Fill(TruthTauA1.Pt() - TauA1_EF.at(AmbiguitySolution).Pt(),1);
	      TauMuPtResolQC.at(t).Fill(TruthTauMu.Pt() -TauMu_EF.at(AmbiguitySolution).Pt(),1);
	    
	      TauA1EResolQC.at(t).Fill(TruthTauA1.E() - TauA1_EF.at(AmbiguitySolution).E(),1);
	      TauMuEResolQC.at(t).Fill(TruthTauMu.E() - TauMu_EF.at(AmbiguitySolution).E(),1);


	      TauA1PzResolQC.at(t).Fill(TruthTauA1.Pz() - TauA1_EF.at(AmbiguitySolution).Pz(),1);
	      TauMuPzResolQC.at(t).Fill(TruthTauMu.Pz() - TauMu_EF.at(AmbiguitySolution).Pz(),1);
	    }
	    ZMass.at(t).Fill(Z_EF.at(AmbiguitySolution).M(),1);

	    
	    if( EventFit_Ok.at(1) && !EventFit_Ok.at(2))AmbiguityEfficiency.at(t).Fill(1,1);
	    if( EventFit_Ok.at(2) && !EventFit_Ok.at(1))AmbiguityEfficiency.at(t).Fill(2,1);
	    if( EventFit_Ok.at(1) &&  EventFit_Ok.at(2)){AmbiguityEfficiency.at(t).Fill(2,1); AmbiguityEfficiency.at(t).Fill(1,1);}


	    EtaTauMuEtaTauA1.at(t).Fill(TauA1_EF.at(AmbiguitySolution).Eta(),TauMu_EF.at(AmbiguitySolution).Eta(),1);


	    ResTauMuPtEtaTauMu.at(t).Fill(TruthTauMu.Pt() -TauMu_EF.at(AmbiguitySolution).Pt(),TauMu_EF.at(AmbiguitySolution).Eta(),1);
	    ResTauMuEEtaTauMu.at(t).Fill(TruthTauMu.E() -TauMu_EF.at(AmbiguitySolution).E(),TauMu_EF.at(AmbiguitySolution).Eta(),1);
	    ResTauMuPtPtTauMu.at(t).Fill(TruthTauMu.Pt() -TauMu_EF.at(AmbiguitySolution).Pt(),TauMu_EF.at(AmbiguitySolution).Pt(),1);


	    ResTauA1PtEtaTauA1.at(t).Fill(TruthTauA1.Pt() -TauA1_EF.at(AmbiguitySolution).Pt(),TauA1_EF.at(AmbiguitySolution).Eta(),1);
	    ResTauA1EEtaTauA1.at(t).Fill(TruthTauA1.E() -TauA1_EF.at(AmbiguitySolution).E(),TauA1_EF.at(AmbiguitySolution).Eta(),1);
	    ResTauA1PtPtTauA1.at(t).Fill(TruthTauA1.Pt() -TauA1_EF.at(AmbiguitySolution).Pt(),TauA1_EF.at(AmbiguitySolution).Pt(),1);

	   
	    ResTauMuPtProb.at(t).Fill(TruthTauMu.Pt() -TauMu_EF.at(AmbiguitySolution).Pt(),Chi2ProbabilityEventFit.at(AmbiguitySolution),1);
	    ResTauMuEProb.at(t).Fill(TruthTauMu.E() -TauMu_EF.at(AmbiguitySolution).E(),Chi2ProbabilityEventFit.at(AmbiguitySolution),1);
	    ResTauA1PtProb.at(t).Fill(TruthTauA1.Pt() -TauA1_EF.at(AmbiguitySolution).Pt(),Chi2ProbabilityEventFit.at(AmbiguitySolution),1);
	    ResTauA1EProb.at(t).Fill(TruthTauA1.E() -TauA1_EF.at(AmbiguitySolution).E(),Chi2ProbabilityEventFit.at(AmbiguitySolution),1);
  



	    
	    if( EventFit_Ok.at(1) && EventFit_Ok.at(2)){
	      
	      int CorrectAmbiguity; int IncorrectAmbiguity;
	      if(fabs(TauA1_EF.at(1).E() - TruthTauA1.E()) < fabs(TauA1_EF.at(2).E() - TruthTauA1.E())){  CorrectAmbiguity=1; IncorrectAmbiguity =2;}
	      if(fabs(TauA1_EF.at(2).E() - TruthTauA1.E()) < fabs(TauA1_EF.at(1).E() - TruthTauA1.E())){  CorrectAmbiguity=2;IncorrectAmbiguity =1; }
	      
	      ProbabilityOfCorrectZPt.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),TZ1.Pt(),1);
	      ProbabilityOfCorrectPtBalance.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),TauA1_EF.at(CorrectAmbiguity).Pt() - TauMu_EF.at(CorrectAmbiguity).Pt(),1);

	      std::cout<<" pt of mu   "<< TauMu_EF.at(CorrectAmbiguity).Pt() <<"   pt of a1 "<<TauA1_EF.at(CorrectAmbiguity).Pt()<<" diff  "<<TauA1_EF.at(CorrectAmbiguity).Pt() - TauMu_EF.at(CorrectAmbiguity).Pt()<<std::endl;

	      csum2Prob.at(t).Fill(csum_GF.at(CorrectAmbiguity),Chi2ProbabilityEventFit.at(CorrectAmbiguity),1);
	      ProbabilityOfCorrect.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),1);

	      TauA1PhiResolCorrectAmbiga.at(t).Fill(TruthTauA1.Phi() - TauA1_EF.at(CorrectAmbiguity).Phi(),1);
	      TauMuPhiResolCorrectAmbiga.at(t).Fill(TruthTauMu.Phi() -TauMu_EF.at(CorrectAmbiguity).Phi(),1);
	    
	      TauA1EtaResolCorrectAmbiga.at(t).Fill(TruthTauA1.Eta() - TauA1_EF.at(CorrectAmbiguity).Eta(),1);
	      TauMuEtaResolCorrectAmbiga.at(t).Fill(TruthTauMu.Eta() -TauMu_EF.at(CorrectAmbiguity).Eta(),1);
		    
	      TauA1PtResolCorrectAmbiga.at(t).Fill(TruthTauA1.Pt() - TauA1_EF.at(CorrectAmbiguity).Pt(),1);
	      TauMuPtResolCorrectAmbiga.at(t).Fill(TruthTauMu.Pt() -TauMu_EF.at(CorrectAmbiguity).Pt(),1);
		    
	      TauA1EResolCorrectAmbiga.at(t).Fill(TruthTauA1.E() - TauA1_EF.at(CorrectAmbiguity).E(),1);
	      TauMuEResolCorrectAmbiga.at(t).Fill(TruthTauMu.E() - TauMu_EF.at(CorrectAmbiguity).E(),1);





	      TauA1PhiResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),TruthTauA1.Phi() - TauA1_EF.at(AmbiguitySolution).Phi(),1);
	      TauMuPhiResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),TruthTauMu.Phi() -TauMu_EF.at(AmbiguitySolution).Phi(),1);
	    
	      TauA1EtaResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),TruthTauA1.Eta() - TauA1_EF.at(AmbiguitySolution).Eta(),1);
	      TauMuEtaResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),TruthTauMu.Eta() -TauMu_EF.at(AmbiguitySolution).Eta(),1);
	    
	      TauA1PtResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),TruthTauA1.Pt() - TauA1_EF.at(AmbiguitySolution).Pt(),1);
	      TauMuPtResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),TruthTauMu.Pt() -TauMu_EF.at(AmbiguitySolution).Pt(),1);
	    
	      TauA1EResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),TruthTauA1.E() - TauA1_EF.at(AmbiguitySolution).E(),1);
	      TauMuEResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),TruthTauMu.E() - TauMu_EF.at(AmbiguitySolution).E(),1);
		






	      MissingEnergyCheckCorrectAmbiga.at(t).Fill( (TauA1Neutrino_TPF.at(CorrectAmbiguity) + TauMuNetrino_EF.at(CorrectAmbiguity)).Et(),Ntp->MET_Uncorr_et(),1);
	      MissingEnergyCheckInCorrectAmbiga.at(t).Fill((TauA1Neutrino_TPF.at(IncorrectAmbiguity) + TauMuNetrino_EF.at(IncorrectAmbiguity)).Et(),Ntp->MET_Uncorr_et(),1);

// 	      std::cout<<" proverka 1 " << TauA1Neutrino_TPF.at(CorrectAmbiguity).Pt() << std::endl;
// 	      std::cout<<" proverka 2 " << TauMuNetrino_EF.at(CorrectAmbiguity).Pt() << std::endl;

	      Ambiguity.at(t).Fill(AmbiguitySolution,1);
	      if(AmbiguitySolution==CorrectAmbiguity){
		AmbiguityIsCorrect.at(t).Fill(1.,1.);
		RightlyTakenAmbiguity.at(t).Fill(AmbiguitySolution,1.);

		deltaAmbigaEnergyForCorreAmb.at(t).Fill(TauA1_EF.at(2).E() - TauA1_EF.at(1).E());



	      }
	      if(AmbiguitySolution!=CorrectAmbiguity){
		deltaAmbigaEnergyForWrongAmb.at(t).Fill(TauA1_EF.at(2).E() - TauA1_EF.at(1).E());
		WronglyTakenAmbiguity.at(t).Fill(AmbiguitySolution,1.);
		AmbiguityIsCorrect.at(t).Fill(0.,1.);
	      }

	      // std::cout<<"TruthAmibiga  " <<TruthAmibiga << " CorrectAmbiguity "<< CorrectAmbiguity<<std::endl;



	    }
	  }
	  
	}
    
	//------------- Fill plots for ambiguity zero 
	if(EventFit_Ok.at(0) &&  (!EventFit_Ok.at(1) && !EventFit_Ok.at(2)) ){
	  isInPhysicalRegion =0;
	  if(Ntp->CheckDecayID(2,5)){
   
	    TLorentzVector TruthA1 = Ntp->GetTruthTauProductLV(5, 20213);
	    TLorentzVector TruthMu = Ntp->GetTruthTauProductLV(2, 13);
	    
	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	    
	    int idToFill(0);
	    if(id ==998)idToFill = 1;
	    if(id ==10230833)idToFill = 2;
	    if(id ==10231433)idToFill = 3;
	    if(id ==10231833)idToFill = 4;
	    idPassGFZero.at(t).Fill(idToFill,1);

	    TauA1PhiResolZero.at(t).Fill(TruthTauA1.Phi() - TauA1_EF.at(0).Phi(),1);
	    TauMuPhiResolZero.at(t).Fill(TruthTauMu.Phi() -TauMu_EF.at(0).Phi(),1);
	    
	    TauA1EtaResolZero.at(t).Fill(TruthTauA1.Eta() - TauA1_EF.at(0).Eta(),1);
	    TauMuEtaResolZero.at(t).Fill(TruthTauMu.Eta() -TauMu_EF.at(0).Eta(),1);
	    
	    TauA1PtResolZero.at(t).Fill(TruthTauA1.Pt() - TauA1_EF.at(0).Pt(),1);
	    TauMuPtResolZero.at(t).Fill(TruthTauMu.Pt() -TauMu_EF.at(0).Pt(),1);
	    
	    TauA1EResolZero.at(t).Fill(TruthTauA1.E() - TauA1_EF.at(0).E(),1);
	    TauMuEResolZero.at(t).Fill(TruthTauMu.E() - TauMu_EF.at(0).E(),1);

	    TauA1PhiZero.at(t).Fill(TauA1_EF.at(0).Phi(),1);
	    TauA1EtaZero.at(t).Fill(TauA1_EF.at(0).Eta(),1);
	    TauA1PtZero.at(t).Fill(TauA1_EF.at(0).Pt(),1);
	    
	    TauMuPhiZero.at(t).Fill(TauMu_EF.at(0).Phi(),1);
	    TauMuEtaZero.at(t).Fill(TauMu_EF.at(0).Eta(),1);
	    TauMuPtZero.at(t).Fill(TauMu_EF.at(0).Pt(),1);
	    

	    csumZero.at(t).Fill(csum_GF.at(0),1);
	  }
	  
	}


	AmbiguitySolverH.at(t).Fill(isInPhysicalRegion,1);
	AmbiguitySolverHvsSign.at(t).Fill(isInPhysicalRegion,significan,1);

    }
  }
}





void  FitTest::Finish(){
  unsigned int t=1;
  //    char name[5];
  for(int iBin = 1; iBin < TauMuEResolVsSignificanceSlices.at(t).GetNbinsX() + 1; iBin++){
    //   sprintf (name, "_%d", iBin);
    TauMuEResolVsSignificanceSlices.at(t).SetBinContent(iBin,TauMuEResolVsSignificance.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    TauMuEResolVsSignificanceSlices.at(t).SetBinError(iBin,TauMuEResolVsSignificance.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
    //   std::cout<<" iBin  RMS  "<< iBin<<" "<<iBin+1<<"   "<<TauMuEResolVsSignificance.at(t).ProjectionX(" ",iBin, iBin)->GetRMS()<<std::endl;
  }


  for(int iBin = 1; iBin < ResTauMuPtEtaTauMuProfile.at(t).GetNbinsX() + 1; iBin++){
    ResTauMuPtEtaTauMuProfile.at(t).SetBinContent(iBin,ResTauMuPtEtaTauMu.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    ResTauMuPtEtaTauMuProfile.at(t).SetBinError(iBin,ResTauMuPtEtaTauMu.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
  }

  for(int iBin = 1; iBin < ResTauMuEEtaTauMuProfile.at(t).GetNbinsX() + 1; iBin++){
    ResTauMuEEtaTauMuProfile.at(t).SetBinContent(iBin,ResTauMuEEtaTauMu.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    ResTauMuEEtaTauMuProfile.at(t).SetBinError(iBin,ResTauMuEEtaTauMu.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
  }

  for(int iBin = 1; iBin < ResTauMuPtPtTauMuProfile.at(t).GetNbinsX() + 1; iBin++){
    ResTauMuPtPtTauMuProfile.at(t).SetBinContent(iBin,ResTauMuPtPtTauMu.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    ResTauMuPtPtTauMuProfile.at(t).SetBinError(iBin,ResTauMuPtPtTauMu.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
  }


  for(int iBin = 1; iBin < ResTauA1PtEtaTauA1Profile.at(t).GetNbinsX() + 1; iBin++){
    ResTauA1PtEtaTauA1Profile.at(t).SetBinContent(iBin,ResTauA1PtEtaTauA1.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    ResTauA1PtEtaTauA1Profile.at(t).SetBinError(iBin,ResTauA1PtEtaTauA1.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
  }


  for(int iBin = 1; iBin < ResTauA1EEtaTauA1Profile.at(t).GetNbinsX() + 1; iBin++){
    ResTauA1EEtaTauA1Profile.at(t).SetBinContent(iBin,ResTauA1EEtaTauA1.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    ResTauA1EEtaTauA1Profile.at(t).SetBinError(iBin,ResTauA1EEtaTauA1.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
  }

  for(int iBin = 1; iBin < ResTauA1PtPtTauA1Profile.at(t).GetNbinsX() + 1; iBin++){
    ResTauA1PtPtTauA1Profile.at(t).SetBinContent(iBin,ResTauA1PtPtTauA1.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    ResTauA1PtPtTauA1Profile.at(t).SetBinError(iBin,ResTauA1PtPtTauA1.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
  }


  for(int iBin = 1; iBin < ResTauMuPtProbProfile.at(t).GetNbinsX() + 1; iBin++){
    ResTauMuPtProbProfile.at(t).SetBinContent(iBin,ResTauMuPtProb.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    ResTauMuPtProbProfile.at(t).SetBinError(iBin,ResTauMuPtProb.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
  }


  for(int iBin = 1; iBin < ResTauMuEProbProfile.at(t).GetNbinsX() + 1; iBin++){
    ResTauMuEProbProfile.at(t).SetBinContent(iBin,ResTauMuEProb.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    ResTauMuEProbProfile.at(t).SetBinError(iBin,ResTauMuEProb.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
  }


  for(int iBin = 1; iBin < ResTauA1PtProbProfile.at(t).GetNbinsX() + 1; iBin++){
    ResTauA1PtProbProfile.at(t).SetBinContent(iBin,ResTauA1PtProb.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    ResTauA1PtProbProfile.at(t).SetBinError(iBin,ResTauA1PtProb.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
  }

  for(int iBin = 1; iBin < ResTauA1EProbProfile.at(t).GetNbinsX() + 1; iBin++){
    ResTauA1EProbProfile.at(t).SetBinContent(iBin,ResTauA1EProb.at(t).ProjectionX(" ",iBin, iBin)->GetRMS());
    ResTauA1EProbProfile.at(t).SetBinError(iBin,ResTauA1EProb.at(t).ProjectionX(" ",iBin, iBin)->GetRMSError());
  }



  if(Nminus0.at(0).at(1).Integral()!=0){ if(HConfig.GetHisto(false,998,t)){AmbiguityEfficiency.at(t).Scale(1/Nminus0.at(0).at(t).Integral());}}
  // if(Nminus0.at(0).at(1).Integral()!=0){ if(HConfig.GetHisto(false,998,t)){CorrectAmbiguityTruth.at(t).Scale(1/Nminus0.at(0).at(t).Integral());}}


  Selection::Finish();
}



