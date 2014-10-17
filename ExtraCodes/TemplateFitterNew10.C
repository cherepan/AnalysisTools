#include "stdlib.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TTree.h"
#include "TComplex.h"
#include "TMath.h"
#include "TMinuit.h"

static int NUMBER_OF_BINS=20;
//static int NUMBER_OF_BINS=10;

static int Min =-1.;
static int Max = 1;




//TFile *templateSignalMinus = TFile::Open("LOCAL_COMBINED_signaltemplatePlus_default.root");
TFile *templateSignalPlus = TFile::Open("LOCAL_COMBINED_templatesmc_default.root");
TFile *data = TFile::Open("LOCAL_COMBINED_templatesdata_default.root");



//       TH1F *template_minus15 = (TH1F *)templateSignalPlus->Get("polarizationtemplates_default_x2minusSignal");
//       TH1F *template_plus15= (TH1F *)templateSignalPlus->Get("polarizationtemplates_default_x2plusSignal");

//       TH1F *h_data15  = (TH1F *)data->Get("basicselection_default_xambData");
//       TH1F *h_bkg2  = (TH1F *)data->Get("basicselection_default_xambMC_WJ1");
//       TH1F *h_bkg1  = (TH1F *)data->Get("basicselection_default_xambMC_QCD");


//   KEY: TH1D     templatesdata_default_RecOmegaA1pData;1 Data
//   KEY: TH1D     templatesdata_default_RecOmegaA1mData;1 Data


//   KEY: TH1D     templatesdata_default_RecOmegaA1zData;1 Data
//   KEY: TH1D     templatesdata_default_RecOmegaMupData;1 Data
//   KEY: TH1D     templatesdata_default_RecOmegaMumData;1 Data


//   KEY: TH1D     templatesdata_default_RecOmegaMuzData;1 Data
//   KEY: TH1D     templatesdata_default_RecOmegaCombinedpData;1   Data
//   KEY: TH1D     templatesdata_default_RecOmegaCombinedmData;1   Data
//   KEY: TH1D     templatesdata_default_RecOmegaCombinedzData;1   Data
 



TH1F *template_minus15 = (TH1F *)templateSignalPlus->Get("templatesmc_default_RecOmegaA1mSignal");
TH1F *template_plus15= (TH1F *)templateSignalPlus->Get("templatesmc_default_RecOmegaA1pSignal");

TH1F *h_data15  = (TH1F *)data->Get("templatesdata_default_RecOmegaA1pData");

TH1F *h_bkg2  = (TH1F *)data->Get("templatesdata_default_RecOmegaA1pMC_WJ1");
TH1F *h_bkg1  = (TH1F *)data->Get("templatesdata_default_RecOmegaA1pQCD");

TH1F *forEff_bkg2  = (TH1F *)data->Get("templatesdata_default_TauA1VisiblePtMC_WJ1");
TH1F *forEff_bkg1  = (TH1F *)data->Get("templatesdata_default_TauA1VisiblePtQCD");


TH1F *h_bkgUnSc2  = (TH1F *)data->Get("templatesdata_default_RecOmegaA1pMC_WJ1");
TH1F *h_bkgUnSc1  = (TH1F *)data->Get("templatesdata_default_RecOmegaA1pQCD");



// template_plus15->Rebin(2);
// template_minus15->Rebin(2);
// h_data15->Rebin(2);
// h_bkg2->Rebin(2);
// h_bkg1->Rebin(2);



//h_data15->Draw();

 TH1F *plus = new TH1F("plus","plus",20,-1.0,1.0);
 TH1F *mins = new TH1F("mins","mins",20,-1.0,1.0);

 TH1F *bkg2 = new TH1F("bkg2","bkg2",20,-1.0,1.0);
 TH1F *bkg1 = new TH1F("bkg1","bkg1",20,-1.0,1.0);


// TH1F *plus = new TH1F("plus","plus",10,0.0,1.5);
// TH1F *mins = new TH1F("mins","mins",10,0.0,1.5);


TH1F *plusUnSc = new TH1F("plusUnSc","plusUnSc",20,-1.0,1.0);
TH1F *minsUnSc = new TH1F("minsUnSc","minsUnSc",20,1.0,1.0);

double NDN;


void PerformHisto(){

  template_plus15->Sumw2();
  template_minus15->Sumw2();

  h_data15->Sumw2();
  h_bkg2->Sumw2();
  h_bkg1->Sumw2();



  double bkg1Int = h_bkg1->Integral();
  double bkg2Int = h_bkg2->Integral();

  double effqcd =h_bkgUnSc1->Integral()/ forEff_bkg1->Integral();
  double effwj  =h_bkgUnSc2->Integral()/ forEff_bkg2->Integral();


//   h_bkg1->Scale(760.1*effqcd/bkg1Int);
//   h_bkg2->Scale(34.397*effwj/bkg2Int);


 
  for(int ib = 1; ib <NUMBER_OF_BINS + 1; ib++ ){
    plus->SetBinContent(ib,0);
    mins->SetBinContent(ib,0);

    
    plus->SetBinContent(ib, template_plus15->GetBinContent(ib));
    mins->SetBinContent(ib, template_minus15->GetBinContent(ib));
    //     plus->SetBinError(ib, sqrt(template_plus->GetBinContent(ib)));
    //     mins->SetBinError(ib, sqrt(template_minus->GetBinContent(ib)));
    plusUnSc->SetBinContent(ib, template_plus15->GetBinContent(ib));
    minsUnSc->SetBinContent(ib, template_minus15->GetBinContent(ib));
    
    bkg2 ->SetBinContent(ib, h_bkg2->GetBinContent(ib));
    bkg1 ->SetBinContent(ib, h_bkg1->GetBinContent(ib));



    //    h_data ->SetBinContent(ib,h_data->GetBinContent(ib) );
          printf("template_qcd->GetBinContent(ib) %d     template_ewk->GetBinContent(ib) %d \n",plus->GetBinContent(ib),mins->GetBinContent(ib));

  }
  double bkg1Int = h_bkg1->Integral();
  double bkg2Int = h_bkg2->Integral();

  double effqcd =h_bkgUnSc1->Integral()/ forEff_bkg1->Integral();
  double effwj  =h_bkgUnSc2->Integral()/ forEff_bkg2->Integral();


//   h_bkg1->Scale(760.1*effqcd/bkg1Int);
//   h_bkg2->Scale(343.97*effwj/bkg2Int);
  
//   bkg1 ->Scale(760.1*effqcd/bkg1Int);
//   bkg2 ->Scale(343.97*effwj/bkg2Int);


   printf("------------------------mins %f   plus %f  \n",plus->Integral(),mins->Integral() );
   plus ->Scale(1/(plus->Integral() ) ); 
   mins ->Scale(1/(mins->Integral() ));
   printf("------------------------mins %f   plus %f  data %f   %f  %f  \n",plus->Integral(),mins->Integral(), h_data15->Integral(),h_bkg1->Integral(), h_bkg2->Integral());
   plus ->Scale(h_data15->Integral() - h_bkg1->Integral() - h_bkg2->Integral()); 
   mins ->Scale(h_data15->Integral() - h_bkg1->Integral() - h_bkg2->Integral());
   NDN = (h_data15->Integral()  - h_bkg1->Integral() - h_bkg2->Integral())/(plusUnSc->Integral() + minsUnSc->Integral());
   //  printf("------------------------mins %f   plus %f  data %f   %f  %f\n",mins->Integral(),plus->Integral(),h_data15->Integral(), minsUnSc->Integral()+plusUnSc->Integral(),NDN);
   //  printf("------------------------mins %f   plus %f  data %f   %f  %f  %f  %f\n",effqcd,effwj,h_data15->Integral(),h_bkg2->Integral(),h_bkgUnSc2->Integral(),bkg2Int);
  
 
}


//  



double N_mc(int iBin, double *par){
  return plus->GetBinContent(iBin[0])*(1+par[0])/2 + mins->GetBinContent(iBin[0])*(1-par[0])/2;
}

 double N1(int iBin, double *par){
   //  printf(" >>>>>>>>>>>>>>>>>>>>>>>>>>>>  %f   %f    %f   %f   %f   NDN %f \n", plus->GetBinContent(iBin[0]) ,mins->GetBinContent(iBin[0]), h_bkg1->GetBinContent(iBin[0]) , h_bkg2->GetBinContent(iBin[0]),(plus->GetBinContent(iBin[0])*(1+par[0])/2 + mins->GetBinContent(iBin[0])*(1-par[0])/2)*NDN + h_bkg1->GetBinContent(iBin[0]) + h_bkg2->GetBinContent(iBin[0]),NDN );

   return (plus->GetBinContent(iBin[0])*(1+par[0])/2 + mins->GetBinContent(iBin[0])*(1-par[0])/2) + h_bkg1->GetBinContent(iBin[0]) + h_bkg2->GetBinContent(iBin[0]) ;

 }

 double NR(int iBin, double *par){
   return plusUnSc->GetBinContent(iBin[0])*log(plus->GetBinContent(iBin[0])) - plus->GetBinContent(iBin[0]);
 }

 double NL(int iBin, double *par){
   return minsUnSc->GetBinContent(iBin[0])*log(mins->GetBinContent(iBin[0])) - mins->GetBinContent(iBin[0]);
 }
 double NB1(int iBin, double *par){
   if( h_bkg1->GetBinContent(iBin[0]) ==0) return 0;
   return h_bkgUnSc1->GetBinContent(iBin[0])*log(h_bkg1->GetBinContent(iBin[0])) - h_bkg1->GetBinContent(iBin[0]);
 }
 double NB2(int iBin, double *par){
   if( h_bkg2->GetBinContent(iBin[0]) ==0) return 0;

   return h_bkgUnSc2->GetBinContent(iBin[0])*log(h_bkg2->GetBinContent(iBin[0])) - h_bkg2->GetBinContent(iBin[0]);

 

 }




double N_mc_check(int iBin, double p){
 
  
  double Nmc_ev = plus->GetBinContent(iBin[0])*(1+p)/2 + mins->GetBinContent(iBin[0])*(1-p)/2;
       return Nmc_ev;
}


int N_data(int iBin){
  return h_data15->GetBinContent(iBin[0]);
  //  return h_data15->GetBinContent(iBin[0]) - h_bkg1->GetBinContent(iBin[0]));
  //return h_data15->GetBinContent(iBin[0]) - h_bkg1->GetBinContent(iBin[0]);

}


double PolarUncer(double p,double m, double dp,double dm){


  double deltaPol;
  deltaPol = (fabs(p - m)/(p+m)/(p+m)/(p+m))*sqrt(dp*dp + dm*dm);
  return deltaPol;
}



// -----------------------  main function ------------------------------------
void TemplateFitterNew10(string type, string ScanOrFit = "f"){


  PerformHisto();
//   mins->Draw();
//   plus->Draw("same");
  gStyle->SetOptStat(0000);

   // Double LLScan[95],ParScan[95];

  Double_t par[1],dpar[1];
  TMinuit *myMin = new TMinuit(3);

  if(type == 'c')   myMin->SetFCN(Chi2); myMin->SetErrorDef(1);
  if(type == 'l')   myMin->SetFCN(Likelihood); myMin->SetErrorDef(0.5);
 
  Double_t arglist[2];
  Int_t ierflg = 1;


 
  
   myMin->mnparm( 0, "polarization",  -0.05,   0.005,   -1, 1., ierflg);
  // myMin->mnparm( 0, "polarization",  -1.,   0.05,   -1.00., 0., ierflg);
 //myMin->FixParameter(0);
   myMin->SetPrintLevel(1);

  arglist[0]=2;
  arglist[1]=2;
  myMin->mnexcm("SET STR", arglist ,1,ierflg);


  myMin->Migrad();

  myMin->mnhess();
  myMin->mnhess();
  myMin->mnhess();
  myMin->mnmnos();
  
  get_parameters(myMin,par,dpar);

//    for(int ib = 1; ib <15 + 1; ib++ ){
        
//      printf("   MCCHECK %f \n",N_mc_check(ib,-0.525000));

//    }



  TH1F *Sum = new TH1F("Sum","Sum",20,-1.0,1.0);
  
  plus->SetLineColor(3);
  mins->SetLineColor(2);
  h_data15->SetLineWidth(2);
  
  plus->Clone("plusClone");
  mins->Clone("minsClone");

  bkg1->Clone("bkg1Clone");
  bkg2->Clone("bkg2Clone");

  
  
  plusClone->SetName("plusClone");
  minsClone->SetName("minsClone");
  
  plusClone->Scale((1+par[0])/2);
  minsClone->Scale((1-par[0])/2);
//   plusClone->Scale(0.4);
//   minsClone->Scale(0.6);

  template_minus15->Scale((1-par[0])/2);
  template_plus15->Scale((1+par[0])/2);
  double  llcheck =0;
  for(int ib = 1; ib <20 + 1; ib++ ){
    
    Sum->SetBinContent(ib, plusClone->GetBinContent(ib) + minsClone->GetBinContent(ib)  +bkg1Clone->GetBinContent(ib) + bkg2Clone->GetBinContent(ib)  );
    //    printf("Sum  %f    MCCHECK %f \n",Sum->GetBinContent(ib),N_mc_check(ib,-0.525000));
    if(Sum->GetBinContent(ib)!=0)  llcheck+= (h_data15->GetBinContent(ib) - Sum->GetBinContent(ib))*(h_data15->GetBinContent(ib) - Sum->GetBinContent(ib))/Sum->GetBinContent(ib);
  }
  
  
  //      myMin->Command("SCAn 0");
  //      TGraph *gr = (TGraph*)myMin->GetPlot();
  //      gr->Draw("alp"); 
  
  printf("Nplus  %f  Nminus  %f   polar  %f   unc %f     llcheck %f\n",template_plus15->Integral(),template_minus15->Integral(), par[0], dpar[0],2*llcheck);
   h_data15->SetStats(0);
   h_data15->SetLineWidth(2);
   h_data15->SetMarkerStyle(20);
   h_data15->SetMarkerSize(1.2);
   
   h_data15->SetTitle("");
   h_data15->SetYTitle("");
   h_data15->SetXTitle("E_{#mu}/E_{#tau}");
   
   Sum->SetLineWidth(2);
   Sum->SetLineColor(4);
   
     Sum->Draw();
     h_data15->Draw("SAMESAMESAME");
     plusClone->Draw("HISTSAME");
     minsClone->Draw("HISTSAME");
     bkg1->SetLineWidth(2);
     //    bkg1->SetLineWidth(2);

     bkg1->Draw("HISTSAME");
     bkg2->Draw("HISTSAME");
   
   TLegend *leg = new TLegend(0.5201149,0.4957627,0.8965517,0.8432203,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(82);
   leg->SetLineColor(0);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("NULL","Data","lpf");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Sum","Fit","f");
   entry->SetFillStyle(1001);
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("plusClone","h_{#tau}  = +1","f");
   entry->SetFillStyle(1001);
   entry->SetLineColor(3);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("minsClone","h_{#tau}  = -1","f");
   entry->SetFillStyle(1001);
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);

   
   leg->Draw();
    TPaveText leg2(0.1,0.925,.9,0.975,"NDC");
   leg2.SetBorderSize(0);
   leg2.SetFillStyle(4000);
   leg2.SetFillColor(0);
   leg2.SetTextSize(0.03);
   leg2.SetMargin(0.15);
   leg2.AddText("CMS, work in progress, 19.7 fb^{-1}, #sqrt{s} = 8 TeV");

   leg2.Draw("SAMESAMESAME");





   if(ScanOrFit == "s"){
   myMin->Command("SCAn 0");
   TGraph *gr = (TGraph*)myMin->GetPlot();
   gr->GetXaxis()->SetTitle("P_{#tau}");
   gr->GetYaxis()->SetTitle("#chi^{2}");

   gr->Draw("alp"); 
   }
   
   
   
   
}

void get_parameters(TMinuit * minuit,double *par, double *par_err ){
  for( int i =0 ; i <1; i++){
    minuit->GetParameter(i,par[i], par_err[i]);
  }
}
//  Minimization function definition




void Likelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double LL=0;
  double LL2=0;
  int bin[1]={0};

  //  printf("---------------- \n");
  
  for ( int i = 1; i< 20 + 1; i++){
    bin[0] = i;
    double  mc = N_mc(i, par);
    double  Sig = N1(i, par);

    // printf("---------------- N1(i, par)  %f  N_MC   %f \n", Sig, mc);
    double  SR =  NR(i, par);
    double  SL =  NL(i, par);

    double  B1 =  NB1(i, par);
    double  B2 =  NB2(i, par);

   
    int scale = 1;
    int scalebk1 = 1;
    int scalebk2 = 1;
    if(Sig==0){ scale =0; Sig =1;}
    else scale = 1;

//     if(B1 <0&& B1==0 && B1 > 0){scalebk1 = 0; B1=0}
//     else scalebk1=1;

//     if(B1 <0&& B2==0 && B2 > 0){scalebk2 = 0; B2=0}
//     else scalebk2=1;

 
    // LL+= (N_data(i)*log(mc)  - scale*mc);
    LL+= (N_data(i)*log(Sig)  - scale*Sig)  + SR  + SL + B1+ B2; 
    //  LL+= (N_data(i)*log(Sig)  - scale*Sig);
      printf("----------------First %f Sig  %f SR  %f  SL  %f  B1   %f  B2  %f    and LL   %f  \n", (N_data(i)*log(Sig)  - scale*Sig),Sig, SR, SL,B1,B2, LL);

  }
 
  f= -LL;
}





void Chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double LL=0;
  int bin[1]={0};
  
  //  printf("---------------- \n");
  int count =0;
  for ( int i = 1; i< 20 + 1; i++){
    //if(i==6) continue;
    bin[0] = i;
    double  mc = N_mc(i, par);
    int scale = 1;
    if(mc==0){ scale =0; mc =1;}
    else scale = 1;
       LL+= scale*(mc -N_data(i))*(mc - N_data(i))/sqrt(mc)/sqrt(mc);
        printf("data %f   mc %f    chi2 %f\n",N_data(i), mc, scale*(mc -N_data(i))*(mc - N_data(i))/sqrt(mc)/sqrt(mc));
       count++;
       
  }
  
     f= LL;    
     //f= LL/2.0;
            printf("Chi2  %f     par %f     \n",f,par[0]);
}


