
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
#include "TRandom.h"

static int NUMBER_OF_BINS=20;
//static int NUMBER_OF_BINS=10;

static int Min =-1.;
static int Max = 1;

static float Pol = -0.25;
static float Lum = 50;


TFile *templateSignalPlus = TFile::Open("LOCAL_COMBINED_templatesmc_default.root");
TFile *data = TFile::Open("LOCAL_COMBINED_templatesdata_default.root");



//        TH1F *template_minus15 = (TH1F *)templateSignalPlus->Get("polarizationtemplates_default_x2minusSignal");
//        TH1F *template_plus15= (TH1F *)templateSignalPlus->Get("polarizationtemplates_default_x2plusSignal");

TH1F *template_minus15 = (TH1F *)templateSignalPlus->Get("templatesmc_default_RecOmegaA1mSignal");
TH1F *template_plus15= (TH1F *)templateSignalPlus->Get("templatesmc_default_RecOmegaA1pSignal");
TH1F *h_data15  = (TH1F *)data->Get("templatesdata_default_RecOmegaA1pData");




TH1F *plus = new TH1F("plus","plus",20,-1.0,1.0);
TH1F *mins = new TH1F("mins","mins",20,-1.0,1.0);


static float Scale = Lum*h_data15->Integral();

void PerformHisto(){

 
  for(int ib = 1; ib <NUMBER_OF_BINS + 1; ib++ ){
    plus->SetBinContent(ib,0);
    mins->SetBinContent(ib,0);
    plus->SetBinContent(ib, template_plus15->GetBinContent(ib));
    mins->SetBinContent(ib, template_minus15->GetBinContent(ib));

  }

  plus ->Scale(1/plus->Integral()); 
  mins ->Scale(1/mins->Integral());


}

 




void ToyMC(){


  PerformHisto();
  TH1F *FakeData = new TH1F("FakeData","FakeData",20,-1,1.);
  TH1F *FakeData2 = new TH1F("FakeData2","FakeData2",20,-1.0,1.);

  FakeData->Add(plus,0.5*(1+Pol));
  FakeData->Add(mins,0.5*(1-Pol));
  FakeData->Scale(Scale);
  //  FakeDataClone =   FakeData->Clone();
  //  FakeDataClone->SetName("FakeDataClone");
  TRandom *rd = new TRandom();
  for(int ibin=1; ibin < FakeData2->GetNbinsX(); ibin++ ){
    FakeData2->SetBinContent(ibin, rd->Poisson(FakeData->GetBinContent(ibin)));

  }
  TFile *out= new TFile("output.root","RECREATE");


//cout<<rd->Poisson(5)<<endl;


  mins->Scale(Scale);
  plus->Scale(Scale);
  plus->Write();
  mins->Write();
  FakeData2->Write();
  out->Write();
  out->Close();
//     FakeData->Draw();



}
