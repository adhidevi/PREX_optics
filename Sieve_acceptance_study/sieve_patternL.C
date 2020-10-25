#include "CollimatorL.C"
#include "sieveL.C"
#include "TString.h"

void sieve_patternL(){

gStyle->SetOptStat(0);

TChain *T = new TChain("T");
for(int i = 1; i < 6; i++){
T->Add(Form("/work/halla/parity/disk1/ryanrich/carbon/septscaleLHRS_PREX_0.0_%i.root",i));
}

TCut specCut = Form("x_vdc_tr!=-333.");
TCut sieveCut = "sieveL(x_zsieve_tr,y_zsieve_tr)!=1";
TCut colCut = "CollimatorL(x_col_tr,y_col_tr)";
//TCut thphCut = "(th_zsieve_tr>-0.02&&th_zsieve_tr<-0.0009&&ph_zsieve_tr>-0.0125&&ph_zsieve_tr<-0.01)";
TCut thphCut = ""; 

TCut totalCut = specCut&&sieveCut&&colCut&&thphCut;
TH2F *hsieve = new TH2F("hsieve","Sieve Pattern;phi;theta",200,-0.05,0.05,200,-0.06,0.06);

  TCanvas *c = new TCanvas();
  T->Draw("th_zsieve_tr:ph_zsieve_tr >> hsieve",totalCut*"rate");

}
