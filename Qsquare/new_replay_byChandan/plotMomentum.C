#include <TROOT.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <fstream>
#include <iostream>
void plotMomentum(int run,TString target){ 
  gROOT->Reset();
  gStyle->SetOptStat(1002211);
  gStyle->SetOptFit(0);
//  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
 
      TChain *T = new TChain("T");
      if(run<10000){
	//LHRS
      T->Add(Form("/chafs1/work1/prex_counting/QsqRootFiles/prexLHRS_%d_-1*.root",run));
      TCut cutL = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.04&&L.gold.dp<-0.002)");
//      TCut cutL = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)");

      gStyle->SetStatX(0.40);
      gStyle->SetStatY(0.40);
      TH1F* momL = new TH1F("momL",Form("LHRS Momentum Distribution No ADC cut(run%d);L.gold.p (GeV/c)^{2};",run),200,0.91,0.96);
      momL->SetLineColor(1);
      TCanvas *c1 = new TCanvas("c1","c1",1000,700);
      gPad->SetLogy(1);
      T->Draw("L.gold.p>>momL",cutL);
      c1->SaveAs(Form("./temp/momentum_%s_run%d.pdf",target.Data(),run));

      }else{
	//RHRS
      T->Add(Form("/chafs1/work1/prex_counting/QsqRootFiles/prexRHRS_%d_-1*.root",run));
      TCut cutR = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.04&&R.gold.dp<-0.002)");
//      TCut cutR = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)");

      gStyle->SetStatX(0.40);
      gStyle->SetStatY(0.40);
      TH1F* momR = new TH1F("momR",Form("RHRS Momentum Distribution No ADC cut(run%d);R.gold.p (GeV/c)^{2};",run),200,0.91,0.96);
      momR->SetLineColor(1);
      TCanvas *c1 = new TCanvas("c1","c1",1000,700);
      gPad->SetLogy(1);
      T->Draw("R.gold.p>>momR",cutR);
      c1->SaveAs(Form("./temp/momentum_%s_run%d.pdf",target.Data(),run));

      }
}
