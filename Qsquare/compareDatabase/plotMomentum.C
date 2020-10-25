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
void plotMomentum(int run, TString target){ 
  gROOT->Reset();
  gStyle->SetOptStat(1002211);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
 
      TChain *TA = new TChain("T");//F==all sieve database
      TChain *TS = new TChain("T");//P==selected sieve database
      if(run<10000){
	//LHRS
      TA->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run));
      TS->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run));

      TCut cutL = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05)");

      gStyle->SetStatX(0.40);
      gStyle->SetStatY(0.40);

      TH1F* momLA = new TH1F("momLA",Form("LHRS Momentum Distribution No ADC cut(run%d);L.gold.p (GeV);",run),200,0.92,0.95);
      TH1F* momLS = new TH1F("momLS",Form("LHRS Momentum Distribution No ADC cut(run%d);L.gold.p (GeV);",run),200,0.92,0.95);

      momLA->SetLineColor(1);
      momLS->SetLineColor(2);

      TCanvas *c1 = new TCanvas("c1","c1",700,500);
      gPad->SetLogy(1);
      TA->Draw("L.gold.p>>momLA",cutL);
      TS->Draw("L.gold.p>>momLS",cutL,"same");

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(1);
      latex.DrawLatex(0.20,0.85,"DB1");
      latex.SetTextColor(2);
      latex.DrawLatex(0.20,0.80,"DB2");

      c1->SaveAs(Form("./temp/momentum_DB1DB2_%s_run%d.pdf",target.Data(),run));

      }else{

	//RHRS
      TA->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",run));
      TS->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/fewHoleDataBase/prexRHRS_%d_-1*.root",run));

      TCut cutR = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05)");

      gStyle->SetStatX(0.40);
      gStyle->SetStatY(0.40);

      TH1F* momRA = new TH1F("momRA",Form("RHRS Momentum Distribution No ADC cut(run%d);R.gold.p (GeV};",run),200,0.92,0.95);
      TH1F* momRS = new TH1F("momRS",Form("RHRS Momentum Distribution No ADC cut(run%d);R.gold.p (GeV);",run),200,0.92,0.95);

      momRA->SetLineColor(1);
      momRS->SetLineColor(2);

      TCanvas *c1 = new TCanvas("c1","c1",700,500);
      gPad->SetLogy(1);
      TA->Draw("R.gold.p>>momRA",cutR);
      TS->Draw("R.gold.p>>momRS",cutR,"same");

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(1);
      latex.DrawLatex(0.20,0.85,"DB1");
      latex.SetTextColor(2);
      latex.DrawLatex(0.20,0.80,"DB2");

      c1->SaveAs(Form("./temp/momentum_DB1DB2_%s_run%d.pdf",target.Data(),run));

      }
}
