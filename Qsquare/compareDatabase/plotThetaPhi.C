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
void plotThetaPhi(int run_num, TString target){ 
  gROOT->Reset();
  gStyle->SetOptStat(1002211);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
 
      TChain *TA = new TChain("T");
      TChain *TS = new TChain("T");
      if(run_num<10000){
	//LHRS
      TA->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run_num));
      TS->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run_num));
      TCut cutL = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05)");

      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.35);
      gStyle->SetStatH(0.08);
      TCanvas *c1 = new TCanvas("c1","c1",800,600);
      c1->Divide(2,1);
      c1->cd(1);
      TH1F* thLA = new TH1F("thLA",Form("L Vertical Transport Angle (run%d);#theta_{tg} (rad);",run_num),200,-0.08,0.08);
      TH1F* thLS = new TH1F("thLS",Form("L Vertical Transport Angle (run%d);#theta_{tg} (rad);",run_num),200,-0.08,0.08);
      thLA->SetLineColor(kBlack);
      thLS->SetLineColor(kRed);
      TA->Draw("L.gold.th>>thLA",cutL);
      TS->Draw("L.gold.th>>thLS",cutL,"same");
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(kBlack);
      latex.DrawLatex(0.55,0.85,"L.gold.th (DB1)");
      latex.SetTextColor(kRed);
      latex.DrawLatex(0.55,0.80,"L.gold.th (DB2)");
      c1->cd(2);
      TH1F* phLA = new TH1F("phLA",Form("L Horizontal Transport Angle (run%d);#phi_{tg} (rad);",run_num),200,-0.04,0.04);
      TH1F* phLS = new TH1F("phLS",Form("L Horizontal Transport Angle (run%d);#phi_{tg} (rad);",run_num),200,-0.04,0.04);
      phLA->SetLineColor(kBlack);
      phLS->SetLineColor(kRed);
      TA->Draw("L.gold.ph>>phLA",cutL);
      TS->Draw("L.gold.ph>>phLS",cutL,"same");
      latex.SetTextColor(kBlack);
      latex.DrawLatex(0.53,0.85,"L.gold.ph (DB1)");
      latex.SetTextColor(kRed);
      latex.DrawLatex(0.53,0.80,"L.gold.ph (DB2)");
      c1->SaveAs(Form("./temp/transport_angle_angle_%s_run%d.pdf",target.Data(),run_num));
      }else{
	//LHRS
      TA->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",run_num));
      TS->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",run_num));
      TCut cutR = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05)");

      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.35);
      gStyle->SetStatH(0.08);
      TCanvas *c1 = new TCanvas("c1","c1",800,600);
      c1->Divide(2,1);
      c1->cd(1);
      TH1F* thRA = new TH1F("thRA",Form("R Vertical Transport Angle (run%d);#theta_{tg} (rad);",run_num),200,-0.08,0.08);
      TH1F* thRS = new TH1F("thRS",Form("R Vertical Transport Angle (run%d);#theta_{tg} (rad);",run_num),200,-0.08,0.08);
      thRA->SetLineColor(kBlack);
      thRS->SetLineColor(kRed);
      TA->Draw("R.gold.th>>thRA",cutR);
      TS->Draw("R.gold.th>>thRS",cutR,"same");
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(kBlack);
      latex.DrawLatex(0.55,0.85,"R.gold.th (DB1)");
      latex.SetTextColor(kRed);
      latex.DrawLatex(0.55,0.80,"R.gold.th (DB2)");
      c1->cd(2);
      TH1F* phRA = new TH1F("phRA",Form("R Horizontal Transport Angle (run%d);#phi_{tg} (rad);",run_num),200,-0.04,0.04);
      TH1F* phRS = new TH1F("phRS",Form("R Horizontal Transport Angle (run%d);#phi_{tg} (rad);",run_num),200,-0.04,0.04);
      phRA->SetLineColor(kBlack);
      phRS->SetLineColor(kRed);
      TA->Draw("R.gold.ph>>phRA",cutR);
      TS->Draw("R.gold.ph>>phRS",cutR,"same");
      latex.SetTextColor(kBlack);
      latex.DrawLatex(0.15,0.85,"R.gold.ph (DB1)");
      latex.SetTextColor(kRed);
      latex.DrawLatex(0.15,0.80,"R.gold.ph (DB2)");
      c1->SaveAs(Form("./temp/transport_angle_angle_%s_run%d.pdf",target.Data(),run_num));

      }
}
