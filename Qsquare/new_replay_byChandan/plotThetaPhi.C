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
 
      TChain *T1 = new TChain("T");
      TChain *T2 = new TChain("T");
      if(run_num<10000){
	//LHRS
      T1->Add(Form("/chafs1/work1/prex_counting/QsqRootFiles/prexLHRS_%d_-1*.root",run_num));
      T2->Add(Form("/chafs1/work1/prex_counting/RootFilesQsq/prexLHRS_%d_-1*.root",run_num));
      TCut cut_wadcL = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05)");
//      TCut cut_wadcL = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)");

      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.35);
      gStyle->SetStatH(0.08);
      TCanvas *c1 = new TCanvas("c1","c1",800,600);
      c1->Divide(2,1);
      c1->cd(1);
      TH1F* thL1 = new TH1F("thL1",Form("L Vertical Transport Angle (run%d);#theta_{tg} (rad);",run_num),200,-0.08,0.08);
      TH1F* thL2 = new TH1F("thL2",Form("L Vertical Transport Angle (run%d);#theta_{tg} (rad);",run_num),200,-0.08,0.08);
      thL1->SetLineColor(kBlack);
      thL2->SetLineColor(kRed);
      T1->Draw("L.gold.th>>thL1",cut_wadcL);
      T2->Draw("L.gold.th>>thL2",cut_wadcL,"same");
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(kBlack);
      latex.DrawLatex(0.15,0.85,"L.gold.th (old DB)");
      latex.SetTextColor(kRed);
      latex.DrawLatex(0.15,0.80,"L.gold.th (new DB)");
      c1->cd(2);
      TH1F* phL1 = new TH1F("phL1",Form("L Horizontal Transport Angle (run%d);#phi_{tg} (rad);",run_num),200,-0.04,0.04);
      TH1F* phL2 = new TH1F("phL2",Form("L Horizontal Transport Angle (run%d);#phi_{tg} (rad);",run_num),200,-0.04,0.04);
      phL1->SetLineColor(kBlack);
      phL2->SetLineColor(kRed);
      T1->Draw("L.gold.ph>>phL1",cut_wadcL);
      T2->Draw("L.gold.ph>>phL2",cut_wadcL,"same");
      latex.SetTextColor(kBlack);
      latex.DrawLatex(0.15,0.85,"L.gold.ph (old DB)");
      latex.SetTextColor(kRed);
      latex.DrawLatex(0.15,0.80,"L.gold.ph (new DB)");
      ofstream outfile("./TextFiles/tg_thL.csv",ios_base::app);
      outfile<<thL1->GetMean()<<"\t"<<thL1->GetRMS()<<"\t"<<phL1->GetMean()<<"\t"<<phL2->GetRMS()<<endl;
      outfile.close();
      c1->SaveAs(Form("./temp/transport_angle_angle_%s_run%d.pdf",target.Data(),run_num));
      }else{
	//LHRS
      T1->Add(Form("/chafs1/work1/prex_counting/QsqRootFiles/prexRHRS_%d_-1*.root",run_num));
      T2->Add(Form("/chafs1/work1/prex_counting/RootFilesQsq/prexRHRS_%d_-1*.root",run_num));
      TCut cut_wadcR = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05)");
//      TCut cut_wadcR = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)");

      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.35);
      gStyle->SetStatH(0.08);
      TCanvas *c1 = new TCanvas("c1","c1",800,600);
      c1->Divide(2,1);
      c1->cd(1);
      TH1F* thR1 = new TH1F("thR1",Form("R Vertical Transport Angle (run%d);#theta_{tg} (rad);",run_num),200,-0.08,0.08);
      TH1F* thR2 = new TH1F("thR2",Form("R Vertical Transport Angle (run%d);#theta_{tg} (rad);",run_num),200,-0.08,0.08);
      thR1->SetLineColor(kBlack);
      thR2->SetLineColor(kRed);
      T1->Draw("R.gold.th>>thR1",cut_wadcR);
      T2->Draw("R.gold.th>>thR2",cut_wadcR,"same");
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(kBlack);
      latex.DrawLatex(0.15,0.85,"R.gold.th (old DB)");
      latex.SetTextColor(kRed);
      latex.DrawLatex(0.15,0.80,"R.gold.th (new DB)");
      c1->cd(2);
      TH1F* phR1 = new TH1F("phR1",Form("R Horizontal Transport Angle (run%d);#phi_{tg} (rad);",run_num),200,-0.04,0.04);
      TH1F* phR2 = new TH1F("phR2",Form("R Horizontal Transport Angle (run%d);#phi_{tg} (rad);",run_num),200,-0.04,0.04);
      phR1->SetLineColor(kBlack);
      phR2->SetLineColor(kRed);
      T1->Draw("R.gold.ph>>phR1",cut_wadcR);
      T2->Draw("R.gold.ph>>phR2",cut_wadcR,"same");
      latex.SetTextColor(kBlack);
      latex.DrawLatex(0.15,0.85,"R.gold.ph (old DB)");
      latex.SetTextColor(kRed);
      latex.DrawLatex(0.15,0.80,"R.gold.ph (new DB)");
      ofstream outfile("./TextFiles/tg_thR.csv",ios_base::app);
      outfile<<thR1->GetMean()<<"\t"<<thR1->GetRMS()<<"\t"<<phR1->GetMean()<<"\t"<<phR2->GetRMS()<<endl;
      outfile.close();
      c1->SaveAs(Form("./temp/transport_angle_angle_%s_run%d.pdf",target.Data(),run_num));

      }
}
