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
void justQsq(){ 
  gROOT->Reset();
  gStyle->SetOptStat(2211);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
 
      TChain *TPb208_L6 = new TChain("T");
      TChain *TPb208_L7 = new TChain("T");
      TChain *TPb208_L8 = new TChain("T");
      TChain *TPb208_L9 = new TChain("T");
      TChain *TPb208_L10 = new TChain("T");
      TPb208_L6->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",2317));
      TPb208_L7->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",2291));
      TPb208_L8->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",2292));
      TPb208_L9->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",2293));
      TPb208_L10->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",2294));
      double upADCcutL = 480;
      TCut cut_wadcL = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&P.upQadcL>%f)",upADCcutL);
//      TCut cut_wadcL = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.006&&L.gold.dp<-0.002&&P.upQadcL>%f)",upADCcutL);
      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.30);
      gStyle->SetStatH(0.08);
      TCanvas *c1 = new TCanvas("c1","c1",1000,600);
      c1->Divide(2,1);
      c1->cd(1);
      TH1F* LQsq_Pb6 = new TH1F("LQsq_Pb6",Form("LHRS Qsq using T1 trigger;Qsq (GeV/c)^{2}"), 100, 0.002, 0.013);
      TH1F* LQsq_Pb7 = new TH1F("LQsq_Pb7",Form("LHRS Qsq using T1 trigger;Qsq (GeV/c)^{2}"), 100, 0.002, 0.013);
      TH1F* LQsq_Pb8 = new TH1F("LQsq_Pb8",Form("LHRS Qsq using T1 trigger;Qsq (GeV/c)^{2}"), 100, 0.002, 0.013);
      TH1F* LQsq_Pb9 = new TH1F("LQsq_Pb9",Form("LHRS Qsq using T1 trigger;Qsq (GeV/c)^{2}"), 100, 0.002, 0.013);
      TH1F* LQsq_Pb10 = new TH1F("LQsq_Pb10",Form("LHRS Qsq using T1 trigger;Qsq (GeV/c)^{2}"), 100, 0.002, 0.013);
      LQsq_Pb6->SetLineColor(1);
      LQsq_Pb7->SetLineColor(2);
      LQsq_Pb8->SetLineColor(3);
      LQsq_Pb9->SetLineColor(4);
      LQsq_Pb10->SetLineColor(6);

      TPb208_L6->Draw("EK_L.Q2>>LQsq_Pb6", cut_wadcL,"hist");
      TPb208_L7->Draw("EK_L.Q2>>LQsq_Pb7", cut_wadcL,"hist sames");
      TPb208_L8->Draw("EK_L.Q2>>LQsq_Pb8", cut_wadcL,"hist sames");
      TPb208_L9->Draw("EK_L.Q2>>LQsq_Pb9", cut_wadcL,"hist sames");
      TPb208_L10->Draw("EK_L.Q2>>LQsq_Pb10", cut_wadcL,"hist sames");
      Double_t scalePb208_6 = (1./(LQsq_Pb6->GetMaximum()));
      Double_t scalePb208_7 = (1./(LQsq_Pb7->GetMaximum()));
      Double_t scalePb208_8 = (1./(LQsq_Pb8->GetMaximum()));
      Double_t scalePb208_9 = (1./(LQsq_Pb9->GetMaximum()));
      Double_t scalePb208_10 = (1./(LQsq_Pb10->GetMaximum()));

      LQsq_Pb6->Scale(scalePb208_6);
      LQsq_Pb7->Scale(scalePb208_7);
      LQsq_Pb8->Scale(scalePb208_8);
      LQsq_Pb9->Scale(scalePb208_9);
      LQsq_Pb10->Scale(scalePb208_10);
      gPad->Update();

//      LQsq_Pb6->GetYaxis()->SetRangeUser(0,85e-6);
      TPaveStats *Pb208_6 = (TPaveStats*)LQsq_Pb6->FindObject("stats");
      TPaveStats *Pb208_7 = (TPaveStats*)LQsq_Pb7->FindObject("stats");
      TPaveStats *Pb208_8 = (TPaveStats*)LQsq_Pb8->FindObject("stats");
      TPaveStats *Pb208_9 = (TPaveStats*)LQsq_Pb9->FindObject("stats");
      TPaveStats *Pb208_10 = (TPaveStats*)LQsq_Pb10->FindObject("stats");
      Pb208_6->SetY2NDC(0.90);
      Pb208_6->SetY1NDC(0.75);
      Pb208_7->SetY2NDC(0.75);
      Pb208_7->SetY1NDC(0.60);
      Pb208_8->SetY2NDC(0.60);
      Pb208_8->SetY1NDC(0.45);
      Pb208_9->SetY2NDC(0.45);
      Pb208_9->SetY1NDC(0.30);
      Pb208_10->SetY2NDC(0.30);
      Pb208_10->SetY1NDC(0.15);
      Pb208_6->SetTextColor(1);
      Pb208_7->SetTextColor(2);
      Pb208_8->SetTextColor(3);
      Pb208_9->SetTextColor(4);
      Pb208_10->SetTextColor(6);
      gPad->Modified();

      TChain *TPb208_R6 = new TChain("T");
      TChain *TPb208_R7 = new TChain("T");
      TChain *TPb208_R8 = new TChain("T");
      TChain *TPb208_R9 = new TChain("T");
      TChain *TPb208_R10 = new TChain("T");
      TPb208_R6->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",21436));
      TPb208_R7->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",21412));
      TPb208_R8->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",21413));
      TPb208_R9->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",21414));
      TPb208_R10->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",21415));
      double upADCcutR = 500;
      TCut cut_wadcR = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&P.upQadcR>%f)",upADCcutR);
//      TCut cut_wadcR = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.02&&R.gold.dp<0&&P.upQadcR>%f)",upADCcutR);
      c1->cd(2);
      TH1F* RQsq_Pb6 = new TH1F("RQsq_Pb6",Form("RHRS Qsq using T1 trigger;Qsq (GeV/c)^{2}"), 100, 0.002, 0.013);
      TH1F* RQsq_Pb7 = new TH1F("RQsq_Pb7",Form("RHRS Qsq using T1 trigger;Qsq (GeV/c)^{2}"), 100, 0.002, 0.013);
      TH1F* RQsq_Pb8 = new TH1F("RQsq_Pb8",Form("RHRS Qsq using T1 trigger;Qsq (GeV/c)^{2}"), 100, 0.002, 0.013);
      TH1F* RQsq_Pb9 = new TH1F("RQsq_Pb9",Form("RHRS Qsq using T1 trigger;Qsq (GeV/c)^{2}"), 100, 0.002, 0.013);
      TH1F* RQsq_Pb10 = new TH1F("RQsq_Pb10",Form("RHRS Qsq using T1 trigger;Qsq (GeV/c)^{2}"), 100, 0.002, 0.013);
      RQsq_Pb6->SetLineColor(1);
      RQsq_Pb7->SetLineColor(2);
      RQsq_Pb8->SetLineColor(3);
      RQsq_Pb9->SetLineColor(4);
      RQsq_Pb10->SetLineColor(6);
      TPb208_R6->Draw("EK_R.Q2>>RQsq_Pb6", cut_wadcR,"hist");
      TPb208_R7->Draw("EK_R.Q2>>RQsq_Pb7", cut_wadcR,"hist sames");
      TPb208_R8->Draw("EK_R.Q2>>RQsq_Pb8", cut_wadcR,"hist sames");
      TPb208_R9->Draw("EK_R.Q2>>RQsq_Pb9", cut_wadcR,"hist sames");
      TPb208_R10->Draw("EK_R.Q2>>RQsq_Pb10", cut_wadcR,"hist sames");

      scalePb208_6 = (1./(RQsq_Pb6->GetMaximum()));
      scalePb208_7 = (1./(RQsq_Pb7->GetMaximum()));
      scalePb208_8 = (1./(RQsq_Pb8->GetMaximum()));
      scalePb208_9 = (1./(RQsq_Pb9->GetMaximum()));
      scalePb208_10 = (1./(RQsq_Pb10->GetMaximum()));

      RQsq_Pb6->Scale(scalePb208_6);
      RQsq_Pb7->Scale(scalePb208_7);
      RQsq_Pb8->Scale(scalePb208_8);
      RQsq_Pb9->Scale(scalePb208_9);
      RQsq_Pb10->Scale(scalePb208_10);
      gPad->Update();
      
//      RQsq_Pb6->GetYaxis()->SetRangeUser(0,85e-6);
      Pb208_6 = (TPaveStats*)RQsq_Pb6->FindObject("stats");
      Pb208_7 = (TPaveStats*)RQsq_Pb7->FindObject("stats");
      Pb208_8 = (TPaveStats*)RQsq_Pb8->FindObject("stats");
      Pb208_9 = (TPaveStats*)RQsq_Pb9->FindObject("stats");
      Pb208_10 = (TPaveStats*)RQsq_Pb10->FindObject("stats");
      Pb208_6->SetY2NDC(0.90);
      Pb208_6->SetY1NDC(0.75);
      Pb208_7->SetY2NDC(0.75);
      Pb208_7->SetY1NDC(0.60);
      Pb208_8->SetY2NDC(0.60);
      Pb208_8->SetY1NDC(0.45);
      Pb208_9->SetY2NDC(0.45);
      Pb208_9->SetY1NDC(0.30);
      Pb208_10->SetY2NDC(0.30);
      Pb208_10->SetY1NDC(0.15);
      Pb208_6->SetTextColor(1);
      Pb208_7->SetTextColor(2);
      Pb208_8->SetTextColor(3);
      Pb208_9->SetTextColor(4);
      Pb208_10->SetTextColor(6);
      gPad->Modified();

      c1->SaveAs(Form("./plots/Qsq_combined.png"));

}
