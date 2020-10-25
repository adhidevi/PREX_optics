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
void plotQsqOnlyUS(int run, TString target, TString date){ 
  gROOT->Reset();
  gStyle->SetOptStat(1002211);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
 
      TChain *TA = new TChain("T");
      TChain *TS = new TChain("T");
  if(run < 10000){
      ofstream outfile("out_fileL.txt",ios_base::app);
      //LHRS
      double upADCcut_approx = 480;
      TA->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run));
      TS->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run));

      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.25);
      gStyle->SetStatH(0.10);

      TCanvas *c1 = new TCanvas("c1","c1",1000,600);
      c1->Divide(2,1);
      c1->cd(1);
      TH1F* huq = new TH1F("huq","P.upQadcL;ADC CH;Event/CH",250,400,650);
      TA->Draw("P.upQadcL>>huq");

      TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
      huqCpy->GetXaxis()->SetRangeUser(upADCcut_approx-10,upADCcut_approx+10);
      cout<<"ADC cut: "<<huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin())<<endl;
      double upADCcut = huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin());

      TLine* line = new TLine(upADCcut,0.0,upADCcut,huq->GetMaximum());
      line->SetLineColor(kMagenta);
      line->Draw(); 
      TF1 *fped=new TF1("fped","gaus",460,470);
      fped->SetLineWidth(2);
      fped->SetLineColor(kRed);
      huq->Fit("fped","R");
      Double_t fpedpar[3];
      fped->GetParameters(&fpedpar[0]);
      char labelped[64];
      sprintf(labelped,"Pedestal = %4.2f",fpedpar[1]);
      TPaveLabel *pt = new TPaveLabel(0.60,0.51,0.75,0.57,labelped,"NDC");
      pt->SetBorderSize(0);
      pt->SetTextColor(kRed);
      pt->SetTextSize(0.75);
      pt->SetFillColor(0);
      pt->Draw();
      char labelcut[64];
      sprintf(labelcut,"ADC_Cut = %3.0f",upADCcut);
      TPaveLabel *ptcut = new TPaveLabel(0.60,0.45,0.75,0.51,labelcut,"NDC");
      ptcut->SetBorderSize(0);
      ptcut->SetTextColor(kMagenta);
      ptcut->SetTextSize(0.75);
      ptcut->SetFillColor(0);
      ptcut->Draw();
      gPad->Update();
      double pedMean = fpedpar[1];
      cout<<"pedMean:"<<pedMean<<endl;

//      TCut cut_wadc = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&P.upQadcL>%f)",upADCcut);
      TCut cut_wadc = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&P.upQadcL>%f)",upADCcut);

      c1->cd(2);
      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.35);
      gStyle->SetStatH(0.08);

      TH1F* hQsqA = new TH1F("hQsqA",Form("LHRS Qsq from T-tree w/ T1 trigger (run%d,%s)",run,target.Data()), 100, 0.002, 0.013);
      TH1F* hQsqS = new TH1F("hQsqS",Form("LHRS Qsq from T-tree w/ T1 trigger (run%d,%s)",run,target.Data()), 100, 0.002, 0.013);

      TA->Draw("EK_L.Q2>>hQsqA", cut_wadc,"goff");
      TS->Draw("EK_L.Q2>>hQsqS", cut_wadc,"goff");
      hQsqA->SetXTitle("Qsq(GeV/c)^{2}");
      hQsqA->SetLineColor(1);
      hQsqS->SetLineColor(2);
      hQsqA->Draw();
      hQsqS->Draw("sames");
      
      gPad->Update();
      TPaveStats* statA = (TPaveStats*)hQsqA->FindObject("stats");
      TPaveStats* statS = (TPaveStats*)hQsqS->FindObject("stats");
      statA->SetY2NDC(0.90);
      statA->SetY1NDC(0.75);
      statS->SetY2NDC(0.75);
      statS->SetY1NDC(0.60);
      statA->SetTextColor(1);
      statS->SetTextColor(2);
      gPad->Modified();

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(1);
      latex.DrawLatex(0.65,0.55,"DB1");
      latex.SetTextColor(2);
      latex.DrawLatex(0.65,0.50,"DB2");

      outfile<<date.Data()<<"\t"<<target.Data()<<"\t"<<run<<"\t"<<hQsqA->GetMean()<<"\t"<<hQsqA->GetMeanError()<<"\t"<<hQsqS->GetMean()<<"\t"<<hQsqS->GetMeanError()<<endl;

      c1->SaveAs(Form("./temp/LHRS_Qsq_%s_run%d.pdf",target.Data(),run));
    }else{
      ofstream outfile("out_fileR.txt",ios_base::app);
      // RHRS
      double upADCcut_approx = 500;
      TA->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",run));
      TS->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/fewHoleDataBase/prexRHRS_%d_-1*.root",run));

      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.25);
      gStyle->SetStatH(0.10);

      TCanvas* c1=new TCanvas("c1","c1",1000,600);
      c1->Divide(2,1);
      c1->cd(1);
      TH1F* huq = new TH1F("huq", "P.upQadcR;ADC CH;Events/CH", 200, 450, 650);
      TA->Draw("P.upQadcR>>huq");

      TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
      huqCpy->GetXaxis()->SetRangeUser(upADCcut_approx-10,upADCcut_approx+10);
      cout<<"ADC cut: "<<huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin())<<endl;
      double upADCcut = huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin());
      TLine* line = new TLine(upADCcut,0.0,upADCcut,huq->GetMaximum());
      line->SetLineColor(kMagenta);
      line->Draw(); 
      
      TF1 *fped=new TF1("fped","gaus",483,496);
      fped->SetLineWidth(2);
      fped->SetLineColor(kRed);
      huq->Fit("fped","R");
      Double_t fpedpar[3];
      fped->GetParameters(&fpedpar[0]);
      char labelped[64];
      sprintf(labelped,"Pedeatal= %4.2f",fpedpar[1]);
      TPaveLabel *pt = new TPaveLabel(0.60,0.46,0.75,0.52,labelped,"NDC");
      pt->SetBorderSize(0);
      pt->SetTextColor(kRed);
      pt->SetTextSize(0.75);
      pt->SetFillColor(0);
      pt->Draw();
      char labelcut[64];
      sprintf(labelcut,"ADC_Cut = %3.0f",upADCcut);
      TPaveLabel *ptcut = new TPaveLabel(0.60,0.40,0.75,0.46,labelcut,"NDC");
      ptcut->SetBorderSize(0);
      ptcut->SetTextColor(kMagenta);
      ptcut->SetTextSize(0.75);
      ptcut->SetFillColor(0);
      ptcut->Draw();
      gPad->Update();
      double pedMean = fpedpar[1];
      cout<<"pedMean:"<<pedMean<<endl;

//      TCut cut_wadc = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&P.upQadcR>%f)",upADCcut);
      TCut cut_wadc = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&P.upQadcR>%f)",upADCcut);

      c1->cd(2);
      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.35);
      gStyle->SetStatH(0.08);

      TH1F* hQsqA = new TH1F("hQsqA", Form("RHRS Qsq from T-tree w/ T1 trigger (run%d,%s)",run,target.Data()), 100, 0.002, 0.013);
      TH1F* hQsqS = new TH1F("hQsqS", Form("RHRS Qsq from T-tree w/ T1 trigger (run%d,%s)",run,target.Data()), 100, 0.002, 0.013);

      TA->Draw("EK_R.Q2>>hQsqA", cut_wadc,"goff");
      TS->Draw("EK_R.Q2>>hQsqS", cut_wadc,"goff");
      hQsqA->SetXTitle("Qsq(GeV/c)^{2}");
      hQsqA->SetLineColor(1);
      hQsqS->SetLineColor(2);
      hQsqA->Draw();
      hQsqS->Draw("sames");

      gPad->Update();
      TPaveStats* statA = (TPaveStats*)hQsqA->FindObject("stats");
      TPaveStats* statS = (TPaveStats*)hQsqS->FindObject("stats");
      statA->SetY2NDC(0.90);
      statA->SetY1NDC(0.75);
      statS->SetY2NDC(0.75);
      statS->SetY1NDC(0.60);
      statA->SetTextColor(1);
      statS->SetTextColor(2);
      gPad->Modified();
   
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(1);
      latex.DrawLatex(0.65,0.55,"DB1");
      latex.SetTextColor(2);
      latex.DrawLatex(0.65,0.50,"DB2");

      outfile<<date.Data()<<"\t"<<target.Data()<<"\t"<<run<<"\t"<<hQsqA->GetMean()<<"\t"<<hQsqA->GetMeanError()<<"\t"<<hQsqS->GetMean()<<"\t"<<hQsqS->GetMeanError()<<endl;

      c1->SaveAs(Form("./temp/RHRS_Qsq_%s_run%d.pdf",target.Data(),run));
    }

}
