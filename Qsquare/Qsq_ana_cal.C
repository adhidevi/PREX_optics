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
void Qsq_ana_cal(int run, TString target, TString date){ 
  gROOT->Reset();
  gStyle->SetOptStat(1002211);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
      
      Float_t Ebeam = 0.9508;//GeV
      Float_t theta0 = 4.74;//degrees
      theta0 = theta0*3.1415926/180;//radians

      Float_t s0 = TMath::Sin(theta0);
      Float_t c0 = TMath::Cos(theta0);

      TChain *T = new TChain("T");
  if(run < 10000){
      //LHRS
      double upADCcut_approx = 480;
      T->Add(Form("/chafs1/work1/prex_counting/chandan/prexLHRS_%d_-1*.root",run));
      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.25);
      gStyle->SetStatH(0.10);
      TCanvas *c1 = new TCanvas("c1","c1",1000,600);
      c1->Divide(2,1);
      c1->cd(1);
      TH1F* huq = new TH1F("huq","P.upQadcL;ADC CH;Event/CH",250,400,650);
      T->Draw("P.upQadcL>>huq");
      
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
      TH1F* hQsq = new TH1F("hQsq",Form("LHRS Qsq w/ T1 trigger (run%d,%s);Qsq (GeV/c)^{2}",run,target.Data()), 100, 0.002, 0.013);
      TH1F* hQsq_cal = new TH1F("hQsq_cal",Form("LHRS Qsq w/ T1 trigger (run%d,%s);Qsq (GeV/c)^{2}",run,target.Data()), 100, 0.002, 0.013);

      hQsq->SetLineColor(kBlue);
      hQsq_cal->SetLineColor(kRed);

      T->Draw("EK_L.Q2>>hQsq", cut_wadc);

      T->Draw(Form("2*%f*L.gold.p*(1-((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))))>>hQsq_cal",Ebeam,c0,s0),cut_wadc,"sames");

      gPad->Update();

      TPaveStats* ana = (TPaveStats*)hQsq->FindObject("stats");
      TPaveStats* cal = (TPaveStats*)hQsq_cal->FindObject("stats");
      ana->SetY1NDC(0.90);
      ana->SetY2NDC(0.75);
      cal->SetY1NDC(0.75);
      cal->SetY2NDC(0.60);
      ana->SetTextColor(kBlue);
      cal->SetTextColor(kRed);
      gPad->Modified();
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(kBlue);
      latex.DrawLatex(0.57,0.55,"From Analyzer");
      latex.SetTextColor(kRed);
      latex.DrawLatex(0.57,0.50,"My Calculation");

      ofstream out_analyzer("./TextFiles/qsq_analyzerL.csv",ios_base::app);
      ofstream out_calculated("./TextFiles/qsq_calculatedL.csv",ios_base::app);
      out_analyzer<<date.Data()<<"\t"<<target.Data()<<"\t"<<run<<"\t"<<hQsq->GetMean()<<"\t"<<hQsq->GetMeanError()<<endl;
      out_calculated<<date.Data()<<"\t"<<target.Data()<<"\t"<<run<<"\t"<<hQsq_cal->GetMean()<<"\t"<<hQsq_cal->GetMeanError()<<endl;
      out_analyzer.close();
      out_calculated.close();

      c1->SaveAs(Form("./temp/LHRS_Qsq_comp_%s_run%d.pdf",target.Data(),run));

    }else{
      // RHRS
      double upADCcut_approx = 502;
      T->Add(Form("/chafs1/work1/prex_counting/chandan/prexRHRS_%d_-1*.root",run));
      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.25);
      gStyle->SetStatH(0.10);
      TCanvas* c1=new TCanvas("c1","c1",1000,600);
      c1->Divide(2,1);
      c1->cd(1);
      TH1F* huq = new TH1F("huq", "P.upQadcR;ADC CH;Events/CH", 200, 450, 650);
      T->Draw("P.upQadcR>>huq");

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

      TH1F* hQsq = new TH1F("hQsq", Form("RHRS Qsq w/ T1 trigger (run%d,%s);Qsq (GeV/c)^{2}",run,target.Data()), 100, 0.002, 0.013);
      TH1F* hQsq_cal = new TH1F("hQsq_cal",Form("RHRS Qsq w/ T1 trigger (run%d,%s);Qsq (GeV/c)^{2}",run,target.Data()), 100, 0.002, 0.013);

      hQsq->SetLineColor(kBlue);
      hQsq_cal->SetLineColor(kRed);

      T->Draw("EK_R.Q2>>hQsq", cut_wadc);

      T->Draw(Form("2*%f*R.gold.p*(1-((%f+R.gold.ph*%f)/(TMath::Sqrt(1.0+R.gold.th*R.gold.th+R.gold.ph*R.gold.ph))))>>hQsq_cal",Ebeam,c0,s0),cut_wadc,"sames");

      gPad->Update();

      TPaveStats* ana = (TPaveStats*)hQsq->FindObject("stats");
      TPaveStats* cal = (TPaveStats*)hQsq_cal->FindObject("stats");

      ana->SetY1NDC(0.90);
      ana->SetY2NDC(0.75);
      cal->SetY1NDC(0.75);
      cal->SetY2NDC(0.60);
      ana->SetTextColor(kBlue);
      cal->SetTextColor(kRed);
      gPad->Modified();
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(kBlue);
      latex.DrawLatex(0.57,0.55,"From Analyzer");
      latex.SetTextColor(kRed);
      latex.DrawLatex(0.57,0.50,"My Calculation");

      ofstream out_analyzer("./TextFiles/qsq_analyzerR.csv",ios_base::app);
      ofstream out_calculated("./TextFiles/qsq_calculatedR.csv",ios_base::app);
      out_analyzer<<date.Data()<<"\t"<<target.Data()<<"\t"<<run<<"\t"<<hQsq->GetMean()<<"\t"<<hQsq->GetMeanError()<<endl;
      out_calculated<<date.Data()<<"\t"<<target.Data()<<"\t"<<run<<"\t"<<hQsq_cal->GetMean()<<"\t"<<hQsq_cal->GetMeanError()<<endl;
      out_analyzer.close();
      out_calculated.close();

      c1->SaveAs(Form("./temp/RHRS_Qsq_comp_%s_run%d.pdf",target.Data(),run));
    }

}
