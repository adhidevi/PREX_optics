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
void calculated_only_qsq(int runL, int runR, TString target, TString date, Float_t Ebeam){ 
  gROOT->Reset();
  gStyle->SetOptStat(1002211);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
      
      TChain *TL = new TChain("T");
      TChain *TR = new TChain("T");
      //LHRS
      Float_t theta0L = 4.801;//degrees
      Float_t theta0R = 4.772;//degrees
      theta0L = theta0L*3.1415926/180;//radians
      theta0R = theta0R*3.1415926/180;//radians

      Float_t s0L = TMath::Sin(theta0L);
      Float_t s0R = TMath::Sin(theta0R);
      Float_t c0L = TMath::Cos(theta0L);
      Float_t c0R = TMath::Cos(theta0R);

      double upADCcutL_approx = 480;
      double upADCcutR_approx = 502;
      TL->Add(Form("/chafs1/work1/prex_counting/QsqRootFiles/prexLHRS_%d_-1*.root",runL));
      TR->Add(Form("/chafs1/work1/prex_counting/QsqRootFiles/prexRHRS_%d_-1*.root",runR));
      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.25);
      gStyle->SetStatH(0.10);
      TCanvas *c1 = new TCanvas("c1","c1",1000,600);
      c1->Divide(2,1);
      c1->cd(1);
      TH1F* huqL = new TH1F("huqL",Form("Upstream ADC (run%d, %d);ADC CH;Event/CH",runL,runR),250,400,650);
      TH1F* huqR = new TH1F("huqR",Form("Upstream ADC (run%d, %d);ADC CH;Event/CH",runL,runR),250,400,650);
      huqL->SetLineColor(1);
      huqR->SetLineColor(2);
      TL->Draw("P.upQadcL>>huqL","","hist");
      TR->Draw("P.upQadcR>>huqR","","sames hist");
      huqL->Scale(1./huqL->Integral());
      huqR->Scale(1./huqR->Integral());
      
      TH1F* huqLCpy = (TH1F*)huqL->Clone("huqLCpy");
      TH1F* huqRCpy = (TH1F*)huqR->Clone("huqRCpy");
      huqLCpy->GetXaxis()->SetRangeUser(upADCcutL_approx-10,upADCcutL_approx+10);
      huqRCpy->GetXaxis()->SetRangeUser(upADCcutR_approx-10,upADCcutR_approx+10);
      double upADCcutL = huqLCpy->GetXaxis()->GetBinCenter(huqLCpy->GetMinimumBin());
      double upADCcutR = huqRCpy->GetXaxis()->GetBinCenter(huqRCpy->GetMinimumBin());
      cout<<"LHRS ADC cut: "<<upADCcutL<<endl;
      cout<<"RHRS ADC cut: "<<upADCcutR<<endl;
      gPad->Update();
      TPaveStats* ADCl = (TPaveStats*)huqL->FindObject("stats");
      TPaveStats* ADCr = (TPaveStats*)huqR->FindObject("stats");
      ADCl->SetY1NDC(0.90);
      ADCl->SetY2NDC(0.75);
      ADCr->SetY1NDC(0.75);
      ADCr->SetY2NDC(0.60);
      ADCl->SetTextColor(1);
      ADCr->SetTextColor(2);
      gPad->Modified();

      TLine* lineL = new TLine(upADCcutL,0.0,upADCcutL,huqL->GetMaximum());
      TLine* lineR = new TLine(upADCcutR,0.0,upADCcutR,huqR->GetMaximum());
      lineL->SetLineColor(4);
      lineR->SetLineColor(6);
      lineL->Draw(); 
      lineR->Draw();
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.04);
      latex.SetTextColor(4);
      latex.DrawLatex(0.55,0.55,Form("PedL Cut = %.1f",upADCcutL));
      latex.SetTextColor(6);
      latex.DrawLatex(0.55,0.50,Form("PedR Cut = %.1f",upADCcutR));
      TCut cut_wadcL = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.04&&L.gold.ph<0.05&&L.gold.dp<-0.002&&P.upQadcL>%f)",upADCcutL);
      TCut cut_wadcR = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.04&&R.gold.ph<0.05&&R.gold.dp<-0.002&&P.upQadcR>%f)",upADCcutR);
      c1->cd(2);
      gStyle->SetStatX(0.90);
      gStyle->SetStatY(0.90);
      gStyle->SetStatW(0.35);
      gStyle->SetStatH(0.08);
      TH1F* hQsqL = new TH1F("hQsqL",Form("Qsq w/ T1 trigger (run%d, %d);Qsq (GeV/c)^{2}",runL,runR), 100, 0.002, 0.013);
      TH1F* hQsqR = new TH1F("hQsqR",Form("Qsq w/ T1 trigger (run%d, %d);Qsq (GeV/c)^{2}",runL,runR), 100, 0.002, 0.013);

      hQsqL->SetLineColor(1);
      hQsqR->SetLineColor(2);


      TL->Draw(Form("2*%f*L.gold.p*(1-((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))))>>hQsqL",Ebeam,c0L,s0L),cut_wadcL,"hist");
      TR->Draw(Form("2*%f*R.gold.p*(1-((%f+R.gold.ph*%f)/(TMath::Sqrt(1.0+R.gold.th*R.gold.th+R.gold.ph*R.gold.ph))))>>hQsqR",Ebeam,c0R,s0R),cut_wadcR,"sames hist");
      hQsqL->Scale(1./hQsqL->Integral());
      hQsqR->Scale(1./hQsqR->Integral());
      gPad->Update();

      TPaveStats* left = (TPaveStats*)hQsqL->FindObject("stats");
      TPaveStats* right = (TPaveStats*)hQsqR->FindObject("stats");
      left->SetY1NDC(0.90);
      left->SetY2NDC(0.75);
      right->SetY1NDC(0.75);
      right->SetY2NDC(0.60);
      left->SetTextColor(1);
      right->SetTextColor(2);
      gPad->Modified();
      latex.SetTextColor(1);
      latex.DrawLatex(0.57,0.55,"LHRS");
      latex.SetTextColor(2);
      latex.DrawLatex(0.57,0.50,"RHRS");

      ofstream outLR("./TextFiles/qsq_LR.csv",ios_base::app);
      outLR<<date.Data()<<"\t"<<target.Data()<<"\t"<<runL<<"\t"<<hQsqL->GetMean()<<"\t"<<hQsqL->GetMeanError()<<"\t"<<runR<<"\t"<<hQsqR->GetMean()<<"\t"<<hQsqR->GetMeanError()<<endl;
      outLR.close();

      c1->SaveAs(Form("./temp/Qsq_comp_LR_%d_run%d.pdf",runL,runR));
}
