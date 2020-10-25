#include "langau.h"
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
void spectraUS(int run_num){ 

  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
   
  // Cut for FBUS
  double upadc_cutL_approx = 485;
  double upadc_cutR_approx = 505;
  // RHRS = R.*
  // LHRS = L.*
//  int run_num = (int)T->GetMaximum("fEvtHdr.fRun");
      TChain *T = new TChain("T");

  if(run_num < 10000){
      //LHRS
      T->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run_num));
      TCut trig_cut = "";
//      TCut trig_cut = "fEvtHdr.fEvtType==1";
//      TCut trig_cut = "(P.evtypebits&64)==64";
      // only one cluster in each VDC plane
      TCut vdc_cut = "L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1";
      //track cut on theta and phi
      TCut tr_cut = "(L.tr.th[0]<0.05&&L.tr.th[0]>-0.2&&L.tr.ph[0]<0.1&&L.tr.ph[0]>-0.1)";
      //track cut on target
      TCut tg_cut = "(L.tr.tg_th[0]<0.055&&L.tr.tg_th[0]>-0.055&&L.tr.tg_ph[0]>-0.018&&L.tr.tg_ph[0]<0.026)";

      TCanvas* cadc = new TCanvas("cadc","cadc",700,500);
      TH1F* huq = new TH1F("huq","P.upQadcL",300,400,700);
      T->Draw("P.upQadcL>>huq");
      huq->SetTitle("P.upQadcL;P.upQadcL raw;Events/CH");
      double ped_peak = huq->GetBinCenter(huq->GetMaximumBin());
      TF1 *fped=new TF1("fped","gaus",ped_peak-5,ped_peak+5);
      fped->SetLineWidth(2);
      fped->SetLineColor(kRed);
      huq->Fit("fped","R");
      Double_t fpedpar[3];
      fped->GetParameters(&fpedpar[0]);
      char labelped[64];
      sprintf(labelped,"RMS/Mean = %4.2f",fpedpar[2]/fpedpar[1]);
      TPaveLabel *pt = new TPaveLabel(0.60,0.46,0.75,0.52,labelped,"NDC");
      pt->SetBorderSize(0);
      pt->SetTextColor(kRed);
      pt->SetTextSize(0.75);
      pt->SetFillColor(0);
      pt->Draw();
      double pedMean = fpedpar[1];
      cout<<"pedMean:"<<pedMean<<endl;
      TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
      huqCpy->GetXaxis()->SetRangeUser(upadc_cutL_approx-10,upadc_cutL_approx+10);
      double upadc_cutL = huqCpy->GetBinCenter(huqCpy->GetMinimumBin());
      TLine* ped_line = new TLine(upadc_cutL,0.0,upadc_cutL,huq->GetMaximum());
      ped_line->SetLineColor(6);
      ped_line->SetLineWidth(2);
      ped_line->Draw();

      // FBUS adccuts
      TCut adc_cut_up = Form("P.upQadcL> %f", upadc_cutL);
      //FBUS <adccuts
      TCut adc_cut_up_ped = Form("P.upQadcL< %f", upadc_cutL);
      //cut on radiative tail
//      TCut x_cut = "(L.tr.x[0]+0.9*L.tr.th[0]) > -0.0665";
      TCut x_cut = "";
      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;
      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;
      
      TCanvas* c0 = new TCanvas("c0","dp LHRS",1000,700);
      TH1F* h;
//      T->Draw("L.gold.dp[0]>>h(500,-0.015,0.0)",cut_wadc);
      T->Draw("L.gold.dp[0]>>h(300,-0.015,0.0)",cut);
      h = (TH1F*)gDirectory->FindObject("h");
      h->SetDirectory(gROOT);
      h->SetTitle(Form("dp/p LHRS (run%d) with us_adcCut;#frac{dp}{p}",run_num));
      h->Draw();
      c0->SaveAs(Form("plots/dp_distribution_run%d_upstream_adcCut.pdf",run_num));

      TCanvas* c2 = new TCanvas("dp LHRS","dp LHRS",1000,700);
      T->Draw("L.tr.x[0]+0.9*L.tr.th[0]:L.gold.dp[0]>>hdp(500,-0.06,0.04,200,-1.8,1.0)",cut,"prof");
      TH2D* hdp = (TH2D*)gDirectory->FindObject("hdp");
      hdp->SetLineColor(kBlue);
      hdp->SetTitle(Form("Projected x on Det plane vs dp/p (run%d);#frac{dp}{p};Projected x on Det plane",run_num));
      TF1* fit = new TF1("fit","pol1",-0.04,-0.0025);
      fit->SetLineWidth(2);
      fit->SetLineColor(kRed);
      hdp->Fit("fit","R");
      Double_t fitpar[2];
      fit->GetParameters(&fitpar[0]);
      Double_t slope = fitpar[1];
//      Double_t slope = 14.35;
      cout<<Form("Dispersive Constant is %f per 1%s dp/p",slope,"%")<<endl;
      c2->SaveAs(Form("plots/det_X_vs_dp_run%d_upstream.pdf",run_num)); 

      TCanvas* c2FP = new TCanvas("dp LHRS FP","dp LHRS FP",1000,700);
      T->Draw("L.tr.x[0]:L.gold.dp[0]>>hdpFP(500,-0.06,0.04,200,-1.8,1.0)",cut,"prof");
      TH2D* hdpFP = (TH2D*)gDirectory->FindObject("hdpFP");
      hdpFP->SetLineColor(kBlue);
      hdpFP->SetTitle(Form("Dispersive x on focal plane vs dp/p (run%d);#frac{dp}{p};Dispersive x on focal plane",run_num));
      TF1* fitFP = new TF1("fitFP","pol1",-0.04,-0.0025);
      fitFP->SetLineWidth(2);
      fitFP->SetLineColor(kRed);
      hdpFP->Fit("fitFP","R");
      Double_t fitFPpar[2];
      fitFP->GetParameters(&fitFPpar[0]);
      Double_t slopeFP = fitFPpar[1];
//      Double_t slopeFP = 14.35;
      cout<<Form("Dispersive Constant at FP is %f per 1%s dp/p",slopeFP,"%")<<endl;
      c2FP->SaveAs(Form("plots/det_X_vs_dp_FP_run%d_upstream.pdf",run_num)); 

      TCanvas* c3 = new TCanvas("dp LHRS Scatter","dp LHRS Scatter",1000,700);
      TH2F* hdps = new TH2F("hdps", "Projected x on Det plane vs dp/p", 500, -0.02, 0.0, 200, -0.13, 0.02);
      T->Project(hdps->GetName(),"L.tr.x[0]+0.9*L.tr.th[0]:L.gold.dp[0]",cut);
      hdps->SetLineColor(kBlue);
      hdps->SetTitle(Form("Projected x on Det plane vs dp/p (run%d);#frac{dp}{p};Projected x on Det plane",run_num));
      hdps->Draw();
      c3->SaveAs(Form("plots/det_X_vs_dp_Scatter_run%d_upstream.pdf",run_num)); 

      TCanvas* c1 = new TCanvas("spectra LHRS","spectra LHRS",1000,700);
      TH1F* hx1 = new TH1F("hx1", Form("Projected x on detector plane (run%d)",run_num), 200, -0.10, 0.07);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc>cut", 200, -0.10, 0.07);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc<cut", 200, -0.10, 0.07);

      T->Project(hx1->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]+0.0335", cut); 
      T->Project(hx2->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]+0.0335", cut_wadc);
      T->Project(hx3->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]+0.0335", cut_wadc_ped);
      hx1->SetXTitle("L.tr.x[0]+0.9*L.tr.th[0] (m)");
      hx1->SetLineColor(kBlack);
      hx2->SetLineColor(kRed);
      hx3->SetLineColor(kBlue);
      hx1->Draw();
      hx2->Draw("same");
      hx3->Draw("same");
      int maxEntry = hx1->GetMaximum();
      TLine *line[5];
      Double_t inep0 = -slope*0.0/9.508/100;
      Double_t inep1 = -slope*2.615/9.508/100;
      Double_t inep2 = -slope*3.198/9.508/100;
      Double_t inep3 = -slope*3.475/9.508/100;
      Double_t inep4 = -slope*3.708/9.508/100;
      line[0] = new TLine(inep1,0.0,inep1,maxEntry+100);
      line[1] = new TLine(inep2,0.0,inep2,maxEntry+100);
      line[2] = new TLine(inep3,0.0,inep3,maxEntry+100);
      line[3] = new TLine(inep4,0.0,inep4,maxEntry+100);
      line[4] = new TLine(inep0,0.0,inep0,maxEntry+100);
      line[0]->SetLineColor(kGreen);
      line[1]->SetLineColor(kGreen+3);
      line[2]->SetLineColor(kCyan);
      line[3]->SetLineColor(kOrange+7);
      line[4]->SetLineColor(kMagenta);
      TPaveLabel* leb[6];
      leb[4] = new TPaveLabel(0.18,0.84,0.25,0.88,"Spectra Pb-208","NDC");
      leb[5] = new TPaveLabel(0.18,0.76,0.25,0.80,"Elastic","NDC");
      leb[0] = new TPaveLabel(0.18,0.72,0.25,0.76,"2.615MeV","NDC");
      leb[1] = new TPaveLabel(0.18,0.68,0.25,0.72,"3.198MeV","NDC");
      leb[2] = new TPaveLabel(0.18,0.64,0.25,0.68,"3.475MeV","NDC");
      leb[3] = new TPaveLabel(0.18,0.60,0.25,0.64,"3.708MeV","NDC");
      leb[0]->SetTextColor(kGreen);
      leb[1]->SetTextColor(kGreen+3);
      leb[2]->SetTextColor(kCyan);
      leb[3]->SetTextColor(kOrange+7);
      leb[4]->SetTextColor(kBlack);
      leb[5]->SetTextColor(kMagenta);
    for(int i=0;i<6;i++){
      leb[i]->SetBorderSize(0);
      leb[i]->SetTextSize(0.75);
      leb[i]->SetFillColor(0);
      leb[i]->Draw();
     }
    for(int i=0;i<5;i++)
      line[i]->Draw();
      c1->SaveAs(Form("plots/spectra_Lead208_run%d_upstream.pdf",run_num));
      
      TCanvas* c11 = new TCanvas("Pb208 spectra LHRS","Pb208 spectra LHRS",1000,700);
      TH1F* hx11 = new TH1F("hx11", Form("Pb208 spectrum on detector plane (run%d)",run_num), 200, 0.940, 0.9525);
      TH1F* hx22 = new TH1F("hx22", "Pb208 spectra on detector plane + adc>cut", 200, 0.940, 0.9525);
      TH1F* hx33 = new TH1F("hx33", "Pb208 spectra on detector plane + adc<cut", 200, 0.940, 0.9525);

      T->Project(hx11->GetName(), Form("0.9508*((L.tr.x[0]+0.9*L.tr.th[0])/%f+1)",slope), cut); 
      T->Project(hx22->GetName(), Form("0.9508*((L.tr.x[0]+0.9*L.tr.th[0])/%f+1)",slope), cut_wadc);
      T->Project(hx33->GetName(), Form("0.9508*((L.tr.x[0]+0.9*L.tr.th[0])/%f+1)",slope), cut_wadc_ped);
      hx11->SetXTitle("Energy(GeV)");
      hx11->SetLineColor(kBlack);
      hx22->SetLineColor(kRed);
      hx33->SetLineColor(kBlue);
      hx11->Draw();
      hx22->Draw("same");
      hx33->Draw("same");
      TLine* ln[6];
      int max = hx11->GetMaximum()+100;
      int maxBin = hx11->GetMaximumBin();
      float elasticPeak = hx11->GetXaxis()->GetBinCenter(maxBin);
      ln[0] = new TLine(0.9508,0.0,0.9508,max);
      ln[1] = new TLine(elasticPeak,0.0,elasticPeak,max);
      ln[2] = new TLine(elasticPeak-2.614/1000.,0.0,elasticPeak-2.614/1000.,max);
      ln[3] = new TLine(elasticPeak-3.198/1000.,0.0,elasticPeak-3.198/1000.,max);
      ln[4] = new TLine(elasticPeak-3.475/1000.,0.0,elasticPeak-3.475/1000.,max);
      ln[5] = new TLine(elasticPeak-3.708/1000.,0.0,elasticPeak-3.708/1000.,max);
      ln[0]->SetLineColor(kViolet+2);
      ln[1]->SetLineColor(kMagenta);
      ln[2]->SetLineColor(kGreen);
      ln[3]->SetLineColor(kGreen+3);
      ln[4]->SetLineColor(kCyan);
      ln[5]->SetLineColor(kOrange+7);
      for(int i=0;i<6;i++)
      ln[i]->Draw();
      TPaveLabel* beamE = new TPaveLabel(0.18,0.80,0.25,0.84,"Beam E","NDC");
      beamE->SetTextColor(kViolet+2);
      beamE->SetBorderSize(0);
      beamE->SetTextSize(0.75);
      beamE->SetFillColor(0);
      for(int i=0;i<6;i++)
      leb[i]->Draw();
      beamE->Draw();
      c11->SaveAs(Form("./plots/Pb208_E_spectra_run%d.pdf",run_num));
    }else{
      // RHRS
      T->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",run_num));
      TCut trig_cut = "";
//      TCut trig_cut = "fEvtHdr.fEvtType==1";
//      TCut trig_cut = "(P.evtypebits&64)==64";
      TCut vdc_cut = "R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1";
      TCut tr_cut = "(R.tr.th[0]<0.05&&R.tr.th[0]>-0.2&&R.tr.ph[0]<0.1&&R.tr.ph[0]>-0.1)";
      TCut tg_cut = "(R.tr.tg_th[0]<0.05&&R.tr.tg_th[0]>-0.055&&R.tr.tg_ph[0]>-0.026&&R.tr.tg_ph[0]<0.016)";

      TCanvas* cadc = new TCanvas("cadc","cadc",700,500);
      TH1F* huq = new TH1F("huq", "P.upQadcR", 300, 400, 700);
      T->Project(huq->GetName(),"P.upQadcR");
      huq->SetTitle("P.upQadcR raw;ADC CH");
      huq->SetYTitle("Events/CH");
      double ped_peak = huq->GetBinCenter(huq->GetMaximumBin());
      TF1 *fped=new TF1("fped","gaus",ped_peak-5,ped_peak+5);
      fped->SetLineWidth(2);
      fped->SetLineColor(kRed);
      huq->Fit("fped","R");
      Double_t fpedpar[3];
      fped->GetParameters(&fpedpar[0]);
      char labelped[64];
      sprintf(labelped,"RMS/Mean = %4.2f",fpedpar[2]/fpedpar[1]);
      TPaveLabel *pt = new TPaveLabel(0.60,0.46,0.75,0.52,labelped,"NDC");
      pt->SetBorderSize(0);
      pt->SetTextColor(kRed);
      pt->SetTextSize(0.75);
      pt->SetFillColor(0);
      pt->Draw();
      double pedMean = fpedpar[1];
      cout<<"pedMean:"<<pedMean<<endl;
      TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
      huqCpy->GetXaxis()->SetRangeUser(upadc_cutR_approx-10,upadc_cutR_approx+10);
      double upadc_cutR = huqCpy->GetBinCenter(huqCpy->GetMinimumBin());
      TLine* ped_line = new TLine(upadc_cutR,0.0,upadc_cutR,huq->GetMaximum());
      ped_line->SetLineColor(6);
      ped_line->SetLineWidth(2);
      ped_line->Draw();


      // FBUS adccuts
      TCut adc_cut_up = Form("P.upQadcR> %f", upadc_cutR);
      //FBUS <adccuts
      TCut adc_cut_up_ped = Form("P.upQadcR< %f", upadc_cutR);

//      TCut x_cut = "(R.tr.x[0]+0.9*R.tr.th[0]) > -0.072";
      TCut x_cut = "";
      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;
      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;

      TCanvas* c0 = new TCanvas("c0","dp RHRS",1000,700);
      TH1F* h;
      T->Draw("R.gold.dp[0]>>h(500,-0.015,0.002)",cut);
      h = (TH1F*)gDirectory->FindObject("h");
      h->SetDirectory(gROOT);
      h->SetTitle(Form("dp/p RHRS (run%d);#frac{dp}{p}",run_num));
      h->Draw();
      c0->SaveAs(Form("plots/dp_distribution_run%d_upstream.pdf",run_num));

      TCanvas* c2 = new TCanvas("dp RHRS","dp RHRS",1000,700);
//      TH2F* hdp = new TH2F("hdp", "Projected x on Det plane vs dp ", 500, -0.015, 0.0, 200, -0.15, 0.0);
      T->Draw("R.tr.x[0]+0.9*R.tr.th[0]:R.gold.dp[0]>>hdp(500,-0.05,0.01,200,-0.15,0.0)",cut,"prof");
      TH2F* hdp = (TH2F*)gDirectory->FindObject("hdp");
      hdp->SetLineColor(kBlue);
      hdp->SetTitle(Form("Projected x on Det plane vs dp/p (run%d);#frac{dp}{p};Projected x on Det plane",run_num));
      TF1* fit = new TF1("fit","pol1",-0.035,-0.002);
      fit->SetLineWidth(2);
      fit->SetLineColor(kRed);
      hdp->Fit("fit","R");
      Double_t fitpar[2];
      fit->GetParameters(&fitpar[0]);
      Double_t slope = fitpar[1];
//      Double_t slope = 14.26;
      cout<<Form("Dispersive Constant is %f per 1%s dp/p",slope,"%")<<endl;
      c2->SaveAs(Form("plots/det_X_vs_dp_run%d_upstream.pdf",run_num)); 

      TCanvas* c3 = new TCanvas("dp RHRS Scatter","dp RHRS Scatter",1000,700);
      TH2F* hdps = new TH2F("hdps", "Projected x on Det plane vs dp/p ", 500, -0.05, 0.01, 200, -0.15, 0.0);
      T->Project(hdps->GetName(),"R.tr.x[0]+0.9*R.tr.th[0]:R.gold.dp[0]",cut);
      hdps->SetLineColor(kBlue);
      hdps->SetTitle(Form("Projected x on Det plane vs dp/p (run%d);#frac{dp}{p};Projected x on Det plane",run_num));
      hdps->Draw();
      c3->SaveAs(Form("plots/det_X_vs_dp_scatter_run%d_upstream.pdf",run_num));

      TCanvas* c1 = new TCanvas("spectra RHRS","spectra RHRS",1000,700);
//      gPad->SetGridx();
      TH1F* hx1 = new TH1F("hx1", Form("Projected x on detector plane (run%d)",run_num), 200, -0.10, 0.07);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc> cut", 200, -0.10, 0.07);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc< cut", 200, -0.10, 0.07);

      T->Project(hx1->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]+0.0425", cut); 
      T->Project(hx2->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]+0.0425", cut_wadc);
      T->Project(hx3->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]+0.0425", cut_wadc_ped);
      hx1->SetXTitle("R.tr.x[0]+0.9*R.tr.th[0] (m)");
      hx1->SetLineColor(kBlack);
      hx2->SetLineColor(kRed);
      hx3->SetLineColor(kBlue);
      hx1->Draw();
      hx2->Draw("same");
      hx3->Draw("same");
      int maxEntry = hx1->GetMaximum();
      TLine *line[5];
      Double_t inep0 = -slope*0.0/9.508/100;
      Double_t inep1 = -slope*2.614/9.508/100;
      Double_t inep2 = -slope*3.198/9.508/100;
      Double_t inep3 = -slope*3.475/9.508/100;
      Double_t inep4 = -slope*3.708/9.508/100;
      line[0] = new TLine(inep1,0.0,inep1,maxEntry+100);
      line[1] = new TLine(inep2,0.0,inep2,maxEntry+100);
      line[2] = new TLine(inep3,0.0,inep3,maxEntry+100);
      line[3] = new TLine(inep4,0.0,inep4,maxEntry+100);
      line[4] = new TLine(inep0,0.0,inep0,maxEntry+100);
      line[0]->SetLineColor(kGreen);
      line[1]->SetLineColor(kGreen+3);
      line[2]->SetLineColor(kCyan);
      line[3]->SetLineColor(kOrange+7);
      line[4]->SetLineColor(kMagenta);
      TPaveLabel* leb[6];
      leb[4] = new TPaveLabel(0.18,0.84,0.25,0.88,"Spectra Pb208","NDC");
      leb[5] = new TPaveLabel(0.18,0.76,0.25,0.80,"Elastic","NDC");
      leb[0] = new TPaveLabel(0.18,0.72,0.25,0.76,"2.614MeV","NDC");
      leb[1] = new TPaveLabel(0.18,0.68,0.25,0.72,"3.198MeV","NDC");
      leb[2] = new TPaveLabel(0.18,0.64,0.25,0.68,"3.475MeV","NDC");
      leb[3] = new TPaveLabel(0.18,0.60,0.25,0.64,"3.708MeV","NDC");
      leb[0]->SetTextColor(kGreen);
      leb[1]->SetTextColor(kGreen+3);
      leb[2]->SetTextColor(kCyan);
      leb[3]->SetTextColor(kOrange+7);
      leb[4]->SetTextColor(kBlack);
    for(int i=0;i<6;i++){
      leb[i]->SetBorderSize(0);
      leb[i]->SetTextSize(0.75);
      leb[i]->SetFillColor(0);
      leb[i]->Draw();
     }
    for(int i=0;i<5;i++)
      line[i]->Draw();

      c1->SaveAs(Form("plots/spectra_Lead208_run%d_upstream.pdf",run_num));

     TCanvas* c11 = new TCanvas("Pb208 spectra RHRS","Pb208 spectra RHRS",1000,700);
      TH1F* hx11 = new TH1F("hx11", Form("Pb208 spectrum on detector plane (run%d)",run_num), 200, 0.940, 0.9525);
      TH1F* hx22 = new TH1F("hx22", "Pb208 spectra on detector plane + adc>cut", 200, 0.940, 0.9525);
      TH1F* hx33 = new TH1F("hx33", "Pb208 spectra on detector plane + adc<cut", 200, 0.940, 0.9525);

      T->Project(hx11->GetName(), Form("0.9508*((R.tr.x[0]+0.9*R.tr.th[0])/%f+1)",slope), cut); 
      T->Project(hx22->GetName(), Form("0.9508*((R.tr.x[0]+0.9*R.tr.th[0])/%f+1)",slope), cut_wadc);
      T->Project(hx33->GetName(), Form("0.9508*((R.tr.x[0]+0.9*R.tr.th[0])/%f+1)",slope), cut_wadc_ped);
      hx11->SetXTitle("Energy(GeV)");
      hx11->SetLineColor(kBlack);
      hx22->SetLineColor(kRed);
      hx33->SetLineColor(kBlue);
      hx11->Draw();
      hx22->Draw("same");
      hx33->Draw("same");
      TLine* ln[6];
      int max = hx11->GetMaximum()+100;
      int maxBin = hx11->GetMaximumBin();
      float elasticPeak = hx11->GetXaxis()->GetBinCenter(maxBin);
      ln[0] = new TLine(0.9508,0.0,0.9508,max);
      ln[1] = new TLine(elasticPeak,0.0,elasticPeak,max);
      ln[2] = new TLine(elasticPeak-2.614/1000.,0.0,elasticPeak-2.614/1000.,max);
      ln[3] = new TLine(elasticPeak-3.198/1000.,0.0,elasticPeak-3.198/1000.,max);
      ln[4] = new TLine(elasticPeak-3.475/1000.,0.0,elasticPeak-3.475/1000.,max);
      ln[5] = new TLine(elasticPeak-3.708/1000.,0.0,elasticPeak-3.708/1000.,max);
      ln[0]->SetLineColor(kViolet+2);
      ln[1]->SetLineColor(kMagenta);
      ln[2]->SetLineColor(kGreen);
      ln[3]->SetLineColor(kGreen+3);
      ln[4]->SetLineColor(kCyan);
      ln[5]->SetLineColor(kOrange+7);
      for(int i=0;i<6;i++)
      ln[i]->Draw();
      TPaveLabel* beamE = new TPaveLabel(0.18,0.80,0.25,0.84,"Beam E","NDC");
      beamE->SetTextColor(kViolet+2);
      beamE->SetBorderSize(0);
      beamE->SetTextSize(0.75);
      beamE->SetFillColor(0);
      for(int i=0;i<6;i++)
      leb[i]->Draw();
      beamE->Draw();
      c11->SaveAs(Form("./plots/Pb208_E_spectra_run%d.pdf",run_num));
}

}
