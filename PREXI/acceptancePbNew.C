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
void acceptancePbNew(int run_num, TString target){ 

  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
   
//  double upadc_cutL = 480;//for run2886
  double upadc_cutL = 482;//for run2977
  double upadc_cutR = 502;

      TChain *T = new TChain("T");

  if(run_num > 10000){
      //LHRS
      T->Add(Form("/w/halla-scifs17exp/parity/disk1/ryanrich/prexIrootfiles/happexsp_%d.root",run_num));
//      T->Add(Form("/w/halla-scifs17exp/parity/disk1/bob/podd2016/analyzer/Afile_%d.root",run_num));
      TCut trig_cut = "";
//      TCut trig_cut = "(P.evtypebits&2)==2";
      // only one cluster in each VDC plane
      TCut vdc_cut = "L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1";
      //track cut on theta and phi
      TCut tr_cut = "(L.tr.th[0]<0.05&&L.tr.th[0]>-0.2&&L.tr.ph[0]<0.1&&L.tr.ph[0]>-0.1)";
      //track cut on target
      TCut tg_cut = "(L.tr.tg_th[0]<0.055&&L.tr.tg_th[0]>-0.055&&L.tr.tg_ph[0]>-0.018&&L.tr.tg_ph[0]<0.026)";

      // FBUS adccuts
      TCut adc_cut_up = Form("P.upQadcL> %f", upadc_cutL);
      //FBUS <adccuts
      TCut adc_cut_up_ped = Form("P.upQadcL< %f", upadc_cutL);
      //cut on radiative tail
      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;
      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;

      TCanvas* c2 = new TCanvas("Momentum LHRS","dp RHRS",1200,700);
      gPad->SetLogy();
      TH1F* hp = new TH1F("hp", Form("Momentum %s (run%d) without adcCut;E_{beam}(MeV)",target.Data(),run_num), 200, 1053, 1065);
      T->Project(hp->GetName(),"1063*(1.+L.tr.tg_dp)",cut);
      hp->SetLineColor(kBlue);
      hp->Draw();

//////fit gaussian on the elastic peak//////////////////////
      auto e_peak = hp->GetBinCenter(hp->GetMaximumBin());
      hp->Fit("gaus","R0","ep",e_peak-0.7,e_peak+0.7);
      double_t gCrystalPar[5];
      TF1* gCrystal = new TF1("gCrystal","crystalball",hp->GetFunction("gaus")->GetParameter(1)-6*hp->GetFunction("gaus")->GetParameter(2),
		hp->GetFunction("gaus")->GetParameter(1)+3*hp->GetFunction("gaus")->GetParameter(2));
      gCrystal->SetParameters(hp->GetFunction("gaus")->GetParameter(0), hp->GetFunction("gaus")->GetParameter(1),
		hp->GetFunction("gaus")->GetParameter(2),1.64,1.1615);
      hp->Fit("gCrystal","R+","ep",gCrystal->GetXmin(),gCrystal->GetXmax());
      gCrystal->GetParameters(gCrystalPar);
      TLine* ground = new TLine(hp->GetFunction("gaus")->GetParameter(1),0.0,hp->GetFunction("gaus")->GetParameter(1),hp->GetMaximum());
      TLine* fstIn = new TLine(hp->GetFunction("gaus")->GetParameter(1)-2.615,0.0,hp->GetFunction("gaus")->GetParameter(1)-2.615,hp->GetMaximum());
      ground->SetLineColor(2);
      fstIn->SetLineColor(6);
      ground->Draw();
      fstIn->Draw();      
}else{
      //RHRS
      T->Add(Form("/w/halla-scifs17exp/parity/disk1/ryanrich/prexIrootfiles/happexsp_%d.dat.0.root",run_num));
//      T->Add(Form("/w/halla-scifs17exp/parity/disk1/bob/podd2016/analyzer/Afile_%d.root",run_num));
      TCut trig_cut = "";
//      TCut trig_cut = "(P.evtypebits&2)==2";
      // only one cluster in each VDC plane
      TCut vdc_cut = "R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1";
      //track cut on theta and phi
      TCut tr_cut = "(R.tr.th[0]<0.005&&R.tr.th[0]>-0.02&&R.tr.ph[0]<0.02&&R.tr.ph[0]>-0.02)";
      //track cut on target
      TCut tg_cut = "(R.tr.tg_th[0]<0.05&&R.tr.tg_th[0]>-0.05&&R.tr.tg_ph[0]>-0.02&&R.tr.tg_ph[0]<0.02)";

      // FBUS adccuts
      TCut adc_cut_up = Form("P.upQadcR> %f", upadc_cutR);
      //FBUS <adccuts
      TCut adc_cut_up_ped = Form("P.upQadcR< %f", upadc_cutR);
      //cut on radiative tail
      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;
      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;

      TCanvas* c2 = new TCanvas("Momentum RHRS","dp RHRS",1200,700);
      gPad->SetLogy();
      TH1F* hp = new TH1F("hp", Form("Momentum %s (run%d) without adcCut;E_{beam}(MeV)",target.Data(),run_num), 200, 1053, 1065);
      T->Project(hp->GetName(),"1063*(1.+R.tr.tg_dp)",cut);
      hp->SetLineColor(kBlue);
      hp->Draw();

//////fit gaussian on the elastic peak//////////////////////
      auto e_peak = hp->GetBinCenter(hp->GetMaximumBin());
      hp->Fit("gaus","R0","ep",e_peak-0.7,e_peak+0.5);
      double_t gCrystalPar[5];
      TF1* gCrystal = new TF1("gCrystal","crystalball",hp->GetFunction("gaus")->GetParameter(1)-4*hp->GetFunction("gaus")->GetParameter(2),
		hp->GetFunction("gaus")->GetParameter(1)+3*hp->GetFunction("gaus")->GetParameter(2));
      gCrystal->SetParameters(hp->GetFunction("gaus")->GetParameter(0), hp->GetFunction("gaus")->GetParameter(1),
		hp->GetFunction("gaus")->GetParameter(2),1.64,1.1615);
      hp->Fit("gCrystal","R+","ep",gCrystal->GetXmin(),gCrystal->GetXmax());
      gCrystal->GetParameters(gCrystalPar);
      
      TLine* ground = new TLine(hp->GetFunction("gaus")->GetParameter(1),0.0,hp->GetFunction("gaus")->GetParameter(1),hp->GetMaximum());
      TLine* fstIn = new TLine(hp->GetFunction("gaus")->GetParameter(1)-2.615,0.0,hp->GetFunction("gaus")->GetParameter(1)-2.615,hp->GetMaximum());
      ground->SetLineColor(2);
      fstIn->SetLineColor(6);
      ground->Draw();
      fstIn->Draw();
}
}
