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
void alignUS(int run_num){ 

  gROOT->Reset();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
   
  const double z_scale = 100;
  // Cut for FBUS
  double upadc_cutL_approx = 481;
  double upadc_cutR_approx = 510;
  double USLPMT = 5370;
  double USLgain = 7.5;//7.5e5
  double USRPMT = 5401;
  double USRgain = 6.46;//6.46e5
  // RHRS = R.*
  // LHRS = L.*
//  int run_num = (int)T->GetMaximum("fEvtHdr.fRun");
      TChain *T = new TChain("T");

  if(run_num < 10000){
      //LHRS
       T->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run_num));
//      TCut trig_cut = "";
      TCut trig_cut = "(P.evtypebits&2)==2";
      // only one cluster in each VDC plane
      TCut vdc_cut = "L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1";
      //track cut on theta and phi
      TCut tr_cut = "(L.tr.th[0]<0.05&&L.tr.th[0]>-0.2&&L.tr.ph[0]<0.1&&L.tr.ph[0]>-0.1)";
      //track cut on target
      TCut tg_cut = "(L.tr.tg_th[0]<0.055&&L.tr.tg_th[0]>-0.055&&L.tr.tg_ph[0]>-0.018&&L.tr.tg_ph[0]<0.026)";

      TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
      c2->Divide(3,3);
      c2->cd(7);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(1);
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* huq = new TH1F("huq", "P.upQadcL", 300,400,700);
      T->Project(huq->GetName(),"P.upQadcL");
      huq->SetXTitle("P.upQadcL raw");
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
      TCut x_cut = "(L.tr.x[0]+0.9*L.tr.th[0]) > -0.0665";
      TCut y_cut = "";
//      TCut y_cut = "((L.tr.y[0]+0.9*L.tr.ph[0]) < -0.013)&&((L.tr.y[0]+0.9*L.tr.ph[0]) > -0.035)";

//      TCut x_cut = "";
      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;
      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;

      c2->cd(1);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(0);
      gPad->SetGridx();
      TH1F* hy1 = new TH1F("hy1", "Projected y on detector plane", 200, -0.06, 0.01);
      TH1F* hy2 = new TH1F("hy2", "Projected y on detector plane + adc>cut", 200, -0.06, 0.01);
      TH1F* hy3 = new TH1F("hy3", "Projected y on detector plane + adc<cut", 200, -0.06, 0.01);

      T->Project(hy1->GetName(), "L.tr.y[0]+0.9*L.tr.ph[0]", cut + x_cut); 
      T->Project(hy2->GetName(), "L.tr.y[0]+0.9*L.tr.ph[0]", cut_wadc + x_cut);
      T->Project(hy3->GetName(), "L.tr.y[0]+0.9*L.tr.ph[0]", cut_wadc_ped + x_cut+y_cut);
      hy1->SetXTitle("L.tr.y[0]+0.9*L.tr.ph[0] (m)");
      hy1->SetLineColor(kBlack);
      hy2->SetLineColor(kRed);
      hy3->SetLineColor(kBlue);
      hy1->Draw();
      hy2->Draw("same");
      hy3->Draw("same");

      double sum_hy1 = hy1->Integral();
      double sum_hy2 = hy2->Integral();
      double frac_y = sum_hy2 / sum_hy1;

//      TLatex latex;
//      latex.SetTextSize(0.05);
//      latex.SetTextAlign(13);
//      latex.DrawLatex(-.045, hy1->GetMaximum()/4., Form("fract passed_adc_cut/all: %.2f", frac_y));

      TH1F* hy_clone = (TH1F*)hy2->Clone();
      hy_clone->SetTitle("ratio (y + adc>cut)/(y + adc<cut)");
      hy_clone->Divide(hy3);

      c2->cd(6);
      hy_clone->Draw();

      c2->cd(2);
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* hx1 = new TH1F("hx1", "Projected x on detector plane", 300, -0.14, 0.03);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc>cut", 300, -0.14, 0.03);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc<cut", 300, -0.14, 0.03);

      T->Project(hx1->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]", cut); 
      T->Project(hx2->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]", cut_wadc);
      T->Project(hx3->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]", cut_wadc_ped+y_cut);
      hx1->SetXTitle("L.tr.x[0]+0.9*L.tr.th[0] (m)");
      hx1->SetLineColor(kBlack);
      hx2->SetLineColor(kRed);
 //     TF1 *fx=new TF1("fx","gaus",-0.05,-0.02);
 //     fx->SetLineWidth(2);
 //     fx->SetLineColor(kRed);
 //     hx2->Fit("fx","R");
 //     Double_t fxpar[3];
 //     fx->GetParameters(&fxpar[0]);
      hx3->SetLineColor(kBlue);
      hx3->SetLineColor(kBlue);
      hx1->Draw();
      hx2->Draw("same");
      hx3->Draw("same");
 //     TString xfit = Form("Mean = %2.3f",fxpar[1]);
 //     TPaveLabel *xpt = new TPaveLabel(0.20,0.6,0.45,0.65,xfit,"NDC");
 //     xpt->SetBorderSize(0);
 //     xpt->SetTextColor(kRed);
 //     xpt->SetTextSize(0.75);
 //     xpt->SetFillColor(0);
 //     xpt->Draw();
      
      c2->cd(3);
      gPad->SetGridx();
      gPad->SetGridy();

      TH2F* h2d1 = new TH2F("h2d1", "Projected x vs y", 150, -0.12, 0.08, 150, -0.15, 0.05);
      T->Project(h2d1->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]:L.tr.y[0]+0.9*L.tr.ph[0]",cut);
      h2d1->SetXTitle("Projected y (m)");
      h2d1->SetYTitle("Projected x (m)");
      // h2d1->GetZaxis()->SetRangeUser(0, z_scale);
      h2d1->Draw("COLZ");

      c2->cd(4);
      gPad->SetGridx();
      gPad->SetGridy();

      TH2F* h2d2 = new TH2F("h2d2", "Projected x vs y w/ adc cut", 150, -0.12, 0.08, 150, -0.15, 0.05);
      float nAccept = T->Project(h2d2->GetName(),"L.tr.x[0]+0.9*L.tr.th[0]:L.tr.y[0]+0.9*L.tr.ph[0]",cut_wadc);
      h2d2->SetXTitle("Projected y (m)");
      h2d2->SetYTitle("Projected x (m)");
      h2d2->Draw("COLZ");

      c2->cd(5);
      gPad->SetGridx();
      gPad->SetGridy();

      TCut ePeak = "L.tr.y[0]+0.9*L.tr.ph[0]>-0.02 && L.tr.x[0]+0.9*L.tr.th[0]>-0.05";
      float nMiss = T->Draw("L.tr.x[0]+0.9*L.tr.th[0]:L.tr.y[0]+0.9*L.tr.ph[0]",cut_wadc_ped + ePeak);
      cout<<"Accept:"<<nAccept<<"\t"<<"Miss:"<<nMiss<<"\t"<<"Ratio:"<<nMiss/nAccept<<endl;
      TH2F* h2d3 = new TH2F("h2d3", "Projected x vs y adc < cut", 150, -0.12, 0.08, 150, -0.15, 0.05);
      T->Project(h2d3->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]:L.tr.y[0]+0.9*L.tr.ph[0]",cut_wadc_ped);
      h2d3->SetXTitle("Projected y (m)");
      h2d3->SetYTitle("Projected x (m)");
      h2d3->Draw("COLZ");

      c2->cd(8);
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogy(1);
      TH1F* huq_cor = new TH1F("huq_cor", "P.upQadcL w/adccut", 800, 0, 800);
      T->Project(huq_cor->GetName(), Form("(P.upQadcL-%f)",pedMean),adc_cut_up);
      huq_cor->SetXTitle("ADC CH");
      huq_cor->SetYTitle("Events/CH");
      huq_cor->Draw();
      double sig_peak = huq_cor->GetBinCenter(huq_cor->GetMaximumBin());
      TF1 *f1 = new TF1("f1","gaus",sig_peak-5,sig_peak+5);
      f1->SetLineWidth(2);
      f1->SetLineColor(kRed);
      huq_cor->Fit("f1","R"); 
      Double_t f1par[3];
      f1->GetParameters(&f1par[0]);
      TF1 *fun = new TF1("fun",langaufun,f1par[1]-f1par[2]*1.5,f1par[1]+f1par[2]*2,4);
      fun->SetParNames("Lwidth","MPV","Integral","GSigma");
      fun->SetParameters(0.5,f1par[1],huq_cor->GetEntries(), f1par[2]/2.);
      fun->SetLineWidth(1);
      fun->SetLineColor(kRed);
      huq_cor->Fit("fun","R");
      Double_t lgpar[4];
      fun->GetParameters(&lgpar[0]);
      char label_corSigma[64];
      Double_t Sigma_cor = sqrt(lgpar[3]*lgpar[3]-fpedpar[2]*fpedpar[2]);
      Double_t Resolution_cor = Sigma_cor/lgpar[1];
      int pes = int(lgpar[1]*lgpar[1]/Sigma_cor/Sigma_cor);
      sprintf(label_corSigma,"#sigma_{cor}= %3.2f",Sigma_cor);
      char labelPE[64];
      sprintf(labelPE,"PE's (from stat)= %d",pes);
      TPaveLabel *ptPE = new TPaveLabel(0.60,0.29, 0.73,0.34,labelPE,"NDC");
      ptPE->SetBorderSize(0);
      ptPE->SetTextColor(kRed);
      ptPE->SetTextSize(1.0);
      ptPE->SetFillColor(0);
      ptPE->Draw();
      char labelRes[64];
      sprintf(labelRes,"#sigma_{cor}/MPV= %5.3f",Resolution_cor);
      TPaveLabel *ptRes = new TPaveLabel(0.60,0.35, 0.73,0.40,labelRes,"NDC");
      ptRes->SetBorderSize(0);
      ptRes->SetTextColor(kRed);
      ptRes->SetTextSize(1.0);
      ptRes->SetFillColor(0);
      ptRes->Draw();
      TPaveLabel *ptSig = new TPaveLabel(0.60,0.41, 0.73,0.46,label_corSigma,"NDC");
      ptSig->SetBorderSize(0);
      ptSig->SetTextColor(kRed);
      ptSig->SetTextSize(1.0);
      ptSig->SetFillColor(0);
      ptSig->Draw();
      char label2[64];
      sprintf(label2,"GSigma/MPV = %4.2f",lgpar[3]/lgpar[1]);
      TPaveLabel *pt2 = new TPaveLabel(0.60,0.46,0.75,0.52,label2,"NDC");
      pt2->SetBorderSize(0);
      pt2->SetTextColor(kRed);
      pt2->SetTextSize(0.75);
      pt2->SetFillColor(0);
      pt2->Draw();

      c2->cd(9);
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogy(1);
      TH1F* huqpe = new TH1F("huqpe", "PE w/adccut", 200, 0, 400);
      T->Project(huqpe->GetName(), Form("5.0*(P.upQadcL-%f)/%f/1.6",pedMean,USLgain),adc_cut_up);
      huqpe->SetXTitle("PE from gain");
      huqpe->SetYTitle("Events/2PE");
      huqpe->Draw();
      double binmax = huqpe->GetMaximumBin(); 
      double PeakPEs = huqpe->GetXaxis()->GetBinCenter(binmax);
      TString labelPMT = Form("PMT ZK%4.0f, Gain = %1.2fE5",USLPMT,USLgain);
      TPaveLabel *ptPMT = new TPaveLabel(0.60,0.50,0.75,0.55,labelPMT,"NDC");
      ptPMT->SetBorderSize(0);
      ptPMT->SetTextColor(kBlue);
      ptPMT->SetTextSize(0.75);
      ptPMT->SetFillColor(0);
      ptPMT->Draw();
      char labelPPE[64];
      sprintf(labelPPE,"PEs (from gain) = %3.0f", PeakPEs);
      TPaveLabel *ptPPE = new TPaveLabel(0.60,0.45,0.75,0.50,labelPPE,"NDC");
      ptPPE->SetBorderSize(0);
      ptPPE->SetTextColor(kBlue);
      ptPPE->SetTextSize(1.0);
      ptPPE->SetFillColor(0);
      ptPPE->Draw();
      c2->SaveAs(Form("plots/USL_alignment_run%d.pdf",run_num));
    }else{
      // RHRS
      T->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",run_num));
//      TCut trig_cut = "";
      TCut trig_cut = "fEvtHdr.fEvtType==1";
//      TCut trig_cut = "(P.evtypebits&2)==2";
      TCut vdc_cut = "R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1";
      TCut tr_cut = "(R.tr.th[0]<0.05&&R.tr.th[0]>-0.2&&R.tr.ph[0]<0.1&&R.tr.ph[0]>-0.1)";
      TCut tg_cut = "(R.tr.tg_th[0]<0.05&&R.tr.tg_th[0]>-0.055&&R.tr.tg_ph[0]>-0.026&&R.tr.tg_ph[0]<0.016)";


      TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
      c2->Divide(3,3);

      c2->cd(7);
      gPad->SetGridx();
      gPad->SetGridy();
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

      TCut x_cut = "(R.tr.x[0]+0.9*R.tr.th[0]) > -0.072";
//      TCut x_cut = "";
      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;

      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;
      c2->cd(1);
      gPad->SetGridx();
      TH1F* hy1 = new TH1F("hy1", "Projected y on detector plane", 200, -0.04, 0.04);
      TH1F* hy2 = new TH1F("hy2", "Projected y on detector plane + adc> cut", 200, -0.04, 0.04);
      TH1F* hy3 = new TH1F("hy3", "Projected y on detector plane + adc< cut", 200, -0.04, 0.04);

      T->Project(hy1->GetName(), "R.tr.y[0]+0.9*R.tr.ph[0]", cut + x_cut); 
      T->Project(hy2->GetName(), "R.tr.y[0]+0.9*R.tr.ph[0]", cut_wadc + x_cut);
      T->Project(hy3->GetName(), "R.tr.y[0]+0.9*R.tr.ph[0]", cut_wadc_ped + x_cut);
      hy1->SetXTitle("R.tr.y[0]+0.9*R.tr.ph[0] (m)");
      hy1->SetLineColor(kBlack);
      hy2->SetLineColor(kRed);
      hy3->SetLineColor(kBlue);
      hy1->Draw();
      hy2->Draw("same");
      hy3->Draw("same");

      double sum_hy1 = hy1->Integral();
      double sum_hy2 = hy2->Integral();
      double frac_y = sum_hy2 / sum_hy1;

      TLatex latex;
      latex.SetTextSize(0.05);
      latex.SetTextAlign(13);
      latex.DrawLatex(-.03, hy1->GetMaximum()/4., Form("fract passed_adc_cut/all: %.2f", frac_y));

      TH1F* hy_clone = (TH1F*)hy2->Clone();
      hy_clone->SetTitle("ratio (y + adc>cut)/(y + adc<cut)");
      hy_clone->Divide(hy3);

      c2->cd(6);
      hy_clone->Draw();

      c2->cd(2);
      gPad->SetGridx();
      TH1F* hx1 = new TH1F("hx1", "Projected x on detector plane", 300, -0.2, 0.1);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc> cut", 300, -0.2, 0.1);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc< cut", 300, -0.2, 0.1);

      T->Project(hx1->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]", cut); 
      T->Project(hx2->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]", cut_wadc);
      T->Project(hx3->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]", cut_wadc_ped);
      hx1->SetXTitle("R.tr.x[0]+0.9*R.tr.th[0] (m)");
      hx1->SetLineColor(kBlack);
      hx2->SetLineColor(kRed);
      TF1 *fx=new TF1("fx","gaus",-0.050,-0.025);
      fx->SetLineWidth(2);
      fx->SetLineColor(kRed);
      hx2->Fit("fx","R");
      Double_t fxpar[3];
      fx->GetParameters(&fxpar[0]);
      hx3->SetLineColor(kBlue);
      hx1->Draw();
      hx2->Draw("same");
      hx3->Draw("same");
      TString xfit = Form("Mean = %2.3f",fxpar[1]);
      TPaveLabel *xpt = new TPaveLabel(0.20,0.6,0.45,0.65,xfit,"NDC");
      xpt->SetBorderSize(0);
      xpt->SetTextColor(kRed);
      xpt->SetTextSize(0.75);
      xpt->SetFillColor(0);
      xpt->Draw();

      c2->cd(3);
      gPad->SetGridx();
      gPad->SetGridy();

      TH2F* h2d1 = new TH2F("h2d1", "Projected x vs y", 150, -0.05, 0.05, 150, -0.10, 0.01);
      T->Project(h2d1->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]:R.tr.y[0]+0.9*R.tr.ph[0]",cut);
      h2d1->SetXTitle("Projected y (m)");
      h2d1->SetYTitle("Projected x (m)");
      h2d1->Draw("COLZ");

      c2->cd(4);
      gPad->SetGridx();
      gPad->SetGridy();

      TH2F* h2d2 = new TH2F("h2d2", "Projected x vs y w/ adc > cut", 150, -0.05, 0.05, 150, -0.10, 0.01);
      float nAccept = T->Project(h2d2->GetName(),"R.tr.x[0]+0.9*R.tr.th[0]:R.tr.y[0]+0.9*R.tr.ph[0]",cut_wadc);
      h2d2->SetXTitle("Projected y (m)");
      h2d2->SetYTitle("Projected x (m)");
      h2d2->Draw("COLZ");

      c2->cd(5);
      gPad->SetGridx();
      gPad->SetGridy();
      TCut ePeak = "R.tr.y[0]+0.9*R.tr.ph[0]<-0.01 && R.tr.x[0]+0.9*R.tr.th[0]>-0.05";
      float nMiss = T->Draw("R.tr.x[0]+0.9*R.tr.th[0]:R.tr.y[0]+0.9*R.tr.ph[0]",cut_wadc_ped + ePeak);
      cout<<"Accept:"<<nAccept<<"\t"<<"Miss:"<<nMiss<<"\t"<<"Ratio:"<<nMiss/nAccept<<endl;

      TH2F* h2d3 = new TH2F("h2d3", "Projected x vs y adc < cut", 150, -0.05, 0.05, 150, -0.10, 0.01);
      T->Project(h2d3->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]:R.tr.y[0]+0.9*R.tr.ph[0]",cut_wadc_ped);
      h2d3->SetXTitle("Projected y (m)");
      h2d3->SetYTitle("Projected x (m)");
      h2d3->Draw("COLZ");

      c2->cd(8);
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogy(1);
      TH1F* huq_cor = new TH1F("huq_cor", "P.upQadcR w/adccut", 400, 0, 400);
      T->Project(huq_cor->GetName(), Form("(P.upQadcR-%f)",pedMean),adc_cut_up);
      huq_cor->SetXTitle("ADC CH");
      huq_cor->SetYTitle("Events/CH");
      huq_cor->Draw();
      double sig_peak = huq_cor->GetBinCenter(huq_cor->GetMaximumBin());
      TF1 *f1 = new TF1("f1","gaus",sig_peak-5,sig_peak+5);
      f1->SetLineWidth(2);
      f1->SetLineColor(kRed);
      huq_cor->Fit("f1","R"); 
      Double_t f1par[3];
      f1->GetParameters(&f1par[0]);
      TF1 *fun = new TF1("fun",langaufun,f1par[1]-f1par[2]*1.5,f1par[1]+f1par[2]*2,4);
      fun->SetParNames("Lwidth","MPV","Integral","GSigma");
      fun->SetParameters(0.5,f1par[1],huq_cor->GetEntries(), f1par[2]/2.);
      fun->SetLineWidth(1);
      fun->SetLineColor(kRed);
      huq_cor->Fit("fun","R");
      Double_t lgpar[4];
      fun->GetParameters(&lgpar[0]);
      char label_corSigma[64];
      Double_t Sigma_cor = sqrt(lgpar[3]*lgpar[3]-fpedpar[2]*fpedpar[2]);
      Double_t Resolution_cor = Sigma_cor/lgpar[1];
      int pes = int(lgpar[1]*lgpar[1]/Sigma_cor/Sigma_cor);
      sprintf(label_corSigma,"#sigma_{cor}= %3.2f",Sigma_cor);
      char labelPE[64];
      sprintf(labelPE,"PE's (from stat)= %d",pes);
      TPaveLabel *ptPE = new TPaveLabel(0.60,0.29, 0.73,0.34,labelPE,"NDC");
      ptPE->SetBorderSize(0);
      ptPE->SetTextColor(kRed);
      ptPE->SetTextSize(1.0);
      ptPE->SetFillColor(0);
      ptPE->Draw();
      char labelRes[64];
      sprintf(labelRes,"#sigma_{cor}/MPV= %5.3f",Resolution_cor);
      TPaveLabel *ptRes = new TPaveLabel(0.60,0.35, 0.73,0.40,labelRes,"NDC");
      ptRes->SetBorderSize(0);
      ptRes->SetTextColor(kRed);
      ptRes->SetTextSize(1.0);
      ptRes->SetFillColor(0);
      ptRes->Draw();
      TPaveLabel *ptSig = new TPaveLabel(0.60,0.41, 0.73,0.46,label_corSigma,"NDC");
      ptSig->SetBorderSize(0);
      ptSig->SetTextColor(kRed);
      ptSig->SetTextSize(1.0);
      ptSig->SetFillColor(0);
      ptSig->Draw();
      char label2[64];
      sprintf(label2,"GSigma/MPV = %4.2f",lgpar[3]/lgpar[1]);
      TPaveLabel *pt2 = new TPaveLabel(0.60,0.46,0.75,0.52,label2,"NDC");
      pt2->SetBorderSize(0);
      pt2->SetTextColor(kRed);
      pt2->SetTextSize(0.75);
      pt2->SetFillColor(0);
      pt2->Draw();

      c2->cd(9);
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogy(1);
      TH1F* huqpe = new TH1F("huqpe", "PE w/adccut", 200, 0, 400);
      T->Project(huqpe->GetName(), Form("5.0*(P.upQadcR-%f)/%f/1.6",pedMean,USRgain),adc_cut_up);
      huqpe->SetXTitle("PE from gain");
      huqpe->SetYTitle("Events/2PE");
      huqpe->Draw();
      double binmax = huqpe->GetMaximumBin(); 
      double PeakPEs = huqpe->GetXaxis()->GetBinCenter(binmax);
      TString labelPMT = Form("PMT ZK%4.0f, Gain = %1.2fE5",USRPMT,USRgain);
      TPaveLabel *ptPMT = new TPaveLabel(0.60,0.50,0.75,0.55,labelPMT,"NDC");
      ptPMT->SetBorderSize(0);
      ptPMT->SetTextColor(kBlue);
      ptPMT->SetTextSize(0.75);
      ptPMT->SetFillColor(0);
      ptPMT->Draw();
      char labelPPE[64];
      sprintf(labelPPE,"PEs (from gain) = %3.0f", PeakPEs);
      TPaveLabel *ptPPE = new TPaveLabel(0.60,0.45,0.75,0.50,labelPPE,"NDC");
      ptPPE->SetBorderSize(0);
      ptPPE->SetTextColor(kBlue);
      ptPPE->SetTextSize(1.0);
      ptPPE->SetFillColor(0);
      ptPPE->Draw();
      c2->SaveAs(Form("plots/USR_alignment_run%d.pdf",run_num));
    }

}
