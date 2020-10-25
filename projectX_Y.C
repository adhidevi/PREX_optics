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
void projectX_Y(int run_num){ 

  gROOT->Reset();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
   
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

      TCanvas *c2 = new TCanvas("c2","c2",800,600);
      c2->Divide(2,2);
      c2->cd(1);
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
//      TCut x_cut = "";
      TCut x_cut90 = "(L.tr.x[0]+0.9*L.tr.th[0]) > -0.0665";
      TCut x_cut140 = "(L.tr.x[0]+1.4*L.tr.th[0]) > -0.0725";
      TCut y_cut = "";
//      TCut y_cut = "((L.tr.y[0]+0.9*L.tr.ph[0]) < -0.013)&&((L.tr.y[0]+0.9*L.tr.ph[0]) > -0.035)";
      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;
      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;

      c2->cd(2);
      gPad->SetGridx();
      TH1F* hx1 = new TH1F("hx1", "Projected x on detector plane", 300, -0.14, 0.03);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc>cut", 300, -0.14, 0.03);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc<cut", 300, -0.14, 0.03);
      TH1F* hxN1 = new TH1F("hxN1", "Projected x on detector plane", 300, -0.14, 0.03);
      TH1F* hxN2 = new TH1F("hxN2", "Projected x on detector plane + adc>cut", 300, -0.14, 0.03);
      TH1F* hxN3 = new TH1F("hxN3", "Projected x on detector plane + adc<cut", 300, -0.14, 0.03);

      T->Project(hx1->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]", cut); 
      T->Project(hx2->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]", cut_wadc);
      T->Project(hx3->GetName(), "L.tr.x[0]+0.9*L.tr.th[0]", cut_wadc_ped+y_cut);
      T->Project(hxN1->GetName(), "L.tr.x[0]+1.4*L.tr.th[0]", cut); 
      T->Project(hxN2->GetName(), "L.tr.x[0]+1.4*L.tr.th[0]", cut_wadc);
      T->Project(hxN3->GetName(), "L.tr.x[0]+1.4*L.tr.th[0]", cut_wadc_ped+y_cut);

      hx1->SetXTitle("Transport x projected to quartz (m)");
      hx1->SetLineColor(1);
      hx2->SetLineColor(2);
      hx3->SetLineColor(4);
      hxN1->SetLineColor(6);
      hxN2->SetLineColor(7);
      hxN3->SetLineColor(8);

      hx1->Draw();
      hx2->Draw("same");
      hx3->Draw("same");
      hxN1->Draw("same");
      hxN2->Draw("same");
      hxN3->Draw("same");
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(1);
      latex.DrawLatex(0.15,0.85, Form("All @ 90 cm"));
      latex.SetTextColor(2);
      latex.DrawLatex(0.15,0.80, Form("Accepted @ 90 cm"));
      latex.SetTextColor(4);
      latex.DrawLatex(0.15,0.75, Form("Missed @ 90 cm"));
      latex.SetTextColor(6);
      latex.DrawLatex(0.15,0.70, Form("All @ 140 cm"));
      latex.SetTextColor(7);
      latex.DrawLatex(0.15,0.65, Form("Accepted @ 140 cm"));
      latex.SetTextColor(8);
      latex.DrawLatex(0.15,0.60, Form("Missed @ 140 cm"));

      
      c2->cd(3);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(0);
      gPad->SetGridx();
      TH1F* hy1 = new TH1F("hy1", "Projected y on detector plane", 200, -0.06, 0.01);
      TH1F* hy2 = new TH1F("hy2", "Projected y on detector plane + adc>cut", 200, -0.06, 0.01);
      TH1F* hy3 = new TH1F("hy3", "Projected y on detector plane + adc<cut", 200, -0.06, 0.01);
      TH1F* hyN1 = new TH1F("hyN1", "Projected y on detector plane", 200, -0.06, 0.01);
      TH1F* hyN2 = new TH1F("hyN2", "Projected y on detector plane + adc>cut", 200, -0.06, 0.01);
      TH1F* hyN3 = new TH1F("hyN3", "Projected y on detector plane + adc<cut", 200, -0.06, 0.01);

      T->Project(hy1->GetName(), "L.tr.y[0]+0.9*L.tr.ph[0]", cut + x_cut90); 
      T->Project(hy2->GetName(), "L.tr.y[0]+0.9*L.tr.ph[0]", cut_wadc + x_cut90);
      T->Project(hy3->GetName(), "L.tr.y[0]+0.9*L.tr.ph[0]", cut_wadc_ped + x_cut90+y_cut);
      T->Project(hyN1->GetName(), "L.tr.y[0]+1.4*L.tr.ph[0]", cut + x_cut140); 
      T->Project(hyN2->GetName(), "L.tr.y[0]+1.4*L.tr.ph[0]", cut_wadc + x_cut140);
      T->Project(hyN3->GetName(), "L.tr.y[0]+1.4*L.tr.ph[0]", cut_wadc_ped + x_cut140+y_cut);
      hy1->SetXTitle("Transport y projected to quartz (m)");
      hy1->SetLineColor(1);
      hy2->SetLineColor(2);
      hy3->SetLineColor(4);
      hyN1->SetLineColor(6);
      hyN2->SetLineColor(7);
      hyN3->SetLineColor(8);
      hy1->Draw();
      hy2->Draw("same");
      hy3->Draw("same");
      hyN1->Draw("same");
      hyN2->Draw("same");
      hyN3->Draw("same");

      c2->SaveAs(Form("plots/projectedL_x_y_run%d.pdf",run_num));
    }else{
      // RHRS
      T->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",run_num));
//      TCut trig_cut = "";
      TCut trig_cut = "fEvtHdr.fEvtType==1";
//      TCut trig_cut = "(P.evtypebits&2)==2";
      TCut vdc_cut = "R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1";
      TCut tr_cut = "(R.tr.th[0]<0.05&&R.tr.th[0]>-0.2&&R.tr.ph[0]<0.1&&R.tr.ph[0]>-0.1)";
      TCut tg_cut = "(R.tr.tg_th[0]<0.05&&R.tr.tg_th[0]>-0.055&&R.tr.tg_ph[0]>-0.026&&R.tr.tg_ph[0]<0.016)";


      TCanvas *c2 = new TCanvas("c2","c2",800,600);
      c2->Divide(2,2);
      c2->cd(1);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(1);
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* huq = new TH1F("huq", "P.upQadcR", 300,400,700);
      T->Project(huq->GetName(),"P.upQadcR");
      huq->SetXTitle("P.upQadcR raw");
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
      //cut on radiative tail
//      TCut x_cut = "";
      TCut x_cut90 = "(R.tr.x[0]+0.9*R.tr.th[0]) > -0.066";
      TCut x_cut140 = "(R.tr.x[0]+1.4*R.tr.th[0]) > -0.077";
      TCut y_cut = "";
//      TCut y_cut = "((R.tr.y[0]+0.9*R.tr.ph[0]) < -0.013)&&((R.tr.y[0]+0.9*R.tr.ph[0]) > -0.035)";
      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;
      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;

      c2->cd(2);
      gPad->SetGridx();
      TH1F* hx1 = new TH1F("hx1", "Projected x on detector plane", 300, -0.14, 0.03);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc>cut", 300, -0.14, 0.03);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc<cut", 300, -0.14, 0.03);
      TH1F* hxN1 = new TH1F("hxN1", "Projected x on detector plane", 300, -0.14, 0.03);
      TH1F* hxN2 = new TH1F("hxN2", "Projected x on detector plane + adc>cut", 300, -0.14, 0.03);
      TH1F* hxN3 = new TH1F("hxN3", "Projected x on detector plane + adc<cut", 300, -0.14, 0.03);

      T->Project(hx1->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]", cut); 
      T->Project(hx2->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]", cut_wadc);
      T->Project(hx3->GetName(), "R.tr.x[0]+0.9*R.tr.th[0]", cut_wadc_ped+y_cut);
      T->Project(hxN1->GetName(), "R.tr.x[0]+1.4*R.tr.th[0]", cut); 
      T->Project(hxN2->GetName(), "R.tr.x[0]+1.4*R.tr.th[0]", cut_wadc);
      T->Project(hxN3->GetName(), "R.tr.x[0]+1.4*R.tr.th[0]", cut_wadc_ped+y_cut);

      hx1->SetXTitle("Transport x projected to quartz (m)");
      hx1->SetLineColor(1);
      hx2->SetLineColor(2);
      hx3->SetLineColor(4);
      hxN1->SetLineColor(6);
      hxN2->SetLineColor(7);
      hxN3->SetLineColor(8);

      hx1->Draw();
      hx2->Draw("same");
      hx3->Draw("same");
      hxN1->Draw("same");
      hxN2->Draw("same");
      hxN3->Draw("same");
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(1);
      latex.DrawLatex(0.15,0.85, Form("All @ 90 cm"));
      latex.SetTextColor(2);
      latex.DrawLatex(0.15,0.80, Form("Accepted @ 90 cm"));
      latex.SetTextColor(4);
      latex.DrawLatex(0.15,0.75, Form("Missed @ 90 cm"));
      latex.SetTextColor(6);
      latex.DrawLatex(0.15,0.70, Form("All @ 140 cm"));
      latex.SetTextColor(7);
      latex.DrawLatex(0.15,0.65, Form("Accepted @ 140 cm"));
      latex.SetTextColor(8);
      latex.DrawLatex(0.15,0.60, Form("Missed @ 140 cm"));

      
      c2->cd(3);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(0);
      gPad->SetGridx();
      TH1F* hy1 = new TH1F("hy1", "Projected y on detector plane", 200, -0.04, 0.04);
      TH1F* hy2 = new TH1F("hy2", "Projected y on detector plane + adc>cut", 200, -0.04, 0.04);
      TH1F* hy3 = new TH1F("hy3", "Projected y on detector plane + adc<cut", 200, -0.04, 0.04);
      TH1F* hyN1 = new TH1F("hyN1", "Projected y on detector plane", 200, -0.04, 0.04);
      TH1F* hyN2 = new TH1F("hyN2", "Projected y on detector plane + adc>cut", 200, -0.04, 0.04);
      TH1F* hyN3 = new TH1F("hyN3", "Projected y on detector plane + adc<cut", 200, -0.04, 0.04);

      T->Project(hy1->GetName(), "R.tr.y[0]+0.9*R.tr.ph[0]", cut + x_cut90); 
      T->Project(hy2->GetName(), "R.tr.y[0]+0.9*R.tr.ph[0]", cut_wadc + x_cut90);
      T->Project(hy3->GetName(), "R.tr.y[0]+0.9*R.tr.ph[0]", cut_wadc_ped + x_cut90+y_cut);
      T->Project(hyN1->GetName(), "R.tr.y[0]+1.4*R.tr.ph[0]", cut + x_cut140); 
      T->Project(hyN2->GetName(), "R.tr.y[0]+1.4*R.tr.ph[0]", cut_wadc + x_cut140);
      T->Project(hyN3->GetName(), "R.tr.y[0]+1.4*R.tr.ph[0]", cut_wadc_ped + x_cut140+y_cut);
      hy1->SetXTitle("Transport y projected to quartz (m)");
      hy1->SetLineColor(1);
      hy2->SetLineColor(2);
      hy3->SetLineColor(4);
      hyN1->SetLineColor(6);
      hyN2->SetLineColor(7);
      hyN3->SetLineColor(8);
      hyN1->Draw();
      hy2->Draw("same");
      hy3->Draw("same");
      hy1->Draw("same");
      hyN2->Draw("same");
      hyN3->Draw("same");

      c2->SaveAs(Form("plots/projectedR_x_y_run%d.pdf",run_num));
    }

}
