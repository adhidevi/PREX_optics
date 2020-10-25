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
void project_different_z(int run, double detZ){ 
//change the x_cut on rad tail based on your prefered Z location
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
   
  // Cut for FBUS
  double upadc_cutL_approx = 481;
  double loadc_cutL_approx = 600;
  double upadc_cutR_approx = 500;
  double USLPMT = 5370;
  double USLgain = 7.5;//7.5e5
  double USRPMT = 5401;
  double USRgain = 6.46;//6.46e5
      TChain *T = new TChain("T");

  if(run < 10000){
      //LHRS
       T->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run));

      TCanvas *c2 = new TCanvas("c2","c2",900,700);
      c2->Divide(2,2);
      c2->cd(1);
      gPad->Modified();gPad->Update();
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* huq = new TH1F("huq", Form("ADC raw (run%d);ADC raw;Events/CH",run), 400,400,800);
      huq->SetLineColor(1);
      T->Project(huq->GetName(),"P.upQadcL");
      huq->Draw();

      TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
      huqCpy->GetXaxis()->SetRangeUser(upadc_cutL_approx-10,upadc_cutL_approx+10);
      double upadc_cutL = huqCpy->GetBinCenter(huqCpy->GetMinimumBin());

      TLine* ped_line = new TLine(upadc_cutL,0.0,upadc_cutL,huq->GetMaximum());
      ped_line->SetLineColor(9);
      ped_line->SetLineWidth(2);
      ped_line->Draw();

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(1);
      latex.DrawLatex(0.70,0.85,"P.upQadcL");
      latex.SetTextColor(9);
      latex.DrawLatex(0.70,0.75,"Ped. cut");

      TCut cut = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)");
      TCut cutADCup = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)&&P.upQadcL> %f",upadc_cutL);
      TCut cutPEDup = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)&&P.upQadcL< %f",upadc_cutL);
      double x_edge;
      if(detZ==90)
             x_edge = -0.0665;// for z=90cm
      else if(detZ==100)
             x_edge = -0.0680;// for z=100cm
      else if(detZ==110)
             x_edge = -0.0695;// for z=110cm
      else if(detZ==120)
             x_edge = -0.0705;// for z=120cm
      else if(detZ==130)
             x_edge = -0.0715;// for z=130cm
      else if(detZ==140)
             x_edge = -0.0725;// for z=140cm
      else if(detZ==150)
             x_edge = -0.0735;// for z=150cm
      
      TCut x_cut = Form("(L.tr.x[0]+%f*L.tr.th[0]) > %f",detZ/100.0,x_edge);

      c2->cd(2);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(0);
      gPad->SetGridx();
      TH1F* hy1 = new TH1F("hy1", "Projected y on detector plane", 200, -0.06, 0.01);
      TH1F* hyup2 = new TH1F("hyup2", "Projected y on detector plane + adc>cut", 200, -0.06, 0.01);
      TH1F* hyup3 = new TH1F("hyup3", "Projected y on detector plane + adc<cut", 200, -0.06, 0.01);

      T->Project(hy1->GetName(), Form("L.tr.y[0]+%f*L.tr.ph[0]",detZ/100.0), cut+x_cut); 
      T->Project(hyup2->GetName(), Form("L.tr.y[0]+%f*L.tr.ph[0]",detZ/100.0), cutADCup+x_cut);
      T->Project(hyup3->GetName(), Form("L.tr.y[0]+%f*L.tr.ph[0]",detZ/100.0), cutPEDup+x_cut);
      hy1->SetXTitle(Form("L.tr.y[0]+%1.1f*L.tr.ph[0] (m)",detZ/100.0));
      hy1->SetLineColor(1);
      hyup2->SetLineColor(2);
      hyup3->SetLineColor(4);
      hy1->Draw();
      hyup2->Draw("same");
      hyup3->Draw("same");
     
      latex.SetTextColor(1);
      latex.DrawLatex(0.20,0.85,"All events");
      latex.SetTextColor(2);
      latex.DrawLatex(0.20,0.80,"US accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.20,0.75,"US missed");

      c2->cd(3);
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* hx1 = new TH1F("hx1", "Projected x on detector plane", 300, -0.14, 0.03);
      TH1F* hxup2 = new TH1F("hxup2", "Projected x on detector plane + adc>cut", 300, -0.14, 0.03);
      TH1F* hxup3 = new TH1F("hxup3", "Projected x on detector plane + adc<cut", 300, -0.14, 0.03);

      T->Project(hx1->GetName(), Form("L.tr.x[0]+%f*L.tr.th[0]",detZ/100.0), cut); 
      T->Project(hxup2->GetName(), Form("L.tr.x[0]+%f*L.tr.th[0]",detZ/100.0), cutADCup);
      T->Project(hxup3->GetName(), Form("L.tr.x[0]+%f*L.tr.th[0]",detZ/100.0), cutPEDup);
      hx1->SetXTitle(Form("L.tr.x[0]+%1.1f*L.tr.th[0] (m)",detZ/100.0));
      hx1->SetLineColor(1);
      hxup2->SetLineColor(2);
      hxup3->SetLineColor(4);
      hx1->Draw();
      hxup2->Draw("same");
      hxup3->Draw("same");
      TF1* fit_x = new TF1("fit_x","gaus",hxup2->GetBinCenter(hxup2->GetMaximumBin())-0.01,hxup2->GetBinCenter(hxup2->GetMaximumBin())+0.01);
      hxup2->Fit(fit_x,"R");
      
      latex.SetTextColor(1);
      latex.DrawLatex(0.20,0.85,"All events");
      latex.SetTextColor(2);
      latex.DrawLatex(0.20,0.80,"US accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.20,0.75,"US missed");
      latex.SetTextColor(6);
      latex.DrawLatex(0.20,0.70,Form("Peak = %1.4f",fit_x->GetParameter(1)));
      latex.DrawLatex(0.20,0.65,Form("Edge = %1.4f",x_edge));

      c2->cd(4);
      gPad->SetGridx();
      gPad->SetGridy();

      TH2F* h2d2 = new TH2F("h2d2", "Projected x vs y w/ up_adc cut", 150, -0.12, 0.08, 150, -0.15, 0.05);
      float nAccept = T->Project(h2d2->GetName(),Form("L.tr.x[0]+%f*L.tr.th[0]:L.tr.y[0]+%f*L.tr.ph[0]",detZ/100.0,detZ/100.0),cutADCup);
      h2d2->SetXTitle("Projected y (m)");
      h2d2->SetYTitle("Projected x (m)");
      h2d2->Draw("COLZ");

      c2->SaveAs(Form("plots/quartz_z_%3.0fcm_run%d.pdf",detZ,run));
    }else{
      // RHRS
      T->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",run));
//      TCut trig_cut = "";
      TCut trig_cut = "fEvtHdr.fEvtType==1";
//      TCut trig_cut = "(P.evtypebits&2)==2";
      TCut vdc_cut = "R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1";
      TCut tr_cut = "(R.tr.th[0]<0.05&&R.tr.th[0]>-0.2&&R.tr.ph[0]<0.1&&R.tr.ph[0]>-0.1)";
      TCut tg_cut = "(R.tr.tg_th[0]<0.05&&R.tr.tg_th[0]>-0.055&&R.tr.tg_ph[0]>-0.026&&R.tr.tg_ph[0]<0.016)";

      TCanvas *c2 = new TCanvas("c2","c2",900,700);
      c2->Divide(2,2);

      c2->cd(1);
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* huq = new TH1F("huq", Form("ADC raw (run%d);P.upQadcR;Events/CH",run), 300, 400, 700);
      T->Project(huq->GetName(),"P.upQadcR");
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

//      TCut x_cut = "";
      TCut x_cut = Form("(R.tr.x[0]+%f*R.tr.th[0]) > %f",detZ/100.0,-0.072);
      //x_cut = -0.072 for z=90cm
      //x_cut = -0.073 for z=100cm
      //x_cut = -0.074 for z=110cm
      //x_cut = -0.075 for z=120cm
      //x_cut = -0.076 for z=130cm
      //x_cut = -0.077 for z=140cm
      //x_cut = -0.078 for z=150cm

      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;

      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;
      c2->cd(2);
      gPad->SetGridx();
      TH1F* hy1 = new TH1F("hy1", "Projected y on detector plane", 200, -0.04, 0.04);
      TH1F* hy2 = new TH1F("hy2", "Projected y on detector plane + adc> cut", 200, -0.04, 0.04);
      TH1F* hy3 = new TH1F("hy3", "Projected y on detector plane + adc< cut", 200, -0.04, 0.04);

      T->Project(hy1->GetName(), Form("R.tr.y[0]+%f*R.tr.ph[0]",detZ/100.0), cut + x_cut); 
      T->Project(hy2->GetName(), Form("R.tr.y[0]+%f*R.tr.ph[0]",detZ/100.0), cut_wadc + x_cut);
      T->Project(hy3->GetName(), Form("R.tr.y[0]+%f*R.tr.ph[0]",detZ/100.0), cut_wadc_ped + x_cut);
      hy1->SetXTitle(Form("R.tr.y[0]+%1.1f*R.tr.ph[0] (m)",detZ/100.0));
      hy1->SetLineColor(kBlack);
      hy2->SetLineColor(kRed);
      hy3->SetLineColor(kBlue);
      hy1->Draw();
      hy2->Draw("same");
      hy3->Draw("same");

      c2->cd(3);
      gPad->SetGridx();
      TH1F* hx1 = new TH1F("hx1", "Projected x on detector plane", 300, -0.2, 0.1);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc> cut", 300, -0.2, 0.1);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc< cut", 300, -0.2, 0.1);

      T->Project(hx1->GetName(), Form("R.tr.x[0]+%f*R.tr.th[0]",detZ/100.0), cut); 
      T->Project(hx2->GetName(), Form("R.tr.x[0]+%f*R.tr.th[0]",detZ/100.0), cut_wadc);
      T->Project(hx3->GetName(), Form("R.tr.x[0]+%f*R.tr.th[0]",detZ/100.0), cut_wadc_ped);
      hx1->SetXTitle(Form("R.tr.x[0]+%1.1f*R.tr.th[0] (m)",detZ/100.0));
      hx1->SetLineColor(kBlack);
      hx2->SetLineColor(kRed);
      TF1 *fx=new TF1("fx","gaus",hx2->GetBinCenter(hx2->GetMaximumBin())-0.01,hx2->GetBinCenter(hx2->GetMaximumBin())+0.01);
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

      c2->cd(4);
      gPad->SetGridx();
      gPad->SetGridy();

      TH2F* h2d2 = new TH2F("h2d2", "Projected x vs y w/ adc > cut", 150, -0.05, 0.05, 150, -0.10, 0.01);
      float nAccept = T->Project(h2d2->GetName(),Form("R.tr.x[0]+%f*R.tr.th[0]:R.tr.y[0]+%f*R.tr.ph[0]",detZ/100.0,detZ/100.0),cut_wadc);
      h2d2->SetXTitle("Projected y (m)");
      h2d2->SetYTitle("Projected x (m)");
      h2d2->Draw("COLZ");

//      TLine* y1 = new TLine(-0.016,-0.1,-0.016,0.01);
//      TLine* y2 = new TLine(-0.016+0.035,-0.1,-0.016+0.035,0.01);
//      y1->SetLineColor(2);
//      y2->SetLineColor(2);
//      y1->Draw();
//      y2->Draw();

      c2->SaveAs(Form("plots/quartz_z_%3.0fcm_run%d.pdf",detZ,run));
    }

}
