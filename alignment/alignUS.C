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
void alignUS(int run){ 

  gROOT->Reset();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
   
  // Cut for FBUS
  double upADCcutL_approx = 480;
  double upADCcutR_approx = 502;
  double USLPMT = 5370;
  double USLgain = 7.5;//7.5e5
  double USRPMT = 5401;
  double USRgain = 6.46;//6.46e5
  // RHRS = R.*
  // LHRS = L.*
      TChain *T = new TChain("T");

  if(run < 10000){
      //LHRS
       T->Add(Form("/chafs1/work1/prex_counting/chandan/prexLHRS_%d_-1*.root",run));


      TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
      c2->Divide(3,3);
      c2->cd(7);
      TH1F* huq = new TH1F("huq","P.upQadcL;ADC CH;Event/CH",250,400,650);
      T->Draw("P.upQadcL>>huq");

      TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
      huqCpy->GetXaxis()->SetRangeUser(upADCcutL_approx-10,upADCcutL_approx+10);
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
      TPaveLabel *pt = new TPaveLabel(0.60,0.48,0.75,0.54,labelped,"NDC");
      pt->SetBorderSize(0);
      pt->SetTextColor(kRed);
      pt->SetTextSize(0.75);
      pt->SetFillColor(0);
      pt->Draw();
      char labelcut[64];
      sprintf(labelcut,"ADC_Cut = %3.0f",upADCcut);
      TPaveLabel *ptcut = new TPaveLabel(0.60,0.42,0.75,0.48,labelcut,"NDC");
      ptcut->SetBorderSize(0);
      ptcut->SetTextColor(kMagenta);
      ptcut->SetTextSize(0.75);
      ptcut->SetFillColor(0);
      ptcut->Draw();
      gPad->Update();
      double pedMean = fpedpar[1];

      TCut cut = Form("(P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04");
      TCut cut_wadc = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&P.upQadcL>%f)",upADCcut);
      TCut cut_wped = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&P.upQadcL<%f)",upADCcut);
 
     double q_edge = -0.0715;
     TCut x_cut = Form("(L.tr.x[0]+1.3*L.tr.th[0]) > %f",q_edge);
     TCut y_cut = "";

      c2->cd(1);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(0);
      gPad->SetGridx();
      TH1F* hy1 = new TH1F("hy1", "Projected y on detector plane", 200, -0.06, 0.01);
      TH1F* hy2 = new TH1F("hy2", "Projected y on detector plane + adc>cut", 200, -0.06, 0.01);
      TH1F* hy3 = new TH1F("hy3", "Projected y on detector plane + adc<cut", 200, -0.06, 0.01);

      T->Project(hy1->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cut + x_cut); 
      T->Project(hy2->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cut_wadc + x_cut);
      T->Project(hy3->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cut_wped + x_cut+y_cut);
      hy1->SetXTitle("L.tr.y[0]+1.3*L.tr.ph[0] (m)");
      hy1->SetLineColor(kBlack);
      hy2->SetLineColor(kRed);
      hy3->SetLineColor(kBlue);
      hy1->Draw();
      hy2->Draw("same");
      hy3->Draw("same");

      double sum_hy1 = hy1->Integral();
      double sum_hy2 = hy2->Integral();
      double frac_y = sum_hy2 / sum_hy1;

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

      T->Project(hx1->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cut); 
      T->Project(hx2->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cut_wadc);
      T->Project(hx3->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cut_wped+y_cut);
      hx1->SetXTitle("L.tr.x[0]+1.3*L.tr.th[0] (m)");
      hx1->SetLineColor(kBlack);
      hx2->SetLineColor(kRed);
      TF1 *fx=new TF1("fx","gaus",hx2->GetBinCenter(hx2->GetMaximumBin())-0.012,hx2->GetBinCenter(hx2->GetMaximumBin())+0.012);
      fx->SetLineWidth(2);
      fx->SetLineColor(kRed);
      hx2->Fit("fx","R");
      Double_t fxpar[3];
      fx->GetParameters(&fxpar[0]);
      hx3->SetLineColor(kBlue);
      hx3->SetLineColor(kBlue);
      hx1->Draw();
      hx2->Draw("same");
      hx3->Draw("same");

      TLine* qedgeLine = new TLine(q_edge,0.0,q_edge,hx1->GetMaximum()/2.0);
      qedgeLine->SetLineColor(6);
      qedgeLine->SetLineWidth(2);
      qedgeLine->Draw();

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(2);
      latex.DrawLatex(0.2,0.8,"Elastic Peak");
      latex.SetTextColor(6);
      latex.DrawLatex(0.2,0.75,"q_edge");

      c2->cd(3);
      gPad->SetGridx();
      gPad->SetGridy();

      TH2F* h2d1 = new TH2F("h2d1", "Projected x vs y", 150, -0.12, 0.08, 150, -0.15, 0.05);
      T->Project(h2d1->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]:L.tr.y[0]+1.3*L.tr.ph[0]",cut);
      h2d1->SetXTitle("Projected y (m)");
      h2d1->SetYTitle("Projected x (m)");
      h2d1->Draw("COLZ");

      c2->cd(4);
      gPad->SetGridx();
      gPad->SetGridy();

      TH2F* h2d2 = new TH2F("h2d2", "Projected x vs y w/ adc cut", 150, -0.12, 0.08, 150, -0.15, 0.05);
      float nAccept = T->Project(h2d2->GetName(),"L.tr.x[0]+1.3*L.tr.th[0]:L.tr.y[0]+1.3*L.tr.ph[0]",cut_wadc);
      h2d2->SetXTitle("Projected y (m)");
      h2d2->SetYTitle("Projected x (m)");
      h2d2->Draw("COLZ");

      c2->cd(5);
      gPad->SetGridx();
      gPad->SetGridy();

      TCut ePeak = "L.tr.y[0]+1.3*L.tr.ph[0]>-0.02 && L.tr.x[0]+1.3*L.tr.th[0]>-0.05";
      float nMiss = T->Draw("L.tr.x[0]+1.3*L.tr.th[0]:L.tr.y[0]+1.3*L.tr.ph[0]",cut_wped + ePeak);
      cout<<"Accept:"<<nAccept<<"\t"<<"Miss:"<<nMiss<<"\t"<<"Ratio:"<<nMiss/nAccept<<endl;
      TH2F* h2d3 = new TH2F("h2d3", "Projected x vs y adc < cut", 150, -0.12, 0.08, 150, -0.15, 0.05);
      T->Project(h2d3->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]:L.tr.y[0]+1.3*L.tr.ph[0]",cut_wped);
      h2d3->SetXTitle("Projected y (m)");
      h2d3->SetYTitle("Projected x (m)");
      h2d3->Draw("COLZ");

      c2->cd(8);
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogy(1);
      TH1F* huq_cor = new TH1F("huq_cor", "P.upQadcL w/adccut", 800, 0, 800);
      T->Project(huq_cor->GetName(), Form("(P.upQadcL-%f)",pedMean),Form("P.upQadcL>%f",upADCcut));
      huq_cor->SetXTitle("ADC CH");
      huq_cor->SetYTitle("Events/CH");
      huq_cor->Draw();
      TF1 *f1 = new TF1("f1","gaus",40,90);
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
      T->Project(huqpe->GetName(), Form("5.0*(P.upQadcL-%f)/%f/1.6",pedMean,USLgain),Form("P.upQadcL>%f",upADCcut));
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
      c2->SaveAs(Form("plots/USL_alignment_run%d.pdf",run));
    }else{
      // RHRS
      T->Add(Form("/chafs1/work1/prex_counting/chandan/prexRHRS_%d_-1*.root",run));


      TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
      c2->Divide(3,3);
      c2->cd(7);
      TH1F* huq = new TH1F("huq", "P.upQadcR;ADC CH;Events/CH", 200, 450, 650);
      T->Draw("P.upQadcR>>huq");

      TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
      huqCpy->GetXaxis()->SetRangeUser(upADCcutR_approx-10,upADCcutR_approx+10);
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

      TCut cut = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04)");
      TCut cut_wadc = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&P.upQadcR>%f)",upADCcut);
      TCut cut_wped = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&P.upQadcR<%f)",upADCcut);
      double q_edge = -0.076;
      TCut x_cut = Form("(R.tr.x[0]+1.3*R.tr.th[0]) > %f",q_edge);
      TCut y_cut = "";

      c2->cd(1);
      gPad->SetGridx();
      TH1F* hy1 = new TH1F("hy1", "Projected y on detector plane", 200, -0.03, 0.04);
      TH1F* hy2 = new TH1F("hy2", "Projected y on detector plane + adc> cut", 200, -0.03, 0.04);
      TH1F* hy3 = new TH1F("hy3", "Projected y on detector plane + adc< cut", 200, -0.03, 0.04);

      T->Project(hy1->GetName(), "R.tr.y[0]+1.3*R.tr.ph[0]", cut + x_cut); 
      T->Project(hy2->GetName(), "R.tr.y[0]+1.3*R.tr.ph[0]", cut_wadc + x_cut);
      T->Project(hy3->GetName(), "R.tr.y[0]+1.3*R.tr.ph[0]", cut_wped + x_cut);
      hy1->SetXTitle("R.tr.y[0]+1.3*R.tr.ph[0] (m)");
      hy1->SetLineColor(kBlack);
      hy2->SetLineColor(kRed);
      hy3->SetLineColor(kBlue);
      hy1->Draw();
      hy2->Draw("same");
      hy3->Draw("same");

      double sum_hy1 = hy1->Integral();
      double sum_hy2 = hy2->Integral();
      double frac_y = sum_hy2 / sum_hy1;

      TH1F* hy_clone = (TH1F*)hy2->Clone();
      hy_clone->SetTitle("ratio (y + adc>cut)/(y + adc<cut)");
      hy_clone->Divide(hy3);

      c2->cd(6);
      hy_clone->Draw();

      c2->cd(2);
      gPad->SetGridx();
      TH1F* hx1 = new TH1F("hx1", "Projected x on detector plane", 300, -0.2, 0.05);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc> cut", 300, -0.2, 0.05);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc< cut", 300, -0.2, 0.05);

      T->Project(hx1->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cut); 
      T->Project(hx2->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cut_wadc);
      T->Project(hx3->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cut_wped);
      hx1->SetXTitle("R.tr.x[0]+1.3*R.tr.th[0] (m)");
      hx1->SetLineColor(kBlack);
      hx2->SetLineColor(kRed);
      TF1 *fx=new TF1("fx","gaus",hx1->GetBinCenter(hx1->GetMaximumBin())-0.012,hx1->GetBinCenter(hx1->GetMaximumBin())+0.012);
      fx->SetLineWidth(2);
      fx->SetLineColor(kRed);
      hx2->Fit("fx","R");
      Double_t fxpar[3];
      fx->GetParameters(&fxpar[0]);
      hx3->SetLineColor(kBlue);
      hx1->Draw();
      hx2->Draw("same");
      hx3->Draw("same");

      TLine* qedge_line = new TLine(q_edge,0.0,q_edge,hx1->GetMaximum()/2.0);
      qedge_line->SetLineColor(6);
      qedge_line->SetLineWidth(2);
      qedge_line->Draw();

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(6);
      latex.DrawLatex(0.20,0.80,Form("q_edge = %.4f m",q_edge));
      latex.SetTextColor(2);
      latex.DrawLatex(0.20,0.85,Form("Elastic Peak = %.3f m",fxpar[1]));

      c2->cd(3);
      gPad->SetGridx();
      gPad->SetGridy();

      TH2F* h2d1 = new TH2F("h2d1", "Projected x vs y", 150, -0.08, 0.08, 150, -0.12, 0.04);
      T->Project(h2d1->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]:R.tr.y[0]+1.3*R.tr.ph[0]",cut);
      h2d1->SetXTitle("Projected y (m)");
      h2d1->SetYTitle("Projected x (m)");
      h2d1->Draw("COLZ");

      c2->cd(4);
      gPad->SetGridx();
      gPad->SetGridy();

      TH2F* h2d2 = new TH2F("h2d2", "Projected x vs y w/ adc > cut", 150, -0.08, 0.08, 150, -0.12, 0.04);
      float nAccept = T->Project(h2d2->GetName(),"R.tr.x[0]+1.3*R.tr.th[0]:R.tr.y[0]+1.3*R.tr.ph[0]",cut_wadc);
      h2d2->SetXTitle("Projected y (m)");
      h2d2->SetYTitle("Projected x (m)");
      h2d2->Draw("COLZ");

      c2->cd(5);
      gPad->SetGridx();
      gPad->SetGridy();
      TCut ePeak = "R.tr.y[0]+1.3*R.tr.ph[0]<-0.01 && R.tr.x[0]+1.3*R.tr.th[0]>-0.05";
      float nMiss = T->Draw("R.tr.x[0]+1.3*R.tr.th[0]:R.tr.y[0]+1.3*R.tr.ph[0]",cut_wped + ePeak);
      cout<<"Accept:"<<nAccept<<"\t"<<"Miss:"<<nMiss<<"\t"<<"Ratio:"<<nMiss/nAccept<<endl;

      TH2F* h2d3 = new TH2F("h2d3", "Projected x vs y adc < cut", 150, -0.08, 0.08, 150, -0.12, 0.04);
      T->Project(h2d3->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]:R.tr.y[0]+1.3*R.tr.ph[0]",cut_wped);
      h2d3->SetXTitle("Projected y (m)");
      h2d3->SetYTitle("Projected x (m)");
      h2d3->Draw("COLZ");

      c2->cd(8);
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogy(1);
      TH1F* huq_cor = new TH1F("huq_cor", "P.upQadcR w/adccut", 400, 0, 400);
      T->Project(huq_cor->GetName(), Form("(P.upQadcR-%f)",pedMean),Form("P.upQadcR>%f",upADCcut));
      huq_cor->SetXTitle("ADC CH");
      huq_cor->SetYTitle("Events/CH");
      huq_cor->Draw();
      TF1 *f1 = new TF1("f1","gaus",huq_cor->GetBinCenter(huq_cor->GetMaximumBin())-10,huq_cor->GetBinCenter(huq_cor->GetMaximumBin())+10);
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
      T->Project(huqpe->GetName(), Form("5.0*(P.upQadcR-%f)/%f/1.6",pedMean,USRgain),Form("P.upQadcR>%f",upADCcut));
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
      c2->SaveAs(Form("plots/USR_alignment_run%d.pdf",run));
    }

}
