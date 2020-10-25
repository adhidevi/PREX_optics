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
void spectraUS(int run){ 

  gROOT->Reset();
  gStyle->SetOptStat(112200);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
   

      TChain *T = new TChain("T");

  if(run < 10000){
      //LHRS
      T->Add(Form("/chafs1/work1/prex_counting/QsqRootFiles/prexLHRS_%d_-1*.root",run));

      double upADCcut_approx = 480;

      TCanvas* cADC = new TCanvas("cADC","LHRS US ADC",800,400);
      cADC->Divide(2,1);
      cADC->cd(1);
      TH1F* huq = new TH1F("huq","P.upQadcL;ADC CH;Event/CH",250,400,650);
      T->Draw("P.upQadcL>>huq");

      TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
      huqCpy->GetXaxis()->SetRangeUser(upADCcut_approx-10,upADCcut_approx+10);
      cout<<"ADC cut: "<<huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin())<<endl;
      double upADCcut = huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin());

      TLine* adc_line = new TLine(upADCcut,0.0,upADCcut,huq->GetMaximum());
      adc_line->SetLineColor(kMagenta);
      adc_line->SetLineWidth(2);
      adc_line->Draw();
      TF1 *fped=new TF1("fped","gaus",460,470);
      fped->SetLineWidth(2);
      fped->SetLineColor(kRed);
      huq->Fit("fped","R");
      Double_t fpedpar[3];
      fped->GetParameters(&fpedpar[0]);
      char labelres[64];
      sprintf(labelres,"#frac{RMS}{Mean} = %4.2f",fpedpar[2]/fpedpar[1]);
      TPaveLabel *pt = new TPaveLabel(0.65,0.38,0.83,0.44,labelres,"NDC");
      pt->SetBorderSize(0);
      pt->SetTextColor(kRed);
      pt->SetTextSize(0.75);
      pt->SetFillColor(0);
      pt->Draw();
      char labelped[64];
      sprintf(labelped,"Pedestal = %4.1f",fpedpar[1]);
      TPaveLabel *ptped = new TPaveLabel(0.65,0.31,0.83,0.37,labelped,"NDC");
      ptped->SetBorderSize(0);
      ptped->SetTextColor(kRed);
      ptped->SetTextSize(0.75);
      ptped->SetFillColor(0);
      ptped->Draw();
      char labelcut[64];
      sprintf(labelcut,"ADC_Cut = %3.1f",upADCcut);
      TPaveLabel *ptcut = new TPaveLabel(0.65,0.25,0.83,0.31,labelcut,"NDC");
      ptcut->SetBorderSize(0);
      ptcut->SetTextColor(kMagenta);
      ptcut->SetTextSize(0.75);
      ptcut->SetFillColor(0);
      ptcut->Draw();
      double pedMean = fpedpar[1];
      cout<<"pedMean:"<<pedMean<<endl;

      cADC->cd(2);
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
      TPaveLabel *ptPE = new TPaveLabel(0.60,0.20, 0.73,0.25,labelPE,"NDC");
      ptPE->SetBorderSize(0);
      ptPE->SetTextColor(kRed);
      ptPE->SetTextSize(1.0);
      ptPE->SetFillColor(0);
      ptPE->Draw();
      char labelRes[64];
      sprintf(labelRes,"#sigma_{cor}/MPV= %5.3f",Resolution_cor);
      TPaveLabel *ptRes = new TPaveLabel(0.60,0.26, 0.73,0.31,labelRes,"NDC");
      ptRes->SetBorderSize(0);
      ptRes->SetTextColor(kRed);
      ptRes->SetTextSize(1.0);
      ptRes->SetFillColor(0);
      ptRes->Draw();
      TPaveLabel *ptSig = new TPaveLabel(0.60,0.32, 0.73,0.37,label_corSigma,"NDC");
      ptSig->SetBorderSize(0);
      ptSig->SetTextColor(kRed);
      ptSig->SetTextSize(1.0);
      ptSig->SetFillColor(0);
      ptSig->Draw();

      cADC->SaveAs(Form("./temp/USL_adc_run%d.pdf",run));

      TCut cut = Form("(P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<-0.002");
      TCut cut_wadc = Form("(P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<-0.002&&P.upQadcL>%f",upADCcut);
      TCut cut_wped = Form("(P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<-0.002&&P.upQadcL<%f",upADCcut);
 
      TCanvas* c0 = new TCanvas("c0","p LHRS",1000,700);
      TH1F* h;
      T->Draw("L.gold.p[0]>>h(300,0.91,0.96)",cut);
      h = (TH1F*)gDirectory->FindObject("h");
      h->SetDirectory(gROOT);
      h->SetTitle(Form("Momentum LHRS (run%d);L.gold.p (GeV/c);",run));
      h->Draw();
      gPad->Update();
      h->SetStats(0);
      double el_peak_approx = h->GetBinCenter(h->GetMaximumBin());
      TF1* fitEP = new TF1("fitEP","gaus",el_peak_approx-0.0003,el_peak_approx+0.0003);
      h->Fit(fitEP,"R0");
      gPad->Modified();
      c0->SaveAs(Form("temp/p_distribution_run%d.pdf",run));
      
      double firstIP = 2.615/1000;
      double secondIP = 3.198/1000;
      double thirdIP = 3.475/1000;
      double fourthIP = 3.708/1000;
      double fstIP = fitEP->GetParameter(1)-firstIP;
      double sndIP = fitEP->GetParameter(1)-secondIP;
      double thrdIP = fitEP->GetParameter(1)-thirdIP;
      double frthIP = fitEP->GetParameter(1)-fourthIP;

      TCanvas* c0a = new TCanvas("c0a","p LHRS",1000,700);
      TH1F* ha = new TH1F("ha",Form("Momentum LHRS (run%d);L.gold.p (GeV/c);",run),300,0.94,0.955);
      TH1F* ha1 = new TH1F("ha1",Form("Momentum LHRS (run%d);L.gold.p (GeV/c);",run),300,0.94,0.955);
      TH1F* ha2 = new TH1F("ha2",Form("Momentum LHRS (run%d);L.gold.p (GeV/c);",run),300,0.94,0.955);
      ha->SetLineColor(1);
      ha1->SetLineColor(2);
      ha2->SetLineColor(4);
      ha->SetStats(0);
      ha1->SetStats(0);
      ha2->SetStats(0);
      T->Draw("L.gold.p[0]>>ha",cut);
      T->Draw("L.gold.p[0]>>ha1",cut_wadc,"same");
      T->Draw("L.gold.p[0]>>ha2",cut_wped,"same");
      TLine* El0 = new TLine(fitEP->GetParameter(1),0.0,fitEP->GetParameter(1),ha->GetMaximum());
      TLine* El1 = new TLine(fstIP,0.0,fstIP,ha->GetMaximum()*2/3.);
      TLine* El2 = new TLine(sndIP,0.0,sndIP,ha->GetMaximum()*2/3.);
      TLine* El3 = new TLine(thrdIP,0.0,thrdIP,ha->GetMaximum()*2/3.);
      TLine* El4 = new TLine(frthIP,0.0,frthIP,ha->GetMaximum()*2/3.);
      El0->SetLineColor(kMagenta);
      El1->SetLineColor(kGreen);
      El2->SetLineColor(kGreen+3);
      El3->SetLineColor(kCyan);
      El4->SetLineColor(kOrange+7);
      El0->SetLineWidth(2);
      El1->SetLineWidth(2);
      El2->SetLineWidth(2);
      El3->SetLineWidth(2);
      El4->SetLineWidth(2);
      El0->Draw();
      El1->Draw();
      El2->Draw();
      El3->Draw();
      El4->Draw();

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(1);
      latex.DrawLatex(0.20,0.85,"Total Flux");
      latex.SetTextColor(2);
      latex.DrawLatex(0.20,0.80,"Accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.20,0.75,"Missed");
      latex.SetTextColor(kMagenta);
      latex.DrawLatex(0.75,0.85,"Elastic");
      latex.SetTextColor(kGreen);
      latex.DrawLatex(0.75,0.80,Form("%.3f MeV",1000*firstIP));
      latex.SetTextColor(kGreen+3);
      latex.DrawLatex(0.75,0.75,Form("%.3f MeV",1000*secondIP));
      latex.SetTextColor(kCyan);
      latex.DrawLatex(0.75,0.70,Form("%.3f MeV",1000*thirdIP));
      latex.SetTextColor(kOrange+7);
      latex.DrawLatex(0.75,0.65,Form("%.3f MeV",1000*fourthIP));
      c0a->SaveAs(Form("temp/p_distribution_truncated_run%d.pdf",run));

      TCanvas* c0l = new TCanvas("c0l","p LHRS",1000,700);
      gPad->SetLogy();
      TH1F* hCpy = (TH1F*)h->Clone("hCpy");
      hCpy->SetTitle(Form("Momentum LHRS (run%d);L.gold.p (GeV/c);",run));
      hCpy->SetStats(0);
      hCpy->Draw();
      c0l->SaveAs(Form("temp/p_distribution_log_run%d.pdf",run));

      TCanvas* c0la = new TCanvas("c0la","p LHRS",1000,700);
      gPad->SetLogy();
      TH1F* haCpy = (TH1F*)ha->Clone("haCpy");
      TH1F* ha1Cpy = (TH1F*)ha1->Clone("ha1Cpy");
      TH1F* ha2Cpy = (TH1F*)ha2->Clone("ha2Cpy");
      haCpy->SetTitle(Form("Momentum LHRS (run%d);L.gold.p (GeV/c);",run));
      haCpy->SetStats(0);
      ha1Cpy->SetStats(0);
      ha2Cpy->SetStats(0);
      haCpy->Draw();
      ha1Cpy->Draw("same");
      ha2Cpy->Draw("same");
      El0->Draw();
      El1->Draw();
      El2->Draw();
      El3->Draw();
      El4->Draw();
      latex.SetTextColor(1);
      latex.DrawLatex(0.20,0.85,"Total Flux");
      latex.SetTextColor(2);
      latex.DrawLatex(0.20,0.80,"Accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.20,0.75,"Missed");
      latex.SetTextColor(kMagenta);
      latex.DrawLatex(0.75,0.85,"Elastic");
      latex.SetTextColor(kGreen);
      latex.DrawLatex(0.75,0.80,Form("%.3f MeV",1000*firstIP));
      latex.SetTextColor(kGreen+3);
      latex.DrawLatex(0.75,0.75,Form("%.3f MeV",1000*secondIP));
      latex.SetTextColor(kCyan);
      latex.DrawLatex(0.75,0.70,Form("%.3f MeV",1000*thirdIP));
      latex.SetTextColor(kOrange+7);
      latex.DrawLatex(0.75,0.65,Form("%.3f MeV",1000*fourthIP));
      c0l->SaveAs(Form("temp/p_distribution_truncated_log_run%d.pdf",run));

      TCanvas* c2 = new TCanvas("dp LHRS","dp LHRS",1000,700);
      T->Draw("L.tr.x[0]+1.3*L.tr.th[0]:L.gold.dp[0]>>hdp(500,-0.06,0.04,200,-1.8,1.0)",cut,"prof");
      TH2D* hdp = (TH2D*)gDirectory->FindObject("hdp");
      hdp->SetLineColor(kBlue);
      hdp->SetTitle(Form("Projected x on Det plane vs dp/p (run%d);#frac{dp}{p};Projected x on Det plane",run));
      TF1* fit = new TF1("fit","pol1",-0.04,-0.0025);
      fit->SetLineWidth(2);
      fit->SetLineColor(kRed);
      hdp->Fit("fit","R");
      Double_t fitpar[2];
      fit->GetParameters(&fitpar[0]);
      Double_t slope = fitpar[1];
//      Double_t slope = 14.35;
      cout<<Form("Dispersive Constant is %f per 1%s dp/p",slope,"%")<<endl;
      c2->SaveAs(Form("temp/det_X_vs_dp_upstream_run%d.pdf",run)); 

      TCanvas* c2FP = new TCanvas("dp LHRS FP","dp LHRS FP",1000,700);
      T->Draw("L.tr.x[0]:L.gold.dp[0]>>hdpFP(500,-0.06,0.04,200,-1.8,1.0)",cut,"prof");
      TH2D* hdpFP = (TH2D*)gDirectory->FindObject("hdpFP");
      hdpFP->SetLineColor(kBlue);
      hdpFP->SetTitle(Form("Dispersive x on focal plane vs dp/p (run%d);#frac{dp}{p};Dispersive x on focal plane",run));
      TF1* fitFP = new TF1("fitFP","pol1",-0.04,-0.0025);
      fitFP->SetLineWidth(2);
      fitFP->SetLineColor(kRed);
      hdpFP->Fit("fitFP","R");
      Double_t fitFPpar[2];
      fitFP->GetParameters(&fitFPpar[0]);
      Double_t slopeFP = fitFPpar[1];
//      Double_t slopeFP = 14.35;
      cout<<Form("Dispersive Constant at FP is %f per 1%s dp/p",slopeFP,"%")<<endl;
      c2FP->SaveAs(Form("temp/det_X_vs_dp_FP_upstream_run%d.pdf",run)); 

      TCanvas* c3 = new TCanvas("dp LHRS Scatter","dp LHRS Scatter",1000,700);
      TH2F* hdps = new TH2F("hdps", "Projected x on Det plane vs dp/p", 500, -0.02, 0.0, 200, -0.13, 0.02);
      T->Project(hdps->GetName(),"L.tr.x[0]+1.3*L.tr.th[0]:L.gold.dp[0]",cut);
      hdps->SetLineColor(kBlue);
      hdps->SetTitle(Form("Projected x on Det plane vs dp/p (run%d);#frac{dp}{p};Projected x on Det plane",run));
      hdps->SetStats(0);
      hdps->Draw();
      c3->SaveAs(Form("temp/det_X_vs_dp_Scatter_upstream_run%d.pdf",run)); 

      TCanvas* c1 = new TCanvas("spectra LHRS","spectra LHRS",1000,700);
      TH1F* hx1 = new TH1F("hx1", Form("Projected x on detector plane (run%d)",run), 200, -0.10, 0.07);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc>cut", 200, -0.10, 0.07);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc<cut", 200, -0.10, 0.07);

      T->Project(hx1->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]+0.0362", cut); 
      T->Project(hx2->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]+0.0362", cut_wadc);
      T->Project(hx3->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]+0.0362", cut_wped);
      hx1->SetXTitle("L.tr.x[0]+1.3*L.tr.th[0] (m)");
      hx1->SetLineColor(kBlack);
      hx2->SetLineColor(kRed);
      hx3->SetLineColor(kBlue);
      hx1->SetStats(0);
      hx2->SetStats(0);
      hx3->SetStats(0);
      hx1->Draw();
      hx2->Draw("same");
      hx3->Draw("same");
      int maxEntry = hx1->GetMaximum();
      TLine* line[5];
      double_t inep0 = -slope*0.0/9.53374/100;
      double_t inep1 = -slope*2.615/9.53374/100;
      double_t inep2 = -slope*3.198/9.53374/100;
      double_t inep3 = -slope*3.475/9.53374/100;
      double_t inep4 = -slope*3.708/9.53374/100;
      line[0] = new TLine(inep1,0.0,inep1,maxEntry/2.0);
      line[1] = new TLine(inep2,0.0,inep2,maxEntry/2.0);
      line[2] = new TLine(inep3,0.0,inep3,maxEntry/2.0);
      line[3] = new TLine(inep4,0.0,inep4,maxEntry/2.0);
      line[4] = new TLine(inep0,0.0,inep0,maxEntry+100);
      line[0]->SetLineColor(kGreen);
      line[1]->SetLineColor(kGreen+3);
      line[2]->SetLineColor(kCyan);
      line[3]->SetLineColor(kOrange+7);
      line[4]->SetLineColor(kMagenta);

      latex.SetTextColor(kBlack);
      latex.DrawLatex(0.20,0.85,"Spectra Pb-208");
      latex.SetTextColor(kMagenta);
      latex.DrawLatex(0.20,0.80,"Elastic");
      latex.SetTextColor(kGreen);
      latex.DrawLatex(0.20,0.75,"2.615MeV");
      latex.SetTextColor(kGreen+3);
      latex.DrawLatex(0.20,0.70,"3.198MeV");
      latex.SetTextColor(kCyan);
      latex.DrawLatex(0.20,0.65,"3.475MeV");
      latex.SetTextColor(kOrange+7);
      latex.DrawLatex(0.20,0.60,"3.708MeV");
    for(int i=0;i<5;i++)
      line[i]->Draw();

      int BinContent[100];
      int EdgeBin;
      int diff[100];
      cout<<"bin"<<"\t"<<"abs(diff)"<<endl;
      for(int i=1;i<=100;i++){
      diff[i] = abs(hx2->GetBinContent(i)-hx3->GetBinContent(i));
      cout<<i<<"\t"<<diff[i]<<endl;
      }
      int mindiff = {diff[1]};
      for(int i=1;i<=100;i++){
      if(mindiff>diff[i]){
      mindiff = diff[i];
      EdgeBin = i;
      }
      }
      double QEdge = hx1->GetBinCenter(EdgeBin);
      cout<<"Edge Bin: "<<EdgeBin<<"\t QEdge: "<<QEdge<<endl;

      TLine* edge = new TLine(QEdge,0.0,QEdge,hx1->GetMaximum()*2/3.);
      edge->SetLineColor(9);
      edge->SetLineWidth(2);
      edge->Draw();
      latex.SetTextColor(1);
      latex.DrawLatex(0.70,0.85,"Total Flux");
      latex.SetTextColor(2);
      latex.DrawLatex(0.70,0.80,"Accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.70,0.75,"Missed");
      latex.SetTextColor(9);
      latex.DrawLatex(0.70,0.70,Form("Quartz edge"));
      c1->SaveAs(Form("temp/spectra_Lead208_upstream_run%d.pdf",run));
      

      TCanvas* c1y = new TCanvas("projY","projY LHRS",1000,700);
      TH1F* hy1 = new TH1F("hy1", Form("Projected y on detector plane (run%d)",run), 200, -0.06, 0.01);
      TH1F* hy2 = new TH1F("hy2", "Projected y on detector plane + adc>cut", 200, -0.06, 0.01);
      TH1F* hy3 = new TH1F("hy3", "Projected y on detector plane + adc<cut", 200, -0.06, 0.01);

      TCut x_cut = Form("L.tr.x[0]+1.3*L.tr.th[0]>%f",QEdge-0.0362);
      T->Project(hy1->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cut+x_cut); 
      T->Project(hy2->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cut_wadc+x_cut);
      T->Project(hy3->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cut_wped+x_cut);
      hy1->SetXTitle("L.tr.y[0]+1.3*L.tr.ph[0] (m)");
      hy1->SetLineColor(kBlack);
      hy2->SetLineColor(kRed);
      hy3->SetLineColor(kBlue);
      hy1->SetStats(0);
      hy2->SetStats(0);
      hy3->SetStats(0);
      hy1->Draw();
      hy2->Draw("same");
      hy3->Draw("same");
      latex.SetTextColor(1);
      latex.DrawLatex(0.70,0.85,"Total Flux");
      latex.SetTextColor(2);
      latex.DrawLatex(0.70,0.80,"Accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.70,0.75,"Missed");

      c1->SaveAs(Form("temp/spectra_Lead208_upstream_run%d.pdf",run));
      TCanvas* cxy = new TCanvas("cxy","cxy",600,600);
      TH2F* hxy = new TH2F("hxy",Form("Projected x vs y on det. plane (run%d);L.tr.y[0]+1.3*L.tr.ph[0] (m);L.tr.x[0]+1.3*L.tr.th[0] (m)",run),200,-0.1,0.1,200,-0.1,0.1);
      hxy->SetStats(0);
      T->Draw("L.tr.x[0]+1.3*L.tr.th[0]:L.tr.y[0]+1.3*L.tr.ph[0]>>hxy",cut_wadc,"colz");

      TCanvas* c11 = new TCanvas("Pb208 spectra LHRS","Pb208 spectra LHRS",1000,700);
      TH1F* hx11 = new TH1F("hx11", Form("Pb208 spectrum on detector plane (run%d)",run), 200, 0.940, 0.955);
      TH1F* hx22 = new TH1F("hx22", "Pb208 spectra on detector plane + adc>cut", 200, 0.940, 0.955);
      TH1F* hx33 = new TH1F("hx33", "Pb208 spectra on detector plane + adc<cut", 200, 0.940, 0.955);

      T->Project(hx11->GetName(), Form("0.953374*((L.tr.x[0]+1.3*L.tr.th[0])/%f+1)",slope), cut); 
      T->Project(hx22->GetName(), Form("0.953374*((L.tr.x[0]+1.3*L.tr.th[0])/%f+1)",slope), cut_wadc);
      T->Project(hx33->GetName(), Form("0.953374*((L.tr.x[0]+1.3*L.tr.th[0])/%f+1)",slope), cut_wped);
      hx11->SetXTitle("Energy(GeV)");
      hx11->SetLineColor(kBlack);
      hx22->SetLineColor(kRed);
      hx33->SetLineColor(kBlue);
      hx11->SetStats(0);
      hx22->SetStats(0);
      hx33->SetStats(0);
      hx11->Draw();
      hx22->Draw("same");
      hx33->Draw("same");
      TLine* ln[6];
      int max = hx11->GetMaximum()+100;
      int maxBin = hx11->GetMaximumBin();
      float elasticPeak = hx11->GetXaxis()->GetBinCenter(maxBin);
      ln[0] = new TLine(0.953374,0.0,0.953374,max);
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
      beamE->Draw();
      c11->SaveAs(Form("./temp/Pb208_E_spectra_run%d.pdf",run));

      gSystem->Exec(Form("pdfunite ./temp/*_run%d.pdf ./plots/spectrum_run%d.pdf",run,run));
      gSystem->Exec(Form("rm -rf ./temp/*_run%d.pdf",run));
    }else{
      // RHRS
        double upADCcut_approx = 502;

      T->Add(Form("/chafs1/work1/prex_counting/QsqRootFiles/prexRHRS_%d_-1*.root",run));

      TCanvas* cADC = new TCanvas("cADC","RHRS US ADC",700,500);
      TH1F* huq = new TH1F("huq", "P.upQadcR;ADC CH;Events/CH", 200, 450, 650);
      T->Draw("P.upQadcR>>huq");

      TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
      huqCpy->GetXaxis()->SetRangeUser(upADCcut_approx-10,upADCcut_approx+10);
      cout<<"ADC cut: "<<huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin())<<endl;
      double upADCcut = huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin());
      TLine* adc_line = new TLine(upADCcut,0.0,upADCcut,huq->GetMaximum());
      adc_line->SetLineColor(kMagenta);
      adc_line->Draw();

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

      cADC->SaveAs(Form("./temp/USR_adc_run%d.pdf",run));

      TCut cut = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04)");
      TCut cut_wadc = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&P.upQadcR>%f)",upADCcut);
      TCut cut_wped = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&P.upQadcR<%f)",upADCcut);

      TCanvas* c0 = new TCanvas("c0","p RHRS",1000,700);
      TH1F* h;
      T->Draw("R.gold.p[0]>>h(300,0.91,0.96)",cut);
      h = (TH1F*)gDirectory->FindObject("h");
      h->SetDirectory(gROOT);
      h->SetTitle(Form("Momentum RHRS (run%d);R.gold.p (GeV/c)",run));
      h->Draw();
      c0->SaveAs(Form("temp/p_distribution_run%d.pdf",run));

      TCanvas* c0l = new TCanvas("c0l","p RHRS",1000,700);
      gPad->SetLogy();
      TH1F* hCpy = (TH1F*)h->Clone("hCpy");
      hCpy->SetTitle(Form("Momentum RHRS (run%d);R.gold.p (GeV/c);",run));
      hCpy->Draw();
      c0l->SaveAs(Form("temp/p_distribution_log_run%d.pdf",run));

      TCanvas* c2 = new TCanvas("dp RHRS","dp RHRS",1000,700);
      T->Draw("R.tr.x[0]+1.3*R.tr.th[0]:R.gold.dp[0]>>hdp(500,-0.05,0.01,200,-0.15,0.0)",cut,"prof");
      TH2F* hdp = (TH2F*)gDirectory->FindObject("hdp");
      hdp->SetLineColor(kBlue);
      hdp->SetTitle(Form("Projected x on Det plane vs dp/p (run%d);#frac{dp}{p};Projected x on Det plane",run));
      TF1* fit = new TF1("fit","pol1",-0.035,-0.002);
      fit->SetLineWidth(2);
      fit->SetLineColor(kRed);
      hdp->Fit("fit","R");
      Double_t fitpar[2];
      fit->GetParameters(&fitpar[0]);
      Double_t slope = fitpar[1];
//      Double_t slope = 14.26;
      cout<<Form("Dispersive Constant is %f per 1%s dp/p",slope,"%")<<endl;
      c2->SaveAs(Form("temp/det_X_vs_dp_upstream_run%d.pdf",run)); 

      TCanvas* c3 = new TCanvas("dp RHRS Scatter","dp RHRS Scatter",1000,700);
      TH2F* hdps = new TH2F("hdps", "Projected x on Det plane vs dp/p ", 500, -0.05, 0.01, 200, -0.15, 0.0);
      T->Project(hdps->GetName(),"R.tr.x[0]+1.3*R.tr.th[0]:R.gold.dp[0]",cut);
      hdps->SetLineColor(kBlue);
      hdps->SetTitle(Form("Projected x on Det plane vs dp/p (run%d);#frac{dp}{p};Projected x on Det plane",run));
      hdps->Draw();
      c3->SaveAs(Form("temp/det_X_vs_dp_scatter_upstream_run%d.pdf",run));

      TCanvas* c1 = new TCanvas("spectra RHRS","spectra RHRS",1000,700);
//      gPad->SetGridx();
      TH1F* hx1 = new TH1F("hx1", Form("Projected x on detector plane (run%d)",run), 200, -0.10, 0.07);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc> cut", 200, -0.10, 0.07);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc< cut", 200, -0.10, 0.07);

      T->Project(hx1->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]+0.0395", cut); 
      T->Project(hx2->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]+0.0395", cut_wadc);
      T->Project(hx3->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]+0.0395", cut_wped);
      hx1->SetXTitle("R.tr.x[0]+1.3*R.tr.th[0] (m)");
      hx1->SetLineColor(kBlack);
      hx2->SetLineColor(kRed);
      hx3->SetLineColor(kBlue);
      hx1->Draw();
      hx2->Draw("same");
      hx3->Draw("same");
      int maxEntry = hx1->GetMaximum();
      TLine* line[5];
      
      Double_t inep0 = -slope*0.0/9.53374/100;
      Double_t inep1 = -slope*2.614/9.53374/100;
      Double_t inep2 = -slope*3.198/9.53374/100;
      Double_t inep3 = -slope*3.475/9.53374/100;
      Double_t inep4 = -slope*3.708/9.53374/100;
      line[0] = new TLine(inep1,0.0,inep1,maxEntry/2.0);
      line[1] = new TLine(inep2,0.0,inep2,maxEntry/2.0);
      line[2] = new TLine(inep3,0.0,inep3,maxEntry/2.0);
      line[3] = new TLine(inep4,0.0,inep4,maxEntry/2.0);
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

      c1->SaveAs(Form("temp/spectra_Lead208_upstream_run%d.pdf",run));

     TCanvas* c11 = new TCanvas("Pb208 spectra RHRS","Pb208 spectra RHRS",1000,700);
      TH1F* hx11 = new TH1F("hx11", Form("Pb208 spectrum on detector plane (run%d)",run), 200, 0.940, 0.9525);
      TH1F* hx22 = new TH1F("hx22", "Pb208 spectra on detector plane + adc>cut", 200, 0.940, 0.9525);
      TH1F* hx33 = new TH1F("hx33", "Pb208 spectra on detector plane + adc<cut", 200, 0.940, 0.9525);

      T->Project(hx11->GetName(), Form("0.953374*((R.tr.x[0]+1.3*R.tr.th[0])/%f+1)",slope), cut); 
      T->Project(hx22->GetName(), Form("0.953374*((R.tr.x[0]+1.3*R.tr.th[0])/%f+1)",slope), cut_wadc);
      T->Project(hx33->GetName(), Form("0.953374*((R.tr.x[0]+1.3*R.tr.th[0])/%f+1)",slope), cut_wped);
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
      ln[0] = new TLine(0.953374,0.0,0.953374,max);
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
      c11->SaveAs(Form("./temp/Pb208_E_spectra_run%d.pdf",run));

      gSystem->Exec(Form("pdfunite ./temp/*_run%d.pdf ./plots/spectrum_run%d.pdf",run,run));
      gSystem->Exec(Form("rm -rf ./temp/*_run%d.pdf",run));
}

}
