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
void projXY_and_findQsq(int run, TString target, TString date, double Ebeam){ 

  gROOT->Reset();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.5);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
   
      TChain *T = new TChain("T");

  if(run < 10000){
      //LHRS
      T->Add(Form("/chafs1/work1/prex_counting/Qsq_Oct22/prexLHRS_%d_-1*.root",run));

      double upADCcut_approx = 480;
      TCut cut = Form("(P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002");

      TCanvas* cSlp = new TCanvas("cSlp","slope",600,400);
      T->Draw("L.tr.x[0]+1.3*L.tr.th[0]:L.gold.dp[0]>>hdp(500,-0.05,0.01,200,-0.15,0.0)",cut,"prof");
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


      TCanvas* c1 = new TCanvas("c1","LHRS US ADC",800,700);
      c1->Divide(2,2);
      c1->cd(1);
      TH1F* huq = new TH1F("huq",Form("P.upQadcL (run%d), %s, %s;ADC CH;Event/CH",run,target.Data(),date.Data()),250,400,650);
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

      TCut cut_wadc = Form("(P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&P.upQadcL>%f",upADCcut);
      TCut cut_wped = Form("(P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&P.upQadcL<%f",upADCcut);
 
      double firstIP = 2.615/1000;
      double secondIP = 3.198/1000;
      double thirdIP = 3.475/1000;
      double fourthIP = 3.708/1000;
 
      c1->cd(2);
      TH1F* hx1 = new TH1F("hx1", Form("Projected x on detector plane (run%d)",run), 200, -0.14, 0.03);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc>cut", 200, -0.14, 0.03);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc<cut", 200, -0.14, 0.03);

      T->Project(hx1->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cut); 
      T->Project(hx2->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cut_wadc);
      T->Project(hx3->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cut_wped);
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
      TF1* fitEP = new TF1("fitEP","gaus",hx2->GetBinCenter(hx2->GetMaximumBin())-0.007,hx2->GetBinCenter(hx2->GetMaximumBin())+0.007);
      hx2->Fit(fitEP,"R0");

      TLine* line[5];
      double_t inep0 = -slope*0.0/Ebeam/1000+fitEP->GetParameter(1);
      double_t inep1 = -slope*2.615/Ebeam/1000+fitEP->GetParameter(1);
      double_t inep2 = -slope*3.198/Ebeam/1000+fitEP->GetParameter(1);
      double_t inep3 = -slope*3.475/Ebeam/1000+fitEP->GetParameter(1);
      double_t inep4 = -slope*3.708/Ebeam/1000+fitEP->GetParameter(1);
      line[0] = new TLine(inep1,0.0,inep1,maxEntry/2.0);
      line[1] = new TLine(inep2,0.0,inep2,maxEntry/2.0);
      line[2] = new TLine(inep3,0.0,inep3,maxEntry/2.0);
      line[3] = new TLine(inep4,0.0,inep4,maxEntry/2.0);
      line[4] = new TLine(inep0,0.0,inep0,maxEntry+100);
      line[0]->SetLineWidth(2);
      line[1]->SetLineWidth(2);
      line[2]->SetLineWidth(2);
      line[3]->SetLineWidth(2);
      line[4]->SetLineWidth(2);
      line[0]->SetLineColor(kGreen);
      line[1]->SetLineColor(kGreen+3);
      line[2]->SetLineColor(kCyan);
      line[3]->SetLineColor(kOrange+7);
      line[4]->SetLineColor(kMagenta);

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
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
      int peakBin = hx2->GetMaximumBin();
      int BinContent[peakBin];
      int EdgeBin;
      int diff[peakBin];
      cout<<"bin"<<"\t"<<"abs(diff)"<<endl;
      for(int i=1;i<=peakBin;i++){
      diff[i] = abs(hx2->GetBinContent(i)-hx3->GetBinContent(i));
      cout<<i<<"\t"<<diff[i]<<endl;
      }
      int mindiff = {diff[1]};
      for(int i=1;i<=peakBin;i++){
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
      latex.DrawLatex(0.70,0.70,Form("Q_{edge}@ %.4f",QEdge));

      c1->cd(3);
      TH1F* hy1 = new TH1F("hy1", Form("Projected y on detector plane (run%d)",run), 200, -0.06, 0.01);
      TH1F* hy2 = new TH1F("hy2", "Projected y on detector plane + adc>cut", 200, -0.06, 0.01);
      TH1F* hy3 = new TH1F("hy3", "Projected y on detector plane + adc<cut", 200, -0.06, 0.01);

      TCut xCut = Form("L.tr.x[0]+1.3*L.tr.th[0]>%f",QEdge);
      T->Project(hy1->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cut+xCut); 
      T->Project(hy2->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cut_wadc+xCut);
      T->Project(hy3->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cut_wped+xCut);
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

      TLine* edge_lo = new TLine(-0.0478,0.0,-0.0478,hy1->GetMaximum()*2/3.);
      edge_lo->SetLineColor(9);
      edge_lo->SetLineWidth(2);
      edge_lo->Draw();

      TLine* edge_hi = new TLine(-0.012,0.0,-0.012,hy1->GetMaximum()*2/3.);
      edge_hi->SetLineColor(9);
      edge_hi->SetLineWidth(2);
      edge_hi->Draw();

      double ylo = -0.0478;
      double yhi = -0.0120;
      TCut ylocut = Form("L.tr.y[0]+1.3*L.tr.ph[0]>%f",ylo);
      TCut yhicut = Form("L.tr.y[0]+1.3*L.tr.ph[0]<%f",yhi);

      latex.SetTextColor(1);
      latex.DrawLatex(0.70,0.85,"Total Flux");
      latex.SetTextColor(2);
      latex.DrawLatex(0.70,0.80,"Accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.70,0.75,"Missed");
      latex.SetTextColor(9);
      latex.DrawLatex(0.70,0.70,"Q_{edge}");

      c1->cd(4);
      TH2F* hxy = new TH2F("hxy",Form("Projected x vs y on det. plane (run%d);L.tr.y[0]+1.3*L.tr.ph[0] (m);L.tr.x[0]+1.3*L.tr.th[0] (m)",run),200,-0.1,0.1,200,-0.1,0.1);
      hxy->SetStats(0);
      T->Draw("L.tr.x[0]+1.3*L.tr.th[0]:L.tr.y[0]+1.3*L.tr.ph[0]>>hxy",cut_wadc,"colz");
      TBox* qbox = new TBox(ylo,QEdge,yhi,QEdge+0.16);
      qbox->SetFillStyle(0);
      qbox->SetFillColor(0);
      qbox->SetLineColor(2);
      qbox->SetLineWidth(2);
      qbox->Draw();
      c1->SaveAs(Form("./temp1/projectXY_run%d.pdf",run));
     
      Float_t theta0 = 4.806;//degrees
      theta0 = theta0*TMath::Pi()/180;//radians

      Float_t s0 = TMath::Sin(theta0);
      Float_t c0 = TMath::Cos(theta0);

      TCanvas* cQsq = new TCanvas("cQsq","cQsq",800,500);
      cQsq->Divide(2,1);
      cQsq->cd(1);
      gPad->Modified();gPad->Update();
      TH1F* hp1 = new TH1F("hp1",Form("LHRS momentum (run%d)",run),200,0.94,0.955);
      TH1F* hp2 = new TH1F("hp2","ADC passed",200,0.94,0.955);
      TH1F* hp3 = new TH1F("hp3","ADC missed",200,0.94,0.955);
      hp1->SetLineColor(1);
      hp2->SetLineColor(2);
      hp3->SetLineColor(4);
      T->Draw("L.gold.p>>hp1",cut);
      T->Draw("L.gold.p>>hp2",cut_wadc,"sames");
      T->Draw("L.gold.p>>hp3",cut_wped,"sames");
      double peak_approx = hp2->GetBinCenter(hp2->GetMaximumBin());
      TF1* momfit = new TF1("momfit","gaus",peak_approx-0.0004,peak_approx+0.0004);
      hp2->Fit(momfit,"R0");

      int pBin = hp2->GetMaximumBin();
      int BinContentP[pBin];
      int EdgeBinP;
      int diffP[pBin];
      cout<<"bin"<<"\t"<<"abs(diff)"<<endl;
      for(int i=1;i<=pBin;i++){
      diffP[i] = abs(hp2->GetBinContent(i)-hp3->GetBinContent(i));
      cout<<i<<"\t"<<diffP[i]<<endl;
      }
      int mindiffP = {diffP[1]};
      for(int i=1;i<=pBin;i++){
      if(mindiffP>diffP[i]){
      mindiffP = diffP[i];
      EdgeBinP = i;
      }
      }
      double peakpos = momfit->GetParameter(1);;
      double efact = (Ebeam-0.001)/peakpos;

      double EdgeP = hp1->GetBinCenter(EdgeBinP);
      cout<<"Edge Bin for P: "<<EdgeBinP<<"\t Edge P: "<<EdgeP<<endl;

      latex.SetTextColor(6);
      latex.DrawLatex(0.13,0.85,Form("p@edge = %1.4f m",EdgeP));
      latex.SetTextColor(1);
      latex.DrawLatex(0.15,0.80,"All events");
      latex.SetTextColor(2);
      latex.DrawLatex(0.15,0.75,"US accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.15,0.70,"US missed");
      latex.SetTextColor(46);
      latex.DrawLatex(0.13,0.65,Form("el_peak = %1.4f m",peakpos));
      hp1->SetStats(0);
      hp2->SetStats(0);
      hp3->SetStats(0);
      TLine* pEdge = new TLine(EdgeP,0.0,EdgeP,hp1->GetMaximum()/2.2);
      pEdge->SetLineColor(6);
      pEdge->SetLineWidth(2);
      pEdge->Draw();
      TLine* ppeak = new TLine(peakpos,0.0,peakpos,hp1->GetMaximum());
      ppeak->SetLineColor(46);
      ppeak->SetLineWidth(2);
      ppeak->Draw();

  TCut pCut = Form("L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&(P.evtypebits&2)==2&&L.gold.p> %f",EdgeP);

      cQsq->cd(2);
      TH1F *Q2adc= new TH1F("Q2adc",Form("Q^{2} (GeV/c)^{2} (run%d), %s, %s;Q^(2) (GeV/c)^{2};",run,target.Data(),date.Data()),150,0,0.015);
      TH1F *Q2xcut= new TH1F("Q2xcut",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
      TH1F *Q2pcut= new TH1F("Q2pcut",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
      TH1F *Q2xylocut= new TH1F("Q2xylocut",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
      TH1F *Q2xyhicut= new TH1F("Q2xyhicut",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
      Q2adc->SetLineColor(1);
      Q2xcut->SetLineColor(2);
      Q2pcut->SetLineColor(4);
      Q2xylocut->SetLineColor(6);
      Q2xyhicut->SetLineColor(8);
      T->Draw(Form("2*%f*%f*L.gold.p*(1-((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))))>>Q2adc",Ebeam,efact,c0,s0),cut_wadc,"hist");
      T->Draw(Form("2*%f*%f*L.gold.p*(1-((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))))>>Q2xcut",Ebeam,efact,c0,s0),cut+xCut,"hist sames");
      T->Draw(Form("2*%f*%f*L.gold.p*(1-((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))))>>Q2pcut",Ebeam,efact,c0,s0),pCut,"hist sames");
      T->Draw(Form("2*%f*%f*L.gold.p*(1-((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))))>>Q2xylocut",Ebeam,efact,c0,s0),cut+xCut+ylocut,"hist sames");
      T->Draw(Form("2*%f*%f*L.gold.p*(1-((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))))>>Q2xyhicut",Ebeam,efact,c0,s0),cut+xCut+yhicut,"hist sames");
      Q2xcut->Scale(Q2adc->Integral()/Q2xcut->Integral());
      Q2pcut->Scale(Q2adc->Integral()/Q2pcut->Integral());
      Q2xylocut->Scale(Q2adc->Integral()/Q2xylocut->Integral());
      Q2xyhicut->Scale(Q2adc->Integral()/Q2xyhicut->Integral());
      gPad->Update();
      TPaveStats* statQ1 = (TPaveStats*)Q2adc->FindObject("stats");
      TPaveStats* statQ2 = (TPaveStats*)Q2xcut->FindObject("stats");
      TPaveStats* statQ3 = (TPaveStats*)Q2pcut->FindObject("stats");
      TPaveStats* statQ4 = (TPaveStats*)Q2xylocut->FindObject("stats");
      TPaveStats* statQ5 = (TPaveStats*)Q2xyhicut->FindObject("stats");
      statQ1->SetY2NDC(0.90);
      statQ1->SetY1NDC(0.75);
      statQ2->SetY2NDC(0.75);
      statQ2->SetY1NDC(0.60);
      statQ3->SetY2NDC(0.60);
      statQ3->SetY1NDC(0.45);
      statQ4->SetY2NDC(0.45);
      statQ4->SetY1NDC(0.30);
      statQ5->SetY2NDC(0.30);
      statQ5->SetY1NDC(0.15);
      statQ1->SetTextColor(1);
      statQ2->SetTextColor(2);
      statQ3->SetTextColor(4);
      statQ4->SetTextColor(6);
      statQ5->SetTextColor(8);
      gPad->Modified();
      double Q2_adc = Q2adc->GetMean();
      double Q2_xcut = Q2xcut->GetMean();
      double Q2_pcut = Q2pcut->GetMean();
      double Q2_xylocut = Q2xylocut->GetMean();
      double Q2_xyhicut = Q2xyhicut->GetMean();
      latex.SetTextSize(0.04);
      latex.SetTextColor(2);
      latex.DrawLatex(0.65,0.70,Form("%.2f%s",Q2_xcut/Q2_adc*100.,"%"));
      latex.SetTextColor(4);
      latex.DrawLatex(0.65,0.55,Form("%.2f%s",Q2_pcut/Q2_adc*100.,"%"));
      latex.SetTextColor(6);
      latex.DrawLatex(0.65,0.40,Form("%.2f%s",Q2_xylocut/Q2_adc*100.,"%"));
      latex.SetTextColor(8);
      latex.DrawLatex(0.65,0.25,Form("%.2f%s",Q2_xyhicut/Q2_adc*100.,"%"));
      
      cQsq->SaveAs(Form("./temp1/qsq_analysis_run%d.pdf",run));
      ofstream outfile("./TextFiles/qsquare_diagnostic_quantitiesL.csv",ios_base::app);
      outfile<<run<<"\t"<<target<<"\t"<<date<<"\t"<<upADCcut<<"\t"<<fitEP->GetParameter(1)<<"\t"<<QEdge<<"\t"<<peakpos<<"\t"<<EdgeP<<"\t"<<Q2adc->GetMean()<<"\t"<<Q2adc->GetRMS()<<"\t"<<Q2xcut->GetMean()<<"\t"<<Q2xcut->GetRMS()<<"\t"<<Q2xylocut->GetMean()<<"\t"<<Q2xylocut->GetRMS()<<"\t"<<Q2xyhicut->GetMean()<<"\t"<<Q2xyhicut->GetRMS()<<endl;
      outfile.close();
    }else{
      //RHRS
      T->Add(Form("/chafs1/work1/prex_counting/Qsq_Oct22/prexRHRS_%d_-1*.root",run));

      double upADCcut_approx = 502;
      TCut cut = Form("(P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<0.002");

      TCanvas* cSlp = new TCanvas("cSlp","slope",600,400);
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


      TCanvas* c1 = new TCanvas("c1","LHRS US ADC",800,700);
      c1->Divide(2,2);
      c1->cd(1);
      TH1F* huq = new TH1F("huq",Form("P.upQadcR (run%d), %s, %s;ADC CH;Event/CH",run,target.Data(),date.Data()),250,400,650);
      T->Draw("P.upQadcR>>huq");

      TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
      huqCpy->GetXaxis()->SetRangeUser(upADCcut_approx-10,upADCcut_approx+10);
      cout<<"ADC cut: "<<huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin())<<endl;
      double upADCcut = huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin());

      TLine* adc_line = new TLine(upADCcut,0.0,upADCcut,huq->GetMaximum());
      adc_line->SetLineColor(kMagenta);
      adc_line->SetLineWidth(2);
      adc_line->Draw();
      TF1 *fped=new TF1("fped","gaus",huq->GetBinCenter(huq->GetMaximumBin())-5,huq->GetBinCenter(huq->GetMaximumBin())+5);
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

      TCut cut_wadc = Form("(P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<0.002&&P.upQadcR>%f",upADCcut);
      TCut cut_wped = Form("(P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<0.002&&P.upQadcR<%f",upADCcut);
 
      double firstIP = 2.615/1000;
      double secondIP = 3.198/1000;
      double thirdIP = 3.475/1000;
      double fourthIP = 3.708/1000;
 
      c1->cd(2);
      TH1F* hx1 = new TH1F("hx1", Form("Projected x on detector plane (run%d)",run), 200, -0.14, 0.03);
      TH1F* hx2 = new TH1F("hx2", "Projected x on detector plane + adc>cut", 200, -0.14, 0.03);
      TH1F* hx3 = new TH1F("hx3", "Projected x on detector plane + adc<cut", 200, -0.14, 0.03);

      T->Project(hx1->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cut); 
      T->Project(hx2->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cut_wadc);
      T->Project(hx3->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cut_wped);
      hx1->SetXTitle("R.tr.x[0]+1.3*R.tr.th[0] (m)");
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
      TF1* fitEP = new TF1("fitEP","gaus",hx2->GetBinCenter(hx2->GetMaximumBin())-0.007,hx2->GetBinCenter(hx2->GetMaximumBin())+0.007);
      hx2->Fit(fitEP,"R0");

      TLine* line[5];
      double_t inep0 = -slope*0.0/Ebeam/1000+fitEP->GetParameter(1);
      double_t inep1 = -slope*2.615/Ebeam/1000+fitEP->GetParameter(1);
      double_t inep2 = -slope*3.198/Ebeam/1000+fitEP->GetParameter(1);
      double_t inep3 = -slope*3.475/Ebeam/1000+fitEP->GetParameter(1);
      double_t inep4 = -slope*3.708/Ebeam/1000+fitEP->GetParameter(1);
      line[0] = new TLine(inep1,0.0,inep1,maxEntry/2.0);
      line[1] = new TLine(inep2,0.0,inep2,maxEntry/2.0);
      line[2] = new TLine(inep3,0.0,inep3,maxEntry/2.0);
      line[3] = new TLine(inep4,0.0,inep4,maxEntry/2.0);
      line[4] = new TLine(inep0,0.0,inep0,maxEntry+100);
      line[0]->SetLineWidth(2);
      line[1]->SetLineWidth(2);
      line[2]->SetLineWidth(2);
      line[3]->SetLineWidth(2);
      line[4]->SetLineWidth(2);
      line[0]->SetLineColor(kGreen);
      line[1]->SetLineColor(kGreen+3);
      line[2]->SetLineColor(kCyan);
      line[3]->SetLineColor(kOrange+7);
      line[4]->SetLineColor(kMagenta);

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
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
      int peakBin = hx2->GetMaximumBin();
      int BinContent[peakBin];
      int EdgeBin;
      int diff[peakBin];
      cout<<"bin"<<"\t"<<"abs(diff)"<<endl;
      for(int i=1;i<=peakBin;i++){
      diff[i] = abs(hx2->GetBinContent(i)-hx3->GetBinContent(i));
      cout<<i<<"\t"<<diff[i]<<endl;
      }
      int mindiff = {diff[1]};
      for(int i=1;i<=peakBin;i++){
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
      latex.DrawLatex(0.70,0.70,Form("Q_{edge}@ %.4f",QEdge));

      c1->cd(3);
      TH1F* hy1 = new TH1F("hy1", Form("Projected y on detector plane (run%d)",run), 200, -0.03, 0.04);
      TH1F* hy2 = new TH1F("hy2", "Projected y on detector plane + adc>cut", 200, -0.03, 0.04);
      TH1F* hy3 = new TH1F("hy3", "Projected y on detector plane + adc<cut", 200, -0.03, 0.04);

      TCut xCut = Form("R.tr.x[0]+1.3*R.tr.th[0]>%f",QEdge);
      T->Project(hy1->GetName(), "R.tr.y[0]+1.3*R.tr.ph[0]", cut+xCut); 
      T->Project(hy2->GetName(), "R.tr.y[0]+1.3*R.tr.ph[0]", cut_wadc+xCut);
      T->Project(hy3->GetName(), "R.tr.y[0]+1.3*R.tr.ph[0]", cut_wped+xCut);
      hy1->SetXTitle("R.tr.y[0]+1.3*R.tr.ph[0] (m)");
      hy1->SetLineColor(kBlack);
      hy2->SetLineColor(kRed);
      hy3->SetLineColor(kBlue);
      hy1->SetStats(0);
      hy2->SetStats(0);
      hy3->SetStats(0);
      hy1->Draw();
      hy2->Draw("same");
      hy3->Draw("same");

      double ylo = -0.012;
      double yhi = 0.023;
      TLine* edge_lo = new TLine(ylo,0.0,ylo,hy1->GetMaximum()*2/3.);
      edge_lo->SetLineColor(9);
      edge_lo->SetLineWidth(2);
      edge_lo->Draw();

      TLine* edge_hi = new TLine(yhi,0.0,yhi,hy1->GetMaximum()*2/3.);
      edge_hi->SetLineColor(9);
      edge_hi->SetLineWidth(2);
      edge_hi->Draw();

      TCut ylocut = Form("R.tr.y[0]+1.3*R.tr.ph[0]>%f",ylo);
      TCut yhicut = Form("R.tr.y[0]+1.3*R.tr.ph[0]<%f",yhi);

      latex.SetTextColor(1);
      latex.DrawLatex(0.70,0.85,"Total Flux");
      latex.SetTextColor(2);
      latex.DrawLatex(0.70,0.80,"Accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.70,0.75,"Missed");
      latex.SetTextColor(9);
      latex.DrawLatex(0.70,0.70,"Q_{edge}");

      c1->cd(4);
      TH2F* hxy = new TH2F("hxy",Form("Projected x vs y on det. plane (run%d);R.tr.y[0]+1.3*R.tr.ph[0] (m);R.tr.x[0]+1.3*R.tr.th[0] (m)",run),200,-0.1,0.1,200,-0.1,0.1);
      hxy->SetStats(0);
      T->Draw("R.tr.x[0]+1.3*R.tr.th[0]:R.tr.y[0]+1.3*R.tr.ph[0]>>hxy",cut_wadc,"colz");
      TBox* qbox = new TBox(ylo,QEdge,yhi,QEdge+0.16);
      qbox->SetFillStyle(0);
      qbox->SetFillColor(0);
      qbox->SetLineColor(2);
      qbox->SetLineWidth(2);
      qbox->Draw();
      c1->SaveAs(Form("./temp1/projectXY_run%d.pdf",run));
     
      Float_t theta0 = 4.764;//degrees
      theta0 = theta0*TMath::Pi()/180;//radians

      Float_t s0 = TMath::Sin(theta0);
      Float_t c0 = TMath::Cos(theta0);

      TCanvas* cQsq = new TCanvas("cQsq","cQsq",800,500);
      cQsq->Divide(2,1);
      cQsq->cd(1);
      gPad->Modified();gPad->Update();
      TH1F* hp1 = new TH1F("hp1",Form("LHRS momentum (run%d)",run),200,0.94,0.955);
      TH1F* hp2 = new TH1F("hp2","ADC passed",200,0.94,0.955);
      TH1F* hp3 = new TH1F("hp3","ADC missed",200,0.94,0.955);
      hp1->SetLineColor(1);
      hp2->SetLineColor(2);
      hp3->SetLineColor(4);
      T->Draw("R.gold.p>>hp1",cut);
      T->Draw("R.gold.p>>hp2",cut_wadc,"sames");
      T->Draw("R.gold.p>>hp3",cut_wped,"sames");
      double peak_approx = hp2->GetBinCenter(hp2->GetMaximumBin());
      TF1* momfit = new TF1("momfit","gaus",peak_approx-0.0004,peak_approx+0.0004);
      hp2->Fit(momfit,"R0");

      int pBin = hp2->GetMaximumBin();
      int BinContentP[pBin];
      int EdgeBinP;
      int diffP[pBin];
      cout<<"bin"<<"\t"<<"abs(diff)"<<endl;
      for(int i=1;i<=pBin;i++){
      diffP[i] = abs(hp2->GetBinContent(i)-hp3->GetBinContent(i));
      cout<<i<<"\t"<<diffP[i]<<endl;
      }
      int mindiffP = {diffP[1]};
      for(int i=1;i<=pBin;i++){
      if(mindiffP>diffP[i]){
      mindiffP = diffP[i];
      EdgeBinP = i;
      }
      }
      double peakpos = momfit->GetParameter(1);;
      double efact = (Ebeam-0.001)/peakpos;

      double EdgeP = hp1->GetBinCenter(EdgeBinP);
      cout<<"Edge Bin for P: "<<EdgeBinP<<"\t Edge P: "<<EdgeP<<endl;

      latex.SetTextColor(6);
      latex.DrawLatex(0.13,0.85,Form("p@edge = %1.4f m",EdgeP));
      latex.SetTextColor(1);
      latex.DrawLatex(0.15,0.80,"All events");
      latex.SetTextColor(2);
      latex.DrawLatex(0.15,0.75,"US accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.15,0.70,"US missed");
      latex.SetTextColor(46);
      latex.DrawLatex(0.13,0.65,Form("el_peak = %1.4f m",peakpos));
      hp1->SetStats(0);
      hp2->SetStats(0);
      hp3->SetStats(0);
      TLine* pEdge = new TLine(EdgeP,0.0,EdgeP,hp1->GetMaximum()/2.2);
      pEdge->SetLineColor(6);
      pEdge->SetLineWidth(2);
      pEdge->Draw();
      TLine* ppeak = new TLine(peakpos,0.0,peakpos,hp1->GetMaximum());
      ppeak->SetLineColor(46);
      ppeak->SetLineWidth(2);
      ppeak->Draw();

  TCut pCut = Form("R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<0.002&&(P.evtypebits&2)==2&&R.gold.p> %f",EdgeP);

      cQsq->cd(2);
      TH1F *Q2adc= new TH1F("Q2adc",Form("Q^{2} (GeV/c)^{2} (run%d), %s, %s;Q^(2) (GeV/c)^{2};",run,target.Data(),date.Data()),150,0,0.015);
      TH1F *Q2xcut= new TH1F("Q2xcut",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
      TH1F *Q2pcut= new TH1F("Q2pcut",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
      TH1F *Q2xylocut= new TH1F("Q2xylocut",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
      TH1F *Q2xyhicut= new TH1F("Q2xyhicut",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
      Q2adc->SetLineColor(1);
      Q2xcut->SetLineColor(2);
      Q2pcut->SetLineColor(4);
      Q2xylocut->SetLineColor(6);
      Q2xyhicut->SetLineColor(8);
      T->Draw(Form("2*%f*%f*R.gold.p*(1-((%f+R.gold.ph*%f)/(TMath::Sqrt(1.0+R.gold.th*R.gold.th+R.gold.ph*R.gold.ph))))>>Q2adc",Ebeam,efact,c0,s0),cut_wadc,"hist");
      T->Draw(Form("2*%f*%f*R.gold.p*(1-((%f+R.gold.ph*%f)/(TMath::Sqrt(1.0+R.gold.th*R.gold.th+R.gold.ph*R.gold.ph))))>>Q2xcut",Ebeam,efact,c0,s0),cut+xCut,"hist sames");
      T->Draw(Form("2*%f*%f*R.gold.p*(1-((%f+R.gold.ph*%f)/(TMath::Sqrt(1.0+R.gold.th*R.gold.th+R.gold.ph*R.gold.ph))))>>Q2pcut",Ebeam,efact,c0,s0),pCut,"hist sames");
      T->Draw(Form("2*%f*%f*R.gold.p*(1-((%f+R.gold.ph*%f)/(TMath::Sqrt(1.0+R.gold.th*R.gold.th+R.gold.ph*R.gold.ph))))>>Q2xylocut",Ebeam,efact,c0,s0),cut+xCut+ylocut,"hist sames");
      T->Draw(Form("2*%f*%f*R.gold.p*(1-((%f+R.gold.ph*%f)/(TMath::Sqrt(1.0+R.gold.th*R.gold.th+R.gold.ph*R.gold.ph))))>>Q2xyhicut",Ebeam,efact,c0,s0),cut+xCut+yhicut,"hist sames");
      Q2xcut->Scale(Q2adc->Integral()/Q2xcut->Integral());
      Q2pcut->Scale(Q2adc->Integral()/Q2pcut->Integral());
      Q2xylocut->Scale(Q2adc->Integral()/Q2xylocut->Integral());
      Q2xyhicut->Scale(Q2adc->Integral()/Q2xyhicut->Integral());
      gPad->Update();
      TPaveStats* statQ1 = (TPaveStats*)Q2adc->FindObject("stats");
      TPaveStats* statQ2 = (TPaveStats*)Q2xcut->FindObject("stats");
      TPaveStats* statQ3 = (TPaveStats*)Q2pcut->FindObject("stats");
      TPaveStats* statQ4 = (TPaveStats*)Q2xylocut->FindObject("stats");
      TPaveStats* statQ5 = (TPaveStats*)Q2xyhicut->FindObject("stats");
      statQ1->SetY2NDC(0.90);
      statQ1->SetY1NDC(0.75);
      statQ2->SetY2NDC(0.75);
      statQ2->SetY1NDC(0.60);
      statQ3->SetY2NDC(0.60);
      statQ3->SetY1NDC(0.45);
      statQ4->SetY2NDC(0.45);
      statQ4->SetY1NDC(0.30);
      statQ5->SetY2NDC(0.30);
      statQ5->SetY1NDC(0.15);
      statQ1->SetTextColor(1);
      statQ2->SetTextColor(2);
      statQ3->SetTextColor(4);
      statQ4->SetTextColor(6);
      statQ5->SetTextColor(8);
      gPad->Modified();
      double Q2_adc = Q2adc->GetMean();
      double Q2_xcut = Q2xcut->GetMean();
      double Q2_pcut = Q2pcut->GetMean();
      double Q2_xylocut = Q2xylocut->GetMean();
      double Q2_xyhicut = Q2xyhicut->GetMean();
      latex.SetTextSize(0.04);
      latex.SetTextColor(2);
      latex.DrawLatex(0.65,0.70,Form("%.2f%s",Q2_xcut/Q2_adc*100.,"%"));
      latex.SetTextColor(4);
      latex.DrawLatex(0.65,0.55,Form("%.2f%s",Q2_pcut/Q2_adc*100.,"%"));
      latex.SetTextColor(6);
      latex.DrawLatex(0.65,0.40,Form("%.2f%s",Q2_xylocut/Q2_adc*100.,"%"));
      latex.SetTextColor(8);
      latex.DrawLatex(0.65,0.25,Form("%.2f%s",Q2_xyhicut/Q2_adc*100.,"%"));
      
      cQsq->SaveAs(Form("./temp1/qsq_analysis_run%d.pdf",run));

      ofstream outfile("./TextFiles/qsquare_diagnostic_quantitiesR.csv",ios_base::app);
      outfile<<run<<"\t"<<target<<"\t"<<date<<"\t"<<upADCcut<<"\t"<<fitEP->GetParameter(1)<<"\t"<<QEdge<<"\t"<<peakpos<<"\t"<<EdgeP<<"\t"<<Q2adc->GetMean()<<"\t"<<Q2adc->GetRMS()<<"\t"<<Q2xcut->GetMean()<<"\t"<<Q2xcut->GetRMS()<<"\t"<<Q2xylocut->GetMean()<<"\t"<<Q2xylocut->GetRMS()<<"\t"<<Q2xyhicut->GetMean()<<"\t"<<Q2xyhicut->GetRMS()<<endl;
      outfile.close();
}
}
