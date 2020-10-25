#include "TChain.h"
#include "TMath.h"

void Qsq_only_pCut(int run, Float_t E){

   TString path = "/chafs1/work1/prex_counting/QsqRootFiles";
 
   TChain *T = new TChain("T");

   if(run<10000){
   T->Add(Form("%s/prexLHRS_%i_-1*.root",path.Data(),run));
  double upadc_cutL_approx = 480;
  TCanvas *c1 = new TCanvas("c1","c1",900,700);
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetPad(0.0,0.5,0.4,1.0);
  gPad->Modified();gPad->Update();
  gPad->SetGridx();
  gPad->SetGridy();
  TH1F* huq = new TH1F("huq", Form("ADC raw (run%d);ADC raw;Events/CH",run), 300,400,700);
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
  gPad->Update();
  TPaveStats* statADC = (TPaveStats*)huq->GetFunction("stats");
  statADC->SetY2NDC(0.90);
  statADC->SetY1NDC(0.75);
  statADC->SetX2NDC(0.90);
  statADC->SetX1NDC(0.60);
  gPad->Modified();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextColor(1);
  latex.DrawLatex(0.53,0.70,"P.upQadcL");
  latex.SetTextColor(9);
  latex.DrawLatex(0.53,0.65,Form("Ped. cut = %3.1f",upadc_cutL));
  TCut cut = Form("L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&fEvtHdr.fEvtType==1");
  TCut cutADCup = Form("L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&fEvtHdr.fEvtType==1&&P.upQadcL> %f",upadc_cutL);
  TCut cutPEDup = Form("L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&fEvtHdr.fEvtType==1&&P.upQadcL< %f",upadc_cutL);

  c1->cd(2);
  gPad->SetPad(0.4,0.5,1.0,1.0);
  gPad->SetBottomMargin(0);
  gPad->Modified();gPad->Update();
  gStyle->SetOptStat(1);
  gPad->SetGridx();
  TH1F* hp1 = new TH1F("hp1", Form("LHRS momentum run%d",run), 200, 0.944, 0.953);
  TH1F* hpup2 = new TH1F("hpup2", "LHRS momentum + adc>cut", 200, 0.944, 0.953);
  TH1F* hpup3 = new TH1F("hpup3", "LHRS momentum + adc<cut", 200, 0.944, 0.953);
  T->Project(hp1->GetName(), "L.gold.p", cut);
  T->Project(hpup2->GetName(), "L.gold.p", cutADCup);
  T->Project(hpup3->GetName(), "L.gold.p", cutPEDup);
  hp1->SetXTitle("L.gold.p (GeV)");
  hp1->SetLineColor(1);
  hpup2->SetLineColor(2);
  hpup3->SetLineColor(4);
  hp1->Draw();
  hpup2->Draw("sames");
  hpup3->Draw("sames");

  TF1* epfit = new TF1("epfit","gaus",hpup2->GetBinCenter(hpup2->GetMaximumBin())-0.0003,hpup2->GetBinCenter(hpup2->GetMaximumBin())+0.0003);
  hpup2->Fit(epfit,"R");

  int BinContent[100];
  int EdgeBin;
  int diff[100];
  cout<<"bin"<<"\t"<<"abs(diff)"<<endl;
  for(int i=1;i<=100;i++){
  diff[i] = abs(hpup2->GetBinContent(i)-hpup3->GetBinContent(i));
  cout<<i<<"\t"<<diff[i]<<endl;
  }
  int mindiff = {diff[1]};
  for(int i=1;i<=100;i++){
  if(mindiff>diff[i]){
  mindiff = diff[i];
  EdgeBin = i;
  }
  }
  double EdgeP = hp1->GetBinCenter(EdgeBin);
  cout<<"Edge Bin: "<<EdgeBin<<"\t Edge p: "<<EdgeP<<endl;

  int nCut = 7;
  double pCut[nCut];
  for(int iCut=-3;iCut<=3;iCut++){
  pCut[iCut+3] = EdgeP+(5*iCut*1e-4);
  cout<<"pCut: "<<pCut[iCut+3]<<endl;
  }

  latex.SetTextColor(6);
  latex.DrawLatex(0.13,0.85,Form("p@edge = %1.4f GeV",EdgeP));
  latex.SetTextColor(1);
  latex.DrawLatex(0.15,0.80,"All events");
  latex.SetTextColor(2);
  latex.DrawLatex(0.15,0.75,"US accepted");
  latex.SetTextColor(4);
  latex.DrawLatex(0.15,0.70,"US missed");
  latex.SetTextColor(46);
  latex.DrawLatex(0.13,0.65,Form("el_peak = %1.4f GeV",epfit->GetParameter(1)));
  gPad->Update();
  TPaveStats* stat1 = (TPaveStats*)hp1->FindObject("stats");
  TPaveStats* stat2 = (TPaveStats*)hpup2->FindObject("stats");
  TPaveStats* stat3 = (TPaveStats*)hpup3->FindObject("stats");
  stat1->SetY2NDC(0.90);
  stat1->SetY1NDC(0.75);
  stat2->SetY2NDC(0.75);
  stat2->SetY1NDC(0.60);
  stat3->SetY2NDC(0.60);
  stat3->SetY1NDC(0.45);
  stat1->SetTextColor(1);
  stat2->SetTextColor(2);
  stat3->SetTextColor(4);
  gPad->Modified();
  TLine* pEdge = new TLine(EdgeP,0.0,EdgeP,hp1->GetMaximum()/2.2);
  pEdge->SetLineColor(6);
  pEdge->SetLineWidth(2);
  pEdge->Draw();

  c1->cd(3);
  gPad->SetPad(0.0,0.0,0.4,0.5);
  gPad->SetRightMargin(0);
  c1->cd(4);
  gPad->SetPad(0.4,0.0,1.0,0.5);
  gPad->SetTopMargin(0);
  TH1F* hpup2Cpy = (TH1F*)hpup2->Clone("hpup2Cpy");
  hpup2Cpy->Divide(hp1);
  hpup2Cpy->SetTitle(";L.gold.p (GeV);");
  hpup2Cpy->SetStats(0);
  hpup2Cpy->SetMarkerColor(1);
  hpup2Cpy->SetLineColor(1);
  hpup2Cpy->Draw("hist");
  pEdge->Draw();
  latex.SetTextColor(1);
  latex.DrawLatex(0.40,0.60,"quartz acceptance (accepted/all)");
  latex.SetTextColor(6);
  latex.DrawLatex(0.40,0.50,Form("acceptance @ q_edge = %.3f",hpup2Cpy->Interpolate(EdgeP)));

  c1->SaveAs("./temp1/plot1.pdf");

  TH1F *Q2[nCut];

 for(int iCut=0;iCut<nCut;iCut++){
  TCut momCut = Form("L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&fEvtHdr.fEvtType==1&&L.gold.p> %f",pCut[iCut]);
 Q2[iCut] = new TH1F(Form("Q2_%.4f",pCut[iCut]),Form("Q^{2} (GeV/c)^{2} with pCut (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
 T->Draw(Form("EK_L.Q2>>Q2_%.4f",pCut[iCut]),momCut,"goff");
 if(iCut==4)
 Q2[iCut]->SetLineColor(42);
 else if(iCut==6)
 Q2[iCut]->SetLineColor(9);
 else
 Q2[iCut]->SetLineColor(iCut+1);
 }
 TCanvas *q = new TCanvas("q","q",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 Q2[iCut]->Draw("sames hist");
 Q2[iCut]->Scale(1.0/Q2[iCut]->GetMaximum());
 }
 gPad->Update();
 TPaveStats* statsQ[nCut];
 for(int iCut=0;iCut<nCut;iCut++){
 statsQ[iCut] = (TPaveStats*)Q2[iCut]->FindObject("stats");
 }
 statsQ[0]->SetY2NDC(0.90);
 statsQ[0]->SetY1NDC(0.75);
 statsQ[1]->SetY2NDC(0.75);
 statsQ[1]->SetY1NDC(0.60);
 statsQ[2]->SetY2NDC(0.60);
 statsQ[2]->SetY1NDC(0.45);
 statsQ[3]->SetY2NDC(0.45);
 statsQ[3]->SetY1NDC(0.30);

 statsQ[0]->SetX2NDC(0.90);
 statsQ[0]->SetX1NDC(0.70);
 statsQ[1]->SetX2NDC(0.90);
 statsQ[1]->SetX1NDC(0.70);
 statsQ[2]->SetX2NDC(0.90);
 statsQ[2]->SetX1NDC(0.70);
 statsQ[3]->SetX2NDC(0.90);
 statsQ[3]->SetX1NDC(0.70);

 statsQ[4]->SetY2NDC(0.90);
 statsQ[4]->SetY1NDC(0.75);
 statsQ[5]->SetY2NDC(0.75);
 statsQ[5]->SetY1NDC(0.60);
 statsQ[6]->SetY2NDC(0.60);
 statsQ[6]->SetY1NDC(0.45);

 statsQ[4]->SetX2NDC(0.33);
 statsQ[4]->SetX1NDC(0.13);
 statsQ[5]->SetX2NDC(0.33);
 statsQ[5]->SetX1NDC(0.13);
 statsQ[6]->SetX2NDC(0.33);
 statsQ[6]->SetX1NDC(0.13);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 statsQ[iCut]->SetTextColor(42);
 else if(iCut==6)
 statsQ[iCut]->SetTextColor(9);
 else
 statsQ[iCut]->SetTextColor(iCut+1);
 }
 gPad->Modified();
 q->SaveAs("./temp1/plot2.pdf");

 gSystem->Exec(Form("pdfunite ./temp1/plot*.pdf ./plots/Qsq_pCut_run%d.pdf",run));
 gSystem->Exec(Form("rm -rf ./temp1/plot*.pdf"));

 gSystem->Exec(Form("rm -rf ./TextFiles/Qsq_pCut_run%d.csv",run));

 ofstream outfile(Form("./TextFiles/Qsq_pCut_run%d.csv",run));
 ofstream edgeCut("./TextFiles/Qsq_edgeCut_allRunsL.csv",ios_base::app);

 for(int iCut=0;iCut<nCut;iCut++){
 if(EdgeP==pCut[iCut])
 edgeCut<<EdgeP<<"\t"<<pCut[iCut]<<"\t"<<Q2[iCut]->GetMean()<<"\t"<<Q2[iCut]->GetRMS()<<"\t"<<epfit->GetParameter(1)-pCut[iCut]<<"\t"<<hpup2Cpy->Interpolate(pCut[iCut])<<endl;
 else
 outfile<<EdgeP<<"\t"<<pCut[iCut]<<"\t"<<Q2[iCut]->GetMean()<<"\t"<<Q2[iCut]->GetRMS()<<"\t"<<epfit->GetParameter(1)-pCut[iCut]<<"\t"<<hpup2Cpy->Interpolate(pCut[iCut])<<endl;
}
 outfile.close();
}else{
   T->Add(Form("%s/prexRHRS_%i_-1*.root",path.Data(),run));
  double upadc_cutR_approx = 502;
  TCanvas *c1 = new TCanvas("c1","c1",900,700);
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetPad(0.0,0.5,0.4,1.0);
  gPad->Modified();gPad->Update();
  gPad->SetGridx();
  gPad->SetGridy();
  TH1F* huq = new TH1F("huq", Form("ADC raw (run%d);ADC raw;Events/CH",run), 200,450,650);
  huq->SetLineColor(1);
  T->Project(huq->GetName(),"P.upQadcR");
  huq->Draw();

  TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
  huqCpy->GetXaxis()->SetRangeUser(upadc_cutR_approx-10,upadc_cutR_approx+10);
  double upadc_cutR = huqCpy->GetBinCenter(huqCpy->GetMinimumBin());

  TLine* ped_line = new TLine(upadc_cutR,0.0,upadc_cutR,huq->GetMaximum());
  ped_line->SetLineColor(9);
  ped_line->SetLineWidth(2);
  ped_line->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextColor(1);
  latex.DrawLatex(0.53,0.70,"P.upQadcR");
  latex.SetTextColor(9);
  latex.DrawLatex(0.53,0.65,Form("Ped. cut = %3.1f",upadc_cutR));
  TCut cut = Form("R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<-0.002&&fEvtHdr.fEvtType==1");
  TCut cutADCup = Form("R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<-0.002&&fEvtHdr.fEvtType==1&&P.upQadcR> %f",upadc_cutR);
  TCut cutPEDup = Form("R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<-0.002&&fEvtHdr.fEvtType==1&&P.upQadcR< %f",upadc_cutR);

  c1->cd(2);
  gPad->SetPad(0.4,0.5,1.0,1.0);
  gPad->SetBottomMargin(0);
  gPad->Modified();gPad->Update();
  gStyle->SetOptStat(1);
  gPad->SetGridx();
  TH1F* hp1 = new TH1F("hp1", Form("RHRS momentum (run%d)",run), 200, 0.942, 0.951);
  TH1F* hpup2 = new TH1F("hpup2", Form("RHRS momentum (run%d) + adc>cut",run), 200, 0.942, 0.951);
  TH1F* hpup3 = new TH1F("hpup3", Form("RHRS momentum (run%d) + adc<cut",run), 200, 0.942, 0.951);
  T->Project(hp1->GetName(), "R.gold.p", cut);
  T->Project(hpup2->GetName(), "R.gold.p", cutADCup);
  T->Project(hpup3->GetName(), "R.gold.p", cutPEDup);
  hp1->SetXTitle("R.gold.p (GeV)");
  hp1->SetLineColor(1);
  hpup2->SetLineColor(2);
  hpup3->SetLineColor(4);
  hp1->Draw();
  hpup2->Draw("sames");
  hpup3->Draw("sames");

  TF1* epfit = new TF1("epfit","gaus",hpup2->GetBinCenter(hpup2->GetMaximumBin())-0.0003,hpup2->GetBinCenter(hpup2->GetMaximumBin())+0.0003);
  hpup2->Fit(epfit,"R");

  int BinContent[100];
  int EdgeBin;
  int diff[100];
  cout<<"bin"<<"\t"<<"abs(diff)"<<endl;
  for(int i=1;i<=100;i++){
  diff[i] = abs(hpup2->GetBinContent(i)-hpup3->GetBinContent(i));
  cout<<i<<"\t"<<diff[i]<<endl;
  }
  int mindiff = {diff[1]};
  for(int i=1;i<=100;i++){
  if(mindiff>diff[i]){
  mindiff = diff[i];
  EdgeBin = i;
  }
  }
  double EdgeP = hp1->GetBinCenter(EdgeBin);
  cout<<"Edge Bin: "<<EdgeBin<<"\t Edge p: "<<EdgeP<<endl;

  int nCut = 7;
  double pCut[nCut];
  for(int iCut=-3;iCut<=3;iCut++){
  pCut[iCut+3] = EdgeP+(5*iCut*1e-4);
  cout<<"pCut: "<<pCut[iCut+3]<<endl;
  }

  latex.SetTextColor(6);
  latex.DrawLatex(0.15,0.85,Form("p@edge = %1.4f GeV",EdgeP));
  latex.SetTextColor(1);
  latex.DrawLatex(0.17,0.80,"All events");
  latex.SetTextColor(2);
  latex.DrawLatex(0.17,0.75,"US accepted");
  latex.SetTextColor(4);
  latex.DrawLatex(0.17,0.70,"US missed");
  latex.SetTextColor(46);
  latex.DrawLatex(0.13,0.65,Form("el_peak = %1.4f GeV",epfit->GetParameter(1)));
  gPad->Update();
  TPaveStats* stat1 = (TPaveStats*)hp1->FindObject("stats");
  TPaveStats* stat2 = (TPaveStats*)hpup2->FindObject("stats");
  TPaveStats* stat3 = (TPaveStats*)hpup3->FindObject("stats");
  stat1->SetY2NDC(0.90);
  stat1->SetY1NDC(0.75);
  stat2->SetY2NDC(0.75);
  stat2->SetY1NDC(0.60);
  stat3->SetY2NDC(0.60);
  stat3->SetY1NDC(0.45);
  stat1->SetTextColor(1);
  stat2->SetTextColor(2);
  stat3->SetTextColor(4);
  gPad->Modified();
  TLine* pEdge = new TLine(EdgeP,0.0,EdgeP,hp1->GetMaximum()/2.2);
  pEdge->SetLineColor(6);
  pEdge->SetLineWidth(2);
  pEdge->Draw();

  c1->cd(3);
  gPad->SetPad(0.0,0.0,0.4,0.5);
  gPad->SetRightMargin(0);
  c1->cd(4);
  gPad->SetPad(0.4,0.0,1.0,0.5);
  gPad->SetTopMargin(0);
  TH1F* hpup2Cpy = (TH1F*)hpup2->Clone("hpup2Cpy");
  hpup2Cpy->Divide(hp1);
  hpup2Cpy->SetTitle(";R.gold.p (GeV);");
  hpup2Cpy->SetStats(0);
  hpup2Cpy->SetMarkerColor(1);
  hpup2Cpy->SetLineColor(1);
  hpup2Cpy->Draw("hist");
  pEdge->Draw();
  latex.SetTextColor(1);
  latex.DrawLatex(0.40,0.60,"quartz acceptance (accepted/all)");
  latex.SetTextColor(6);
  latex.DrawLatex(0.40,0.50,Form("acceptance @ q_edge = %.3f",hpup2Cpy->Interpolate(EdgeP)));

  c1->SaveAs("./temp1/plot1.pdf");

 TH1F* Q2[nCut];

 for(int iCut=0;iCut<nCut;iCut++){
 TCut momCut = Form("R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<-0.002&&fEvtHdr.fEvtType==1&&R.gold.p> %f",pCut[iCut]);
 Q2[iCut] = new TH1F(Form("Q2_%.4f",pCut[iCut]),Form("Q^{2} (GeV/c)^{2} with pCut (run%d);Q^{2} (GeV/c)^{2}",run),150,0,0.015);
 T->Draw(Form("EK_R.Q2>>Q2_%.4f",pCut[iCut]),momCut,"goff");
 if(iCut==4)
 Q2[iCut]->SetLineColor(42);
 else if(iCut==6)
 Q2[iCut]->SetLineColor(9);
 else
 Q2[iCut]->SetLineColor(iCut+1);
 }
 TCanvas *q = new TCanvas("q","q",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 Q2[iCut]->Draw("sames hist");
 Q2[iCut]->Scale(1.0/Q2[iCut]->GetMaximum());
 }
 gPad->Update();
 TPaveStats* statsQ[nCut];
 for(int iCut=0;iCut<nCut;iCut++){
 statsQ[iCut] = (TPaveStats*)Q2[iCut]->FindObject("stats");
 }
 statsQ[0]->SetY2NDC(0.90);
 statsQ[0]->SetY2NDC(0.90);
 statsQ[0]->SetY1NDC(0.75);
 statsQ[1]->SetY2NDC(0.75);
 statsQ[1]->SetY1NDC(0.60);
 statsQ[2]->SetY2NDC(0.60);
 statsQ[2]->SetY1NDC(0.45);
 statsQ[3]->SetY2NDC(0.45);
 statsQ[3]->SetY1NDC(0.30);

 statsQ[0]->SetX2NDC(0.90);
 statsQ[0]->SetX1NDC(0.70);
 statsQ[1]->SetX2NDC(0.90);
 statsQ[1]->SetX1NDC(0.70);
 statsQ[2]->SetX2NDC(0.90);
 statsQ[2]->SetX1NDC(0.70);
 statsQ[3]->SetX2NDC(0.90);
 statsQ[3]->SetX1NDC(0.70);

 statsQ[4]->SetY2NDC(0.90);
 statsQ[4]->SetY1NDC(0.75);
 statsQ[5]->SetY2NDC(0.75);
 statsQ[5]->SetY1NDC(0.60);
 statsQ[6]->SetY2NDC(0.60);
 statsQ[6]->SetY1NDC(0.45);

 statsQ[4]->SetX2NDC(0.35);
 statsQ[4]->SetX1NDC(0.15);
 statsQ[5]->SetX2NDC(0.35);
 statsQ[5]->SetX1NDC(0.15);
 statsQ[6]->SetX2NDC(0.35);
 statsQ[6]->SetX1NDC(0.15);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 statsQ[iCut]->SetTextColor(42);
 else if(iCut==6)
 statsQ[iCut]->SetTextColor(9);
 else
 statsQ[iCut]->SetTextColor(iCut+1);
 }
 gPad->Modified();
 q->SaveAs("./temp1/plot2.pdf");
 
 gSystem->Exec(Form("pdfunite ./temp1/plot*.pdf ./plots/Qsq_pCut_run%d.pdf",run));
 gSystem->Exec(Form("rm -rf ./temp1/plot*.pdf"));
 
 gSystem->Exec(Form("rm -rf ./TextFiles/Qsq_pCut_run%d.csv",run));
 ofstream outfile(Form("./TextFiles/Qsq_pCut_run%d.csv",run));
 ofstream edgeCut("./TextFiles/Qsq_edgeCut_allRunsR.csv",ios_base::app);

 for(int iCut=0;iCut<nCut;iCut++){
 if(EdgeP==pCut[iCut])
 edgeCut<<EdgeP<<"\t"<<pCut[iCut]<<"\t"<<Q2[iCut]->GetMean()<<"\t"<<Q2[iCut]->GetRMS()<<"\t"<<epfit->GetParameter(1)-pCut[iCut]<<"\t"<<hpup2Cpy->Interpolate(pCut[iCut])<<endl;
 else
 outfile<<EdgeP<<"\t"<<pCut[iCut]<<"\t"<<Q2[iCut]->GetMean()<<"\t"<<Q2[iCut]->GetRMS()<<"\t"<<epfit->GetParameter(1)-pCut[iCut]<<"\t"<<hpup2Cpy->Interpolate(pCut[iCut])<<endl;
}
 outfile.close();
 edgeCut.close();
}
}
