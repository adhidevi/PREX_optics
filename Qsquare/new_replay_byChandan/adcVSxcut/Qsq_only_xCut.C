#include "TChain.h"
#include "TMath.h"

void Qsq_only_xCut(int run, double Ebeam, TString target, TString date){

   TString path = "/chafs1/work1/prex_counting/RootFilesQsq";
 
   TChain *T = new TChain("T");

   if(run<10000){

  Float_t theta0 = 4.806;//degrees
  theta0 = theta0*3.1415926/180;//radians

  Float_t s0 = TMath::Sin(theta0);
  Float_t c0 = TMath::Cos(theta0);

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
  TCut cut_wadc = Form("L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&fEvtHdr.fEvtType==1&&P.upQadcL> %f",upadc_cutL);
  TCut cut_wped = Form("L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&fEvtHdr.fEvtType==1&&P.upQadcL< %f",upadc_cutL);

  c1->cd(2);
  gPad->SetPad(0.4,0.5,1.0,1.0);
//  gPad->SetBottomMargin(0);
  gPad->Modified();gPad->Update();
  gStyle->SetOptStat(1);
  gPad->SetGridx();
  TH1F* hx1 = new TH1F("hx1", Form("LHRS projected x run%d",run), 200, -0.14, 0.03);
  TH1F* hxup2 = new TH1F("hxup2", "LHRS projected x + adc>cut", 200, -0.14, 0.03);
  TH1F* hxup3 = new TH1F("hxup3", "LHRS projected x + adc<cut", 200, -0.14, 0.03);
  T->Project(hx1->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cut);
  T->Project(hxup2->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cut_wadc);
  T->Project(hxup3->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cut_wped);
  hx1->SetXTitle("L.tr.x[0]+1.3*L.tr.th[0] (m)");
  hx1->SetLineColor(1);
  hxup2->SetLineColor(2);
  hxup3->SetLineColor(4);
  hx1->Draw();
  hxup2->Draw("sames");
  hxup3->Draw("sames");

  TF1* epfit = new TF1("epfit","gaus",hxup2->GetBinCenter(hxup2->GetMaximumBin())-0.007,hxup2->GetBinCenter(hxup2->GetMaximumBin())+0.007);
  hxup2->Fit(epfit,"R");

  int BinContent[100];
  int EdgeBin;
  int diff[100];
  cout<<"bin"<<"\t"<<"abs(diff)"<<endl;
  for(int i=1;i<=100;i++){
  diff[i] = abs(hxup2->GetBinContent(i)-hxup3->GetBinContent(i));
  cout<<i<<"\t"<<diff[i]<<endl;
  }
  int mindiff = {diff[1]};
  for(int i=1;i<=100;i++){
  if(mindiff>diff[i]){
  mindiff = diff[i];
  EdgeBin = i;
  }
  }
  double EdgeX = hx1->GetBinCenter(EdgeBin);
  double el_peak = epfit->GetParameter(1);
  cout<<"Edge Bin: "<<EdgeBin<<"\t Edge x: "<<EdgeX<<endl;

  latex.SetTextColor(6);
  latex.DrawLatex(0.13,0.85,Form("x@edge = %1.4f m",EdgeX));
  latex.SetTextColor(1);
  latex.DrawLatex(0.15,0.80,"All events");
  latex.SetTextColor(2);
  latex.DrawLatex(0.15,0.75,"US accepted");
  latex.SetTextColor(4);
  latex.DrawLatex(0.15,0.70,"US missed");
  latex.SetTextColor(46);
  latex.DrawLatex(0.13,0.65,Form("el_peak = %1.4f m",epfit->GetParameter(1)));
  gPad->Update();
  TPaveStats* stat1 = (TPaveStats*)hx1->FindObject("stats");
  TPaveStats* stat2 = (TPaveStats*)hxup2->FindObject("stats");
  TPaveStats* stat3 = (TPaveStats*)hxup3->FindObject("stats");
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
  TLine* xEdge = new TLine(EdgeX,0.0,EdgeX,hx1->GetMaximum()/2.2);
  xEdge->SetLineColor(6);
  xEdge->SetLineWidth(2);
  xEdge->Draw();
  TLine* xpeak = new TLine(el_peak,0.0,el_peak,hx1->GetMaximum());
  xpeak->SetLineColor(46);
  xpeak->SetLineWidth(2);
  xpeak->Draw();

  TCut xCut = Form("L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&fEvtHdr.fEvtType==1&&L.tr.x[0]+1.3*L.tr.th[0]> %f",EdgeX);

  c1->cd(3);
  gPad->SetPad(0.0,0.0,0.4,0.5);
//  gPad->SetRightMargin(0);
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
  TF1* momfit = new TF1("momfit","gaus",peak_approx-0.0003,peak_approx+0.0003);
  hp2->Fit(momfit,"R");

  int BinContentP[100];
  int EdgeBinP;
  int diffP[100];
  cout<<"bin"<<"\t"<<"abs(diff)"<<endl;
  for(int i=1;i<=100;i++){
  diffP[i] = abs(hp2->GetBinContent(i)-hp3->GetBinContent(i));
  cout<<i<<"\t"<<diffP[i]<<endl;
  }
  int mindiffP = {diffP[1]};
  for(int i=1;i<=100;i++){
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
  gPad->Update();
  TPaveStats* statp1 = (TPaveStats*)hp1->FindObject("stats");
  TPaveStats* statp2 = (TPaveStats*)hp2->FindObject("stats");
  TPaveStats* statp3 = (TPaveStats*)hp3->FindObject("stats");
  statp1->SetY2NDC(0.90);
  statp1->SetY1NDC(0.75);
  statp2->SetY2NDC(0.75);
  statp2->SetY1NDC(0.60);
  statp3->SetY2NDC(0.60);
  statp3->SetY1NDC(0.45);
  statp1->SetTextColor(1);
  statp2->SetTextColor(2);
  statp3->SetTextColor(4);
  gPad->Modified();
  TLine* pEdge = new TLine(EdgeP,0.0,EdgeP,hp1->GetMaximum()/2.2);
  pEdge->SetLineColor(6);
  pEdge->SetLineWidth(2);
  pEdge->Draw();
  TLine* ppeak = new TLine(peakpos,0.0,peakpos,hp1->GetMaximum());
  ppeak->SetLineColor(46);
  ppeak->SetLineWidth(2);
  ppeak->Draw();
  TCut pCut = Form("L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&L.gold.dp>-0.04&&L.gold.dp<0.002&&fEvtHdr.fEvtType==1&&L.gold.p> %f",EdgeP);

  c1->cd(4);
  gPad->SetPad(0.4,0.0,1.0,0.5);
  gPad->Modified();gPad->Update();
//  gPad->SetTopMargin(0);
  TH1F *Q2adc= new TH1F("Q2adc",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
  TH1F *Q2xcut= new TH1F("Q2xcut",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
  TH1F *Q2pcut= new TH1F("Q2pcut",Form("Q^{2} (GeV/c)^{2} (run%d);Q^(2) (GeV/c)^{2};",run),150,0,0.015);
  Q2adc->SetLineColor(1);
  Q2xcut->SetLineColor(2);
  Q2pcut->SetLineColor(4);
  T->Draw(Form("2*%f*%f*L.gold.p*(1-((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))))>>Q2adc",Ebeam,efact,c0,s0),cut_wadc,"hist");
  T->Draw(Form("2*%f*%f*L.gold.p*(1-((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))))>>Q2xcut",Ebeam,efact,c0,s0),xCut,"hist sames");
  T->Draw(Form("2*%f*%f*L.gold.p*(1-((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))))>>Q2pcut",Ebeam,efact,c0,s0),pCut,"hist sames");
  Q2xcut->Scale(Q2adc->Integral()/Q2xcut->Integral());
  Q2pcut->Scale(Q2adc->Integral()/Q2pcut->Integral());
  gPad->Update();
  TPaveStats* statQ1 = (TPaveStats*)Q2adc->FindObject("stats");
  TPaveStats* statQ2 = (TPaveStats*)Q2xcut->FindObject("stats");
  TPaveStats* statQ3 = (TPaveStats*)Q2pcut->FindObject("stats");
  statQ1->SetY2NDC(0.90);
  statQ1->SetY1NDC(0.75);
  statQ2->SetY2NDC(0.75);
  statQ2->SetY1NDC(0.60);
  statQ3->SetY2NDC(0.60);
  statQ3->SetY1NDC(0.45);
  statQ1->SetTextColor(1);
  statQ2->SetTextColor(2);
  statQ3->SetTextColor(4);
  gPad->Modified();
  c1->SaveAs(Form("./temp1/qsq_analysis_run%d.pdf",run));

  ofstream outfile("./TextFiles/adc_x_p_cut_cpmpare.csv",ios_base::app);
  outfile<<run<<"\t"<<Ebeam<<"\t"<<target.Data()<<"\t"<<date.Data()<<"\t"<<uqadc_cutL<<"\t"<<EdgeX<<"\t"<<EdgeP<<"\t"<<Q2adc->GetMean()<<"\t"<<Q2adc->GetRMS()<<"\t"<<Q2xcut->GetMean()<<"\t"<<Q2xcut->GetRMS()<<"\t"<<Q2pcut->GetMean()<<"\t"<<Q2pcut->GetRMS()<<endl;
  outfile.close();

 gSystem->Exec(Form("pdfunite ./temp1/plot*.pdf ./plots/Qsq_xCut_run%d.pdf",run));
 gSystem->Exec(Form("rm -rf ./temp1/plot*.pdf"));

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
  TCut cut_wadc = Form("R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<-0.002&&fEvtHdr.fEvtType==1&&P.upQadcR> %f",upadc_cutR);
  TCut cut_wped = Form("R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<-0.002&&fEvtHdr.fEvtType==1&&P.upQadcR< %f",upadc_cutR);

  c1->cd(2);
  gPad->SetPad(0.4,0.5,1.0,1.0);
  gPad->SetBottomMargin(0);
  gPad->Modified();gPad->Update();
  gStyle->SetOptStat(1);
  gPad->SetGridx();
  TH1F* hx1 = new TH1F("hx1", Form("RHRS projected x (run%d)",run), 200, -0.14, 0.03);
  TH1F* hxup2 = new TH1F("hxup2", Form("RHRS projected x (run%d) + adc>cut",run), 200, -0.14, 0.03);
  TH1F* hxup3 = new TH1F("hxup3", Form("RHRS projected x (run%d) + adc<cut",run), 200, -0.14, 0.03);
  T->Project(hx1->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cut);
  T->Project(hxup2->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cut_wadc);
  T->Project(hxup3->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cut_wped);
  hx1->SetXTitle("R.tr.x[0]+1.3*R.tr.th[0] (m)");
  hx1->SetLineColor(1);
  hxup2->SetLineColor(2);
  hxup3->SetLineColor(4);
  hx1->Draw();
  hxup2->Draw("sames");
  hxup3->Draw("sames");

  TF1* epfit = new TF1("epfit","gaus",hxup2->GetBinCenter(hxup2->GetMaximumBin())-0.007,hxup2->GetBinCenter(hxup2->GetMaximumBin())+0.007);
  hxup2->Fit(epfit,"R");

  int BinContent[100];
  int EdgeBin;
  int diff[100];
  cout<<"bin"<<"\t"<<"abs(diff)"<<endl;
  for(int i=1;i<=100;i++){
  diff[i] = abs(hxup2->GetBinContent(i)-hxup3->GetBinContent(i));
  cout<<i<<"\t"<<diff[i]<<endl;
  }
  int mindiff = {diff[1]};
  for(int i=1;i<=100;i++){
  if(mindiff>diff[i]){
  mindiff = diff[i];
  EdgeBin = i;
  }
  }
  double EdgeX = hx1->GetBinCenter(EdgeBin);
  cout<<"Edge Bin: "<<EdgeBin<<"\t Edge x: "<<EdgeX<<endl;

  latex.SetTextColor(6);
  latex.DrawLatex(0.15,0.85,Form("x@edge = %1.4f m",EdgeX));
  latex.SetTextColor(1);
  latex.DrawLatex(0.17,0.80,"All events");
  latex.SetTextColor(2);
  latex.DrawLatex(0.17,0.75,"US accepted");
  latex.SetTextColor(4);
  latex.DrawLatex(0.17,0.70,"US missed");
  latex.SetTextColor(46);
  latex.DrawLatex(0.13,0.65,Form("el_peak = %1.4f m",epfit->GetParameter(1)));
  gPad->Update();
  TPaveStats* stat1 = (TPaveStats*)hx1->FindObject("stats");
  TPaveStats* stat2 = (TPaveStats*)hxup2->FindObject("stats");
  TPaveStats* stat3 = (TPaveStats*)hxup3->FindObject("stats");
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
  TLine* pEdge = new TLine(EdgeX,0.0,EdgeX,hx1->GetMaximum()/2.2);
  pEdge->SetLineColor(6);
  pEdge->SetLineWidth(2);
  pEdge->Draw();

  c1->cd(3);
  gPad->SetPad(0.0,0.0,0.4,0.5);
  gPad->SetRightMargin(0);
  c1->cd(4);
  gPad->SetPad(0.4,0.0,1.0,0.5);
  gPad->SetTopMargin(0);
  TH1F* hxup2Cpy = (TH1F*)hxup2->Clone("hxup2Cpy");
  hxup2Cpy->Divide(hx1);
  hxup2Cpy->SetTitle(";R.tr.x[0]+1.3*R.tr.th[0] (m);");
  hxup2Cpy->SetStats(0);
  hxup2Cpy->SetMarkerColor(1);
  hxup2Cpy->SetLineColor(1);
  hxup2Cpy->Draw("hist");
  pEdge->Draw();
  latex.SetTextColor(1);
  latex.DrawLatex(0.40,0.60,"quartz acceptance (accepted/all)");
  latex.SetTextColor(6);
  latex.DrawLatex(0.40,0.50,Form("acceptance @ q_edge = %.3f",hxup2Cpy->Interpolate(EdgeX)));

  c1->SaveAs("./temp1/plot1.pdf");

 TCut xCut = Form("R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&R.gold.dp>-0.04&&R.gold.dp<-0.002&&fEvtHdr.fEvtType==1&&R.tr.x[0]+1.3*R.tr.th[0]> %f",EdgeX);
 TH1F* Q2 = new TH1F("Q2",Form("Q^{2} (GeV/c)^{2} (run%d);Q^{2} (GeV/c)^{2}",run),150,0,0.015);
 TCanvas *q = new TCanvas("q","q",600,500);
 T->Draw("EK_R.Q2>>Q2",xCut);
 q->SaveAs("./temp1/plot2.pdf");
 
 gSystem->Exec(Form("pdfunite ./temp1/plot*.pdf ./plots/Qsq_xCut_run%d.pdf",run));
 gSystem->Exec(Form("rm -rf ./temp1/plot*.pdf"));
 
}
}
