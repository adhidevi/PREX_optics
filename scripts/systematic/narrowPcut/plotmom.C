#include "CREXdata.h"
#include "TChain.h"
#include "TMath.h"

void plotmom(int run, TString date, TString target){

  //Loads the Horowitz tables
  LoadTable("horpb.dat",0); //unstretched 
  LoadTable("horpb1.dat",1); //stretched

  //For the interpolation
  int j = 0;
  
  //ROOT 6 memory management makes me use vectors
  vector <double> Qsq, Angle, ASYM, ASYM_ST, Sens; 
  
  TString path = "/chafs1/work1/prex_counting/chandan";
  TChain *T = new TChain("T");
     
  if(run<10000){
  T->Add(Form("%s/prexLHRS_%i_-1*.root",path.Data(),run));
  double upadc_cutL_approx = 480;
  TH1F* huq = new TH1F(Form("huq"), Form("ADC raw (run%d), %s, %s;ADC raw;Events/CH",run,date.Data(),target.Data()), 300,400,700);
  huq->SetLineColor(1);
  T->Project(huq->GetName(),"P.upQadcL","","goff");
  huq->Draw();

  TH1F* huqCpy = (TH1F*)huq->Clone(Form("huqCpy"));
  huqCpy->GetXaxis()->SetRangeUser(upadc_cutL_approx-10,upadc_cutL_approx+10);
  double upadc_cutL = huqCpy->GetBinCenter(huqCpy->GetMinimumBin());

  TCut cut = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)");
  TCut cutADCup = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)&&P.upQadcL> %f",upadc_cutL);
  TCut cutPEDup = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)&&P.upQadcL< %f",upadc_cutL);

  TH1F* hp1 = new TH1F(Form("hp1"), Form("LHRS momentum (run%d)",run), 200, 0.944, 0.953);
  TH1F* hpup2 = new TH1F(Form("hpup2"), "LHRS momentum + adc>cut", 200, 0.944, 0.953);
  TH1F* hpup3 = new TH1F(Form("hpup3"), "LHRS momentum + adc<cut", 200, 0.944, 0.953);
  T->Project(hp1->GetName(), "L.gold.p", cut);
  T->Project(hpup2->GetName(), "L.gold.p", cutADCup);
  T->Project(hpup3->GetName(), "L.gold.p", cutPEDup);
  hp1->SetXTitle("L.gold.p (GeV)");
  hp1->SetLineColor(1);
  hpup2->SetLineColor(2);
  hpup3->SetLineColor(4);
  TCanvas* c1 = new TCanvas(Form("c1"),"c1",900,700);
  c1->Divide(2,2);

  c1->cd(1);
  gPad->SetPad(0.0,0.5,0.4,1.0);
  gPad->Modified();gPad->Update();
  gPad->SetGridx();
  gPad->SetGridy();
  huq->Draw();
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextColor(1);
  latex.DrawLatex(0.53,0.70,"P.upQadcL");
  latex.SetTextColor(9);
  latex.DrawLatex(0.53,0.65,Form("Ped Cut = %.1f",upadc_cutL));
  TLine* pedL = new TLine(upadc_cutL,0.0,upadc_cutL,huq->GetMaximum());
  pedL->SetLineColor(9);
  pedL->SetLineWidth(2);
  pedL->Draw();

  c1->cd(2);
  gPad->SetPad(0.4,0.5,1.0,1.0);
  gPad->SetBottomMargin(0);
  gPad->Modified();gPad->Update();
  gStyle->SetOptStat(1);
  gPad->SetGridx();
  hp1->Draw();
  hpup2->Draw("sames hist");
  hpup3->Draw("sames hist");
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

  TLine* q_edge = new TLine(EdgeP,0.0,EdgeP,hp1->GetMaximum()/2.0);
  q_edge->SetLineColor(6);
  q_edge->SetLineWidth(2);
  q_edge->Draw();

  latex.SetTextColor(6);
  latex.DrawLatex(0.15,0.85,Form("q_edge: %.4f GeV",EdgeP));
  latex.SetTextColor(1);
  latex.DrawLatex(0.15,0.80,"All events");
  latex.SetTextColor(2);
  latex.DrawLatex(0.15,0.75,"US accepted");
  latex.SetTextColor(4);
  latex.DrawLatex(0.15,0.70,"US missed");
  c1->cd(3);
  gPad->SetPad(0.0,0.0,0.4,0.5);
  gPad->SetRightMargin(0);
  c1->cd(4);
  gPad->SetLogy();
  gPad->SetPad(0.4,0.0,1.0,0.5);
  gPad->SetTopMargin(0);
  TH1F* hpup2Cpy = (TH1F*)hpup2->Clone("hpup2Cpy");
  hpup2Cpy->Divide(hp1);
  hpup2Cpy->SetTitle(";L.gold.p (GeV);");
  hpup2Cpy->SetStats(0);
  hpup2Cpy->SetMarkerColor(1);
  hpup2Cpy->SetLineColor(1);
  hpup2Cpy->Draw("hist");
  q_edge->Draw();
  latex.SetTextColor(1);
  latex.DrawLatex(0.40,0.60,"quartz acceptance (accepted/all)");
  double accept = hpup2Cpy->Interpolate(EdgeP);
  latex.SetTextColor(6);
  latex.DrawLatex(0.40,0.50,Form("acceptance @ q_edge = %.3f",accept));

  c1->SaveAs(Form("./temp/plot1.pdf"));

  double thisThtg, thisPhtg, thisP, thisDet;
  double thisu1, thisu2, thisv1, thisv2;
  double thisXvdc[10], thisThvdc[10], thisYvdc[10], thisPhvdc[10];
  double thisAsym, thisAsymSt;
  T->SetMakeClass(1);//This is required to read the sub-branches
  Int_t EvtHdr, Event_Branch;
  Int_t thisTrig, fEvtType;
 
  double xq,yq;//For quartz projection
  double thisAngle; //for reconstructed angle
  double thisCosAng; // for cos angle  
  double thisQsq;//for Qsq
  double E = 0.95;//beam energy

  //Central Angle
  double th0 = 4.74*TMath::Pi()/180;
  double cth0 = TMath::Cos(th0);
  double sth0 = TMath::Sin(th0);

  T->SetBranchAddress("L.gold.th",&thisThtg);
  T->SetBranchAddress("L.gold.ph",&thisPhtg);
  T->SetBranchAddress("L.gold.p",&thisP);
  T->SetBranchAddress("P.upQadcL",&thisDet);
  T->SetBranchAddress("L.vdc.u1.nclust",&thisu1);
  T->SetBranchAddress("L.vdc.v1.nclust",&thisv1);
  T->SetBranchAddress("L.vdc.u2.nclust",&thisu2);
  T->SetBranchAddress("L.vdc.v2.nclust",&thisv2);
  T->SetBranchAddress("Event_Branch",&Event_Branch);
  T->SetBranchAddress("fEvtHdr",&EvtHdr);
  T->SetBranchAddress("fEvtHdr.fEvtType",&thisTrig);
  
  // --- This segfaults when I loop it for some reason
  T->SetBranchAddress("L.tr.x",&thisXvdc[0]);
  T->SetBranchAddress("L.tr.y",&thisYvdc[0]);
  T->SetBranchAddress("L.tr.th",&thisThvdc[0]);
  T->SetBranchAddress("L.tr.ph",&thisPhvdc[0]);
 
  //Asymmetry table is in degrees so we need to convert from radians to degrees
  double radtodeg = 180/TMath::Pi();

  long n = T->GetEntries();

  for(long i = 0; i < n; i++){ 

  T->GetEntry(i);
  if( thisu1==1 && thisv1==1 && thisu2==1 && thisv2==1 && thisP > EdgeP && thisTrig==1){
  thisCosAng = (cth0 - thisPhtg*sth0)/TMath::Sqrt(1+thisThtg*thisThtg+thisPhtg*thisPhtg);
  thisAngle = radtodeg*TMath::ACos(thisCosAng);
  thisQsq = 2*E*thisP*(1-thisCosAng);
  Qsq.push_back(thisQsq); 
  Angle.push_back(thisAngle);
  //Energy table in MeV so need to be careful here
  thisAsym = 1e6*Interpolate(thisP*1000,thisAngle,0,1);//A (ppm)
  thisAsymSt = 1e6*Interpolate(thisP*1000,thisAngle,1,1); // Stretched A (ppm); 

  ASYM.push_back(thisAsym);ASYM_ST.push_back(thisAsymSt);//Stretched A (ppm)
  Sens.push_back(fabs(thisAsym-thisAsymSt)/thisAsym);//Sensitivity - Why not?
  }
 }

 //Histograms 
 TH1F *Q2 = new TH1F("Q2","Q^{2} (GeV/c)^{2}",150,0,0.015);
 TH1F *Theta = new TH1F("Theta","#theta_{lab} (deg)",150,2,8);
 TH1F *Asym = new TH1F("Asym","Asymmetry (ppm)",150,0,1);
 TH1F *SAsym = new TH1F("SAsym","Stretched Asymmetry (ppm)",150,0,1);
 TH1F *sens = new TH1F("sens","Sensitivity",150,0,0.05);

 Q2->SetTitle(Form("Q^{2} (GeV/c)^{2}, EdgeP = %1.4f GeV",E));
 Theta->SetTitle(Form("#theta_{lab} (deg), EdgeP = %1.4f GeV",EdgeP));
 Asym->SetTitle(Form("Asymmetry (ppm), EdgeP = %1.4f GeV",EdgeP));
 SAsym->SetTitle(Form("Stretched Asym. (ppm),  EdgeP = %1.4f GeV",EdgeP));
 sens->SetTitle(Form("Sensitivity, EdgeP = %1.4f GeV",EdgeP));

 Q2->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
 Theta->GetXaxis()->SetTitle("#theta_{lab} (deg)");
 Asym->GetXaxis()->SetTitle("Asymmetry (ppm)");
 SAsym->GetXaxis()->SetTitle("Stretched Asym. (ppm)");

 for(int k = 0; k < Qsq.size(); k++){

 Q2->Fill(Qsq[k]);
 Theta->Fill(Angle[k]);
 Asym->Fill(ASYM[k]);
 SAsym->Fill(ASYM_ST[k]);
 sens->Fill(Sens[k]);
 }

 TCanvas *t = new TCanvas("t","t",600,500);
 Theta->Draw();
 t->SaveAs("./temp/plot2.pdf");

 TCanvas *au = new TCanvas("au","au",600,500);
 Asym->Draw();
 au->SaveAs("./temp/plot3.pdf");
 
 TCanvas *as = new TCanvas("as","as",600,500);
 SAsym->Draw();
 as->SaveAs("./temp/plot4.pdf");

 TCanvas *q = new TCanvas("q","q",600,500);
 Q2->Draw();
 q->SaveAs("./temp/plot5.pdf");

 TCanvas *s = new TCanvas("s","s",600,500);
 sens->Draw();
 s->SaveAs("./temp/plot6.pdf");
 
 gSystem->Exec(Form("pdfunite ./temp/plot*.pdf ./temp/all_EdgeP_L_run%d.pdf",run));
 gSystem->Exec(Form("rm -rf ./temp/plot*.pdf"));

  ofstream outfile("./TextFiles/ped_EdgeP_L.csv",ios_base::app);
  outfile<<run<<"\t"<<upadc_cutL<<"\t"<<EdgeP<<"\t"<<epfit->GetParameter(1)<<"\t"<<Theta->GetMean()<<"\t"<<Theta->GetRMS()<<"\t"<<Asym->GetMean()<<"\t"<<Asym->GetRMS()<<"\t"<<SAsym->GetMean()<<"\t"<<SAsym->GetRMS()<<"\t"<<Q2->GetMean()<<"\t"<<Q2->GetRMS()<<"\t"<<sens->GetMean()<<"\t"<<sens->GetRMS()<<"\t"<<epfit->GetParameter(1)-EdgeP<<"\t"<<accept<<endl;
 outfile.close();

}else{
  T->Add(Form("%s/prexRHRS_%i_-1*.root",path.Data(),run));
  double upadc_cutR_approx = 502;

  TH1F* huq = new TH1F(Form("huq"), Form("ADC raw (run%d);ADC raw;Events/CH",run), 300,400,700);
  huq->SetLineColor(1);
  T->Project(huq->GetName(),"P.upQadcR","","goff");
  huq->Draw();

  TH1F* huqCpy = (TH1F*)huq->Clone(Form("huqCpy"));
  huqCpy->GetXaxis()->SetRangeUser(upadc_cutR_approx-10,upadc_cutR_approx+10);
  double upadc_cutR = huqCpy->GetBinCenter(huqCpy->GetMinimumBin());

  TCut cut = Form("(R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)");
  TCut cutADCup = Form("(R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)&&P.upQadcR> %f",upadc_cutR);
  TCut cutPEDup = Form("(R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)&&P.upQadcR< %f",upadc_cutR);

  TH1F* hp1 = new TH1F(Form("hp1"), Form("RHRS momentum (run%d)",run), 200, 0.942, 0.951);
  TH1F* hpup2 = new TH1F(Form("hpup2"), "RHRS momentum + adc>cut", 200, 0.942, 0.951);
  TH1F* hpup3 = new TH1F(Form("hpup3"), "RHRS momentum + adc<cut", 200, 0.942, 0.951);
  T->Project(hp1->GetName(), "R.gold.p", cut);
  T->Project(hpup2->GetName(), "R.gold.p", cutADCup);
  T->Project(hpup3->GetName(), "R.gold.p", cutPEDup);
  hp1->SetXTitle("R.gold.p (GeV)");
  hp1->SetLineColor(1);
  hpup2->SetLineColor(2);
  hpup3->SetLineColor(4);
  TF1* epfit = new TF1("epfit","gaus",hpup2->GetBinCenter(hpup2->GetMaximumBin())-0.0003,hpup2->GetBinCenter(hpup2->GetMaximumBin())+0.0003);
  hpup2->Fit(epfit,"R");
  TCanvas* c1 = new TCanvas(Form("c1"),"c1",900,700);
  c1->Divide(2,2);

  c1->cd(1);
  gPad->SetPad(0.0,0.5,0.4,1.0);
  gPad->Modified();gPad->Update();
  gPad->SetGridx();
  gPad->SetGridy();
  huq->Draw();
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextColor(1);
  latex.DrawLatex(0.53,0.70,"P.upQadcR");
  latex.SetTextColor(9);
  latex.DrawLatex(0.53,0.65,Form("Ped Cut = %.1f",upadc_cutR));
  TLine* pedL = new TLine(upadc_cutR,0.0,upadc_cutR,huq->GetMaximum());
  pedL->SetLineColor(9);
  pedL->SetLineWidth(2);
  pedL->Draw();

  c1->cd(2);
  gPad->SetPad(0.4,0.5,1.0,1.0);
  gPad->SetBottomMargin(0);
  gPad->Modified();gPad->Update();
  gStyle->SetOptStat(1);
  gPad->SetGridx();
  hp1->Draw();
  hpup2->Draw("sames hist");
  hpup3->Draw("sames hist");

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

  TLine* q_edge = new TLine(EdgeP,0.0,EdgeP,hp1->GetMaximum()/2.0);
  q_edge->SetLineColor(6);
  q_edge->SetLineWidth(2);
  q_edge->Draw();

  latex.SetTextColor(6);
  latex.DrawLatex(0.15,0.85,Form("q_edge: %.4f GeV",EdgeP));
  latex.SetTextColor(1);
  latex.DrawLatex(0.15,0.80,"All events");
  latex.SetTextColor(2);
  latex.DrawLatex(0.15,0.75,"US accepted");
  latex.SetTextColor(4);
  latex.DrawLatex(0.15,0.70,"US missed");

  c1->cd(3);
  gPad->SetPad(0.0,0.0,0.4,0.5);
  gPad->SetRightMargin(0);
  c1->cd(4);
  gPad->SetLogy();
  gPad->SetPad(0.4,0.0,1.0,0.5);
  gPad->SetTopMargin(0);
  TH1F* hpup2Cpy = (TH1F*)hpup2->Clone("hpup2Cpy");
  hpup2Cpy->Divide(hp1);
  hpup2Cpy->SetTitle(";R.gold.p (GeV);");
  hpup2Cpy->SetStats(0);
  hpup2Cpy->SetMarkerColor(1);
  hpup2Cpy->SetLineColor(1);
  hpup2Cpy->Draw("hist");
  q_edge->Draw();
  latex.SetTextColor(1);
  latex.DrawLatex(0.40,0.60,"quartz acceptance (accepted/all)");
  double accept = hpup2Cpy->Interpolate(EdgeP);
  latex.SetTextColor(6);
  latex.DrawLatex(0.40,0.50,Form("acceptance @ q_edge = %.3f",accept));
  c1->SaveAs(Form("./temp/plot1.pdf"));

  double thisThtg, thisPhtg, thisP, thisDet;
  double thisu1, thisu2, thisv1, thisv2;
  double thisXvdc[10], thisThvdc[10], thisYvdc[10], thisPhvdc[10];
  double thisAsym, thisAsymSt;
  T->SetMakeClass(1);//This is required to reead sub-branches
  Int_t EvtHdr, Event_Branch;
  Int_t thisTrig, fEvtType;
 
  double xq,yq;//For quartz projection
  double thisAngle; //for reconstructed angle
  double thisCosAng; // for cos angle  
  double thisQsq;//for Qsq
  double E = 0.95;//beam energy

  //Central Angle
  double th0 = 4.74*TMath::Pi()/180;
  double cth0 = TMath::Cos(th0);
  double sth0 = TMath::Sin(th0);

  T->SetBranchAddress("R.gold.th",&thisThtg);
  T->SetBranchAddress("R.gold.ph",&thisPhtg);
  T->SetBranchAddress("R.gold.p",&thisP);
  T->SetBranchAddress("P.upQadcR",&thisDet);
  T->SetBranchAddress("R.vdc.u1.nclust",&thisu1);
  T->SetBranchAddress("R.vdc.v1.nclust",&thisv1);
  T->SetBranchAddress("R.vdc.u2.nclust",&thisu2);
  T->SetBranchAddress("R.vdc.v2.nclust",&thisv2);
  T->SetBranchAddress("Event_Branch",&Event_Branch);
  T->SetBranchAddress("fEvtHdr",&EvtHdr);
  T->SetBranchAddress("fEvtHdr.fEvtType",&thisTrig);
  
  // --- This segfaults when I loop it for some reason
  T->SetBranchAddress("R.tr.x",&thisXvdc[0]);
  T->SetBranchAddress("R.tr.y",&thisYvdc[0]);
  T->SetBranchAddress("R.tr.th",&thisThvdc[0]);
  T->SetBranchAddress("R.tr.ph",&thisPhvdc[0]);
 
  //Asymmetry table is in degrees so we need to convert from radians to degrees
  double radtodeg = 180/TMath::Pi();

  //run21121 doesn't like full entries (memory issue). Make n=1000000 for this run.
  long n;
  if(run==21121)
  n=1e6;
  else
  n = T->GetEntries();

  for(long i = 0; i < n; i++){ 

    T->GetEntry(i);
    //VDC cut and P cut
     if( thisu1==1 && thisv1==1 && thisu2==1 && thisv2==1 && thisP > EdgeP && thisTrig==1){
 
     thisCosAng = (cth0 + thisPhtg*sth0)/TMath::Sqrt(1+thisThtg*thisThtg+thisPhtg*thisPhtg);
     thisAngle = radtodeg*TMath::ACos(thisCosAng);
    
     thisQsq = 2*E*thisP*(1-thisCosAng);

      Qsq.push_back(thisQsq); 
      Angle.push_back(thisAngle);
      //Energy table in MeV so need to be careful here
      thisAsym = 1e6*Interpolate(thisP*1000,thisAngle,0,1);//A (ppm)
      thisAsymSt = 1e6*Interpolate(thisP*1000,thisAngle,1,1); // Stretched A (ppm); 

      ASYM.push_back(thisAsym);ASYM_ST.push_back(thisAsymSt);//Stretched A (ppm)
      Sens.push_back(fabs(thisAsym-thisAsymSt)/thisAsym);//Sensitivity - Why not?
  }
 }

 //Histograms 
 TH1F *Q2 = new TH1F("Q2","Q^{2} (GeV/c)^{2}",150,0,0.015);
 TH1F *Theta = new TH1F("Theta","#theta_{lab} (deg)",150,2,8);
 TH1F *Asym = new TH1F("Asym","Asymmetry (ppm)",150,0,1);
 TH1F *SAsym = new TH1F("SAsym","Stretched Asymmetry (ppm)",150,0,1);
 TH1F *sens = new TH1F("sens","Sensitivity",150,0,0.05);

 Q2->SetTitle(Form("Q^{2} (GeV/c)^{2}, EdgeP = %1.4f GeV",EdgeP));
 Theta->SetTitle(Form("#theta_{lab} (deg), EdgeP = %1.4f GeV",EdgeP));
 Asym->SetTitle(Form("Asymmetry (ppm), EdgeP = %1.4f GeV",EdgeP));
 SAsym->SetTitle(Form("Stretched Asym. (ppm),  EdgeP = %1.4f GeV",EdgeP));
 sens->SetTitle(Form("Sensitivity, EdgeP = %1.4f GeV",EdgeP));

 Q2->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
 Theta->GetXaxis()->SetTitle("#theta_{lab} (deg)");
 Asym->GetXaxis()->SetTitle("Asymmetry (ppm)");
 SAsym->GetXaxis()->SetTitle("Stretched Asym. (ppm)");

 for(int k = 0; k < Qsq.size(); k++){

 Q2->Fill(Qsq[k]);
 Theta->Fill(Angle[k]);
 Asym->Fill(ASYM[k]);
 SAsym->Fill(ASYM_ST[k]);
 sens->Fill(Sens[k]);
 }

 TCanvas *t = new TCanvas("t","t",600,500);
 Theta->Draw();
 t->SaveAs("./temp/plot2.pdf");

 TCanvas *au = new TCanvas("au","au",600,500);
 Asym->Draw();
 au->SaveAs("./temp/plot3.pdf");
 
 TCanvas *as = new TCanvas("as","as",600,500);
 SAsym->Draw();
 as->SaveAs("./temp/plot4.pdf");

 TCanvas *q = new TCanvas("q","q",600,500);
 Q2->Draw();
 q->SaveAs("./temp/plot5.pdf");

 TCanvas *s = new TCanvas("s","s",600,500);
 sens->Draw();
 s->SaveAs("./temp/plot6.pdf");
 
 gSystem->Exec(Form("pdfunite ./temp/plot*.pdf ./temp/all_EdgeP_R_run%d.pdf",run));
 gSystem->Exec(Form("rm -rf ./temp/plot*.pdf"));

  ofstream outfile("./TextFiles/ped_EdgeP_R.csv",ios_base::app);
  outfile<<run<<"\t"<<upadc_cutR<<"\t"<<EdgeP<<"\t"<<epfit->GetParameter(1)<<"\t"<<Theta->GetMean()<<"\t"<<Theta->GetRMS()<<"\t"<<Asym->GetMean()<<"\t"<<Asym->GetRMS()<<"\t"<<SAsym->GetMean()<<"\t"<<SAsym->GetRMS()<<"\t"<<Q2->GetMean()<<"\t"<<Q2->GetRMS()<<"\t"<<sens->GetMean()<<"\t"<<sens->GetRMS()<<"\t"<<epfit->GetParameter(1)-EdgeP<<"\t"<<accept<<endl;
 outfile.close();
}
}
