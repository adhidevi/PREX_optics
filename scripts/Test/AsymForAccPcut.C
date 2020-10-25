#include "CREXdata.h"
#include "TChain.h"
#include "TMath.h"

void AsymForAccPcut(int run, double pCut){

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
    T->Add(Form("%s/prexLHRS_%i_-1.root",path.Data(),run));
   //plot procected x, y, and adc
      double upadc_cutL_approx = 480;
      TCanvas *c2 = new TCanvas("c2","c2",700,500);
      c2->Divide(2,1);
      c2->cd(1);
      gPad->SetPad(0.0,0.0,0.4,1.0);
      gPad->Modified();gPad->Update();
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* huq = new TH1F("huq", Form("ADC raw (run%d, detZ = 1.3 m);ADC raw;Events/CH",run), 300,400,700);
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
      TCut cut = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)");
      TCut cutADCup = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)&&P.upQadcL> %f",upadc_cutL);
      TCut cutPEDup = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)&&P.upQadcL< %f",upadc_cutL);
      double p_edge = 0.9474;

      TCut p_cut = Form("L.gold.p > %f",p_edge);

      c2->cd(2);
      gPad->SetPad(0.4,0.0,1.0,1.0);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(1);
      gPad->SetGridx();
      TH1F* hp1 = new TH1F("hp1", Form("LHRS momentum run%d",run), 200, 0.91, 0.955);
      TH1F* hpup2 = new TH1F("hpup2", "LHRS momentum + adc>cut", 200, 0.91, 0.955);
      TH1F* hpup3 = new TH1F("hpup3", "LHRS momentum + adc<cut", 200, 0.91, 0.955);
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

      latex.SetTextColor(1);
      latex.DrawLatex(0.20,0.85,"All events");
      latex.SetTextColor(2);
      latex.DrawLatex(0.20,0.80,"US accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.20,0.75,"US missed");
      latex.SetTextColor(6);
      latex.DrawLatex(0.20,0.70,Form("p@edge = %1.4f GeV",p_edge));
      latex.SetTextColor(46);
      latex.DrawLatex(0.20,0.48,Form("pCut = %1.4f GeV",pCut));
      gPad->Update();
      TPaveStats* stat1 = (TPaveStats*)hp1->FindObject("stats");
      TPaveStats* stat2 = (TPaveStats*)hpup2->FindObject("stats");
      TPaveStats* stat3 = (TPaveStats*)hpup3->FindObject("stats");
      stat1->SetY2NDC(0.67);
      stat1->SetY1NDC(0.52);
      stat2->SetY2NDC(0.67);
      stat2->SetY1NDC(0.52);
      stat3->SetY2NDC(0.67);
      stat3->SetY1NDC(0.52);
      stat1->SetX2NDC(0.32);
      stat1->SetX1NDC(0.12);
      stat2->SetX2NDC(0.52);
      stat2->SetX1NDC(0.32);
      stat3->SetX2NDC(0.72);
      stat3->SetX1NDC(0.52);
      stat1->SetTextColor(1);
      stat2->SetTextColor(2);
      stat3->SetTextColor(4);
      gPad->Modified();
      TLine* pEdge = new TLine(p_edge,0.0,p_edge,hp1->GetMaximum()/2.2);
      pEdge->SetLineColor(6);
      pEdge->SetLineWidth(2);
      pEdge->Draw();
      TLine* pcut = new TLine(pCut,0.0,pCut,hp1->GetMaximum()/2.2);
      pcut->SetLineColor(46);
      pcut->SetLineWidth(2);
      pcut->Draw();

      c2->SaveAs("./temp1/plot1.pdf");

  double thisThtg, thisPhtg, thisP, thisDet;
  double thisu1, thisu2, thisv1, thisv2;
  double thisXvdc[10], thisThvdc[10], thisYvdc[10], thisPhvdc[10];
  double thisAsym, thisAsymSt;
 
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
  
     xq = thisXvdc[0]+1.3*thisThvdc[0]; yq = thisYvdc[0]+1.3*thisPhvdc[0];

    //VDC cut and x cut
     if( thisu1==1 && thisv1==1 && thisu2==1 && thisv2==1 && thisP > pCut){
 
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

 Q2->SetTitle(Form("Q^{2} (GeV/c)^{2}, pCut = %1.3f GeV",pCut));
 Theta->SetTitle(Form("#theta_{lab} (deg), pCut = %1.3f GeV",pCut));
 Asym->SetTitle(Form("Asymmetry (ppm), pCut = %1.3f GeV",pCut));
 SAsym->SetTitle(Form("Stretched Asym. (ppm),  pCut = %1.3f GeV",pCut));
 sens->SetTitle(Form("Sensitivity, pCut = %1.3f GeV",pCut));

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
 t->SaveAs("./temp1/plot2.pdf");

 TCanvas *au = new TCanvas("au","au",600,500);
 Asym->Draw();
 au->SaveAs("./temp1/plot3.pdf");
 
 TCanvas *as = new TCanvas("as","as",600,500);
 SAsym->Draw();
 as->SaveAs("./temp1/plot4.pdf");

 TCanvas *q = new TCanvas("q","q",600,500);
 Q2->Draw();
 q->SaveAs("./temp1/plot5.pdf");

 TCanvas *s = new TCanvas("s","s",600,500);
 sens->Draw();
 s->SaveAs("./temp1/plot6.pdf");
 
 gSystem->Exec(Form("pdfunite ./temp1/*.pdf ./plots/all_pCut_%1.3f_run%d.pdf",pCut,run));
 gSystem->Exec(Form("rm -rf ./temp1/*.pdf"));

 ofstream outfile(Form("./TextFiles/output_pCut_run%d.csv",run),ios_base::app);
 outfile<<p_edge<<"\t"<<pCut<<"\t"<<Theta->GetMean()<<"\t"<<Asym->GetMean()<<"\t"<<SAsym->GetMean()<<"\t"<<Q2->GetMean()<<"\t"<<sens->GetMean()<<endl;
 outfile.close();

 std::cout << "Theta " << Theta->GetMean() << "  " << "<A> " << Asym->GetMean() << " " << "Stretched <A>" << SAsym->GetMean() << "   " << "Q2 " << Q2->GetMean() << "  " << sens->GetMean() << std::endl;

}else{
    T->Add(Form("%s/prexRHRS_%i_-1.root",path.Data(),run));
   //plot procected x, y, and adc
      double upadc_cutR_approx = 502;
      TCanvas *c2 = new TCanvas("c2","c2",800,500);
      c2->Divide(2,1);
      c2->cd(1);
      gPad->SetPad(0.0,0.0,0.4,1.0);
      gPad->Modified();gPad->Update();
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* huq = new TH1F("huq", Form("ADC raw (run%d, detZ = 1.3 m);ADC raw;Events/CH",run), 200,450,650);
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
      TCut cut = Form("(R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)");
      TCut cutADCup = Form("(R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)&&P.upQadcR> %f",upadc_cutR);
      TCut cutPEDup = Form("(R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)&&P.upQadcR< %f",upadc_cutR);
      double p_edge = 0.9454;

      TCut p_cut = Form("R.vdc.p > %f",p_edge);

      c2->cd(2);
      gPad->SetPad(0.4,0.0,1.0,1.0);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(1);
      gPad->SetGridx();
      TH1F* hp1 = new TH1F("hp1", Form("RHRS momentum (run%d)",run), 200, 0.91, 0.955);
      TH1F* hpup2 = new TH1F("hpup2", Form("RHRS momentum (run%d) + adc>cut",run), 200, 0.91, 0.955);
      TH1F* hpup3 = new TH1F("hpup3", Form("RHRS momentum (run%d) + adc<cut",run), 200, 0.91, 0.955);
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

      latex.SetTextColor(1);
      latex.DrawLatex(0.20,0.85,"All events");
      latex.SetTextColor(2);
      latex.DrawLatex(0.20,0.80,"US accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.20,0.75,"US missed");
      latex.SetTextColor(6);
      latex.DrawLatex(0.20,0.70,Form("p@edge = %1.4f GeV",p_edge));
      latex.SetTextColor(46);
      latex.DrawLatex(0.20,0.48,Form("pCut = %1.4f GeV",pCut));
      gPad->Update();
      TPaveStats* stat1 = (TPaveStats*)hp1->FindObject("stats");
      TPaveStats* stat2 = (TPaveStats*)hpup2->FindObject("stats");
      TPaveStats* stat3 = (TPaveStats*)hpup3->FindObject("stats");
      stat1->SetY2NDC(0.67);
      stat1->SetY1NDC(0.52);
      stat2->SetY2NDC(0.67);
      stat2->SetY1NDC(0.52);
      stat3->SetY2NDC(0.67);
      stat3->SetY1NDC(0.52);
      stat1->SetX2NDC(0.32);
      stat1->SetX1NDC(0.12);
      stat2->SetX2NDC(0.52);
      stat2->SetX1NDC(0.32);
      stat3->SetX2NDC(0.72);
      stat3->SetX1NDC(0.52);
      stat1->SetTextColor(1);
      stat2->SetTextColor(2);
      stat3->SetTextColor(4);
      gPad->Modified();
      TLine* pEdge = new TLine(p_edge,0.0,p_edge,hp1->GetMaximum()/2.2);
      pEdge->SetLineColor(6);
      pEdge->SetLineWidth(2);
      pEdge->Draw();
      TLine* pcut = new TLine(pCut,0.0,pCut,hp1->GetMaximum()/2.2);
      pcut->SetLineColor(46);
      pcut->SetLineWidth(2);
      pcut->Draw();

      c2->SaveAs("./temp1/plot1.pdf");

  double thisThtg, thisPhtg, thisP, thisDet;
  double thisu1, thisu2, thisv1, thisv2;
  double thisXvdc[10], thisThvdc[10], thisYvdc[10], thisPhvdc[10];
  double thisAsym, thisAsymSt;
 
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
  
  // --- This segfaults when I loop it for some reason
  T->SetBranchAddress("R.tr.x",&thisXvdc[0]);
  T->SetBranchAddress("R.tr.y",&thisYvdc[0]);
  T->SetBranchAddress("R.tr.th",&thisThvdc[0]);
  T->SetBranchAddress("R.tr.ph",&thisPhvdc[0]);
 
  //Asymmetry table is in degrees so we need to convert from radians to degrees
  double radtodeg = 180/TMath::Pi();

  long n = T->GetEntries();

  for(long i = 0; i < n; i++){ 

    T->GetEntry(i);
  
    xq = thisXvdc[0]+1.3*thisThvdc[0]; yq = thisYvdc[0]+1.3*thisPhvdc[0];

    //VDC cut and x cut
     if( thisu1==1 && thisv1==1 && thisu2==1 && thisv2==1 && thisP > pCut){
 
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

 Q2->SetTitle(Form("Q^{2} (GeV/c)^{2}, pCut = %1.3f GeV",pCut));
 Theta->SetTitle(Form("#theta_{lab} (deg), pCut = %1.3f GeV",pCut));
 Asym->SetTitle(Form("Asymmetry (ppm), pCut = %1.3f GeV",pCut));
 SAsym->SetTitle(Form("Stretched Asym. (ppm),  pCut = %1.3f GeV",pCut));
 sens->SetTitle(Form("Sensitivity, pCut = %1.3f GeV",pCut));

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
 t->SaveAs("./temp1/plot2.pdf");

 TCanvas *au = new TCanvas("au","au",600,500);
 Asym->Draw();
 au->SaveAs("./temp1/plot3.pdf");
 
 TCanvas *as = new TCanvas("as","as",600,500);
 SAsym->Draw();
 as->SaveAs("./temp1/plot4.pdf");

 TCanvas *q = new TCanvas("q","q",600,500);
 Q2->Draw();
 q->SaveAs("./temp1/plot5.pdf");

 TCanvas *s = new TCanvas("s","s",600,500);
 sens->Draw();
 s->SaveAs("./temp1/plot6.pdf");
 
 gSystem->Exec(Form("pdfunite ./temp1/*.pdf ./plots/all_pCut_%1.3f_run%d.pdf",pCut,run));
 gSystem->Exec(Form("rm -rf ./temp1/*.pdf"));

 ofstream outfile(Form("./TextFiles/output_pCut_run%d.csv",run),ios_base::app);
 outfile<<p_edge<<"\t"<<pCut<<"\t"<<Theta->GetMean()<<"\t"<<Asym->GetMean()<<"\t"<<SAsym->GetMean()<<"\t"<<Q2->GetMean()<<"\t"<<sens->GetMean()<<endl;
 outfile.close();

 std::cout << "Theta " << Theta->GetMean() << "  " << "<A> " << Asym->GetMean() << " " << "Stretched <A>" << SAsym->GetMean() << "   " << "Q2 " << Q2->GetMean() << "  " << sens->GetMean() << std::endl;
}
}
