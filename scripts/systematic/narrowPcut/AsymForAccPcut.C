#include "CREXdata.h"
#include "TChain.h"
#include "TMath.h"

void AsymForAccPcut(int run, TString date, TString target){

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
      TH1F* huq = new TH1F("huq", Form("ADC raw (run%d), %s, %s;ADC raw;Events/CH",run,date.Data(),target.Data()), 300,400,700);
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
      double p_edge;
      if(run==1983)
      p_edge = 0.946812;
      else if(run==1996)
      p_edge = 0.947217;
      else if(run==2052)
      p_edge = 0.947353;
      else if(run==2199)
      p_edge = 0.947758;
      else if(run==2291)
      p_edge = 0.947353;
      else if(run==2292)
      p_edge = 0.947353;
      else if(run==2293)
      p_edge = 0.947398;
      else if(run==2294)
      p_edge = 0.947353;
      else if(run==2316)
      p_edge = 0.947217;
      else if(run==2317)
      p_edge = 0.947307;
      else if(run==2322)
      p_edge = 0.947398;

      TCut p_cut = Form("L.gold.p > %f",p_edge);

      ifstream pCutfile("./TextFiles/pCutL.list");
      if(pCutfile==NULL){
      cout<<"pCut file doesn't exist! Quiting..."<<endl;
      exit(0);
      }
      double p_Cut;
      vector<double>pCut;
      while(pCutfile>>p_Cut)
      pCut.push_back(p_Cut);
      pCutfile.close();
      int nCut = pCut.size();
      for(int iCut=0;iCut<nCut;iCut++)
      cout<<"Plotting for the cuts: "<<pCut[iCut]<<endl;

      c2->cd(2);
      gPad->SetPad(0.4,0.0,1.0,1.0);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(1);
      gPad->SetGridx();
      TH1F* hp1 = new TH1F("hp1", Form("LHRS momentum run%d",run), 200, 0.94, 0.954);
      TH1F* hpup2 = new TH1F("hpup2", "LHRS momentum + adc>cut", 200, 0.94, 0.954);
      TH1F* hpup3 = new TH1F("hpup3", "LHRS momentum + adc<cut", 200, 0.94, 0.954);
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

      latex.SetTextColor(6);
      latex.DrawLatex(0.13,0.85,Form("p@edge = %1.4f GeV",p_edge));
      latex.SetTextColor(1);
      latex.DrawLatex(0.15,0.80,"All events");
      latex.SetTextColor(2);
      latex.DrawLatex(0.15,0.75,"US accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.15,0.70,"US missed");
      latex.SetTextColor(46);
 //     latex.DrawLatex(0.20,0.48,Form("pCut = %1.4f GeV",pCut));
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
      TLine* pEdge = new TLine(p_edge,0.0,p_edge,hp1->GetMaximum()/2.2);
      pEdge->SetLineColor(6);
      pEdge->SetLineWidth(2);
      pEdge->Draw();
      c2->SaveAs("./temp1/plot1.pdf");

  double thisThtg, thisPhtg, thisP, thisDet;
  double thisu1, thisu2, thisv1, thisv2;
  double thisXvdc[10], thisThvdc[10], thisYvdc[10], thisPhvdc[10];
  double thisAsym, thisAsymSt;
  T->SetMakeClass(1);//This is required to read sub-branches
  Int_t EvtHdr, Event_Branch; 
  Int_t thisTrig, fEvtType;
 
  double xq,yq;//For quartz projection
  double thisAngle; //for reconstructed angle
  double thisCosAng; // for cos angle  
  double thisQsq;//for Qsq
  double E = 0.951;//beam energy

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

  TH1F *Q2[nCut];
  TH1F *Theta[nCut];
  TH1F *Asym[nCut];
  TH1F *SAsym[nCut];
  TH1F *sens[nCut];

  for(int iCut=0;iCut<nCut;iCut++){
  for(long i = 0; i < n; i++){ 

    T->GetEntry(i);
  
    //VDC cut and pCut
     if( thisu1==1 && thisv1==1 && thisu2==1 && thisv2==1 && thisP > pCut[iCut] && thisTrig==1){
 
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
 Q2[iCut] = new TH1F(Form("Q2_%.4ff",pCut[iCut]),"Q^{2} (GeV/c)^{2}",150,0,0.015);
 Theta[iCut] = new TH1F(Form("Theta_%.4f",pCut[iCut]),"#theta_{lab} (deg)",150,2,8);
 Asym[iCut] = new TH1F(Form("Asym_%.4f",pCut[iCut]),"Asymmetry (ppm)",150,0,1);
 SAsym[iCut] = new TH1F(Form("SAsym_%.4f",pCut[iCut]),"Stretched Asymmetry (ppm)",150,0,1);
 sens[iCut] = new TH1F(Form("sens_%.4f",pCut[iCut]),"Sensitivity",150,0,0.05);

 for(int k = 0; k < Qsq.size(); k++){
 Q2[iCut]->Fill(Qsq[k]);
 Theta[iCut]->Fill(Angle[k]);
 Asym[iCut]->Fill(ASYM[k]);
 SAsym[iCut]->Fill(ASYM_ST[k]);
 sens[iCut]->Fill(Sens[k]);
 }
 cout<<Form("Done with pCut: %.4f",pCut[iCut])<<endl;
 Qsq.clear();
 Angle.clear();
 ASYM.clear();
 ASYM_ST.clear();
 Sens.clear();
 }

 Q2[0]->SetTitle(Form("Q^{2} (GeV/c)^{2} with pCut (run%d);Q^{2} (GeV/c)^{2}",run));
 Theta[0]->SetTitle(Form("#theta_{lab} (deg) with pCut (run%d);#theta_{lab} (deg)",run));
 Asym[0]->SetTitle(Form("Asymmetry (ppm) with pCut (run%d);Asymmetry (ppm)",run));
 SAsym[0]->SetTitle(Form("Stretched Asym. (ppm) with pCut (run%d);Stretched Asym. (ppm)",run));
 sens[0]->SetTitle(Form("Sensitivity with pCut (run%d);Sensitivity",run));


 TCanvas *t = new TCanvas("t","t",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 Theta[iCut]->SetLineColor(42);
 else if(iCut==6)
 Theta[iCut]->SetLineColor(9);
 else
 Theta[iCut]->SetLineColor(iCut+1);
 Theta[iCut]->Draw("sames histo");
 Theta[iCut]->Scale(1.0/Theta[iCut]->GetMaximum());
 }

 gPad->Update();
 TPaveStats* statsTh[nCut];
 for(int iCut=0;iCut<nCut;iCut++){
 statsTh[iCut] = (TPaveStats*)Theta[iCut]->FindObject("stats");
 }
 statsTh[0]->SetY2NDC(0.90);
 statsTh[0]->SetY1NDC(0.75);
 statsTh[1]->SetY2NDC(0.75);
 statsTh[1]->SetY1NDC(0.60);
 statsTh[2]->SetY2NDC(0.60);
 statsTh[2]->SetY1NDC(0.45);
 statsTh[3]->SetY2NDC(0.45);
 statsTh[3]->SetY1NDC(0.30);

 statsTh[0]->SetX2NDC(0.90);
 statsTh[0]->SetX1NDC(0.70);
 statsTh[1]->SetX2NDC(0.90);
 statsTh[1]->SetX1NDC(0.70);
 statsTh[2]->SetX2NDC(0.90);
 statsTh[2]->SetX1NDC(0.70);
 statsTh[3]->SetX2NDC(0.90);
 statsTh[3]->SetX1NDC(0.70);

 statsTh[4]->SetY2NDC(0.90);
 statsTh[4]->SetY1NDC(0.75);
 statsTh[5]->SetY2NDC(0.75);
 statsTh[5]->SetY1NDC(0.60);
 statsTh[6]->SetY2NDC(0.60);
 statsTh[6]->SetY1NDC(0.45);

 statsTh[4]->SetX2NDC(0.35);
 statsTh[4]->SetX1NDC(0.15);
 statsTh[5]->SetX2NDC(0.35);
 statsTh[5]->SetX1NDC(0.15);
 statsTh[6]->SetX2NDC(0.35);
 statsTh[6]->SetX1NDC(0.15);

 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 statsTh[iCut]->SetTextColor(42);
 else if(iCut==6)
 statsTh[iCut]->SetTextColor(9);
 else
 statsTh[iCut]->SetTextColor(iCut+1);
 }

 gPad->Modified();
 t->SaveAs("./temp1/plot2.pdf");

 TCanvas *au = new TCanvas("au","au",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 Asym[iCut]->SetLineColor(42);
 else if(iCut==6)
 Asym[iCut]->SetLineColor(9);
 else
 Asym[iCut]->SetLineColor(iCut+1);
 Asym[iCut]->Draw("sames histo");
 Asym[iCut]->Scale(1.0/Asym[iCut]->GetMaximum());
 }
 gPad->Update();
 TPaveStats* statsA[nCut];
 for(int iCut=0;iCut<nCut;iCut++){
 statsA[iCut] = (TPaveStats*)Asym[iCut]->FindObject("stats");
 }
 statsA[0]->SetY2NDC(0.90);
 statsA[0]->SetY1NDC(0.75);
 statsA[1]->SetY2NDC(0.75);
 statsA[1]->SetY1NDC(0.60);
 statsA[2]->SetY2NDC(0.60);
 statsA[2]->SetY1NDC(0.45);
 statsA[3]->SetY2NDC(0.45);
 statsA[3]->SetY1NDC(0.30);

 statsA[0]->SetX2NDC(0.90);
 statsA[0]->SetX1NDC(0.70);
 statsA[1]->SetX2NDC(0.90);
 statsA[1]->SetX1NDC(0.70);
 statsA[2]->SetX2NDC(0.90);
 statsA[2]->SetX1NDC(0.70);
 statsA[3]->SetX2NDC(0.90);
 statsA[3]->SetX1NDC(0.70);

 statsA[4]->SetY2NDC(0.90);
 statsA[4]->SetY1NDC(0.75);
 statsA[5]->SetY2NDC(0.75);
 statsA[5]->SetY1NDC(0.60);
 statsA[6]->SetY2NDC(0.60);
 statsA[6]->SetY1NDC(0.45);

 statsA[4]->SetX2NDC(0.35);
 statsA[4]->SetX1NDC(0.15);
 statsA[5]->SetX2NDC(0.35);
 statsA[5]->SetX1NDC(0.15);
 statsA[6]->SetX2NDC(0.35);
 statsA[6]->SetX1NDC(0.15);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 statsA[iCut]->SetTextColor(42);
 else if(iCut==6)
 statsA[iCut]->SetTextColor(9);
 else
 statsA[iCut]->SetTextColor(iCut+1);
 }
 gPad->Modified();
 au->SaveAs("./temp1/plot3.pdf");

 
 TCanvas *as = new TCanvas("as","as",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 SAsym[iCut]->SetLineColor(42);
 else if(iCut==6)
 SAsym[iCut]->SetLineColor(9);
 else
 SAsym[iCut]->SetLineColor(iCut+1);
 SAsym[iCut]->Draw("sames histo");
 SAsym[iCut]->Scale(1.0/SAsym[iCut]->GetMaximum());
 }
 gPad->Update();
 TPaveStats* statsSA[nCut];
 for(int iCut=0;iCut<nCut;iCut++){
 statsSA[iCut] = (TPaveStats*)SAsym[iCut]->FindObject("stats");
 }
 statsSA[0]->SetY2NDC(0.90);
 statsSA[0]->SetY1NDC(0.75);
 statsSA[1]->SetY2NDC(0.75);
 statsSA[1]->SetY1NDC(0.60);
 statsSA[2]->SetY2NDC(0.60);
 statsSA[2]->SetY1NDC(0.45);
 statsSA[3]->SetY2NDC(0.45);
 statsSA[3]->SetY1NDC(0.30);

 statsSA[0]->SetX2NDC(0.90);
 statsSA[0]->SetX1NDC(0.70);
 statsSA[1]->SetX2NDC(0.90);
 statsSA[1]->SetX1NDC(0.70);
 statsSA[2]->SetX2NDC(0.90);
 statsSA[2]->SetX1NDC(0.70);
 statsSA[3]->SetX2NDC(0.90);
 statsSA[3]->SetX1NDC(0.70);

 statsSA[4]->SetY2NDC(0.90);
 statsSA[4]->SetY1NDC(0.75);
 statsSA[5]->SetY2NDC(0.75);
 statsSA[5]->SetY1NDC(0.60);
 statsSA[6]->SetY2NDC(0.60);
 statsSA[6]->SetY1NDC(0.45);

 statsSA[4]->SetX2NDC(0.35);
 statsSA[4]->SetX1NDC(0.15);
 statsSA[5]->SetX2NDC(0.35);
 statsSA[5]->SetX1NDC(0.15);
 statsSA[6]->SetX2NDC(0.35);
 statsSA[6]->SetX1NDC(0.15);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 statsSA[iCut]->SetTextColor(42);
 else if(iCut==6)
 statsSA[iCut]->SetTextColor(9);
 else
 statsSA[iCut]->SetTextColor(iCut+1);
 }
 gPad->Modified();
 as->SaveAs("./temp1/plot4.pdf");

 TCanvas *q = new TCanvas("q","q",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 Q2[iCut]->SetLineColor(42);
 else if(iCut==6)
 Q2[iCut]->SetLineColor(9);
 else
 Q2[iCut]->SetLineColor(iCut+1);
 Q2[iCut]->Draw("sames histo");
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
 q->SaveAs("./temp1/plot5.pdf");

 TCanvas *s = new TCanvas("s","s",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 sens[iCut]->SetLineColor(42);
 else if(iCut==6)
 sens[iCut]->SetLineColor(9);
 else
 sens[iCut]->SetLineColor(iCut+1);
 sens[iCut]->Draw("sames histo");
 sens[iCut]->Scale(1.0/sens[iCut]->GetMaximum());
 }
 gPad->Update();
 TPaveStats* statsS[nCut];
 for(int iCut=0;iCut<nCut;iCut++){
 statsS[iCut] = (TPaveStats*)sens[iCut]->FindObject("stats");
 }
 statsS[0]->SetY2NDC(0.90);
 statsS[0]->SetY1NDC(0.75);
 statsS[1]->SetY2NDC(0.75);
 statsS[1]->SetY1NDC(0.60);
 statsS[2]->SetY2NDC(0.60);
 statsS[2]->SetY1NDC(0.45);
 statsS[3]->SetY2NDC(0.45);
 statsS[3]->SetY1NDC(0.30);

 statsS[0]->SetX2NDC(0.90);
 statsS[0]->SetX1NDC(0.70);
 statsS[1]->SetX2NDC(0.90);
 statsS[1]->SetX1NDC(0.70);
 statsS[2]->SetX2NDC(0.90);
 statsS[2]->SetX1NDC(0.70);
 statsS[3]->SetX2NDC(0.90);
 statsS[3]->SetX1NDC(0.70);

 statsS[4]->SetY2NDC(0.90);
 statsS[4]->SetY1NDC(0.75);
 statsS[5]->SetY2NDC(0.75);
 statsS[5]->SetY1NDC(0.60);
 statsS[6]->SetY2NDC(0.60);
 statsS[6]->SetY1NDC(0.45);

 statsS[4]->SetX2NDC(0.35);
 statsS[4]->SetX1NDC(0.15);
 statsS[5]->SetX2NDC(0.35);
 statsS[5]->SetX1NDC(0.15);
 statsS[6]->SetX2NDC(0.35);
 statsS[6]->SetX1NDC(0.15);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 statsS[iCut]->SetTextColor(42);
 else if(iCut==6)
 statsS[iCut]->SetTextColor(9);
 else
 statsS[iCut]->SetTextColor(iCut+1);
 }
 gPad->Modified();
 s->SaveAs("./temp1/plot6.pdf");
 
 gSystem->Exec(Form("pdfunite ./temp1/plot*.pdf ./plots/all_pCutL_run%d.pdf",run));
 gSystem->Exec(Form("rm -rf ./temp1/plot*.pdf"));

 gSystem->Exec(Form("rm -rf ./TextFiles/output_pCut_run%d.csv",run));

 ofstream outfile(Form("./TextFiles/output_pCut_run%d.csv",run));
 for(int iCut=0;iCut<nCut;iCut++){
 outfile<<p_edge<<"\t"<<pCut[iCut]<<"\t"<<Theta[iCut]->GetMean()<<"\t"<<Theta[iCut]->GetRMS()<<"\t"<<Asym[iCut]->GetMean()<<"\t"<<Asym[iCut]->GetRMS()<<"\t"<<SAsym[iCut]->GetMean()<<"\t"<<SAsym[iCut]->GetRMS()<<"\t"<<Q2[iCut]->GetMean()<<"\t"<<Q2[iCut]->GetRMS()<<"\t"<<sens[iCut]->GetMean()<<"\t"<<sens[iCut]->GetRMS()<<endl;

 std::cout << "Theta " << Theta[iCut]->GetMean() << "  " << "<A> " << Asym[iCut]->GetMean() << " " << "Stretched <A>" << SAsym[iCut]->GetMean() << "   " << "Q2 " << Q2[iCut]->GetMean() << "  " << sens[iCut]->GetMean() << std::endl;
}
 outfile.close();
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
      TH1F* huq = new TH1F("huq", Form("ADC raw (run%d), %s, %s;ADC raw;Events/CH",run,date.Data(),target.Data()), 200,450,650);
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

      double p_edge;
      if(run==21108)
      p_edge = 0.945488;
      else if(run==21121)
      p_edge = 0.945398;
      else if(run==21185)
      p_edge = 0.945398;
      else if(run==21344)
      p_edge = 0.945893;
      else if(run==21412)
      p_edge = 0.945758;
      else if(run==21413)
      p_edge = 0.945577;
      else if(run==21414)
      p_edge = 0.945712;
      else if(run==21415)
      p_edge = 0.945667;
      else if(run==21435)
      p_edge = 0.945442;
      else if(run==21436)
      p_edge = 0.945488;
      else if(run==21438)
      p_edge = 0.945577;
      TCut p_cut = Form("R.vdc.p > %f",p_edge);

      ifstream pCutfile("./TextFiles/pCutR.list");
      if(pCutfile==NULL){
      cout<<"pCut file doesn't exist! Quiting..."<<endl;
      exit(0);
      }
      double p_Cut;
      vector<double>pCut;
      while(pCutfile>>p_Cut)
      pCut.push_back(p_Cut);
      pCutfile.close();
      int nCut = pCut.size();
      for(int iCut=0;iCut<nCut;iCut++)
      cout<<"Plotting for the cuts: "<<pCut[iCut]<<endl;

      c2->cd(2);
      gPad->SetPad(0.4,0.0,1.0,1.0);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(1);
      gPad->SetGridx();
      TH1F* hp1 = new TH1F("hp1", Form("RHRS momentum (run%d)",run), 200, 0.938, 0.952);
      TH1F* hpup2 = new TH1F("hpup2", Form("RHRS momentum (run%d) + adc>cut",run), 200, 0.938, 0.952);
      TH1F* hpup3 = new TH1F("hpup3", Form("RHRS momentum (run%d) + adc<cut",run), 200, 0.938, 0.952);
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

      latex.SetTextColor(6);
      latex.DrawLatex(0.15,0.85,Form("p@edge = %1.4f GeV",p_edge));
      latex.SetTextColor(1);
      latex.DrawLatex(0.17,0.80,"All events");
      latex.SetTextColor(2);
      latex.DrawLatex(0.17,0.75,"US accepted");
      latex.SetTextColor(4);
      latex.DrawLatex(0.17,0.70,"US missed");
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
      TLine* pEdge = new TLine(p_edge,0.0,p_edge,hp1->GetMaximum()/2.2);
      pEdge->SetLineColor(6);
      pEdge->SetLineWidth(2);
      pEdge->Draw();

      c2->SaveAs("./temp1/plot1.pdf");

  double thisThtg, thisPhtg, thisP, thisDet;
  double thisu1, thisu2, thisv1, thisv2;
  double thisXvdc[10], thisThvdc[10], thisYvdc[10], thisPhvdc[10];
  double thisAsym, thisAsymSt;
  T->SetMakeClass(1);//This is required to read sub-branches
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

  long n = T->GetEntries();

  TH1F *Q2[nCut];
  TH1F *Theta[nCut];
  TH1F *Asym[nCut];
  TH1F *SAsym[nCut];
  TH1F *sens[nCut];

  for(int iCut=0;iCut<nCut;iCut++){
  for(long i = 0; i < n; i++){ 

    T->GetEntry(i);
  
    //VDC cut and pCut
     if( thisu1==1 && thisv1==1 && thisu2==1 && thisv2==1 && thisP > pCut[iCut] && thisTrig==1){
 
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
 Q2[iCut] = new TH1F(Form("Q2_%.4ff",pCut[iCut]),"Q^{2} (GeV/c)^{2}",150,0,0.015);
 Theta[iCut] = new TH1F(Form("Theta_%.4f",pCut[iCut]),"#theta_{lab} (deg)",150,2,8);
 Asym[iCut] = new TH1F(Form("Asym_%.4f",pCut[iCut]),"Asymmetry (ppm)",150,0,1);
 SAsym[iCut] = new TH1F(Form("SAsym_%.4f",pCut[iCut]),"Stretched Asymmetry (ppm)",150,0,1);
 sens[iCut] = new TH1F(Form("sens_%.4f",pCut[iCut]),"Sensitivity",150,0,0.05);

 for(int k = 0; k < Qsq.size(); k++){
 Q2[iCut]->Fill(Qsq[k]);
 Theta[iCut]->Fill(Angle[k]);
 Asym[iCut]->Fill(ASYM[k]);
 SAsym[iCut]->Fill(ASYM_ST[k]);
 sens[iCut]->Fill(Sens[k]);
 }
 cout<<Form("Done with pCut: %.4f",pCut[iCut])<<endl;
 Qsq.clear();
 Angle.clear();
 ASYM.clear();
 ASYM_ST.clear();
 Sens.clear();
 }

 Q2[0]->SetTitle(Form("Q^{2} (GeV/c)^{2} with pCut (run%d);Q^{2} (GeV/c)^{2}",run));
 Theta[0]->SetTitle(Form("#theta_{lab} (deg) with pCut (run%d);#theta_{lab} (deg)",run));
 Asym[0]->SetTitle(Form("Asymmetry (ppm) with pCut (run%d);Asymmetry (ppm)",run));
 SAsym[0]->SetTitle(Form("Stretched Asym. (ppm) with pCut (run%d);Stretched Asym. (ppm)",run));
 sens[0]->SetTitle(Form("Sensitivity with pCut (run%d);Sensitivity",run));


 TCanvas *t = new TCanvas("t","t",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 Theta[iCut]->SetLineColor(42);
 else if(iCut==6)
 Theta[iCut]->SetLineColor(9);
 else
 Theta[iCut]->SetLineColor(iCut+1);
 Theta[iCut]->Draw("sames histo");
 Theta[iCut]->Scale(1.0/Theta[iCut]->GetMaximum());
 }

 gPad->Update();
 TPaveStats* statsTh[nCut];
 for(int iCut=0;iCut<nCut;iCut++){
 statsTh[iCut] = (TPaveStats*)Theta[iCut]->FindObject("stats");
 }
 statsTh[0]->SetY2NDC(0.90);
 statsTh[0]->SetY1NDC(0.75);
 statsTh[1]->SetY2NDC(0.75);
 statsTh[1]->SetY1NDC(0.60);
 statsTh[2]->SetY2NDC(0.60);
 statsTh[2]->SetY1NDC(0.45);
 statsTh[3]->SetY2NDC(0.45);
 statsTh[3]->SetY1NDC(0.30);

 statsTh[0]->SetX2NDC(0.90);
 statsTh[0]->SetX1NDC(0.70);
 statsTh[1]->SetX2NDC(0.90);
 statsTh[1]->SetX1NDC(0.70);
 statsTh[2]->SetX2NDC(0.90);
 statsTh[2]->SetX1NDC(0.70);
 statsTh[3]->SetX2NDC(0.90);
 statsTh[3]->SetX1NDC(0.70);

 statsTh[4]->SetY2NDC(0.90);
 statsTh[4]->SetY1NDC(0.75);
 statsTh[5]->SetY2NDC(0.75);
 statsTh[5]->SetY1NDC(0.60);
 statsTh[6]->SetY2NDC(0.60);
 statsTh[6]->SetY1NDC(0.45);

 statsTh[4]->SetX2NDC(0.35);
 statsTh[4]->SetX1NDC(0.15);
 statsTh[5]->SetX2NDC(0.35);
 statsTh[5]->SetX1NDC(0.15);
 statsTh[6]->SetX2NDC(0.35);
 statsTh[6]->SetX1NDC(0.15);

 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 statsTh[iCut]->SetTextColor(42);
 else if(iCut==6)
 statsTh[iCut]->SetTextColor(9);
 else
 statsTh[iCut]->SetTextColor(iCut+1);
 }

 gPad->Modified();
 t->SaveAs("./temp1/plot2.pdf");

 TCanvas *au = new TCanvas("au","au",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 Asym[iCut]->SetLineColor(42);
 else if(iCut==6)
 Asym[iCut]->SetLineColor(9);
 else
 Asym[iCut]->SetLineColor(iCut+1);
 Asym[iCut]->Draw("sames histo");
 Asym[iCut]->Scale(1.0/Asym[iCut]->GetMaximum());
 }
 gPad->Update();
 TPaveStats* statsA[nCut];
 for(int iCut=0;iCut<nCut;iCut++){
 statsA[iCut] = (TPaveStats*)Asym[iCut]->FindObject("stats");
 }
 statsA[0]->SetY2NDC(0.90);
 statsA[0]->SetY1NDC(0.75);
 statsA[1]->SetY2NDC(0.75);
 statsA[1]->SetY1NDC(0.60);
 statsA[2]->SetY2NDC(0.60);
 statsA[2]->SetY1NDC(0.45);
 statsA[3]->SetY2NDC(0.45);
 statsA[3]->SetY1NDC(0.30);

 statsA[0]->SetX2NDC(0.90);
 statsA[0]->SetX1NDC(0.70);
 statsA[1]->SetX2NDC(0.90);
 statsA[1]->SetX1NDC(0.70);
 statsA[2]->SetX2NDC(0.90);
 statsA[2]->SetX1NDC(0.70);
 statsA[3]->SetX2NDC(0.90);
 statsA[3]->SetX1NDC(0.70);

 statsA[4]->SetY2NDC(0.90);
 statsA[4]->SetY1NDC(0.75);
 statsA[5]->SetY2NDC(0.75);
 statsA[5]->SetY1NDC(0.60);
 statsA[6]->SetY2NDC(0.60);
 statsA[6]->SetY1NDC(0.45);

 statsA[4]->SetX2NDC(0.35);
 statsA[4]->SetX1NDC(0.15);
 statsA[5]->SetX2NDC(0.35);
 statsA[5]->SetX1NDC(0.15);
 statsA[6]->SetX2NDC(0.35);
 statsA[6]->SetX1NDC(0.15);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 statsA[iCut]->SetTextColor(42);
 else if(iCut==6)
 statsA[iCut]->SetTextColor(9);
 else
 statsA[iCut]->SetTextColor(iCut+1);
 }
 gPad->Modified();
 au->SaveAs("./temp1/plot3.pdf");

 
 TCanvas *as = new TCanvas("as","as",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 SAsym[iCut]->SetLineColor(42);
 else if(iCut==6)
 SAsym[iCut]->SetLineColor(9);
 else
 SAsym[iCut]->SetLineColor(iCut+1);
 SAsym[iCut]->Draw("sames histo");
 SAsym[iCut]->Scale(1.0/SAsym[iCut]->GetMaximum());
 }
 gPad->Update();
 TPaveStats* statsSA[nCut];
 for(int iCut=0;iCut<nCut;iCut++){
 statsSA[iCut] = (TPaveStats*)SAsym[iCut]->FindObject("stats");
 }
 statsSA[0]->SetY2NDC(0.90);
 statsSA[0]->SetY1NDC(0.75);
 statsSA[1]->SetY2NDC(0.75);
 statsSA[1]->SetY1NDC(0.60);
 statsSA[2]->SetY2NDC(0.60);
 statsSA[2]->SetY1NDC(0.45);
 statsSA[3]->SetY2NDC(0.45);
 statsSA[3]->SetY1NDC(0.30);

 statsSA[0]->SetX2NDC(0.90);
 statsSA[0]->SetX1NDC(0.70);
 statsSA[1]->SetX2NDC(0.90);
 statsSA[1]->SetX1NDC(0.70);
 statsSA[2]->SetX2NDC(0.90);
 statsSA[2]->SetX1NDC(0.70);
 statsSA[3]->SetX2NDC(0.90);
 statsSA[3]->SetX1NDC(0.70);

 statsSA[4]->SetY2NDC(0.90);
 statsSA[4]->SetY1NDC(0.75);
 statsSA[5]->SetY2NDC(0.75);
 statsSA[5]->SetY1NDC(0.60);
 statsSA[6]->SetY2NDC(0.60);
 statsSA[6]->SetY1NDC(0.45);

 statsSA[4]->SetX2NDC(0.35);
 statsSA[4]->SetX1NDC(0.15);
 statsSA[5]->SetX2NDC(0.35);
 statsSA[5]->SetX1NDC(0.15);
 statsSA[6]->SetX2NDC(0.35);
 statsSA[6]->SetX1NDC(0.15);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 statsSA[iCut]->SetTextColor(42);
 else if(iCut==6)
 statsSA[iCut]->SetTextColor(9);
 else
 statsSA[iCut]->SetTextColor(iCut+1);
 }
 gPad->Modified();
 as->SaveAs("./temp1/plot4.pdf");

 TCanvas *q = new TCanvas("q","q",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 Q2[iCut]->SetLineColor(42);
 else if(iCut==6)
 Q2[iCut]->SetLineColor(9);
 else
 Q2[iCut]->SetLineColor(iCut+1);
 Q2[iCut]->Draw("sames histo");
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
 q->SaveAs("./temp1/plot5.pdf");

 TCanvas *s = new TCanvas("s","s",600,500);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 sens[iCut]->SetLineColor(42);
 else if(iCut==6)
 sens[iCut]->SetLineColor(9);
 else
 sens[iCut]->SetLineColor(iCut+1);
 sens[iCut]->Draw("sames histo");
 sens[iCut]->Scale(1.0/sens[iCut]->GetMaximum());
 }
 gPad->Update();
 TPaveStats* statsS[nCut];
 for(int iCut=0;iCut<nCut;iCut++){
 statsS[iCut] = (TPaveStats*)sens[iCut]->FindObject("stats");
 }
 statsS[0]->SetY2NDC(0.90);
 statsS[0]->SetY1NDC(0.75);
 statsS[1]->SetY2NDC(0.75);
 statsS[1]->SetY1NDC(0.60);
 statsS[2]->SetY2NDC(0.60);
 statsS[2]->SetY1NDC(0.45);
 statsS[3]->SetY2NDC(0.45);
 statsS[3]->SetY1NDC(0.30);

 statsS[0]->SetX2NDC(0.90);
 statsS[0]->SetX1NDC(0.70);
 statsS[1]->SetX2NDC(0.90);
 statsS[1]->SetX1NDC(0.70);
 statsS[2]->SetX2NDC(0.90);
 statsS[2]->SetX1NDC(0.70);
 statsS[3]->SetX2NDC(0.90);
 statsS[3]->SetX1NDC(0.70);

 statsS[4]->SetY2NDC(0.90);
 statsS[4]->SetY1NDC(0.75);
 statsS[5]->SetY2NDC(0.75);
 statsS[5]->SetY1NDC(0.60);
 statsS[6]->SetY2NDC(0.60);
 statsS[6]->SetY1NDC(0.45);

 statsS[4]->SetX2NDC(0.35);
 statsS[4]->SetX1NDC(0.15);
 statsS[5]->SetX2NDC(0.35);
 statsS[5]->SetX1NDC(0.15);
 statsS[6]->SetX2NDC(0.35);
 statsS[6]->SetX1NDC(0.15);
 for(int iCut=0;iCut<nCut;iCut++){
 if(iCut==4)
 statsS[iCut]->SetTextColor(42);
 else if(iCut==6)
 statsS[iCut]->SetTextColor(9);
 else
 statsS[iCut]->SetTextColor(iCut+1);
 }
 gPad->Modified();
 s->SaveAs("./temp1/plot6.pdf");
 
 gSystem->Exec(Form("pdfunite ./temp1/plot*.pdf ./plots/all_pCutR_run%d.pdf",run));
 gSystem->Exec(Form("rm -rf ./temp1/plot*.pdf"));
 
 gSystem->Exec(Form("rm -rf ./TextFiles/output_pCut_run%d.csv",run));
 ofstream outfile(Form("./TextFiles/output_pCut_run%d.csv",run));
 for(int iCut=0;iCut<nCut;iCut++){
 outfile<<p_edge<<"\t"<<pCut[iCut]<<"\t"<<Theta[iCut]->GetMean()<<"\t"<<Theta[iCut]->GetRMS()<<"\t"<<Asym[iCut]->GetMean()<<"\t"<<Asym[iCut]->GetRMS()<<"\t"<<SAsym[iCut]->GetMean()<<"\t"<<SAsym[iCut]->GetRMS()<<"\t"<<Q2[iCut]->GetMean()<<"\t"<<Q2[iCut]->GetRMS()<<"\t"<<sens[iCut]->GetMean()<<"\t"<<sens[iCut]->GetRMS()<<endl;

 std::cout << "Theta " << Theta[iCut]->GetMean() << "  " << "<A> " << Asym[iCut]->GetMean() << " " << "Stretched <A>" << SAsym[iCut]->GetMean() << "   " << "Q2 " << Q2[iCut]->GetMean() << "  " << sens[iCut]->GetMean() << std::endl;
}
 outfile.close();
}
}
