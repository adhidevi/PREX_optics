#include "CREXdata.h"
#include "TChain.h"
#include "TMath.h"

//Macro for plotting Q^2 and the Asymmetry for the data
//Want to see how sensitive the asymmetry is to a cut on the radiative tail... for acceptance function
//This has an explicit cut on radiative tail.. Should run ProjQuartz to see ADC plot and quartz edge...
//Ideally wanted to calibrate ADC with dispersive position

void AsymForAccYcutLo(int run, double yloCut){

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
   //plot procected x, y, and adc
      double upadc_cutL_approx = 481;
      TCanvas *c2 = new TCanvas("c2","c2",900,700);
      c2->Divide(2,2);
      c2->cd(1);
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

      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.05);
      latex.SetTextColor(1);
      latex.DrawLatex(0.60,0.70,"P.upQadcL");
      latex.SetTextColor(9);
      latex.DrawLatex(0.60,0.65,"Ped. cut");
      TCut cut = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)");
      TCut cutADCup = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)&&P.upQadcL> %f",upadc_cutL);
      TCut cutPEDup = Form("(L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1)&&P.upQadcL< %f",upadc_cutL);
      double x_edge = -0.0715;// for z=130cm

      TCut x_cut = Form("(L.tr.x[0]+1.3*L.tr.th[0]) > %f",x_edge);

      c2->cd(2);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(1);
      gPad->SetGridx();
      TH1F* hy1 = new TH1F("hy1", "Projected y on detector plane", 200, -0.07, 0.01);
      TH1F* hyup2 = new TH1F("hyup2", "Projected y on detector plane + adc>cut", 200, -0.07, 0.01);
      TH1F* hyup3 = new TH1F("hyup3", "Projected y on detector plane + adc<cut", 200, -0.07, 0.01);
      T->Project(hy1->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cut+x_cut);
      T->Project(hyup2->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cutADCup+x_cut);
      T->Project(hyup3->GetName(), "L.tr.y[0]+1.3*L.tr.ph[0]", cutPEDup+x_cut);
      hy1->SetXTitle("L.tr.y[0]+1.3*L.tr.ph[0] (m)");
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
      latex.SetTextColor(49);
      latex.DrawLatex(0.20,0.70,Form("yloCut = %1.3f",yloCut));     
      TLine* yloCut_line = new TLine(yloCut,0.0,yloCut,hy1->GetMaximum()/3.0*2.0);
      yloCut_line->SetLineColor(49);
      yloCut_line->SetLineWidth(2);
      yloCut_line->Draw();


      c2->cd(3);
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* hx1 = new TH1F("hx1", "Projected x on detector plane", 300, -0.14, 0.03);
      TH1F* hxup2 = new TH1F("hxup2", "Projected x on detector plane + adc>cut", 300, -0.14, 0.03);
      TH1F* hxup3 = new TH1F("hxup3", "Projected x on detector plane + adc<cut", 300, -0.14, 0.03);

      T->Project(hx1->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cut);
      T->Project(hxup2->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cutADCup);
      T->Project(hxup3->GetName(), "L.tr.x[0]+1.3*L.tr.th[0]", cutPEDup);
      hx1->SetXTitle("L.tr.x[0]+1.3*L.tr.th[0] (m)");
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
      float nAccept = T->Project(h2d2->GetName(),"L.tr.x[0]+1.3*L.tr.th[0]:L.tr.y[0]+1.3*L.tr.ph[0]",cutADCup);
      h2d2->SetXTitle("Projected y (m)");
      h2d2->SetYTitle("Projected x (m)");
      h2d2->Draw("COLZ");

      c2->SaveAs("./temp/plot1.pdf");

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
     if( thisu1==1 && thisv1==1 && thisu2==1 && thisv2==1 && xq>-0.0715 && yq>yloCut){
 
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

 Q2->SetTitle(Form("Q^{2} (GeV/c)^{2}, yloCut = %1.3f m",yloCut));
 Theta->SetTitle(Form("#theta_{lab} (deg), yloCut = %1.3f m",yloCut));
 Asym->SetTitle(Form("Asymmetry (ppm), yloCut = %1.3f m",yloCut));
 SAsym->SetTitle(Form("Stretched Asym. (ppm),  yloCut = %1.3f m",yloCut));
 sens->SetTitle(Form("Sensitivity, yloCut = %1.3f m",yloCut));

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
 
 gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/all_yloCut_%1.3f_run%d.pdf",yloCut,run));
 gSystem->Exec(Form("rm -rf ./temp/*.pdf"));

 ofstream outfile(Form("./TextFiles/output_yloCut_run%d.csv",run),ios_base::app);
 outfile<<x_edge<<"\t"<<yloCut<<"\t"<<Theta->GetMean()<<"\t"<<Asym->GetMean()<<"\t"<<SAsym->GetMean()<<"\t"<<Q2->GetMean()<<"\t"<<sens->GetMean()<<endl;
 outfile.close();

 std::cout << "Theta " << Theta->GetMean() << "  " << "<A> " << Asym->GetMean() << " " << "Stretched <A>" << SAsym->GetMean() << "   " << "Q2 " << Q2->GetMean() << "  " << sens->GetMean() << std::endl;
}else{
    T->Add(Form("%s/prexRHRS_%i_-1*.root",path.Data(),run));
   //plot procected x, y, and adc
      double upadc_cutR_approx = 502;
      TCanvas *c2 = new TCanvas("c2","c2",900,700);
      c2->Divide(2,2);
      c2->cd(1);
      gPad->Modified();gPad->Update();
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* huq = new TH1F("huq", Form("ADC raw (run%d, detZ = 1.3 m);ADC raw;Events/CH",run), 300,400,700);
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
      latex.DrawLatex(0.60,0.70,"P.upQadcR");
      latex.SetTextColor(9);
      latex.DrawLatex(0.60,0.65,"Ped. cut");
      TCut cut = Form("(R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)");
      TCut cutADCup = Form("(R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)&&P.upQadcR> %f",upadc_cutR);
      TCut cutPEDup = Form("(R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)&&P.upQadcR< %f",upadc_cutR);
      double x_edge = -0.076;// for z=130cm

      TCut x_cut = Form("(R.tr.x[0]+1.3*R.tr.th[0]) > %f",x_edge);

      c2->cd(2);
      gPad->Modified();gPad->Update();
      gStyle->SetOptStat(1);
      gPad->SetGridx();
      TH1F* hy1 = new TH1F("hy1", "Projected y on detector plane", 200, -0.04, 0.04);
      TH1F* hyup2 = new TH1F("hyup2", "Projected y on detector plane + adc>cut", 200, -0.04, 0.04);
      TH1F* hyup3 = new TH1F("hyup3", "Projected y on detector plane + adc<cut", 200, -0.04, 0.04);
      T->Project(hy1->GetName(), "R.tr.y[0]+1.3*R.tr.ph[0]", cut+x_cut);
      T->Project(hyup2->GetName(), "R.tr.y[0]+1.3*R.tr.ph[0]", cutADCup+x_cut);
      T->Project(hyup3->GetName(), "R.tr.y[0]+1.3*R.tr.ph[0]", cutPEDup+x_cut);
      hy1->SetXTitle("R.tr.y[0]+1.3*R.tr.ph[0] (m)");
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
      latex.SetTextColor(49);
      latex.DrawLatex(0.20,0.70,Form("yloCut = %1.3f",yloCut));     

      TLine* yloCut_line = new TLine(yloCut,0.0,yloCut,hy1->GetMaximum()/3.0*2.0);
      yloCut_line->SetLineColor(49);
      yloCut_line->SetLineWidth(2);
      yloCut_line->Draw();

      c2->cd(3);
      gPad->SetGridx();
      gPad->SetGridy();
      TH1F* hx1 = new TH1F("hx1", "Projected x on detector plane", 300, -0.14, 0.03);
      TH1F* hxup2 = new TH1F("hxup2", "Projected x on detector plane + adc>cut", 300, -0.14, 0.03);
      TH1F* hxup3 = new TH1F("hxup3", "Projected x on detector plane + adc<cut", 300, -0.14, 0.03);

      T->Project(hx1->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cut);
      T->Project(hxup2->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cutADCup);
      T->Project(hxup3->GetName(), "R.tr.x[0]+1.3*R.tr.th[0]", cutPEDup);
      hx1->SetXTitle("R.tr.x[0]+1.3*R.tr.th[0] (m)");
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
      float nAccept = T->Project(h2d2->GetName(),"R.tr.x[0]+1.3*R.tr.th[0]:R.tr.y[0]+1.3*R.tr.ph[0]",cutADCup);
      h2d2->SetXTitle("Projected y (m)");
      h2d2->SetYTitle("Projected x (m)");
      h2d2->Draw("COLZ");

      c2->SaveAs("./temp/plot1.pdf");

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
     if( thisu1==1 && thisv1==1 && thisu2==1 && thisv2==1 && xq>-0.076 && yq>yloCut){
 
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

 Q2->SetTitle(Form("Q^{2} (GeV/c)^{2}, yloCut = %1.3f m",yloCut));
 Theta->SetTitle(Form("#theta_{lab} (deg), yloCut = %1.3f m",yloCut));
 Asym->SetTitle(Form("Asymmetry (ppm), yloCut = %1.3f m",yloCut));
 SAsym->SetTitle(Form("Stretched Asym. (ppm),  yloCut = %1.3f m",yloCut));
 sens->SetTitle(Form("Sensitivity, yloCut = %1.3f m",yloCut));

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
 
 gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/all_yloCut_%1.3f_run%d.pdf",yloCut,run));
 gSystem->Exec(Form("rm -rf ./temp/*.pdf"));

 ofstream outfile(Form("./TextFiles/output_yloCut_run%d.csv",run),ios_base::app);
 outfile<<x_edge<<"\t"<<yloCut<<"\t"<<Theta->GetMean()<<"\t"<<Asym->GetMean()<<"\t"<<SAsym->GetMean()<<"\t"<<Q2->GetMean()<<"\t"<<sens->GetMean()<<endl;
 outfile.close();

 std::cout << "Theta " << Theta->GetMean() << "  " << "<A> " << Asym->GetMean() << " " << "Stretched <A>" << SAsym->GetMean() << "   " << "Q2 " << Q2->GetMean() << "  " << sens->GetMean() << std::endl;
}
}
