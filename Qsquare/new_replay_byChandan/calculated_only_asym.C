#include "CREXdata.h"
#include "TChain.h"
#include "TMath.h"

void calculated_only_asym(int runL, int runR, Float_t E){

  //Loads the Horowitz tables
  LoadTable("horpb.dat",0); //unstretched 
  LoadTable("horpb1.dat",1); //stretched

  //For the interpolation
   int j = 0;
  
  //ROOT 6 memory management makes me use vectors
  vector <double> QsqL, AngleL, ASYML, ASYM_STL, SensL; 
  vector <double> QsqR, AngleR, ASYMR, ASYM_STR, SensR; 
  
   TString path = "/chafs1/work1/prex_counting/QsqRootFiles";
 
   TChain *TL = new TChain("T");
   TChain *TR = new TChain("T");

   TL->Add(Form("%s/prexLHRS_%i_-1*.root",path.Data(),runL));
   TR->Add(Form("%s/prexRHRS_%i_-1*.root",path.Data(),runR));

  double upADCcutL_approx = 480;
  double upADCcutR_approx = 500;
  TCanvas *c1 = new TCanvas("c1","c1",900,700);
  c1->Divide(2,2);
  c1->cd(1);

      TH1F* huqL = new TH1F("huqL",Form("Upstream ADC (run%d, %d);ADC CH;Event/CH",runL,runR),250,400,650);
      TH1F* huqR = new TH1F("huqR",Form("Upstream ADC (run%d, %d);ADC CH;Event/CH",runL,runR),250,400,650);
      huqL->SetLineColor(1);
      huqR->SetLineColor(2);
      TL->Draw("P.upQadcL>>huqL","","hist");
      TR->Draw("P.upQadcR>>huqR","","sames hist");
      huqL->Scale(1./huqL->Integral());
      huqR->Scale(1./huqR->Integral());

      TH1F* huqLCpy = (TH1F*)huqL->Clone("huqLCpy");
      TH1F* huqRCpy = (TH1F*)huqR->Clone("huqRCpy");
      huqLCpy->GetXaxis()->SetRangeUser(upADCcutL_approx-10,upADCcutL_approx+10);
      huqRCpy->GetXaxis()->SetRangeUser(upADCcutR_approx-10,upADCcutR_approx+10);
      double upADCcutL = huqLCpy->GetXaxis()->GetBinCenter(huqLCpy->GetMinimumBin());
      double upADCcutR = huqRCpy->GetXaxis()->GetBinCenter(huqRCpy->GetMinimumBin());
      cout<<"LHRS ADC cut: "<<upADCcutL<<endl;
      cout<<"RHRS ADC cut: "<<upADCcutR<<endl;
      gPad->Update();
      TPaveStats* ADCl = (TPaveStats*)huqL->FindObject("stats");
      TPaveStats* ADCr = (TPaveStats*)huqR->FindObject("stats");
      ADCl->SetY1NDC(0.90);
      ADCl->SetY2NDC(0.75);
      ADCr->SetY1NDC(0.75);
      ADCr->SetY2NDC(0.60);
      ADCl->SetTextColor(1);
      ADCr->SetTextColor(2);
      gPad->Modified();
      TLine* lineL = new TLine(upADCcutL,0.0,upADCcutL,huqL->GetMaximum());
      TLine* lineR = new TLine(upADCcutR,0.0,upADCcutR,huqR->GetMaximum());
      lineL->SetLineColor(4);
      lineR->SetLineColor(6);
      lineL->Draw();
      lineR->Draw();
      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.04);
      latex.SetTextColor(4);
      latex.DrawLatex(0.55,0.55,Form("PedL Cut = %.1f",upADCcutL));
      latex.SetTextColor(6);
      latex.DrawLatex(0.55,0.50,Form("PedR Cut = %.1f",upADCcutR));
      TCut cut_wadcL = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.04&&L.gold.ph<0.05&&L.gold.dp<-0.002&&P.upQadcL>%f)",upADCcutL);
      TCut cut_wadcR = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.04&&R.gold.ph<0.05&&R.gold.dp<-0.002&&P.upQadcR>%f)",upADCcutR);

  double thisThtgL, thisPhtgL, thisPL, thisdPL, thisDetL;
  double thisu1L, thisu2L, thisv1L, thisv2L;
  double thisXvdcL[10], thisThvdcL[10], thisYvdcL[10], thisPhvdcL[10];
  double thisAsymL, thisAsymStL;
  TL->SetMakeClass(1);//This is required to read sub-branches
  Int_t EvtHdrL, Event_BranchL; 
  Int_t thisTrigL, fEvtTypeL;
 
  double thisAngleL; //for reconstructed angle
  double thisCosAngL; // for cos angle  
  double thisQsqL;//for Qsq

  //Central Angle
  double th0L = 4.801*TMath::Pi()/180;
  double cth0L = TMath::Cos(th0L);
  double sth0L = TMath::Sin(th0L);

  double thisThtgR, thisPhtgR, thisPR, thisdPR, thisDetR;
  double thisu1R, thisu2R, thisv1R, thisv2R;
  double thisXvdcR[10], thisThvdcR[10], thisYvdcR[10], thisPhvdcR[10];
  double thisAsymR, thisAsymStR;
  TR->SetMakeClass(1);//This is required to read sub-branches
  Int_t EvtHdrR, Event_BranchR; 
  Int_t thisTrigR, fEvtTypeR;
 
  double thisAngleR; //for reconstructed angle
  double thisCosAngR; // for cos angle  
  double thisQsqR;//for Qsq

  //Central Angle
  double th0R = 4.772*TMath::Pi()/180;
  double cth0R = TMath::Cos(th0R);
  double sth0R = TMath::Sin(th0R);

  TL->SetBranchAddress("L.gold.th",&thisThtgL);
  TL->SetBranchAddress("L.gold.ph",&thisPhtgL);
  TL->SetBranchAddress("L.gold.p",&thisPL);
  TL->SetBranchAddress("L.gold.dp",&thisdPL);
  TL->SetBranchAddress("P.upQadcL",&thisDetL);
  TL->SetBranchAddress("L.vdc.u1.nclust",&thisu1L);
  TL->SetBranchAddress("L.vdc.v1.nclust",&thisv1L);
  TL->SetBranchAddress("L.vdc.u2.nclust",&thisu2L);
  TL->SetBranchAddress("L.vdc.v2.nclust",&thisv2L);
  TL->SetBranchAddress("Event_Branch",&Event_BranchL);
  TL->SetBranchAddress("fEvtHdr",&EvtHdrL);
  TL->SetBranchAddress("fEvtHdr.fEvtType",&thisTrigL);
  
  TR->SetBranchAddress("R.gold.th",&thisThtgR);
  TR->SetBranchAddress("R.gold.ph",&thisPhtgR);
  TR->SetBranchAddress("R.gold.p",&thisPR);
  TR->SetBranchAddress("R.gold.dp",&thisdPR);
  TR->SetBranchAddress("P.upQadcL",&thisDetR);
  TR->SetBranchAddress("R.vdc.u1.nclust",&thisu1R);
  TR->SetBranchAddress("R.vdc.v1.nclust",&thisv1R);
  TR->SetBranchAddress("R.vdc.u2.nclust",&thisu2R);
  TR->SetBranchAddress("R.vdc.v2.nclust",&thisv2R);
  TR->SetBranchAddress("Event_Branch",&Event_BranchR);
  TR->SetBranchAddress("fEvtHdr",&EvtHdrR);
  TR->SetBranchAddress("fEvtHdr.fEvtType",&thisTrigR);
  
  //Asymmetry table is in degrees so we need to convert from radians to degrees
  double radtodeg = 180/TMath::Pi();

  long nL = TL->GetEntries();
  long nR = TR->GetEntries();
  for(long i = 0; i < nL; i++){ 

   TL->GetEntry(i);
  
  if( thisu1L==1 && thisv1L==1 && thisu2L==1 && thisv2L==1 && thisThtgL>-0.08 && thisThtgL<0.08 && thisPhtgL>-0.04 && thisPhtgL<0.05 && thisdPL <-0.002 && thisTrigL==1){
 
  thisCosAngL = (cth0L - thisPhtgL*sth0L)/TMath::Sqrt(1+thisThtgL*thisThtgL+thisPhtgL*thisPhtgL);
  thisAngleL = radtodeg*TMath::ACos(thisCosAngL);
   
  thisQsqL = 2*E*thisPL*(1-thisCosAngL);

  QsqL.push_back(thisQsqL); 
  AngleL.push_back(thisAngleL);
  //Energy table in MeV so need to be careful here
  thisAsymL = 1e6*Interpolate(thisPL*1000,thisAngleL,0,1);//A (ppm)
  thisAsymStL = 1e6*Interpolate(thisPL*1000,thisAngleL,1,1); // Stretched A (ppm); 

  ASYML.push_back(thisAsymL);ASYM_STL.push_back(thisAsymStL);//Stretched A (ppm)
  SensL.push_back(fabs(thisAsymL-thisAsymStL)/thisAsymL);//Sensitivity - Why not?
  }
 }

  for(long i = 0; i < nR; i++){ 

   TR->GetEntry(i);
  
  if( thisu1R==1 && thisv1R==1 && thisu2R==1 && thisv2R==1 && thisThtgR>-0.08 && thisThtgR<0.08 && thisPhtgR>-0.04 && thisPhtgR<0.05 && thisdPR <-0.002 && thisTrigR==1){
 
  thisCosAngR = (cth0R + thisPhtgR*sth0R)/TMath::Sqrt(1+thisThtgR*thisThtgR+thisPhtgR*thisPhtgR);
  thisAngleR = radtodeg*TMath::ACos(thisCosAngR);
   
  thisQsqR = 2*E*thisPR*(1-thisCosAngR);

  QsqR.push_back(thisQsqR); 
  AngleR.push_back(thisAngleR);
  //Energy table in MeV so need to be careful here
  thisAsymR = 1e6*Interpolate(thisPR*1000,thisAngleR,0,1);//A (ppm)
  thisAsymStR = 1e6*Interpolate(thisPR*1000,thisAngleR,1,1); // Stretched A (ppm); 

  ASYMR.push_back(thisAsymR);ASYM_STR.push_back(thisAsymStR);//Stretched A (ppm)
  SensR.push_back(fabs(thisAsymR-thisAsymStR)/thisAsymR);//Sensitivity - Why not?
  }
 }
 //Histograms 
 TH1F* Q2L = new TH1F("Q2L",Form("Q^{2} (GeV/c)^{2} (run%d, %d)",runL,runR),150,0,0.015);
 TH1F* ThetaL = new TH1F("ThetaL",Form("#theta_{lab} (deg) (run%d, run%d)",runL,runR),150,2,8);
 TH1F* AsymL = new TH1F("AsymL",Form("Asymmetry (ppm) (run%d, %d)",runL,runR),150,0,1);
 TH1F* SAsymL = new TH1F("SAsymL",Form("Stretched Asymmetry (ppm) (run%d, %d)",runL,runR),150,0,1);
 TH1F* sensL = new TH1F("sensL",Form("Sensitivity (run%d, %d)",runL,runR),150,0,0.05);
 TH1F* Q2R = new TH1F("Q2R",Form("Q^{2} (GeV/c)^{2} (run%d, %d)",runL,runR),150,0,0.015);
 TH1F* ThetaR = new TH1F("ThetaR",Form("#theta_{lab} (deg) (run%d, run%d)",runL,runR),150,2,8);
 TH1F* AsymR = new TH1F("AsymR",Form("Asymmetry (ppm) (run%d, %d)",runL,runR),150,0,1);
 TH1F* SAsymR = new TH1F("SAsymR",Form("Stretched Asymmetry (ppm) (run%d, %d)",runL,runR),150,0,1);
 TH1F* sensR = new TH1F("sensR",Form("Sensitivity (run%d, %d)",runL,runR),150,0,0.05);

 for(int k = 0; k < QsqL.size(); k++){
 cout<<ASYML[k]<<"\t"<<ASYM_STL[k]<<"\t"<<QsqL[k]<<SensL[k]<<endl;
 Q2L->Fill(QsqL[k]);
 ThetaL->Fill(AngleL[k]);
 AsymL->Fill(ASYML[k]);
 SAsymL->Fill(ASYM_STL[k]);
 sensL->Fill(SensL[k]);
 }
 for(int k = 0; k < QsqR.size(); k++){
 Q2R->Fill(QsqR[k]);
 ThetaR->Fill(AngleR[k]);
 AsymR->Fill(ASYMR[k]);
 SAsymR->Fill(ASYM_STR[k]);
 sensR->Fill(SensR[k]);
 }

 c1->cd(2);
 AsymL->Draw("hist");
 AsymR->Draw("sames hist");
 AsymL->Scale(1./AsymL->Integral());
 AsymR->Scale(1./AsymR->Integral());
 gPad->Update();
 TPaveStats* statsAsL;
 TPaveStats* statsAsR;
 statsAsL = (TPaveStats*)ThetaL->FindObject("stats");
 statsAsR = (TPaveStats*)ThetaR->FindObject("stats");
 statsAsL->SetY2NDC(0.90);
 statsAsL->SetY1NDC(0.75);
 statsAsR->SetY2NDC(0.75);
 statsAsR->SetY1NDC(0.60);
 statsAsL->SetTextColor(1);
 statsAsR->SetTextColor(2);
 gPad->Modified();
 
 c1->cd(3);
 ThetaL->Draw("hist");
 ThetaR->Draw("sames hist");
 ThetaL->Scale(1./ThetaL->Integral());
 ThetaR->Scale(1./ThetaR->Integral());
 gPad->Update();
 TPaveStats* statsThL;
 TPaveStats* statsThR;
 statsThL = (TPaveStats*)ThetaL->FindObject("stats");
 statsThR = (TPaveStats*)ThetaR->FindObject("stats");
 statsThL->SetY2NDC(0.90);
 statsThL->SetY1NDC(0.75);
 statsThR->SetY2NDC(0.75);
 statsThR->SetY1NDC(0.60);
 statsThL->SetTextColor(1);
 statsThR->SetTextColor(2);
 gPad->Modified();

 c1->cd(4);
 Q2L->Draw("hist");
 Q2R->Draw("sames hist");
 Q2L->Scale(1./Q2L->Integral());
 Q2R->Scale(1./Q2R->Integral());
 gPad->Update();
 TPaveStats* statsQ2L;
 TPaveStats* statsQ2R;
 statsQ2L = (TPaveStats*)Q2L->FindObject("stats");
 statsQ2R = (TPaveStats*)Q2R->FindObject("stats");
 statsQ2L->SetY2NDC(0.90);
 statsQ2L->SetY1NDC(0.75);
 statsQ2R->SetY2NDC(0.75);
 statsQ2R->SetY1NDC(0.60);

 statsQ2L->SetTextColor(1);
 statsQ2R->SetTextColor(2);

 gPad->Modified();
 ofstream outfile("./TextFiles/calculated_qty.csv",ios_base::app);
 outfile<<runL<<AsymL->GetMean()<<"\t"<<AsymL->GetMeanError()<<"\t"<<ThetaL->GetMean()<<"\t"<<ThetaL->GetMeanError()<<"\t"<<Q2L->GetMean()<<"\t"<<Q2L->GetMeanError()<<runR<<AsymR->GetMean()<<"\t"<<AsymR->GetMeanError()<<"\t"<<ThetaR->GetMean()<<"\t"<<ThetaR->GetMeanError()<<"\t"<<Q2R->GetMean()<<"\t"<<Q2R->GetMeanError()<<endl;
 cout<<runL<<AsymL->GetMean()<<"\t"<<AsymL->GetMeanError()<<"\t"<<ThetaL->GetMean()<<"\t"<<ThetaL->GetMeanError()<<"\t"<<Q2L->GetMean()<<"\t"<<Q2L->GetMeanError()<<runR<<AsymR->GetMean()<<"\t"<<AsymR->GetMeanError()<<"\t"<<ThetaR->GetMean()<<"\t"<<ThetaR->GetMeanError()<<"\t"<<Q2R->GetMean()<<"\t"<<Q2R->GetMeanError()<<endl;

 c1->SaveAs("./temp1/calculated_qty_run%d_%d.pdf");

}
