void theta(Int_t run, TString target){
   gROOT->Reset();
   gStyle->SetOptStat(1);
   gStyle->SetStatH(0.3);
   gStyle->SetStatW(0.3);
   gStyle->SetTitleH(0.086);
   gStyle->SetTitleW(0.45);
   gStyle->SetLabelSize(0.04,"x");
   gStyle->SetLabelSize(0.04,"y");
   gROOT->ForceStyle();
   TGaxis::SetMaxDigits(3);

  // Central angle
  Float_t theta0 = 4.74;  // degrees

  theta0 = theta0 * 3.1415926/180;  // radians

  Float_t s0 = TMath::Sin(theta0);
  Float_t c0 = TMath::Cos(theta0);

  char lang[400], rang[400];
  char lcut[400], rcut[400];

  TChain* TA = new TChain("T"); 
  TChain* TS = new TChain("T"); 
  if (run<10000){
  double upADCcut_approx = 480;
  TA->Add(Form("/lustre/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run));
  TS->Add(Form("/lustre/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1*.root",run));
     TH1F* huq = new TH1F("huq","P.upQadcL;ADC CH;Event/CH",250,400,650);
     TA->Draw("P.upQadcL>>huq","","goff");

     TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
     huqCpy->GetXaxis()->SetRangeUser(upADCcut_approx-10,upADCcut_approx+10);
     cout<<"ADC cut: "<<huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin())<<endl;
     double upADCcut = huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin());

     TString lang = Form("57.2958*TMath::ACos((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph)))",c0,s0);

     TCut lcut = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.05&&P.upQadcL>%f)",upADCcut);

     cout << "Left Plot command :"<<endl;
     cout << lang << endl;
     cout << "Left cut "<<lcut<<endl;

     TH1F *hlA = new TH1F("hlA",Form("Scat. Angle (degrees) L-HRS (run%d);#theta (degrees)",run),200,2,8);
     TH1F *hlS = new TH1F("hlS",Form("Scat. Angle (degrees) L-HRS (run%d);#theta (degrees)",run),200,2,8);
     hlA->SetLineColor(1);
     hlS->SetLineColor(2);

     TCanvas* c1 = new TCanvas("c1","c1",1000,550);
     c1->Divide(2,1);;
     c1->cd(1);
     TA->Draw(Form("%s>>hlA",lang.Data()), lcut);
     TS->Draw(Form("%s>>hlS",lang.Data()), lcut, "sames");

     gPad->Update();
     TPaveStats* statA = (TPaveStats*)hlA->FindObject("stats");
     TPaveStats* statS = (TPaveStats*)hlS->FindObject("stats");
     statA->SetY2NDC(0.90);
     statA->SetY1NDC(0.75);
     statS->SetY2NDC(0.75);
     statS->SetY1NDC(0.60);
     statA->SetTextColor(1);
     statS->SetTextColor(2);
     gPad->Modified();
     
     TLatex latex;
     latex.SetNDC();
     latex.SetTextSize(0.05);
     latex.SetTextColor(1);
     latex.DrawLatex(0.70,0.55,"DB1");
     latex.SetTextColor(2);
     latex.DrawLatex(0.70,0.50,"DB2");

     c1->cd(2);
     TH2F *h2A = new TH2F("h2A","tg_th vs tg_ph for L-HRS",100,-0.05,0.05,100,-0.08,0.08);
     TH2F *h2S = new TH2F("h2S","tg_th vs tg_ph for L-HRS",100,-0.05,0.05,100,-0.08,0.08);
     h2A->SetMarkerColor(1);
     h2S->SetMarkerColor(2);
     TA->Draw("L.gold.th:L.gold.ph>>h2A",Form("P.upQadcL>%f",upADCcut));
     TS->Draw("L.gold.th:L.gold.ph>>h2S",Form("P.upQadcL>%f",upADCcut),"same");

     c1->SaveAs(Form("./temp/theta_run%d.pdf",run));
  }else{
  double upADCcut_approx = 500;
  TA->Add(Form("/lustre/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",run));
  TS->Add(Form("/lustre/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_%d_-1*.root",run));

     TH1F* huq = new TH1F("huq", "P.upQadcR;ADC CH;Events/CH", 200, 450, 650);
     TA->Draw("P.upQadcR>>huq");
     TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
     huqCpy->GetXaxis()->SetRangeUser(upADCcut_approx-10,upADCcut_approx+10);
     cout<<"ADC cut: "<<huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin())<<endl;
     double upADCcut = huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin());

    // Note sign of R.gold.ph is, relatively, negative here, in comparison with L-HRS.

     TString rang = Form("57.2958*TMath::ACos((%f+R.gold.ph*%f)/(TMath::Sqrt(1.0+R.gold.th*R.gold.th+R.gold.ph*R.gold.ph)))",c0,s0);

     TCut rcut = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.05&&P.upQadcR>%f)",upADCcut);

     cout << "Right Plot command :"<<endl;
     cout << rang << endl;
     cout << "Right cut "<<rcut<<endl;

     TH1F *hrA = new TH1F("hrA",Form("Scat. Angle (degrees) R-HRS (run%d);#theta (degrees)",run),200,2,8);
     TH1F *hrS = new TH1F("hrS",Form("Scat. Angle (degrees) R-HRS (run%d);#theta (degrees)",run),200,2,8);
     hrA->SetLineColor(1);
     hrS->SetLineColor(2);

     TCanvas* c1 = new TCanvas("c1","c1",1000,550);
     c1->Divide(2,1);
     c1->cd(1);
     TA->Draw(Form("%s>>hrA",rang.Data()), rcut);
     TS->Draw(Form("%s>>hrS",rang.Data()), rcut,"sames");

     gPad->Update();
     TPaveStats* statA = (TPaveStats*)hrA->FindObject("stats");
     TPaveStats* statS = (TPaveStats*)hrS->FindObject("stats");
     statA->SetY2NDC(0.90);
     statA->SetY1NDC(0.75);
     statS->SetY2NDC(0.75);
     statS->SetY1NDC(0.60);
     statA->SetTextColor(1);
     statS->SetTextColor(2);
     gPad->Modified();

     TLatex latex;
     latex.SetNDC();
     latex.SetTextSize(0.05);
     latex.SetTextColor(1);
     latex.DrawLatex(0.70,0.55,"DB1");
     latex.SetTextColor(2);
     latex.DrawLatex(0.70,0.50,"DB2");

     c1->cd(2);
    TH2F *h2A = new TH2F("h2A","tg_th vs tg_ph for R-HRS",100,-0.05,0.05,100,-0.08,0.08);
    TH2F *h2S = new TH2F("h2S","tg_th vs tg_ph for R-HRS",100,-0.05,0.05,100,-0.08,0.08);
    h2A->SetMarkerColor(1);
    h2S->SetMarkerColor(2);
    TA->Draw("R.gold.th:R.gold.ph>>h2A",Form("P.upQadcR>%f",upADCcut));
    TS->Draw("R.gold.th:R.gold.ph>>h2S",Form("P.upQadcR>%f",upADCcut),"same");
     c1->SaveAs(Form("./temp/theta_run%d.pdf",run));
  }
}


