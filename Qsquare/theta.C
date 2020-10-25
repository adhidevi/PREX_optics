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

  TChain* T = new TChain("T"); 
  if (run<10000){
  double upADCcut_approx = 480;
  T->Add(Form("/chafs1/work1/prex_counting/chandan/prexLHRS_%d_-1*.root",run));
     TH1F* huq = new TH1F("huq","P.upQadcL;ADC CH;Event/CH",250,400,650);
     T->Draw("P.upQadcL>>huq","","goff");

     TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
     huqCpy->GetXaxis()->SetRangeUser(upADCcut_approx-10,upADCcut_approx+10);
     cout<<"ADC cut: "<<huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin())<<endl;
     double upADCcut = huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin());

     sprintf(lang,"57.2958*TMath::ACos((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph)))>>hl",c0,s0);

     TCut lcut = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&L.gold.th>-0.08&&L.gold.th<0.08&&L.gold.ph>-0.05&&L.gold.ph<0.04&&L.gold.dp<-0.002&&P.upQadcL>%f)",upADCcut);
//     TCut lcut = Form("((P.evtypebits&2)==2&&L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1&&P.upQadcL>%f)",upADCcut);

     cout << "Left Plot command :"<<endl;
     cout << lang << endl;
     cout << "Left cut "<<lcut<<endl;

     TH1F *hl = new TH1F("hl",Form("Scat. Angle (degrees) L-HRS (run%d);#theta (degrees)",run),200,2,8);

     TCanvas* c1 = new TCanvas("c1","c1",1000,550);
     c1->Divide(2,1);;
     c1->cd(1);
     T->Draw(lang, lcut);
     c1->cd(2);

     TH2F *h2 = new TH2F("h2","tg_th vs tg_ph for L-HRS",100,-0.05,0.05,100,-0.08,0.08);
     T->Draw("L.gold.th:L.gold.ph>>h2","P.upQadcL>480","colz");

     c1->SaveAs(Form("./temp/theta_run%d.pdf",run));
  }else{
  double upADCcut_approx = 500;
  T->Add(Form("/chafs1/work1/prex_counting/chandan/prexRHRS_%d_-1*.root",run));

     TH1F* huq = new TH1F("huq", "P.upQadcR;ADC CH;Events/CH", 200, 450, 650);
     T->Draw("P.upQadcR>>huq");
     TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
     huqCpy->GetXaxis()->SetRangeUser(upADCcut_approx-10,upADCcut_approx+10);
     cout<<"ADC cut: "<<huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin())<<endl;
     double upADCcut = huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin());

    // Note sign of R.gold.ph is, relatively, negative here, in comparison with L-HRS.

     sprintf(rang,"57.2958*TMath::ACos((%f+R.gold.ph*%f)/(TMath::Sqrt(1.0+R.gold.th*R.gold.th+R.gold.ph*R.gold.ph)))>>hr",c0,s0);

     TCut rcut = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&R.gold.th>-0.08&&R.gold.th<0.08&&R.gold.ph>-0.05&&R.gold.ph<0.04&&L.gold.dp<-0.002&&P.upQadcR>%f)",upADCcut);
//     TCut rcut = Form("((P.evtypebits&2)==2&&R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1&&P.upQadcR>%f)",upADCcut);

     cout << "Right Plot command :"<<endl;
     cout << rang << endl;
     cout << "Right cut "<<rcut<<endl;

     TH1F *hr = new TH1F("hr",Form("Scat. Angle (degrees) R-HRS (run%d);#theta (degrees)",run),200,2,8);

     TCanvas* c1 = new TCanvas("c1","c1",1000,550);
     c1->Divide(2,1);
     c1->cd(1);;
     T->Draw(rang, rcut);
     c1->cd(2);
    TH2F *h2 = new TH2F("h2","tg_th vs tg_ph for R-HRS",100,-0.05,0.05,100,-0.08,0.08);
    T->Draw("R.gold.th:R.gold.ph>>h2","P.upQadcR>500","colz");
     c1->SaveAs(Form("./temp/theta_run%d.pdf",run));
  }
}


