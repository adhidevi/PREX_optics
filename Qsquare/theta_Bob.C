void theta_Bob(Int_t run, Int_t which){
   gROOT->Reset();
   gStyle->SetOptStat(1);
   gStyle->SetStatH(0.3);
   gStyle->SetStatW(0.3);
   gStyle->SetTitleH(0.086);
   gStyle->SetTitleW(0.45);
   gStyle->SetLabelSize(0.05,"x");
   gStyle->SetLabelSize(0.05,"y");
   gROOT->ForceStyle();

  TFile* file = new TFile(Form("/lustre/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-1.root",run));
  TTree* T = (TTree*)file->Get("T");
  // Central angle
  Float_t theta0 = 4.74;  // degrees

  theta0 = theta0 * 3.1415926/180;  // radians

  Float_t s0 = TMath::Sin(theta0);
  Float_t c0 = TMath::Cos(theta0);

//  Int_t which = 1;   // 1 = Left,  2 = Right,  

  char lang[400], rang[400];
  char lcut[400], rcut[400];


  if (which == 1) {

     sprintf(lang,"57.2958*TMath::ACos((%f-L.gold.ph*%f)/(TMath::Sqrt(1.0+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph)))>>hl",c0,s0);

     sprintf(lcut,"abs(L.gold.th)<0.06&&abs(L.gold.ph)<0.03&&P.upQadcL>480");

     cout << "Left Plot command :"<<endl;
     cout << lang << endl;
     cout << "Left cut "<<lcut<<endl;

     //     TH1F *hl = new TH1F("hl","Scattering Angle (radians) L-HRS",200,0,0.2);
     TH1F *hl = new TH1F("hl","Scattering Angle (degrees) L-HRS",200,2,8);

     TCanvas c1;

     T->Draw(lang, lcut);
     c1.SaveAs(Form("./plots/theta_run%d.pdf",run));
  }

  if (which == 2) {
 
    // Note sign of R.gold.ph is, relatively, negative here, in comparison with L-HRS.

      sprintf(rang,"57.2958*TMath::ACos((%f+R.gold.ph*%f)/(TMath::Sqrt(1.0+R.gold.th*R.gold.th+R.gold.ph*R.gold.ph)))>>hr",c0,s0);

      sprintf(rcut,"abs(R.gold.th)<0.06&&abs(R.gold.ph)<0.03&&P.upQadcR>500");

      cout << "Right Plot command :"<<endl;
      cout << rang << endl;
      cout << "Right cut "<<rcut<<endl;

      //      TH1F *hr = new TH1F("hr","Angimuthal angle (radians) R-HRS",200,0,0.2);
     TH1F *hr = new TH1F("hr","Scattering Angle (degrees) R-HRS",200,2,8);

      TCanvas c1;

      T->Draw(rang, rcut);

  }
 
  if (which == 3) {
    TH2F *h2 = new TH2F("h2","tg_th vs tg_ph for L-HRS",100,-0.05,0.05,100,-0.08,0.08);
    T->Draw("L.gold.th:L.gold.ph>>h2","P.upQadcL>480");
  }

  if (which == 4) {
    TH2F *h2 = new TH2F("h2","tg_th vs tg_ph for R-HRS",100,-0.05,0.05,100,-0.08,0.08);
    T->Draw("R.gold.th:R.gold.ph>>h2","P.upQadcR>500");
  }

}


