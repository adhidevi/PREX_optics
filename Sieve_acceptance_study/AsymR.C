#include "CollimatorL.C"
#include "UpPlane.C"
#include "DownPlane.C"
#include "CREXdata.h"

//Author Ryan Richards 
//Shift polar angle in sim to match data and compute Q^2, Asym

//Run number, sept scan(0.0 - nominal), ADC, central angle for data, dp offset offset in degrees(1 clearly no shift)
void AsymR(int run,double sept,double ADC,double th0, double dp, double off){

LoadTable("horpb.dat",0);//Lead table - sufficient

double d2r = TMath::Pi()/180;
double r2d = 1/d2r;

//th0 is for data 
double cth0 = TMath::Cos(th0*d2r); double sth0 = TMath::Sin(th0*d2r);

double cthsim = TMath::Cos((4.74/off)*d2r);
double sthsim = TMath::Sin((4.74/off)*d2r);


TH1F *lab[2],*qsq[2], *asym[2];
int color[2] = { 2, 4 };

for(int i = 0; i < 2; i++){ 

   if(i == 0){ qsq[i] = new TH1F(Form("qsq[%i]",i),"RHRS Apparent Q^{2} Data",150,0.0,0.015);
                lab[i] = new TH1F(Form("lab[%i]",i),"RHRS Apparent #theta_{lab} Data",150,1.0,9.0);
                asym[i] = new TH1F(Form("asym[%i]",i),"RHRS Apparent Asymmetry Data",150,0,1);
    } else{
    qsq[i] = new TH1F(Form("qsq[%i]",i),Form("RHRS Apparent Q^{2} Septum Detuned %.1f",sept),150,0.0,0.015);
    lab[i] = new TH1F(Form("lab[%i]",i),Form("RHRS Apparent #theta_{lab} Septum Detuned %.1f",sept),150,1.0,9.0);
    asym[i] = new TH1F(Form("asym[%i]",i),Form("RHRS Apparent Asymmetry Septum Detuend %.1f",sept),150,0,1);

    }

   asym[i]->SetLineColor(color[i]); qsq[i]->SetLineColor(color[i]); lab[i]->SetLineColor(color[i]); 

}

//For simulation
vector <double> Rate,Qsq,Asym,Lab;
double thisRate, thisTh, thisPh, thisQsq, thisXcol, thisYcol, thisXvdc, thisAsym, thisAng,thisCosAng;
double thisE,thisP, thismom;
double thisxu1,thisxu2,thisyu1,thisyu2;
double thisxd1,thisxd2,thisxd3,thisxd4,thisxd5,thisxd6,thisxd7,thisxd8,thisxd9;
double thisyd1,thisyd2,thisyd3,thisyd4,thisyd5,thisyd6,thisyd7,thisyd8,thisyd9;

//For data, need cluster, adc,th,ph,p
vector <double> Angle,Qsqdat, Asymdat;
double thisu1,thisv1,thisu2,thisv2,thisPdat,thisThdat,thisPhdat,thisQsqdat,thisAsymdat,thisCosAngdat,thisAngdat,thisADC;
double thisdPdat, thisEvt;

  TChain *T = new TChain("T");
  T->Add(Form("../prex2Rootfiles/prexRHRS_%i_-1.root",run));


  TChain *T1 = new TChain("T");
  for(int i = 1; i < 6; i++){ T1->Add(Form("../sandwich/SandwichRHRS_PREX_%.1f_%i.root",sept,i)); }

  T->SetBranchAddress("R.vdc.u1.nclust",&thisu1); T->SetBranchAddress("R.vdc.v1.nclust",&thisv1);
  T->SetBranchAddress("R.vdc.u2.nclust",&thisu2); T->SetBranchAddress("R.vdc.v2.nclust",&thisv2);
  T->SetBranchAddress("R.gold.th",&thisThdat); T->SetBranchAddress("R.gold.ph",&thisPhdat);
  T->SetBranchAddress("R.gold.p",&thisPdat); T->SetBranchAddress("P.upQadcR",&thisADC);
  T->SetBranchAddress("R.gold.dp",&thisdPdat);T->SetBranchAddress("P.evtypebits",&thisEvt);

  long n = T->GetEntries();
  
  //Looping over tree in data
  for(long i = 0; i < n; i++){
     T->GetEntry(i);
   
     int thisevent = (int) thisEvt;
      
     if(thisu1 == 1 && thisv1 == 1 && thisu2 == 1 && thisv2 == 1 && thisADC > ADC && thisThdat > -0.08 && thisThdat<0.08 && thisPhdat > -0.05 && thisPhdat < 0.05 && thisdPdat > -0.04 && thisdPdat < -0.002 && ((thisevent&2)==2)  ){ 

       thisCosAngdat = (cth0+thisPhdat*sth0)/(TMath::Sqrt(1+thisThdat*thisThdat+thisPhdat*thisPhdat));
       thisAngdat = r2d*TMath::ACos(thisCosAngdat);
       //Using a hardcoded beam energy -- run dependent, yes?
       thisQsq = 2*0.951*thisPdat*(1-thisCosAngdat);
       thisAsymdat = 1e6*Interpolate(thisPdat*1000,thisAngdat,0,1);

       Angle.push_back(thisAngdat);
       Qsqdat.push_back(thisQsq);     
       Asymdat.push_back(thisAsymdat);
   
    
   }

  } 

 
 T1->SetBranchAddress("rate",&thisRate); T1->SetBranchAddress("x_vdc_tr",&thisXvdc);
 T1->SetBranchAddress("x_col_tr",&thisXcol); T1->SetBranchAddress("y_col_tr",&thisYcol);
 T1->SetBranchAddress("th_ztarg_tr",&thisTh); T1->SetBranchAddress("ph_ztarg_tr",&thisPh);
 T1->SetBranchAddress("ev.beamp",&thisP); T1->SetBranchAddress("ev.ep",&thisE);
 T1->SetBranchAddress("p_ztarg_tr",&thismom);

 T1->SetBranchAddress("x_zup1",&thisxu1); T1->SetBranchAddress("y_zup1",&thisyu1);
 T1->SetBranchAddress("x_zup2",&thisxu2); T1->SetBranchAddress("y_zup2",&thisyu2);
 T1->SetBranchAddress("x_zdown1",&thisxd1); T1->SetBranchAddress("y_zdown1",&thisyd1);
 T1->SetBranchAddress("x_zdown2",&thisxd2); T1->SetBranchAddress("y_zdown2",&thisyd2);
 T1->SetBranchAddress("x_zdown3",&thisxd3); T1->SetBranchAddress("y_zdown3",&thisyd3);
 T1->SetBranchAddress("x_zdown4",&thisxd4); T1->SetBranchAddress("y_zdown4",&thisyd4);
 T1->SetBranchAddress("x_zdown5",&thisxd5); T1->SetBranchAddress("y_zdown5",&thisyd5);
 T1->SetBranchAddress("x_zdown6",&thisxd6); T1->SetBranchAddress("y_zdown6",&thisyd6);
 T1->SetBranchAddress("x_zdown7",&thisxd7); T1->SetBranchAddress("y_zdown7",&thisyd7);
 T1->SetBranchAddress("x_zdown8",&thisxd8); T1->SetBranchAddress("y_zdown8",&thisyd8);
 T1->SetBranchAddress("x_zdown9",&thisxd9); T1->SetBranchAddress("y_zdown9",&thisyd9);
 

 int m = T1->GetEntries();

 //Now loop over sim tree
 for(int j = 0; j < m; j++){ 
 T1->GetEntry(j);


  //This has momentum cut -- needs to be evaluated
  if( thisXvdc!=-333. && CollimatorR(thisXcol,thisYcol) && UpPlane(thisxu1,thisyu1,thisxu2,thisyu2,0) &&
  DownPlane(thisxd1,thisyd1,thisxd2,thisyd2,thisxd3,thisyd3,thisxd4,thisyd4,thisxd5,thisyd5,thisxd6,thisyd6,thisxd7,thisyd7,thisxd8,thisyd8,thisxd9,thisyd9,0)&& ((thismom-50*thisTh) > dp) ){
  thisCosAng = (cthsim+thisPh*sthsim)/TMath::Sqrt(1+thisPh*thisPh+thisTh*thisTh);  
  
  //We will have to be careful here  
  thisAng = r2d*TMath::ACos(thisCosAng);
  //This is the Angle in degrees, now adding offset in degrees - difference in central value btwn hists, add/subtract event by event
  thisAng = thisAng + off;
  //Now with this, compute Q^2 and Asymmetry
  thisQsq = 2*thisP*thisE*(1-TMath::Cos(thisAng*d2r));
  thisAsym = 1e6*Interpolate(thisP*1000,thisAng,0,1);

  Rate.push_back(thisRate); Qsq.push_back(thisQsq);
  Lab.push_back(thisAng); Asym.push_back(thisAsym);

  }



 } 

 //Now fill histograms
 //data
 for(int l = 0; l < Qsqdat.size(); l++){ qsq[0]->Fill(Qsqdat[l]); lab[0]->Fill(Angle[l]); asym[0]->Fill(Asymdat[l]); }
 //simulation
 for(int k = 0; k < Qsq.size(); k++){ qsq[1]->Fill(Qsq[k],Rate[k]); lab[1]->Fill(Lab[k],Rate[k]); asym[1]->Fill(Asym[k],Rate[k]); }

  auto leg = new TLegend(0.1,0.75,0.4,0.9);
  leg->AddEntry(lab[0],"Data","l");
  leg->AddEntry(lab[1],Form("Sim Sept Detuned %.1f Per",sept),"l");
 



 TCanvas *c1 = new TCanvas();
 lab[0]->Draw();
 lab[1]->Draw("HIST same");
 lab[1]->Scale(lab[0]->Integral()/lab[1]->Integral());
 leg->Draw();

 TCanvas *c2 = new TCanvas();
 qsq[0]->Draw(); 
 qsq[1]->Draw("HIST same");
 qsq[1]->Scale(qsq[0]->Integral()/qsq[1]->Integral());
 leg->Draw();

 TCanvas *c3 = new TCanvas();
 asym[0]->Draw();
 asym[1]->Draw("HIST same");
 asym[1]->Scale(asym[0]->Integral()/asym[1]->Integral());
 leg->Draw();


  std::cout << "Data Lab theta mean: " << lab[0]->GetMean() << "  " << "Data Lab theta rms: " << lab[0]->GetRMS() << std::endl;
  std::cout << "Sim Lab theta mean: " << lab[1]->GetMean() << "   " << "Sim Lab theta rms: " << lab[1]->GetRMS() << std::endl;
  std::cout << "Data  Qsq mean: " << qsq[0]->GetMean() << "  " << "Data  Qsq rms: " << qsq[0]->GetRMS() << std::endl;
  std::cout << "Sim  Qsq mean: " << qsq[1]->GetMean() << "   " << "Sim  Qsq rms: " << qsq[1]->GetRMS() << std::endl;
  std::cout << "Data Asymmetry  mean: " << asym[0]->GetMean() << "  " << "Data Asymmetry  rms: " << asym[0]->GetRMS() << std::endl;
  std::cout << "Sim Asymmetry  mean: " << asym[1]->GetMean() << "   " << "Sim Asymmetry  rms: " << asym[1]->GetRMS() << std::endl;








 }
