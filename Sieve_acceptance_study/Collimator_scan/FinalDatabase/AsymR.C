#include "CollimatorR.C"
#include "UpPlane.C"
#include "DownPlane.C"
#include "CREXdata.h"

//Author Ryan Richards 
//Shift polar angle in sim to match data and compute Q^2, Asym

//Run number, sept scan(0.0 - nominal), ADC, central angle for data, dp cut off,  offset in degrees(1 clearly no shift)

void AsymR(int run,double sept,double ADC,double th0, double peakpos, double simMomCut, double Ebeam, double YColl_offset){

double efact = (Ebeam-0.001)/peakpos;
//double efactsim = Ebeam*1000/951.0;
double efactsim = 1;
char targpos[3][20] = { "-5mm","0mm","+5mm" }; 

LoadTable("horpb.dat",0);//Lead table - sufficient
int ngdat = 0;
int ngsim = 0;

double d2r = TMath::Pi()/180;
double r2d = 1/d2r;

//th0 is for data 
double cth0 = TMath::Cos(th0*d2r); double sth0 = TMath::Sin(th0*d2r);

double cthsim = TMath::Cos((4.74)*d2r);
double sthsim = TMath::Sin((4.74)*d2r);

TH1F *lab[2],*qsq[2], *asym[2], *hmom[2];
int color[2] = { 2, 4 };

for(int i = 0; i < 2; i++){ 

   if(i == 0){ qsq[i] = new TH1F(Form("qsq[%i]",i),"RHRS Apparent Q^{2} Data",150,0.0,0.015);
                lab[i] = new TH1F(Form("lab[%i]",i),"RHRS Apparent #theta_{lab} Data",150,1.0,9.0);
                asym[i] = new TH1F(Form("asym[%i]",i),"RHRS Apparent Asymmetry Data",150,0.0,1.0);
                hmom[i] = new TH1F(Form("hmom[%i]",i),"RHRS Apparent Asymmetry Data",150,0.938,0.953);
    } else{
    qsq[i] = new TH1F(Form("qsq[%i]",i),Form("RHRS Apparent Q^{2} Septum Detuned %.1f",sept),150,0.0,0.015);
    lab[i] = new TH1F(Form("lab[%i]",i),Form("RHRS Apparent #theta_{lab} Septum Detuned %.1f",sept),150,1.0,9.0);
    asym[i] = new TH1F(Form("asym[%i]",i),Form("RHRS Apparent Asymmetry Septum Detuned %.1f",sept),150,0.0,1.0);
    hmom[i] = new TH1F(Form("hmom[%i]",i),Form("RHRS Apparent Momentum Septum Detuned %.1f",sept),150,0.938,0.953);
    }

   asym[i]->SetLineColor(color[i]); qsq[i]->SetLineColor(color[i]); lab[i]->SetLineColor(color[i]); hmom[i]->SetLineColor(color[i]);

}
// simulation vertex distribution

TH1F *qsq_v = new TH1F("qsq_v",Form("RHRS sim vertex Q^{2} Septum Detuned %.1f",sept),150,0.0,0.015);
TH1F *lab_v = new TH1F("lab_v",Form("RHRS sim vertex #theta_{lab} Septum Detuned %.1f",sept),150,1.0,9.0);
TH1F *asym_v = new TH1F("asym_v",Form("RHRS sim vertex Asymmetry Septum Detuend %.1f",sept),150,0.0,1.0);

//For simulation
vector <double> Rate,Qsq,Asym,Lab, Mom;
vector <double> Qsq_v, Asym_v, Lab_v, Rate_v;
double thisRate, thisTh, thisPh, thisQsq, thisXcol, thisYcol, thisXvdc, thisAsym, thisAng,thisCosAng;
double thisE,thisP, thismom;
double thisxu1,thisxu2,thisyu1,thisyu2;
double thisxd1,thisxd2,thisxd3,thisxd4,thisxd5,thisxd6,thisxd7,thisxd8,thisxd9;
double thisyd1,thisyd2,thisyd3,thisyd4,thisyd5,thisyd6,thisyd7,thisyd8,thisyd9;
double thisAng1;
double thisQsq_v, thisAsym_v, thisLab_v;
int thisnucleus;
//For data, need cluster, adc,th,ph,p
vector <double> Angle,Qsqdat, Asymdat, Momdat;
double thisu1,thisv1,thisu2,thisv2,thisPdat, thisdPdat, thisThdat,thisPhdat,thisQsqdat,thisAsymdat,thisCosAngdat,thisAngdat,thisADC;
double thisEvt;

  TChain *T = new TChain("T");
  T->Add(Form("/chafs1/work1/prex_counting/Qsq_Oct22/prexRHRS_%i_-1*.root",run));

  TChain *T1 = new TChain("T");
  for(int i = 1; i < 6; i++){ T1->Add(Form("/lustre/expphy/volatile/halla/parity/ryanrich/Ebeam953_4by6/TargNominal/Zero1_SandwichRHRS_PREX_%.1f_%i.root",sept,i)); }

  T->SetBranchAddress("R.vdc.u1.nclust",&thisu1); T->SetBranchAddress("R.vdc.v1.nclust",&thisv1);
  T->SetBranchAddress("R.vdc.u2.nclust",&thisu2); T->SetBranchAddress("R.vdc.v2.nclust",&thisv2);
  T->SetBranchAddress("R.gold.th",&thisThdat); T->SetBranchAddress("R.gold.ph",&thisPhdat);
  T->SetBranchAddress("R.gold.p",&thisPdat); T->SetBranchAddress("P.upQadcR",&thisADC);
  T->SetBranchAddress("R.gold.dp",&thisdPdat); T->SetBranchAddress("P.evtypebits",&thisEvt);
  
  long n = T->GetEntries();
  
  //Looping over tree in data
  for(long i = 0; i < n; i++){
     T->GetEntry(i);
   thisPdat *=efact;  

   int thisevent = (int) thisEvt;
  
   if(thisu1 == 1 && thisv1 == 1 && thisu2 == 1 && thisv2 == 1 && thisADC > ADC && thisThdat > -0.08 && thisThdat<0.08 && thisPhdat > -0.05 && thisPhdat < 0.05 && thisdPdat > -0.04 && thisdPdat < 0.002 && ((thisevent&2)==2)  ){
  
      thisCosAngdat = (cth0+thisPhdat*sth0)/(TMath::Sqrt(1+thisThdat*thisThdat+thisPhdat*thisPhdat));
      thisAngdat = r2d*TMath::ACos(thisCosAngdat);
     
      //Using a hardcoded beam energy -- run dependent, yes?
       thisQsqdat = 2*Ebeam*thisPdat*(1-thisCosAngdat);
       thisAsymdat = 1e6*Interpolate(thisPdat*1000,thisAngdat,0,1);
  
   //    std::cout <<  thisAngdat << "  " << thisQsqdat << "  " << thisAsymdat << std::endl; 

       Angle.push_back(thisAngdat);
       Qsqdat.push_back(thisQsqdat);     
       Asymdat.push_back(thisAsymdat);
       Momdat.push_back(thisPdat);

  } }
 
 T1->SetBranchAddress("rate",&thisRate); T1->SetBranchAddress("x_vdc_tr",&thisXvdc);
 T1->SetBranchAddress("x_col_tr",&thisXcol); T1->SetBranchAddress("y_col_tr",&thisYcol);
 T1->SetBranchAddress("th_ztarg_tr",&thisTh); T1->SetBranchAddress("ph_ztarg_tr",&thisPh);
 T1->SetBranchAddress("ev.beamp",&thisP); T1->SetBranchAddress("ev.ep",&thisE);
 T1->SetBranchAddress("p_ztarg_tr",&thismom); T1->SetBranchAddress("ev.nuclA",&thisnucleus);

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

  //vertex asym, Q2 and theta
  T1->SetBranchAddress("ev.A",&thisAsym_v);
  T1->SetBranchAddress("ev.Q2",&thisQsq_v);
  T1->SetBranchAddress("ev.Th",&thisLab_v); //deg

 int m = T1->GetEntries();

 //Now loop over sim tree
 for(int j = 0; j < m; j++){ 
 T1->GetEntry(j);

  thismom *=efactsim;
  //This has momentum cut -- needs to be evaluated
  if( thisXvdc!=-333. && CollimatorR(thisXcol,thisYcol+YColl_offset) && UpPlane(thisxu1,thisyu1,thisxu2,thisyu2,0) &&
  DownPlane(thisxd1,thisyd1,thisxd2,thisyd2,thisxd3,thisyd3,thisxd4,thisyd4,thisxd5,thisyd5,thisxd6,thisyd6,thisxd7,thisyd7,thisxd8,thisyd8,thisxd9,thisyd9,0)&& ((thismom -50*thisTh) > ((Ebeam - 0.001 - simMomCut)*1000))){

 thisCosAng = (cthsim+thisPh*sthsim)/TMath::Sqrt(1+thisPh*thisPh+thisTh*thisTh);  
  
  //We will have to be careful here  
  thisAng = r2d*TMath::ACos(thisCosAng);
  //This is the Angle in degrees, now adding offset in degrees - difference in central value btwn hists, add/subtract event by event
  thisAng1 = thisAng;

  //Now with this, compute Q^2 and Asymmetry
  thisQsq = 2*thismom/1000.*Ebeam*(1-TMath::Cos(thisAng1*d2r));
  thisAsym = 1e6*Interpolate(thismom,thisAng1,0,1);

  Rate.push_back(thisRate); Qsq.push_back(thisQsq);
  Lab.push_back(thisAng1); Asym.push_back(thisAsym);
  Mom.push_back(thismom);
  if(thisnucleus==208){
  Qsq_v.push_back(thisQsq_v); Asym_v.push_back(thisAsym_v/1000);
  Lab_v.push_back(thisLab_v); Rate_v.push_back(thisRate);
  }
  }
 } 

 //Now fill histograms
 //data
 for(int l = 0; l < Qsqdat.size(); l++){ qsq[0]->Fill(Qsqdat[l]); lab[0]->Fill(Angle[l]); asym[0]->Fill(Asymdat[l]); hmom[0]->Fill(Momdat[l]); }
 //simulation
 for(int k = 0; k < Qsq.size(); k++){ qsq[1]->Fill(Qsq[k],Rate[k]); lab[1]->Fill(Lab[k],Rate[k]); asym[1]->Fill(Asym[k],Rate[k]); hmom[1]->Fill(Mom[k],Rate[k]);
  }
 for(int vi = 0; vi < Qsq_v.size(); vi++){
 qsq_v->Fill(Qsq_v[vi],Rate_v[vi]); asym_v->Fill(Asym_v[vi],Rate_v[vi]); lab_v->Fill(Lab_v[vi],Rate_v[vi]);
  }

  auto leg = new TLegend(0.1,0.75,0.4,0.9);
  leg->AddEntry(lab[0],Form("Data (run%d)",run),"l");
  leg->AddEntry(lab[1],Form("Sim, Sept %.1f Per, TargPos %s",sept,targpos[1]),"l");

 TCanvas *c1 = new TCanvas();
 lab[0]->Draw();
 lab[1]->Draw("HIST sames");
 lab[1]->Scale(lab[0]->Integral()/lab[1]->Integral());
 gPad->Update();
 TPaveStats* statDa = (TPaveStats*)lab[0]->FindObject("stats");
 TPaveStats* statSa = (TPaveStats*)lab[1]->FindObject("stats");
 statDa->SetY1NDC(0.90);
 statDa->SetY2NDC(0.75);
 statSa->SetY1NDC(0.75);
 statSa->SetY2NDC(0.60);
 statDa->SetTextColor(kRed);
 statSa->SetTextColor(kBlue);
 gPad->Modified();
 leg->Draw();
 c1->SaveAs(Form("./temp/%s_PolarAngle_%.1fper.pdf",targpos[1],sept));

 TCanvas *c2 = new TCanvas();
 qsq[0]->Draw(); 
 qsq[1]->Draw("HIST sames");
 qsq[1]->Scale(qsq[0]->Integral()/qsq[1]->Integral());
 gPad->Update();
 TPaveStats* statDq = (TPaveStats*)qsq[0]->FindObject("stats");
 TPaveStats* statSq = (TPaveStats*)qsq[1]->FindObject("stats");
 statDq->SetY1NDC(0.90);
 statDq->SetY2NDC(0.75);
 statSq->SetY1NDC(0.75);
 statSq->SetY2NDC(0.60);
 statDq->SetTextColor(kRed);
 statSq->SetTextColor(kBlue);
 gPad->Modified(); 
 leg->Draw();
 c2->SaveAs(Form("./temp/%s_Qsq_%.1fper.pdf",targpos[1],sept));

 TCanvas *c3 = new TCanvas();
 asym[0]->Draw();
 asym[1]->Draw("HIST sames");
 asym[1]->Scale(asym[0]->Integral()/asym[1]->Integral());
 gPad->Update();
 TPaveStats* statDas = (TPaveStats*)asym[0]->FindObject("stats");
 TPaveStats* statSas = (TPaveStats*)asym[1]->FindObject("stats");
 statDas->SetY1NDC(0.90);
 statDas->SetY2NDC(0.75);
 statSas->SetY1NDC(0.75);
 statSas->SetY2NDC(0.60);
 statDas->SetTextColor(kRed);
 statSas->SetTextColor(kBlue);
 gPad->Modified();
 leg->Draw();
 c3->SaveAs(Form("./temp/%s_Asym_%.1fper.pdf",targpos[1],sept));

 //You can look at this for different target positions
 TCanvas *c4 = new TCanvas();
 c4->Divide(2,2);
 c4->cd(1);
 lab_v->Draw("HIST");
 c4->cd(2);
 qsq_v->Draw("HIST");
 c4->cd(3);
 asym_v->Draw("HIST");
 c4->SaveAs(Form("./temp/%s_Vertex_%.1fper.pdf",targpos[1],sept)); 

  std::cout << "Data Lab theta mean: " << lab[0]->GetMean() << "  " << "Data Lab theta rms: " << lab[0]->GetRMS() << std::endl;
  std::cout << "Sim Lab theta mean: " << lab[1]->GetMean() << "   " << "Sim Lab theta rms: " << lab[1]->GetRMS() << std::endl;
  std::cout << "Data  Qsq mean: " << qsq[0]->GetMean() << "  " << "Data  Qsq rms: " << qsq[0]->GetRMS() << std::endl;
  std::cout << "Sim  Qsq mean: " << qsq[1]->GetMean() << "   " << "Sim  Qsq rms: " << qsq[1]->GetRMS() << std::endl;
  std::cout << "Data Asymmetry  mean: " << asym[0]->GetMean() << "  " << "Data Asymmetry  rms: " << asym[0]->GetRMS() << std::endl;
  std::cout << "Sim Asymmetry  mean: " << asym[1]->GetMean() << "   " << "Sim Asymmetry  rms: " << asym[1]->GetRMS() << std::endl;

  std::cout<<"-------- vertex ----------"<<std::endl;
  std::cout << "Sim Lab theta mean: " << lab_v->GetMean() << "   " << "Sim Lab theta rms: " << lab_v->GetRMS() << std::endl;
  std::cout << "Sim  Qsq mean: " << qsq_v->GetMean() << "   " << "Sim  Qsq rms: " << qsq_v->GetRMS() << std::endl;
  std::cout << "Sim Asymmetry  mean: " << asym_v->GetMean() << "   " << "Sim Asymmetry  rms: " << asym_v->GetRMS() << std::endl;

 ofstream outfile(Form("./TextFiles/finalDB_953MeVsimulation_app_vs_data_sept_scan_Ycoll_scan_run%d.csv",run),ios_base::app);
// ofstream outfile(Form("./TextFiles/scattering_angle_scan_run%d.csv",run),ios_base::app);

 outfile<<sept<<"\t"<<YColl_offset<<"\t"<<lab[0]->GetMean()<<"\t"<<lab[1]->GetMean()<<"\t"<<lab[1]->GetMean()/lab[0]->GetMean()*100<<"\t"<<qsq[0]->GetMean()<<"\t"<<qsq[1]->GetMean()<<"\t"<<qsq[1]->GetMean()/qsq[0]->GetMean()*100<<"\t"<<asym[0]->GetMean()<<"\t"<<asym[1]->GetMean()<<"\t"<<asym[1]->GetMean()/asym[0]->GetMean()*100<<"\t"<<lab[0]->GetRMS()<<"\t"<<lab[1]->GetRMS()<<"\t"<<lab[1]->GetRMS()/lab[0]->GetRMS()*100<<"\t"<<qsq[0]->GetRMS()<<"\t"<<qsq[1]->GetRMS()<<"\t"<<qsq[1]->GetRMS()/qsq[0]->GetRMS()*100<<"\t"<<asym[0]->GetRMS()<<"\t"<<asym[1]->GetRMS()<<"\t"<<asym[1]->GetRMS()/asym[0]->GetRMS()*100<<endl;
 outfile.close();

 ofstream outfile_v(Form("./TextFiles/finalDB_953MeVsimulation_ver_vs_app_sept_scan_Ycoll_scan_run%d.csv",run),ios_base::app);
// ofstream outfile_v(Form("./TextFiles/vertex_scattering_angle_scan_run%d.csv",run),ios_base::app);

 outfile_v<<sept<<"\t"<<YColl_offset<<"\t"<<lab_v->GetMean()<<"\t"<<lab[1]->GetMean()<<"\t"<<lab_v->GetMean()/lab[1]->GetMean()*100.<<"\t"<<qsq_v->GetMean()<<"\t"<<qsq[1]->GetMean()<<"\t"<<qsq_v->GetMean()/qsq[1]->GetMean()*100.<<"\t"<<asym_v->GetMean()<<"\t"<<asym[1]->GetMean()<<"\t"<<asym_v->GetMean()/asym[1]->GetMean()*100.<<"\t"<<lab_v->GetRMS()<<"\t"<<lab[1]->GetRMS()<<"\t"<<lab_v->GetRMS()/lab[1]->GetRMS()*100<<"\t"<<qsq_v->GetRMS()<<"\t"<<qsq[1]->GetRMS()<<"\t"<<qsq_v->GetRMS()/qsq[1]->GetRMS()*100<<"\t"<<asym_v->GetRMS()<<"\t"<<asym[1]->GetRMS()<<"\t"<<asym_v->GetRMS()/asym[1]->GetRMS()*100<<endl;

 outfile_v.close();

 gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/finalDBreplay_953MeVsimulation_septum_%.1f_Ycoll_%.3f_run%d.pdf",sept,YColl_offset,run));
// gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/vertex_angle_scan_septum_%.1f_Ycoll_%.3f.pdf",sept,YColl_offset));
 gSystem->Exec(Form("rm -rf ./temp/*.pdf"));

 }
