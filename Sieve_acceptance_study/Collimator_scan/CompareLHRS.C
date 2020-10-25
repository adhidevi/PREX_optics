#include "CollimatorL.C"
//#include "UpPlane.C"
#define tencm

//Author Ryan Richards, compares theta_tg, phi_tg, polar angle and Q^2 for Apparent quantities  to see how they change with septum current

//Arguments run number, sept scan(0.0 - Nominal, -2.0 is -2% etc), ADC cut value, central angle for data
void CompareLHRS(int run, double sept, double ADC, double th0){ 


  

  double d2r = TMath::Pi()/180; 
  double r2d = 1/d2r;

  TH1F *thtg[2],*phtg[2], *qsq[2], *lab[2];

  int color[2] = { 2, 4 };

  for(int i = 0; i < 2; i++){ 
    thtg[i] = new TH1F(Form("thtg[%i]",i),"LHRS Apparent #theta_{tg}",150,-0.06,0.06); 
    phtg[i] = new TH1F(Form("phtg[%i]",i),"LHRS Apparent #phi_{tg}",150,-0.03,0.03);

    thtg[i]->SetLineColor(color[i]); phtg[i]->SetLineColor(color[i]);
    thtg[i]->GetXaxis()->SetTitle("#theta_{tg} (rad)"); phtg[i]->GetXaxis()->SetTitle("#phi_{tg} (rad)");

    //Qsq and lab histograms will have separate names to compare central values
    if(i == 0){ qsq[i] = new TH1F(Form("qsq[%i]",i),"LHRS Apparent Q^{2} Data",150,0.0,0.015);
                lab[i] = new TH1F(Form("lab[%i]",i),"LHRS Apparent #theta_{lab} Data",150,1.0,9.0);
    } else{ 
    qsq[i] = new TH1F(Form("qsq[%i]",i),Form("LHRS Apparent Q^{2} Septum Detuned %.1f",sept),150,0.0,0.015);
    lab[i] = new TH1F(Form("lab[%i]",i),Form("LHRS Apparent #theta_{lab} Septum Detuned %.1f",sept),150,1.0,9.0);
   
    } 

    qsq[i]->SetLineColor(color[i]); lab[i]->SetLineColor(color[i]);
    lab[i]->GetXaxis()->SetTitle("#theta_{lab} (deg)"); qsq[i]->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");

  }

  TChain *T = new TChain("T");
  T->Add(Form("../prex2Rootfiles/prexLHRS_%i_-1.root",run)); 
  
  TChain *T1 = new TChain("T");
  for(int i = 1; i < 6; i++){ T1->Add(Form("../sandwich/SandwichLHRS_PREX_%.1f_%i.root",sept,i)); }


  //Data cuts
  TCut vdc_cut = "L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1";
  TCut adc_cut = Form("P.upQadcL > %f",ADC);
  TCut cut_dat = vdc_cut+adc_cut;


  //Simulation cuts
  TCut x_vdc = "x_vdc_tr!=-333.";
  TCut sept_cut;// You can comment out if you like...
  TCut col_cut = "CollimatorL(x_col_tr,y_col_tr)";
  //Needs to be thought about 2.5 MeV cut from data
  TCut mom_cut = "(p_ztarg_tr-951)/951 > -0.002628812";//This is 2.5/951 
  TCut cut_sim;
  TCut radCut1, radCut2;
  TCut downCut1, downCut2, downCut3, downCut4, downCut5, downCut6, downCut7, downCut8, downCut9;
 


  double tol = 0.;//Don't need this but I have been too lazy to remove it :)

  #ifdef tencm 
  double xmin[500] = { 0.0, 0.0, 0.082, 0.0862, 0.0989, 0.1084, 0.1195, 0.1236, 0.1326, 0.1958, 0.2149 };
  double xmax[500] = {0.1, 0.1, 0.1256, 0.1318, 0.1505, 0.1643, 0.1798, 0.1853, 0.1975, 0.2777, 0.301 };
  double ymin[500] = {-0.05, -0.05, -0.0447, -0.0469, -0.053, -0.0571, -0.0612, -0.0624, -0.065, -0.0786, -0.0818 };
  double ymax[500] = {0.05, 0.05,0.0445, 0.0466, 0.0529, 0.057, 0.0611, 0.0624, 0.065, 0.0785, 0.0818 };

  radCut1 = Form("x_zup1 > %f && x_zup1 < %f && y_zup1 > %f && y_zup1 < %f",0.0448,0.0707,-0.0258,0.0255);
  radCut2 = Form("x_zup2 > %f && x_zup2 < %f && y_zup2 > %f && y_zup2 < %f",0.0482,0.0757,-0.0275,0.0272);
  downCut1 = Form("x_zdown1 > (%f+%f) && x_zdown1 < (%f-%f) && y_zdown1 > (%f+1.5*%f) && y_zdown1 < (%f-1.5*%f)",xmin[2],tol,xmax[2],tol,ymin[2],tol,ymax[2],tol);
  downCut2 = Form("x_zdown2 > (%f+%f) && x_zdown2 < (%f-%f) && y_zdown2 > (%f+1.5*%f) && y_zdown2 < (%f-1.5*%f)",xmin[3],tol,xmax[3],tol,ymin[3],tol,ymax[3],tol);
  downCut3 = Form("x_zdown3 > (%f+%f) && x_zdown3 < (%f-%f) && y_zdown3 > (%f+1.5*%f) && y_zdown3 < (%f-1.5*%f)",xmin[4],tol,xmax[4],tol,ymin[4],tol,ymax[4],tol);
  downCut4 = Form("x_zdown4 > (%f+%f) && x_zdown4 < (%f-%f) && y_zdown4 > (%f+1.5*%f) && y_zdown4 < (%f-1.5*%f)",xmin[5],tol,xmax[5],tol,ymin[5],tol,ymax[5],tol);
  downCut5 = Form("x_zdown5 > (%f+%f) && x_zdown5 < (%f-%f) && y_zdown5 > (%f+1.5*%f) && y_zdown5 < (%f-1.5*%f)",xmin[6],tol,xmax[6],tol,ymin[6],tol,ymax[6],tol);
  downCut6 = Form("x_zdown6 > (%f+%f) && x_zdown6 < (%f-%f) && y_zdown6 > (%f+1.5*%f) && y_zdown6 < (%f-1.5*%f)",xmin[7],tol,xmax[7],tol,ymin[7],tol,ymax[7],tol);
  downCut7 = Form("x_zdown7 > (%f+%f) && x_zdown7 < (%f-%f) && y_zdown7 > (%f+1.5*%f) && y_zdown7 < (%f-1.5*%f)",xmin[8],tol,xmax[8],tol,ymin[8],tol,ymax[8],tol);
  downCut8 = Form("x_zdown8 > (%f+%f) && x_zdown8 < (%f-%f) && y_zdown8 > (%f+1.5*%f) && y_zdown8 < (%f-1.5*%f)",xmin[9],tol,xmax[9],tol,ymin[9],tol,ymax[9],tol);
  downCut9 = Form("x_zdown9 > (%f+%f) && x_zdown9 < (%f-%f) && y_zdown9 > (%f+1.5*%f) && y_zdown9 < (%f-1.5*%f)",xmin[10],tol,xmax[10],tol,ymin[10],tol,ymax[10],tol);
 
  sept_cut = radCut1+radCut2+downCut1+downCut2+downCut3+downCut4+downCut5+downCut6+downCut7+downCut8+downCut9;
  #endif  

  cut_sim = x_vdc+col_cut+sept_cut+mom_cut;

  //th0 should be 4.77
  double cth0 = TMath::Cos(th0*d2r);
  double sth0 = TMath::Sin(th0*d2r);



  double Off = th0*d2r; //for data

   
  double angsim = 4.74*d2r;
  double cang = TMath::Cos(angsim);
  double sang = TMath::Sin(angsim);

  //Should this have an offset applied too since I compare to simulation? I won't include it for now
  T->SetAlias("th",Form("%f*TMath::ACos((%f-%f*L.gold.ph)/TMath::Sqrt(1+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph))",r2d,cth0,sth0));
  T->SetAlias("Qsq",Form("2*0.951*L.gold.p*(1 -((%f-%f*L.gold.ph)/TMath::Sqrt(1+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph)))",cth0,sth0));

  //Alias for simulation
  T1->SetAlias("thsim",Form("%f*TMath::ACos((%f-%f*ph_ztarg_tr)/TMath::Sqrt(1+th_ztarg_tr*th_ztarg_tr+ph_ztarg_tr*ph_ztarg_tr))",r2d,cang,sang));
  T1->SetAlias("qsim",Form("2*ev.beamp*ev.ep*(1-((%f-%f*ph_ztarg_tr)/TMath::Sqrt(1+th_ztarg_tr*th_ztarg_tr+ph_ztarg_tr*ph_ztarg_tr)))",cang,sang));


  //Now for filling of histograms -- Data
  T->Project(thtg[0]->GetName(),"L.gold.th[0]",cut_dat);
  T->Project(phtg[0]->GetName(),Form("L.gold.ph[0]+%f",Off),cut_dat);//for the data comparison
  T->Project(lab[0]->GetName(),"th",cut_dat);
  T->Project(qsq[0]->GetName(), "Qsq",cut_dat);

  //Simulation -- using apparent quantities
  T1->Project(thtg[1]->GetName(),"th_ztarg_tr",cut_sim*"rate");
  T1->Project(phtg[1]->GetName(),Form("ph_ztarg_tr+%f",angsim),cut_sim*"rate");
  T1->Project(lab[1]->GetName(),"thsim",cut_sim*"rate");
  T1->Project(qsq[1]->GetName(),"qsim",cut_sim*"rate");

  //I find positioning the legend to be annoying...
  auto leg = new TLegend(0.1,0.75,0.4,0.9);
  leg->AddEntry(thtg[0],"Data","l");
  leg->AddEntry(thtg[1],Form("Sim Sept Detuned %.1f Per",sept),"l");




  TCanvas *c1 = new TCanvas();
  thtg[0]->Draw();
  thtg[1]->Draw("HIST same");
  thtg[1]->Scale(thtg[0]->Integral()/thtg[1]->Integral());
  leg->Draw();

 
  TCanvas *c2 = new TCanvas();
  phtg[0]->Draw();
  phtg[1]->Draw("HIST same");
  //phtg[1]->Scale(phtg[0]->Integral()/phtg[1]->Integral());
 // leg->Draw();

  TCanvas *c3 = new TCanvas();
  lab[0]->Draw();
  lab[1]->Draw("HIST same");
  lab[1]->Scale(lab[0]->Integral()/lab[1]->Integral());
    

  TCanvas *c4 = new TCanvas(); 
  qsq[0]->Draw();
  qsq[1]->Draw("HIST same");
  qsq[1]->Scale(qsq[0]->Integral()/qsq[1]->Integral());


  std::cout << "Data Lab theta mean: " << lab[0]->GetMean() << "  " << "Data Lab theta rms" << lab[0]->GetRMS() << std::endl;
  std::cout << "Sim Lab theta mean: " << lab[1]->GetMean() << "   " << "Sim Lab theta rms" << lab[1]->GetRMS() << std::endl;

  std::cout << "Data  Qsq mean: " << qsq[0]->GetMean() << "  " << "Data  Qsq rms" << qsq[0]->GetRMS() << std::endl;
  std::cout << "Sim  Qsq mean: " << qsq[1]->GetMean() << "   " << "Sim  Qsq rms" << qsq[1]->GetRMS() << std::endl;



  





}
