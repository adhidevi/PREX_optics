#include "CollimatorL.C"
#include "CollimatorR.C"
#include "sieveL.C"
#include "sieveR.C"
#include "TString.h"

#define tencm

//Checking sieve holes through the spectrometer

void SieveSpecL(int n, TString pl){

gStyle->SetOptStat(0);


TChain *T = new TChain("T");
//for(int i = 1; i < 6; i++){ T->Add(Form("../septmistune/Rootfiles/TestseptscaleLHRS_PREX_0.0_%i.root",i)); }
for(int i = 1; i < 6; i++){ T->Add(Form("/work/halla/parity/disk1/ryanrich/sandwich/SandwichLHRS_PREX_0.0_%i.root",i)); }


TCut totalCut;
TCut specCut = Form("x_%s_tr!=-333.",pl.Data());
TCut colCut = "CollimatorL(x_col_tr,y_col_tr)";
TCut sieveCut = "sieveL(x_zsieve_tr,y_zsieve_tr)!=1";
TCut septCut;

double tol = 0.0;

#ifdef tencm
double xmin[500] = { 0.0, 0.0, 0.082, 0.0862, 0.0989, 0.1084, 0.1195, 0.1236, 0.1326, 0.1958, 0.2149 };
double xmax[500] = {0.1, 0.1, 0.1256, 0.1318, 0.1505, 0.1643, 0.1798, 0.1853, 0.1975, 0.2777, 0.301 };
double ymin[500] = {-0.05, -0.05, -0.0447, -0.0469, -0.053, -0.0571, -0.0612, -0.0624, -0.065, -0.0786, -0.0818 };
double ymax[500] = {0.05, 0.05,0.0445, 0.0466, 0.0529, 0.057, 0.0611, 0.0624, 0.065, 0.0785, 0.0818 };

TCut radCut1 = Form("x_zup1 > %f && x_zup1 < %f && y_zup1 > %f && y_zup1 < %f",0.0448,0.0707,-0.0258,0.0255);
TCut radCut2 = Form("x_zup2 > %f && x_zup2 < %f && y_zup2 > %f && y_zup2 < %f",0.0482,0.0757,-0.0275,0.0272);
TCut downCut1 = Form("x_zdown1 > (%f+%f) && x_zdown1 < (%f-%f) && y_zdown1 > (%f+1.5*%f) && y_zdown1 < (%f-1.5*%f)",xmin[2],tol,xmax[2],tol,ymin[2],tol,ymax[2]);
TCut downCut2 = Form("x_zdown2 > (%f+%f) && x_zdown2 < (%f-%f) && y_zdown2 > (%f+1.5*%f) && y_zdown2 < (%f-1.5*%f)",xmin[3],tol,xmax[3],tol,ymin[3],tol,ymax[3],tol);
TCut downCut3 = Form("x_zdown3 > (%f+%f) && x_zdown3 < (%f-%f) && y_zdown3 > (%f+1.5*%f) && y_zdown3 < (%f-1.5*%f)",xmin[4],tol,xmax[4],tol,ymin[4],tol,ymax[4],tol);
TCut downCut4 = Form("x_zdown4 > (%f+%f) && x_zdown4 < (%f-%f) && y_zdown4 > (%f+1.5*%f) && y_zdown4 < (%f-1.5*%f)",xmin[5],tol,xmax[5],tol,ymin[5],tol,ymax[5],tol);
TCut downCut5 = Form("x_zdown5 > (%f+%f) && x_zdown5 < (%f-%f) && y_zdown5 > (%f+1.5*%f) && y_zdown5 < (%f-1.5*%f)",xmin[6],tol,xmax[6],tol,ymin[6],tol,ymax[6],tol);
TCut downCut6 = Form("x_zdown6 > (%f+%f) && x_zdown6 < (%f-%f) && y_zdown6 > (%f+1.5*%f) && y_zdown6 < (%f-1.5*%f)",xmin[7],tol,xmax[7],tol,ymin[7],tol,ymax[7],tol);
TCut downCut7 = Form("x_zdown7 > (%f+%f) && x_zdown7 < (%f-%f) && y_zdown7 > (%f+1.5*%f) && y_zdown7 < (%f-1.5*%f)",xmin[8],tol,xmax[8],tol,ymin[8],tol,ymax[8],tol);
TCut downCut8 = Form("x_zdown8 > (%f+%f) && x_zdown8 < (%f-%f) && y_zdown8 > (%f+1.5*%f) && y_zdown8 < (%f-1.5*%f)",xmin[9],tol,xmax[9],tol,ymin[9],tol,ymax[9],tol);
TCut downCut9 = Form("x_zdown9 > (%f+%f) && x_zdown9 < (%f-%f) && y_zdown9 > (%f+1.5*%f) && y_zdown9 < (%f-1.5*%f)",xmin[10],tol,xmax[10],tol,ymin[10],tol,ymax[10],tol);
septCut = radCut1&&radCut2&&downCut1&&downCut2&&downCut3&&downCut4&&downCut5&&downCut6&&downCut7&&downCut8&&downCut9;
#endif

//totalCut = colCut&&sieveCut&&septCut&&specCut;
totalCut = sieveCut&&colCut&&specCut;//&&specCut;

TH2F *hplanes = new TH2F("hplanes",Form("Sieve Pattern Hits %s",pl.Data()),n,-0.05,0.05,n,-0.06,0.06);


 hplanes->GetXaxis()->SetTitle("phi");
 hplanes->GetYaxis()->SetTitle("theta");


  TCanvas *c = new TCanvas();
  T->Draw("th_zsieve_tr:ph_zsieve_tr >> hplanes",totalCut*"rate");
  //c->SaveAs(Form("SieveHolesPassCollTest_%s.png",pl.Data()));
  



//delete c;
//delete hplanes;
//delete T;


}
