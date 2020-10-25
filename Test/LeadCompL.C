#include <iostream>
#include <fstream>
#include "CollimatorL.C"

#define tencm

//For LHRS-- I am separating them to avoid long code
void LeadCompL(int run){

gStyle->SetOptStat(0);


TChain *T = new TChain("T");
T->Add(Form("prexLHRS_%d_-1.root",run));

TChain *T1 = new TChain("T");
 for(int i = 1; i < 6; i++){
   //T1->Add(Form("../septmistune/Rootfiles/septscaleLHRS_PREX_0.0_%i.root",i));  
  T1->Add(Form("Pb208_%i.root",i));

 }


 TH1F *dpdat, *dpsim, *phdat, *phsim;

 //Index 0 is for data, index 1 is for simulation
 TH1F *x[2][4];
 char vars[4][10] = { "x","th","y","ph" };

 for(int i = 0; i < 2; i++){ 
  for(int j = 0; j < 4; j++){ 

   if( j == 0 ){ x[i][j] = new TH1F(Form("x[%i][%i]",i,j),Form("LHRS %s at Detector Plane",vars[j]),150,-0.8,0.1); }
   else if( j == 1 ) { x[i][j] = new TH1F(Form("x[%i][%i]",i,j),Form("LHRS %s at Detector Plane",vars[j]),150,-0.05,0.05); }
   else if( j == 2 ) { x[i][j] = new TH1F(Form("x[%i][%i]",i,j), Form("LHRS %s at Detector Plane",vars[j]),150,-0.1,0.1); }
   else { x[i][j] = new TH1F(Form("x[%i][%i]",i,j), Form("LHRS %s at Detector Plane",vars[j]),150,-0.05,0.05); }

   if( i == 0 ) { x[i][j]->SetLineColor(kBlue); } else { x[i][j]->SetLineColor(kRed); }

  }
 }


 dpdat = new TH1F("dpdat","Dp Distribution (LHRS)",150,-0.01,0.005);
 dpsim = new TH1F("dpsim","Dp Distribution",150,-0.01,0.005);

 phdat = new TH1F("phdat","#phi_{tg} (LHRS)",150,-0.04,0.04);
 phsim = new TH1F("phsim","#phi_{tg}",150,-0.04,0.04);

 dpdat->SetLineColor(kBlue); phdat->SetLineColor(kBlue);
 dpsim->SetLineColor(kRed);  phsim->SetLineColor(kRed);

 TCut cut_dat = "L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1";

 //Full cut for simulation
 TCut cut_sim;
 
 //Extra simulation cuts
 TCut xfp = "x_fp_tr!=-333.";
 TCut colcut = "CollimatorL(x_col_tr,y_col_tr)";

 double tol = 0.;


#ifdef tencm
double xmin[500] = { 0.0, 0.0, 0.082, 0.0862, 0.0989, 0.1084, 0.1195, 0.1236, 0.1326, 0.1958, 0.2149 };
double xmax[500] = {0.1, 0.1, 0.1256, 0.1318, 0.1505, 0.1643, 0.1798, 0.1853, 0.1975, 0.2777, 0.301 };
double ymin[500] = {-0.05, -0.05, -0.0447, -0.0469, -0.053, -0.0571, -0.0612, -0.0624, -0.065, -0.0786, -0.0818 };
double ymax[500] = {0.05, 0.05,0.0445, 0.0466, 0.0529, 0.057, 0.0611, 0.0624, 0.065, 0.0785, 0.0818 };

TCut radCut1 = Form("x_zup1 > %f && x_zup1 < %f && y_zup1 > %f && y_zup1 < %f",0.0448,0.0707,-0.0258,0.0255);
TCut radCut2 = Form("x_zup2 > %f && x_zup2 < %f && y_zup2 > %f && y_zup2 < %f",0.0482,0.0757,-0.0275,0.0272);
TCut downCut1 = Form("x_zdown1 > (%f+%f) && x_zdown1 < (%f-%f) && y_zdown1 > (%f+1.5*%f) && y_zdown1 < (%f-1.5*%f)",xmin[2],tol,xmax[2],tol,ymin[2],tol,ymax[2],tol);
TCut downCut2 = Form("x_zdown2 > (%f+%f) && x_zdown2 < (%f-%f) && y_zdown2 > (%f+1.5*%f) && y_zdown2 < (%f-1.5*%f)",xmin[3],tol,xmax[3],tol,ymin[3],tol,ymax[3],tol);
TCut downCut3 = Form("x_zdown3 > (%f+%f) && x_zdown3 < (%f-%f) && y_zdown3 > (%f+1.5*%f) && y_zdown3 < (%f-1.5*%f)",xmin[4],tol,xmax[4],tol,ymin[4],tol,ymax[4],tol);
TCut downCut4 = Form("x_zdown4 > (%f+%f) && x_zdown4 < (%f-%f) && y_zdown4 > (%f+1.5*%f) && y_zdown4 < (%f-1.5*%f)",xmin[5],tol,xmax[5],tol,ymin[5],tol,ymax[5],tol);
TCut downCut5 = Form("x_zdown5 > (%f+%f) && x_zdown5 < (%f-%f) && y_zdown5 > (%f+1.5*%f) && y_zdown5 < (%f-1.5*%f)",xmin[6],tol,xmax[6],tol,ymin[6],tol,ymax[6],tol);
TCut downCut6 = Form("x_zdown6 > (%f+%f) && x_zdown6 < (%f-%f) && y_zdown6 > (%f+1.5*%f) && y_zdown6 < (%f-1.5*%f)",xmin[7],tol,xmax[7],tol,ymin[7],tol,ymax[7],tol);
TCut downCut7 = Form("x_zdown7 > (%f+%f) && x_zdown7 < (%f-%f) && y_zdown7 > (%f+1.5*%f) && y_zdown7 < (%f-1.5*%f)",xmin[8],tol,xmax[8],tol,ymin[8],tol,ymax[8],tol);
TCut downCut8 = Form("x_zdown8 > (%f+%f) && x_zdown8 < (%f-%f) && y_zdown8 > (%f+1.5*%f) && y_zdown8 < (%f-1.5*%f)",xmin[9],tol,xmax[9],tol,ymin[9],tol,ymax[9],tol);
TCut downCut9 = Form("x_zdown9 > (%f+%f) && x_zdown9 < (%f-%f) && y_zdown9 > (%f+1.5*%f) && y_zdown9 < (%f-1.5*%f)",xmin[10],tol,xmax[10],tol,ymin[10],tol,ymax[10],tol);
cut_sim = xfp&&colcut&&radCut1&&radCut2&&downCut1&&downCut2&&downCut3&&downCut4&&downCut5&&downCut6&&downCut7&&downCut8&&downCut9;
#endif

T->SetAlias("xq","L.tr.x[0]+0.9*L.tr.th[0]");
T->SetAlias("yq","L.tr.y[0]+0.9*L.tr.ph[0]");

x[0][1]->GetYaxis()->SetRangeUser(0,10000);
 
T->Project(dpdat->GetName(),"L.gold.dp[0]+0.003",cut_dat);
T->Project(phdat->GetName(),"L.gold.ph[0]",cut_dat);

T1->Project(dpsim->GetName(),"(p_ztarg_tr-950)/950",cut_sim*"rate");
T1->Project(phsim->GetName(),"ph_tg_tr",cut_sim*"rate");


T->Project(x[0][0]->GetName(),"xq",cut_dat);
T->Project(x[0][1]->GetName(),"L.tr.th[0]",cut_dat);
T->Project(x[0][2]->GetName(),"yq",cut_dat);
T->Project(x[0][3]->GetName(),"L.tr.ph[0]",cut_dat);


for(int i = 0; i < 4; i++){ 
  T1->Project(x[1][i]->GetName(),Form("%s_fp_tr",vars[i]),cut_sim*"rate");
}

auto leg = new TLegend(0.1,0.75,0.4,0.9);
leg->AddEntry(dpdat,"Data","l");
leg->AddEntry(dpsim,"Simulation","l");



TCanvas *c = new TCanvas();
dpdat->Draw();
dpdat->Scale(1.0/dpdat->GetBinContent(dpdat->GetMaximumBin()));
dpsim->Draw("same");
dpsim->Scale(1.0/dpsim->GetBinContent(dpsim->GetMaximumBin()));
leg->Draw();
c->SaveAs(Form("dp_Run%d.jpg",run));


/*
TCanvas *c1 = new TCanvas();
phdat->Draw();
phsim->Draw("same");
phsim->Scale(phdat->Integral()/phsim->Integral());
//leg->Draw();
c1->SaveAs(Form("ph_Run%d.jpg",run));



TCanvas *c2 = new TCanvas();
c2->Divide(2,2);
c2->cd(1);
x[0][0]->Draw();
x[1][0]->Draw("same");
x[1][0]->Scale(x[0][0]->Integral()/x[1][0]->Integral());
c2->cd(2);
x[0][1]->Draw();
x[1][1]->Draw("same");
x[1][1]->Scale(x[0][1]->Integral()/x[1][1]->Integral());
c2->cd(3);
x[0][2]->Draw();
x[1][2]->Draw("same");
x[1][2]->Scale(x[0][2]->Integral()/x[1][2]->Integral());
c2->cd(4);
x[0][3]->Draw();
x[1][3]->Draw("same");
x[1][3]->Scale(x[0][3]->Integral()/x[1][3]->Integral());

c2->SaveAs(Form("fp_Run%d.jpg",run));
*/


}








































