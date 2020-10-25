#include "CREXdata.h"
#include "TMath.h"
#include "TChain.h"
#include "CollimatorL.C"
#include "UpPlaneL.C"
#include "DownPlaneL.C"

//Author -- Ryan Richards (Sept 24, 2020)
//Macro for plotting Q^2 and Asymmetry for Simulation
//Here we want to compare at the vertex and post vertex


void SimAsym(){

 //Loads the Horowitz tables
 LoadTable("horpb.dat",0);
 LoadTable("horpb1.dat",1);

 //Index 0 will be post vertex, 1 at vertex
 TH1F *Qsq[2], *Ang[2], *Aexpt[2], *SAexpt[2];

 char place[2][20] = {"Vertex","Post-Vertex"};  
 int color[2] = {2,4};


 //Use vectors instead
 vector <double> qsq,ang,ASYM,ASYM_ST,Sens,Rate; // This is for vertex 
 vector <double> qsq1,ang1,ASYM1,ASYM_ST1,Sens1; // Post vertex

 for(int i = 0; i < 2; i++){ 

   Qsq[i] = new TH1F(Form("Qsq[%i]",i),Form("Q^{2} %s (GeV/c)^{2}",place[i]),150,0,0.015); 
   Ang[i] = new TH1F(Form("Ang[%i]",i),Form("#theta_{lab} %s (deg)",place[i]),150,2.0,8.0);
   Aexpt[i] = new TH1F(Form("Aexpt[%i]",i),Form("Asymmetry %s (ppm)",place[i]),150,0,1);
   SAexpt[i] = new TH1F(Form("SAexpt[%i]",i),Form("Stretched Asymmetry %s (ppm)",place[i]),150,0,1);

   //I will make the histograms different colors 
   Qsq[i]->SetLineColor(color[i]); Ang[i]->SetLineColor(color[i]); Aexpt[i]->SetLineColor(color[i]); SAexpt[i]->SetLineColor(color[i]); 

}


 TChain *T = new TChain("T");
 for(int i = 1; i < 6; i++){ T->Add(Form("Pb208_%i.root",i)); } 

 double th0 = 4.74*TMath::Pi()/180;
 double cth0 = TMath::Cos(th0);
 double sth0 = TMath::Sin(th0);

 double thisQsq, thisQsq1, thisAngle, thisAngle1;
 double thisCosAng, thisCosAng1;
 double thisR;
 double thisXcol, thisYcol; 
 double thisThtg, thisPhtg, thisP;
 double thisThz, thisPhz, thisPz;
 double this_xzup1, this_xzup2, this_yzup1, this_yzup2;
 double this_xzdown1, this_xzdown2, this_xzdown3, this_xzdown4, this_xzdown5, this_xzdown6, this_xzdown7, this_xzdown8, this_xzdown9;
 double this_yzdown1, this_yzdown2, this_yzdown3, this_yzdown4, this_yzdown5, this_yzdown6, this_yzdown7, this_yzdown8, this_yzdown9;
 double thisAsym, thisAsym1, thisSAsym, thisSAsym1; 
 double thisSen, thisSen1;
 double beamE, thisEf, thisEf1;


 T->SetBranchAddress("x_col_tr",&thisXcol);
 T->SetBranchAddress("y_col_tr",&thisYcol);
 T->SetBranchAddress("rate",&thisR);

 //Vertex variables
 T->SetBranchAddress("th_tg_tr",&thisThtg);
 T->SetBranchAddress("ph_tg_tr",&thisPhtg);
 T->SetBranchAddress("p_tg_tr",&thisP);
 
 //Post Vertex -- choosing a specific part 
 T->SetBranchAddress("th_ztarg_tr",&thisThz);
 T->SetBranchAddress("ph_ztarg_tr",&thisPhz);
 T->SetBranchAddress("p_ztarg_tr",&thisPz);
 
 T->SetBranchAddress("ev.beamp",&beamE);

 //Additional planes -- annoying long 
 T->SetBranchAddress("x_zup1",&this_xzup1); T->SetBranchAddress("y_zup1",&this_yzup1);
 T->SetBranchAddress("x_zup2",&this_xzup2); T->SetBranchAddress("y_zup2",&this_yzup2);
 T->SetBranchAddress("x_zdown1",&this_xzdown1); T->SetBranchAddress("y_zdown1",&this_yzdown1);
 T->SetBranchAddress("x_zdown2",&this_xzdown2); T->SetBranchAddress("y_zdown2",&this_yzdown2);
 T->SetBranchAddress("x_zdown3",&this_xzdown3); T->SetBranchAddress("y_zdown3",&this_yzdown3);
 T->SetBranchAddress("x_zdown4",&this_xzdown4); T->SetBranchAddress("y_zdown4",&this_yzdown4);
 T->SetBranchAddress("x_zdown5",&this_xzdown5); T->SetBranchAddress("y_zdown5",&this_yzdown5);
 T->SetBranchAddress("x_zdown6",&this_xzdown6); T->SetBranchAddress("y_zdown6",&this_yzdown6);
 T->SetBranchAddress("x_zdown7",&this_xzdown7); T->SetBranchAddress("y_zdown7",&this_yzdown7);
 T->SetBranchAddress("x_zdown8",&this_xzdown8); T->SetBranchAddress("y_zdown8",&this_yzdown8);
 T->SetBranchAddress("x_zdown9",&this_xzdown9); T->SetBranchAddress("y_zdown9",&this_yzdown9);

 double theta; 
 T->SetBranchAddress("ev.Q2",&theta);
 vector<double> Theta;
 
 double radtodeg = 180/TMath::Pi();

 double pb208 = 193.729; //mass of lead 208 in GeV


 long n = T->GetEntries();


 for(int i = 0; i < n; i++){ 
   
  T->GetEntry(i);
  
   if( CollimatorL(thisXcol,thisYcol) && thisXcol !=-333. && UpPlaneL(this_xzup1,this_yzup1,this_xzup2,this_yzup2) && DownPlaneL(this_xzdown1,this_yzdown1,this_xzdown2, this_yzdown2,this_xzdown3,this_yzdown3,this_xzdown4,this_yzdown4,this_xzdown5,this_yzdown5,this_xzdown6,this_yzdown6,this_xzdown7,this_yzdown7, this_xzdown8,this_yzdown8,this_xzdown9,this_yzdown9) ) {
    //vertex 
    thisCosAng = (cth0 - thisPhtg*sth0)/TMath::Sqrt(1+thisThtg*thisThtg+thisPhtg*thisPhtg);
    thisAngle = radtodeg*TMath::ACos(thisCosAng);
    thisEf = pb208*beamE/(pb208+beamE*(1-thisCosAng));
    thisQsq = 2*beamE*thisEf*(1-thisCosAng);

    Theta.push_back(theta);

   //post-vertex
    thisCosAng1 = (cth0 -thisPhz*sth0)/TMath::Sqrt(1+thisThz*thisThz+thisPhz*thisPhz);
    thisAngle1 = radtodeg*TMath::ACos(thisCosAng1);
    thisEf1 = pb208*beamE/(pb208+beamE*(1-thisCosAng1));
    thisQsq1 = 2*beamE*thisEf1*(1-thisCosAng1);
    
    Rate.push_back(thisR);
    ang.push_back(thisAngle); ang1.push_back(thisAngle1);
    thisAsym = 1e6*Interpolate(beamE*1000,thisAngle,0,1);
    thisAsym1 = 1e6*Interpolate(beamE*1000,thisAngle1,0,1);
    thisSAsym = 1e6*Interpolate(beamE*1000,thisAngle,1,1);
    thisSAsym1 = 1e6*Interpolate(beamE*1000,thisAngle1,1,1); 
    thisSen = fabs(thisAsym-thisSAsym)/thisAsym; 
    thisSen1 = fabs(thisAsym1-thisSAsym1)/thisAsym1;

    ASYM.push_back(thisAsym); ASYM_ST.push_back(thisAsym1);
    ASYM1.push_back(thisSAsym); ASYM_ST1.push_back(thisSAsym1);

    Sens.push_back(thisSen); Sens1.push_back(thisSen1);

    qsq.push_back(thisQsq); qsq1.push_back(thisQsq1);

 
   }

 

 }   



 for(int k = 0; k < ang.size(); k++){

    Ang[0]->Fill(ang[k],Rate[k]);
    Ang[1]->Fill(ang1[k],Rate[k]);
    Aexpt[0]->Fill(ASYM[k],Rate[k]);
    Aexpt[1]->Fill(ASYM1[k],Rate[k]);
    SAexpt[0]->Fill(ASYM_ST[k],Rate[k]);
    SAexpt[1]->Fill(ASYM_ST1[k],Rate[k]);
    Qsq[0]->Fill(qsq[k],Rate[k]);
//    Qsq[1]->Fill(Theta[k],Rate[k]);
    Qsq[1]->Fill(qsq1[k],Rate[k]);

 } 




TCanvas *c = new TCanvas();
Ang[0]->Draw();
c->SaveAs(Form("Ang%s.png",place[0]));


TCanvas *c1 = new TCanvas();
Ang[1]->Draw();
c1->SaveAs(Form("Ang%s.png",place[1]));

TCanvas *c2 = new TCanvas();
Aexpt[0]->Draw();
c2->SaveAs(Form("Aexpt%s.png",place[0]));

TCanvas *c3 = new TCanvas();
Aexpt[1]->Draw();
c3->SaveAs(Form("Aexpt%s.png",place[1]));



TCanvas *c4 = new TCanvas();
Qsq[0]->Draw();
c4->SaveAs(Form("Qsq%s.png",place[0]));

TCanvas *c5 = new TCanvas();
Qsq[1]->Draw();
c5->SaveAs(Form("Qsq%s.png",place[1]));



}
