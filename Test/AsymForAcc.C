#include "CREXdata.h"
#include "TChain.h"
#include "TMath.h"

//Author -- Ryan Richards
//Macro for plotting Q^2 and the Asymmetry for the data
//Want to see how sensitive the asymmetry is to a cut on the radiative tail... for acceptance function
//This has an explicit cut on radiative tail.. Should run ProjQuartz to see ADC plot and quartz edge...
//Ideally wanted to calibrate ADC with dispersive position

void AsymForAcc(int run, double val){

  //Loads the Horowitz tables
  LoadTable("horpb.dat",0); //unstretched 
  LoadTable("horpb1.dat",1); //stretched

  //For the interpolation
    int j = 0;
  
   //ROOT 6 memory management makes me use vectors
   vector <double> Qsq, Angle, ASYM, ASYM_ST, Sens; 
  
//    char *path = getenv("PREXROOT");
 
    TChain *T = new TChain("T");
    T->Add(Form("/lustre/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%i_-1*.root",run));
  

  double thisThtg, thisPhtg, thisP, thisDet;
  double thisu1, thisu2, thisv1, thisv2;
  double thisXvdc[10], thisThvdc[10], thisYvdc[10], thisPhvdc[10];
  double thisAsym, thisAsymSt;
 
  double xq,yq;//For quartz projection
  double thisAngle; //for reconstructed angle
  double thisCosAng; // for cos angle  
  double thisQsq;//for Qsq
  double E = 0.95;//beam energy

  //Central Angle
  double th0 = 4.74*TMath::Pi()/180;
  double cth0 = TMath::Cos(th0);
  double sth0 = TMath::Sin(th0);

  T->SetBranchAddress("L.gold.th",&thisThtg);
  T->SetBranchAddress("L.gold.ph",&thisPhtg);
  T->SetBranchAddress("L.gold.p",&thisP);
  T->SetBranchAddress("P.upQadcL",&thisDet);
  T->SetBranchAddress("L.vdc.u1.nclust",&thisu1);
  T->SetBranchAddress("L.vdc.v1.nclust",&thisv1);
  T->SetBranchAddress("L.vdc.u2.nclust",&thisu2);
  T->SetBranchAddress("L.vdc.v2.nclust",&thisv2);
  
  // --- This segfaults when I loop it for some reason
  T->SetBranchAddress("L.tr.x",&thisXvdc[0]);
  T->SetBranchAddress("L.tr.y",&thisYvdc[0]);
  T->SetBranchAddress("L.tr.th",&thisThvdc[0]);
  T->SetBranchAddress("L.tr.ph",&thisPhvdc[0]);
 
  //Asymmetry table is in degrees so we need to convert from radians to degrees
  double radtodeg = 180/TMath::Pi();



  long n = T->GetEntries();

  for(long i = 0; i < n; i++){ 

    T->GetEntry(i);
  
    xq = thisXvdc[0]+1.3*thisThvdc[0]; yq = thisYvdc[0]+1.3*thisPhvdc[0];


    //VDC cut and x cut
     if( thisu1==1 && thisv1==1 && thisu2==1 && thisv2==1 && xq>val && xq<0.1 && yq>-0.07 && yq<0.00){
 
     thisCosAng = (cth0 - thisPhtg*sth0)/TMath::Sqrt(1+thisThtg*thisThtg+thisPhtg*thisPhtg);
     thisAngle = radtodeg*TMath::ACos(thisCosAng);
    
     thisQsq = 2*E*thisP*(1-thisCosAng);

      Qsq.push_back(thisQsq); 
      Angle.push_back(thisAngle);
      //Energy table in MeV so need to be careful here
      thisAsym = 1e6*Interpolate(thisP*1000,thisAngle,0,1);//A (ppm)
      thisAsymSt = 1e6*Interpolate(thisP*1000,thisAngle,1,1); // Stretched A (ppm); 

      ASYM.push_back(thisAsym);ASYM_ST.push_back(thisAsymSt);//Stretched A (ppm)
      Sens.push_back(fabs(thisAsym-thisAsymSt)/thisAsym);//Sensitivity - Why not?
  }

 }


 //Histograms 
 TH1F *Q2 = new TH1F("Q2","Q^{2} (GeV/c)^{2}",150,0,0.015);
 TH1F *Theta = new TH1F("Theta","#theta_{lab} (deg)",150,2,8);
 TH1F *Asym = new TH1F("Asym","Asymmetry (ppm)",150,0,1);
 TH1F *SAsym = new TH1F("SAsym","Stretched Asymmetry (ppm)",150,0,1);
 TH1F *sens = new TH1F("sens","Sensitivity",150,0,0.05);

 Q2->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
 Theta->GetXaxis()->SetTitle("#theta_{lab} (deg)");
 Asym->GetXaxis()->SetTitle("Asymmetry (ppm)");
 SAsym->GetXaxis()->SetTitle("Stretched Asym. (ppm)");

 for(int k = 0; k < Qsq.size(); k++){

 Q2->Fill(Qsq[k]);
 Theta->Fill(Angle[k]);
 Asym->Fill(ASYM[k]);
 SAsym->Fill(ASYM_ST[k]);
 sens->Fill(Sens[k]);
 }

 TCanvas *t = new TCanvas();
 Theta->Draw();
 t->SaveAs("./temp/plot1.pdf");

 TCanvas *au = new TCanvas();
 Asym->Draw();
 au->SaveAs("./temp/plot2.pdf");

 TCanvas *as = new TCanvas();
 SAsym->Draw();
 as->SaveAs("./temp/plot3.pdf");

 TCanvas *q = new TCanvas();
 Q2->Draw();
 q->SaveAs("./temp/plot4.pdf");

 TCanvas *s = new TCanvas();
 sens->Draw();
 s->SaveAs("./temp/plot5.pdf");

 gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/plots_radCut_%f.pdf",val));
 gSystem->Exec(Form("rm -rf ./temp/*.pdf"));

 ofstream outfile("./TextFiles/output.csv",ios_base::app);
 outfile<<-0.0715<<"\t"<<val<<"\t"<<Theta->GetMean()<<"\t"<<Asym->GetMean()<<"\t"<<SAsym->GetMean()<<"\t"<<Q2->GetMean()<<"\t"<<sens->GetMean()<<endl;
 outfile.close();
 
 std::cout << "Theta " << Theta->GetMean() << "  " << "<A> " << Asym->GetMean() << " " << "Stretched <A>" << SAsym->GetMean() << "   " << "Q2 " << Q2->GetMean() << std::endl;


}
