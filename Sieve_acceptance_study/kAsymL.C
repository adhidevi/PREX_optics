#include "CollimatorL.C"
#include "UpPlane.C"
#include "DownPlane.C"
#include "CREXdata.h"

//Author Ryan Richards - 
//Shift polar angle in sim to match data and compute Q^2, Asym

//Run number, sept scan(0.0 - nominal), ADC cut, central angle for data, gold.p peak position, momentum cut for th_tg=0 in GeV below the peak, and a "pinch" squeeze parameter. 
void kAsymL(int run,double sept,double ADC,double th0, double peakpos, double simMomCut,double E0,double YColl_offset){
  double pinch=0.0;
  int iprint = 0;
  // expected eloss in lead target about 1 MeV
  double efact = (E0-0.001)/peakpos;
  // Sim E offset  // 953.4 - 1
  double efactsim = E0*1000/951.0;
  
  LoadTable("horpb.dat",0);//Lead table - sufficient
  int ngdat =0;
  int ngsim =0;
  double d2r = TMath::Pi()/180;
  double r2d = 1/d2r;
  
  //th0 is for data 
  double cth0 = TMath::Cos(th0*d2r); double sth0 = TMath::Sin(th0*d2r);
  
  double simth0 = 4.74;
  double cthsim = TMath::Cos(simth0*d2r);
  double sthsim = TMath::Sin(simth0*d2r);
  
//  gStyle->SetOptStat(0);
  TH1F *lab[2],*qsq[2], *asym[2];
  TH1F *hmom[2];
  TH1F *hphtg[2], *hthtg[2];
  int color[2] = { 2, 4 };
  
  for(int i = 0; i < 2; i++){ 
    
    if (i == 0) {
      qsq[i] = new TH1F(Form("qsq[%i]",i),"LHRS Apparent Q^{2} Data",150,0.0,0.012);
      lab[i] = new TH1F(Form("lab[%i]",i),"LHRS Apparent #theta_{lab} Data",150,2.0,8.0);
      asym[i] = new TH1F(Form("asym[%i]",i),"LHRS Apparent Asymmetry Data",200,-100,800);
      hmom[i] = new TH1F(Form("hmom[%i]",i),"LHRS Apparent Momentum Data",200,0.940,0.960);
      hphtg[i] = new TH1F(Form("hphtg[%i]",i),"LHRS Apparent #phi_{tg} Data",200,2,8);
      hthtg[i] = new TH1F(Form("hthtg[%i]",i),"LHRS Apparent #theta_{tg} Data",200,-0.08,0.08);
    } else {
      qsq[i] = new TH1F(Form("qsq[%i]",i),Form("G4HRS-L Apparent Q^{2} Septum %.1f",sept),150,0.0,0.012);
      lab[i] = new TH1F(Form("lab[%i]",i),Form("G4HRS-L Apparent #theta_{lab} Septum %.1f",sept),150,2.0,8.0);
      asym[i] = new TH1F(Form("asym[%i]",i),Form("G4HRS-L Apparent Asymmetry Septum %.1f",sept),200,-100,800);
      hmom[i] = new TH1F(Form("hmom[%i]",i),Form("G4HRS-L Apparent Momentum Septum %.1f",sept),200,0.940,0.960);
      hphtg[i] = new TH1F(Form("hphtg[%i]",i),Form("G4HRS-L Apparent #phi_{tg} Septum %.1f",sept),200,2,8);
      hthtg[i] = new TH1F(Form("hthtg[%i]",i),Form("G4HRS-L Apparent #theta_{tg} Septum %.1f",sept),200,-0.08,0.08);
    }
    
    asym[i]->SetLineColor(color[i]); qsq[i]->SetLineColor(color[i]); lab[i]->SetLineColor(color[i]); 
    hmom[i]->SetLineColor(color[i]); hphtg[i]->SetLineColor(color[i]); hthtg[i]->SetLineColor(color[i]);
  }

  // simulation vertex distribution
  
  TH1F *hqsq_v = new TH1F("hqsq_v",Form("G4HRL-L Vertex Q^{2} Septum %.1f",sept),150,0.0,0.015);
  TH1F *hlab_v = new TH1F("hlab_v",Form("G4HRL-L Vertex #theta_{lab} Septum  %.1f",sept),150,1.0,9.0);
  TH1F *hasym_v = new TH1F("hasym_v",Form("G4HRL-L Vertex Asymmetry Septum %.1f",sept),200,-100,800);
  hlab_v->SetLineColor(6);
  hqsq_v->SetLineColor(6);
  hasym_v->SetLineColor(6);
  
  //For simulation
  vector <double> Rate,Qsq,Asym,Lab,Psim, Phtgsim, Thtgsim, Lab_v, Qsq_v, Asym_v;
  double thisRate, thisTh, thisPh, thisQsq, thisXcol, thisYcol, thisXfp, thisAsym, thisAng,thisCosAng;
  double thisQsq_v, thisAsym_v, thisLab_v;
  double thisE,thisP, thismom;
  double thisxu1,thisxu2,thisyu1,thisyu2;
  double thisxd1,thisxd2,thisxd3,thisxd4,thisxd5,thisxd6,thisxd7,thisxd8,thisxd9;
  double thisyd1,thisyd2,thisyd3,thisyd4,thisyd5,thisyd6,thisyd7,thisyd8,thisyd9;
  
  //For data, need cluster, adc,th,ph,p
  vector <double> Angle,Qsqdat, Asymdat, Pdat, Phtgdat, Thtgdat;
  double thisu1,thisv1,thisu2,thisv2,thisPdat,thisThdat,thisPhdat,thisQsqdat,thisAsymdat,thisCosAngdat,thisAngdat,thisADC;
  double thisTbits;
  
  TChain *T = new TChain("T");
  T->Add(Form("/chafs1/work1/prex_counting/RootFilesQsq/prexLHRS_%i_-1*.root",run));
  
  TChain *T1 = new TChain("T");
  for(int i = 1; i < 6; i++){ T1->Add(Form("/w/halla-scifs17exp/parity/disk1/ryanrich/sandwich/SandwichLHRS_PREX_%.1f_%i.root",sept,i)); }

  T->SetBranchAddress("P.evtypebits",&thisTbits);
  T->SetBranchAddress("L.vdc.u1.nclust",&thisu1); T->SetBranchAddress("L.vdc.v1.nclust",&thisv1);
  T->SetBranchAddress("L.vdc.u2.nclust",&thisu2); T->SetBranchAddress("L.vdc.v2.nclust",&thisv2);
  T->SetBranchAddress("L.gold.th",&thisThdat); T->SetBranchAddress("L.gold.ph",&thisPhdat);
  T->SetBranchAddress("L.gold.p",&thisPdat); T->SetBranchAddress("P.upQadcL",&thisADC);
  //Probably need momentum cut

  int n = T->GetEntries();

  int nufasymdat=0;
  //Looping over tree in data
  bool dcbit;
  bool dcmult;
  bool dcadc;
  bool dcwild;
  for(int i = 0; i < n; i++){
     T->GetEntry(i);

     thisPdat *= efact; // added on 
     dcbit = (((int) thisTbits)&2)==2;
     dcwild = abs(thisThdat)<0.08 && abs(thisPhdat)<0.05 && thisPdat>E0*0.96 && thisPdat<E0*1.002;
     dcadc = thisADC > ADC;
     dcmult = thisu1 == 1 && thisv1 == 1 && thisu2 == 1 && thisv2 == 1;
     //Need momentum cut for this, I think
     //     if(thisu1 == 1 && thisv1 == 1 && thisu2 == 1 && thisv2 == 1 && thisADC > ADC){
     if(dcbit && dcmult && dcadc && dcwild) {

       thisCosAngdat = (cth0-thisPhdat*sth0)/(TMath::Sqrt(1+thisThdat*thisThdat+thisPhdat*thisPhdat));
       thisAngdat = r2d*TMath::ACos(thisCosAngdat);
       //       if (thisAngdat>2 && thisAngdat<9) {
	 ngdat++;
	 //Using a hardcoded beam energy -- run dependent, yes?
	 // Beam Energy 
	 thisQsqdat = 2*E0*thisPdat*(1-thisCosAngdat);
	 thisAsymdat = 1e9*Interpolate(thisPdat*1000,thisAngdat,0,1);
	 
	 Angle.push_back(thisAngdat);
	 Qsqdat.push_back(thisQsqdat);     
	 Asymdat.push_back(thisAsymdat);
	 Pdat.push_back(thisPdat);
	 Phtgdat.push_back(r2d*thisPhdat+th0);
	 Thtgdat.push_back(thisThdat);
	 
	 if(thisAsymdat<100) {
	   nufasymdat++;
	   // cout << ngdat << "    " << thisAsymdat << "   " << thisPhdat*r2d << "   " << thisThdat<< "   " << thisAngdat << "   " << thisQsqdat << "    " << thisPdat << endl;
	 }
	 //std::cout << thisAsymdat << std::endl; 
	 
	 //       }
     }
     
  } 
  
  
  T1->SetBranchAddress("rate",&thisRate); T1->SetBranchAddress("x_vdc_tr",&thisXfp);
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
  
  
  // vertex asym, Q2, theta
  T1->SetBranchAddress("ev.A",&thisAsym_v);
  T1->SetBranchAddress("ev.Q2",&thisQsq_v);
  T1->SetBranchAddress("ev.Th",&thisLab_v);  // deg
  
  int m = T1->GetEntries();
 
  int nufasym =0;
  int nufasym_v =0;

  bool scmom;
  bool scfp;
  bool sccol;
  bool scsvup;
  bool scsvdn;
  bool sctot;

  double pinchscan=1/1000.*(pinch);
  //Now loop over sim tree
  for(int j = 0; j < m; j++){ 
    T1->GetEntry(j);
    
    scfp = thisXfp!=-333.;
    sccol = CollimatorL(thisXcol,thisYcol+YColl_offset);
    scsvup = UpPlane(thisxu1,thisyu1,thisxu2,thisyu2,1);
    scsvdn = DownPlane(thisxd1,thisyd1,thisxd2,thisyd2,
		       thisxd3-pinchscan,thisyd3,thisxd4-pinchscan,thisyd4,
		       thisxd5,thisyd5,
		       thisxd6,thisyd6,thisxd7,thisyd7,thisxd8,thisyd8,thisxd9,thisyd9,1);
    

    thismom*=efactsim;

    //*** The momentum cutneeds to be evaluated and set for each run ***
    //   cut th_tg = 0 at threshold below the expected elastic peak
    scmom = (thismom - 50*thisTh) > ((E0 - 0.001 - simMomCut)*1000);
//    scmom = (thismom - 50*thisTh) > -0.003;
        
    sctot = scfp && sccol && scsvup && scsvdn && scmom;
    
    // if( thisXfp!=-333. && CollimatorL(thisXcol,thisYcol) && UpPlane(thisxu1,thisyu1,thisxu2,thisyu2,1) &&
    // DownPlane(thisxd1,thisyd1,thisxd2,thisyd2,thisxd3,thisyd3,thisxd4,thisyd4,thisxd5,thisyd5, thisxd6,thisyd6,thisxd7,thisyd7,thisxd8,thisyd8,thisxd9,thisyd9,1)&&((thismom-951)/951 > -(2.5/951)) ){
    //  if( thisXfp!=-333. && CollimatorL(thisXcol,thisYcol) &&((thismom-951)/951 > -(2.5/951)) ){
    if (sctot) {
      
      ngsim++;
      thisCosAng = (cthsim-thisPh*sthsim)/TMath::Sqrt(1+thisPh*thisPh+thisTh*thisTh);  
      // 
      //We will have to be careful here  
      thisAng = r2d*TMath::ACos(thisCosAng);
      //This is the Angle in degrees
      thisAng = thisAng;
      //Now with this, compute Q^2 and Asymmetry
      //  thisQsq = 2*thisP*thisE*(1-TMath::Cos(thisAng*d2r));
      //  thisAsym = 1e6*Interpolate(thisP*1000,thisAng,0,1);
      thisQsq = 2*thismom/1000.*E0*(1-TMath::Cos(thisAng*d2r));
      thisAsym = 1e9*Interpolate(thismom,thisAng,0,1);
      
      if(thisAsym<100) {
	nufasym++;
	//	cout << ngsim << "    " << thisAsym << "   " << thisAng << "   " << thisQsq << endl;
      }
      
      Rate.push_back(thisRate); Qsq.push_back(thisQsq);
      Lab.push_back(thisAng); Asym.push_back(thisAsym);
      Psim.push_back(thismom/1000.);
      Phtgsim.push_back(thisPh*r2d +simth0);
      Thtgsim.push_back(thisTh);
      
      //std::cout << thisAsym << std::endl;  
      if(thisAsym_v<100) {
	nufasym_v++;
	//	cout << ngsim << "    " << thisAsym << "   " << thisAng << "   " << thisQsq << endl;
      }
      
      Qsq_v.push_back(thisQsq_v);
      Asym_v.push_back(thisAsym_v);
      Lab_v.push_back(thisLab_v);
      
    }
    
  } 
  
  //Now fill histograms
  //data
  for(int l = 0; l < Qsqdat.size(); l++){
    qsq[0]->Fill(Qsqdat[l]);
    lab[0]->Fill(Angle[l]);
    asym[0]->Fill(Asymdat[l]);
    hmom[0]->Fill(Pdat[l]);
    hphtg[0]->Fill(Phtgdat[l]);
    hthtg[0]->Fill(Thtgdat[l]);
  }
  
  //simulation
  for(int k = 0; k < Qsq.size(); k++){
    qsq[1]->Fill(Qsq[k],Rate[k]);
    lab[1]->Fill(Lab[k],Rate[k]);
    asym[1]->Fill(Asym[k],Rate[k]);
    hmom[1]->Fill(Psim[k],Rate[k]);   
    hphtg[1]->Fill(Phtgsim[k],Rate[k]);   
    hthtg[1]->Fill(Thtgsim[k],Rate[k]);   
    hqsq_v->Fill(Qsq_v[k],Rate[k]);
    hlab_v->Fill(Lab_v[k],Rate[k]);
    hasym_v->Fill(Asym_v[k],Rate[k]);
  }
  
  auto leg = new TLegend(0.1,0.75,0.3,0.9);
  leg->AddEntry(lab[0],Form("Data (run%d)",run),"l");
//  leg->AddEntry(lab[1],Form("Sim"),"l");
    leg->AddEntry(lab[1],Form("Sim Sept Detuned %.1f Per",sept),"l");

  auto leg2 = new TLegend(0.1,0.75,0.3,0.9);
  leg2->AddEntry(hlab_v,"Vertex","l");
  leg2->AddEntry(lab[1],Form("Apparent"),"l");
  //  leg->AddEntry(lab[1],Form("Sim Sept Detuned %.1f Per",sept),"l");

  TCanvas *c1 = new TCanvas();
  lab[1]->Draw("hist");
  lab[0]->Scale(lab[1]->Integral()/lab[0]->Integral());
  lab[0]->Draw("HIST sames");
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
  if(iprint) {
    if (pinch!= 0.0) {
      c1->Print(Form("./temp/Run%4i_sept%.1f_Theta_pinch%.1f.pdf",run,sept,pinch));
    }else{
      c1->Print(Form("./temp/Run%4i_sept%.1f_Theta.pdf",run,sept));
    }
  }	    

  TCanvas *c2 = new TCanvas();
  qsq[1]->Draw("hist"); 
  qsq[0]->Scale(qsq[1]->Integral()/qsq[0]->Integral());
  qsq[0]->Draw("HIST sames");
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

  if(iprint) {
    if (pinch!= 0.0) {
      c2->Print(Form("./temp/Run%4i_sept%.1f_Q2_pinch%.1f.pdf",run,sept,pinch));
    }else{
      c2->Print(Form("./temp/Run%4i_sept%.1f_Q2.pdf",run,sept));
    }
  }	    
  
  TCanvas *c3 = new TCanvas();
  asym[1]->Draw("hist");
  asym[0]->Scale(asym[1]->Integral()/asym[0]->Integral());
  asym[0]->Draw("HIST sames");
  cout << "Underflow in asym in data:    " << nufasymdat << endl;
  cout << "Underflow in asym:    " << nufasym << endl;
  cout << "Underflow in asym_v:    " << nufasym_v << endl;
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
  if(iprint) {
    if (pinch!= 0.0) {
      c3->Print(Form("./temp/Run%4i_sept%.1f_Asym_pinch%.1f.pdf",run,sept,pinch));
    }else{
      c3->Print(Form("./temp/Run%4i_sept%.1f_Asym.pdf",run,sept));
    }
  }	    
  
  TCanvas *c5 = new TCanvas();
  c5->Divide(2,2);
  c5->cd(1);
  hlab_v->Draw("hist");
  lab[1]->Scale(hlab_v->Integral()/lab[1]->Integral());
  lab[1]->Draw("HIST sames");
  leg2->Draw();
  
  c5->cd(3);
  hqsq_v->Draw("hist");
  qsq[1]->Scale(hqsq_v->Integral()/qsq[1]->Integral());
  qsq[1]->Draw("HIST sames");
  leg2->Draw();
  
  c5->cd(4);
  asym[1]->Draw("hist");
  hasym_v->Scale(asym[1]->Integral()/hasym_v->Integral());
  hasym_v->Draw("HIST sames");
  leg2->Draw();
  if(iprint) {
    if (pinch!= 0.0) {
      c5->Print(Form("./temp/Run%4i_sept%.1f_vertex_pinch%.1f.pdf",run,sept,pinch));
    }else{
      c5->Print(Form("./temp/Run%4i_sept%.1f_vertex.pdf",run,sept));
    }
  }	    

  TCanvas *c4 = new TCanvas();
  c4->Divide(2,2);
  c4->cd(1);
  hmom[1]->Draw("hist");
  hmom[0]->Scale(hmom[1]->Integral()/hmom[0]->Integral());
  hmom[0]->Draw("HIST sames");
  leg->Draw();
  
  c4->cd(3);
  hphtg[1]->Draw("hist");
  hphtg[0]->Scale(hphtg[1]->Integral()/hphtg[0]->Integral());
  hphtg[0]->Draw("HIST sames");
  leg->Draw();
  
  c4->cd(4);
  hthtg[1]->Draw("hist");
  hthtg[0]->Scale(hthtg[1]->Integral()/hthtg[0]->Integral());
  hthtg[0]->Draw("HIST sames");
  leg->Draw();
  if(iprint) {
    if (pinch!= 0.0) {
      c4->Print(Form("./temp/Run%4i_sept%.1f_kin_pinch%.1f.pdf",run,sept,pinch));
    }else{
      c4->Print(Form("./temp/Run%4i_sept%.1f_kin.pdf",run,sept));
    }
  }	    
  
  TCanvas *c6 = new TCanvas();
  hmom[1]->Draw("hist");
  hmom[0]->Scale(hmom[1]->Integral()/hmom[0]->Integral());
  hmom[0]->Draw("HIST sames");
  leg->Draw();
  if(iprint) {
    if (pinch!= 0.0) {
      c6->Print(Form("./temp/Run%4i_sept%.1f_mom_pinch%.1f.pdf",run,sept,pinch));
    }else{
      c6->Print(Form("./temp/Run%4i_sept%.1f_mom.pdf",run,sept));
    }
  }	    
  
  
  cout<< Form("ngood data: %i       ngood sim: %i",ngdat,ngsim) << endl;
  std::cout << "Data Lab theta mean: " << lab[0]->GetMean() << "  " << "Data Lab theta rms: " << lab[0]->GetRMS() << std::endl;
  std::cout << "Sim  Lab theta mean: " << lab[1]->GetMean() << "   " << "Sim Lab theta rms: " << lab[1]->GetRMS() << std::endl;
  std::cout << "Data  Qsq mean: " << qsq[0]->GetMean() << "  " << "Data  Qsq rms: " << qsq[0]->GetRMS() << std::endl;
  std::cout << "Sim   Qsq mean: " << qsq[1]->GetMean() << "   " << "Sim  Qsq rms: " << qsq[1]->GetRMS() << std::endl;
  std::cout << "Data Asymmetry  mean: " << asym[0]->GetMean() << "  " << "Data Asymmetry  rms: " << asym[0]->GetRMS() << std::endl;
  std::cout << "Sim  Asymmetry  mean: " << asym[1]->GetMean() << "   " << "Sim Asymmetry  rms: " << asym[1]->GetRMS() << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Vertex Lab theta mean: " << hlab_v->GetMean() << "   " << "Vertex Lab theta rms: " << hlab_v->GetRMS() << std::endl;
  std::cout << "Vertex Qsq mean: " << hqsq_v->GetMean() << "   " << "Vertex  Qsq rms: " << hqsq_v->GetRMS() << std::endl;
  std::cout << "Vertex Asymmetry  mean: " << hasym_v->GetMean() << "   " << "Vertex Asymmetry  rms: " << hasym_v->GetRMS() << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Data Momentum  mean: " << hmom[0]->GetMean() << "   " << "Data Momentum  rms: " << hmom[0]->GetRMS() << std::endl;
  std::cout << "Sim Momentum  mean: " << hmom[1]->GetMean() << "   " << "Sim Momentum  rms: " << hmom[1]->GetRMS() << std::endl;
  std::cout << "Data Phi_tg  mean: " << hphtg[0]->GetMean() << "   " << "Data Phi_tg  rms: " << hphtg[0]->GetRMS() << std::endl;
  std::cout << "Sim Phi_tg  mean: " << hphtg[1]->GetMean() << "   " << "Sim Phi_tg  rms: " << hphtg[1]->GetRMS() << std::endl;
  std::cout << "Data Th_tg  mean: " << hthtg[0]->GetMean() << "   " << "Data Th_tg  rms: " << hthtg[0]->GetRMS() << std::endl;
  std::cout << "Sim Th_tg  mean: " << hthtg[1]->GetMean() << "   " << "Sim Th_tg  rms: " << hthtg[1]->GetRMS() << std::endl;


  cout << endl;
  cout << Form("%i  %3.1f  %3.1f",run,sept,pinch) << endl;
  // ngdata ngsim  0Asym data   sim   sim_v
  cout<< Form("%i   %i   %i   %i  %i",ngdat,ngsim,nufasymdat, nufasym,nufasym_v) << endl; 
  // lab_d lab_s qsq_d Qsq_s asym_d asym_s   labrms_d labrms_s qsqrms_d Qsqrms_s asymrms_d asymrms_s 
  cout<< Form("%6.4f  %9.7f  %7.3f  %6.4f  %9.7f  %6.3f",
	      lab[0]->GetMean(), qsq[0]->GetMean(), asym[0]->GetMean(), 
	      lab[0]->GetRMS() , qsq[0]->GetRMS(),  asym[0]->GetRMS() ) << endl;
  cout<< Form("%6.4f  %9.7f  %7.3f  %6.4f  %9.7f  %6.3f",
	      lab[1]->GetMean(), qsq[1]->GetMean(), asym[1]->GetMean(), 
	      lab[1]->GetRMS() , qsq[1]->GetRMS(),  asym[1]->GetRMS() ) << endl;
  cout<< Form("%6.4f  %9.7f  %7.3f  %6.4f  %9.7f  %6.3f",
	      hlab_v->GetMean(), hqsq_v->GetMean(), hasym_v->GetMean(), 
	      hlab_v->GetRMS() , hqsq_v->GetRMS(),  hasym_v->GetRMS() ) << endl;
  cout<< Form("%7.5f  %6.4f  %6.4f  %7.5f  %8.6f  %6.4f",
	      hmom[0]->GetMean(), hphtg[0]->GetMean(), hthtg[0]->GetMean(), 
	      hmom[0]->GetRMS() , hphtg[0]->GetRMS(),  hthtg[0]->GetRMS() ) << endl;
  cout<< Form("%7.5f  %6.4f  %6.4f  %7.5f  %6.4f  %6.4f",
	      hmom[1]->GetMean(), hphtg[1]->GetMean(), hthtg[1]->GetMean(), 
	      hmom[1]->GetRMS() , hphtg[1]->GetRMS(),  hthtg[1]->GetRMS() ) << endl;

 ofstream outfile(Form("./TextFiles/septum_scan_Ycoll_scanNewReplay_run%d.csv",run),ios_base::app);

 outfile<<sept<<"\t"<<YColl_offset<<"\t"<<lab[0]->GetMean()<<"\t"<<lab[1]->GetMean()<<"\t"<<lab[1]->GetMean()/lab[0]->GetMean()*100<<"\t"<<lab[0]->GetRMS()<<"\t"<<lab[1]->GetRMS()<<"\t"<<lab[1]->GetRMS()/lab[0]->GetRMS()*100<<"\t"<<qsq[0]->GetMean()<<"\t"<<qsq[1]->GetMean()<<"\t"<<qsq[1]->GetMean()/qsq[0]->GetMean()*100<<"\t"<<qsq[0]->GetRMS()<<"\t"<<qsq[1]->GetRMS()<<"\t"<<qsq[1]->GetRMS()/qsq[0]->GetRMS()*100<<"\t"<<asym[0]->GetMean()<<"\t"<<asym[1]->GetMean()<<"\t"<<asym[1]->GetMean()/asym[0]->GetMean()*100<<"\t"<<asym[0]->GetRMS()<<"\t"<<asym[1]->GetRMS()<<"\t"<<asym[1]->GetRMS()/asym[0]->GetRMS()*100<<endl;
 outfile.close();
 gSystem->Exec(Form("pdfunite ./temp/*.pdf ./plots/septum_%.1f_Ycoll_%.3f_NewReplay.pdf",sept,YColl_offset));
 gSystem->Exec(Form("rm -rf ./temp/*.pdf"));

 }
