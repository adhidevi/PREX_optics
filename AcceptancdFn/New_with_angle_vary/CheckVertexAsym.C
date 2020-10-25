#include "LoadAsym.h"
#include "LoadACC.h"
#include "SetCut.h"

//acc_shape = "shift";  // option: original, smear, box, gaus, shift
void CheckVertexAsym(TString acc_shape, double value){
     double ppb = 1e9;
     double amu_c2 = 0.931494028; // GeV

     TChain *T = new TChain("T");
     T->Add("/lustre/expphy/volatile/halla/parity/ryanrich/Ebeam953_4by6/TargNominal/Zero1_SandwichLHRS_PREX_0.0_1.root");
     T->Add("/lustre/expphy/volatile/halla/parity/ryanrich/Ebeam953_4by6/TargNominal/Zero1_SandwichLHRS_PREX_0.0_2.root");
     T->Add("/lustre/expphy/volatile/halla/parity/ryanrich/Ebeam953_4by6/TargNominal/Zero1_SandwichLHRS_PREX_0.0_3.root");
     T->Add("/lustre/expphy/volatile/halla/parity/ryanrich/Ebeam953_4by6/TargNominal/Zero1_SandwichLHRS_PREX_0.0_4.root");
     T->Add("/lustre/expphy/volatile/halla/parity/ryanrich/Ebeam953_4by6/TargNominal/Zero1_SandwichLHRS_PREX_0.0_5.root");

     TH1F *hpztarg = new TH1F("hpztarg","p_ztarg distribution with cut",200,945,954);
     T->Draw("p_ztarg>>hpztarg",(colCut+isPb+XCUT)*"rate");
     hpztarg->Draw("HIST");

     Int_t pbin = hpztarg->GetMaximumBin();
     Double_t p_peak = hpztarg->GetBinCenter(pbin);
     cout<<"p peak:  "<<p_peak<<endl;

     // dp cut for accepted events
     TString dpcut = Form("(%f-p_ztarg)<2.2",p_peak);
     TCut DP = Form("%s",dpcut.Data());

     // dp cut for the incident events
     double delta_p_percent = 2./100.;
     double delta_p = p_peak*delta_p_percent;

     // final cuts
     TCut ACC = DP+isPb+colCut+XCUT;

     // load acceptance table
     double accp_angle[100]={0}, accp[100]={0}, accp_err[100]={0};
//     TString acc_shape = "shift";  // option: original, smear, box, gaus, shift

     int status = LoadACC("accfunction.csv",accp_angle, accp, accp_err);
     if(status==0) exit(0);  // asymmetry table doesn't exist

     Double_t mean0=0, max0=0, rms0=0;
     GetACCshape("accfunction.csv",mean0,rms0,max0);

     // gaus shape
     TF1 *f_gaus = new TF1("f_gaus","gaus",3,8);
     f_gaus->SetParameters(max0, mean0+value, rms0);

     // box shape
     double th_1 = mean0 - sqrt(3.)*rms0;
     double th_2 = mean0 + sqrt(3.)*rms0;
     TF1 *f_box = new TF1("f_box",Form("%f*((x>=%f && x<=%f)? 1 : 0)",max0,th_1+value,th_2+value),3,8);

     // smearing: original + rms_percent
     double rms_percent = value/100.;//value can be 1, 2 or 3. Put this in argument
     double_t accp_angle_p[100]={0};
     for(int ii=0; ii<100; ii++){
	accp_angle_p[ii] = mean0+(accp_angle[ii]-mean0)*(1.+rms_percent); // +rms_percent
     }

     // shift
     double shift_th = value; //Value is in deg. 0.1, 0.2, 0.3, .... Put this in argument
     double_t accp_angle_sf[100]={0};
     for(int ii=0; ii<100; ii++){
	accp_angle_sf[ii] = accp_angle[ii]+shift_th; // +rms_percent
     }

     bool check_ACC = false; // check mean and rms of acceptance function
     TH1F *hACC_check = new TH1F("hACC_check","acceptance function distribution",100,3,8);
     if(check_ACC){
        int check_nbin_th = 100000;
        double dtheta = (8.0-3.0)/(check_nbin_th*1.0);  // delta theta in radius
        for(int ii=0; ii<check_nbin_th; ii++){
          double thisAngle = 3.0 + dtheta*ii;
          double thisACC = 0;
	  if(acc_shape == "original" ) thisACC = FindACC(thisAngle,accp_angle,accp,accp_err,100);  // original acceptance function
	  if(acc_shape == "box" ) thisACC = f_box->Eval(thisAngle);	
	  if(acc_shape == "gaus" ) thisACC = f_gaus->Eval(thisAngle);	
	  if(acc_shape == "smear" ) thisACC = FindACC(thisAngle,accp_angle_p,accp,accp_err,100);  
	  if(acc_shape == "shift" ) thisACC = FindACC(thisAngle,accp_angle_sf,accp,accp_err,100);

          if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

          hACC_check->Fill(thisAngle,thisACC/1000.);
        }

        Double_t mean = hACC_check->GetMean();
        Double_t RMS = hACC_check->GetRMS();
        Double_t mean_bin = hACC_check->FindBin(mean);
        Double_t peak = hACC_check->GetBinContent(mean_bin);
        if(acc_shape=="smear") cout<<Form("smear %.3f",rms_percent)<<" "<<"ACC function peak, mean, rms"<<"  "<<peak<<"  "<<mean<<"  "<<RMS<<endl;
        else {
	  if(acc_shape=="shift") cout<<Form("shift %.3f",shift_th)<<" "<<"ACC function peak, mean, rms"<<"  "<<peak<<"  "<<mean<<"  "<<RMS<<endl;
          else cout<<acc_shape<<" "<<"ACC function peak, mean, rms"<<"  "<<peak<<"  "<<mean<<"  "<<RMS<<endl;
	}
     }

      
     // load asymmetry table
     LoadTable("horpb.dat", 0);

     int nbins=200;
     TH1F *hasym_origin = new TH1F("hasym_origin","vertex asymmetry in g4hrs with ACC cut",nbins,200,1000);
     TH1F *hasym_acc1 = new TH1F("hasym_acc1","vertex asymmetry calculated from acceptance function with dp cut",nbins,200,1000);
     TH1F *hasym_e_acc = new TH1F("hasym_e_acc","elemental asymmetry with acceptance function applied",nbins,200,1000);

     TH1F *hth_origin = new TH1F("hth_origin","vertex theta in g4hrs with ACC cut",100,3,8);
     TH1F *hth_e_acc = new TH1F("hth_e_acc","elemental theta with acceptance function applied",100,3,8);

     TH1F *hq2_origin = new TH1F("hq2_origin","vertex Q2 in g4hrs with ACC cut",100,0,0.015);
     TH1F *hq2_e_acc = new TH1F("hq2_e_acc","elemental Q2 with acceptance function applied",100,0,0.015);

     T->Draw("ev.A>>hasym_origin",ACC*"rate");   // vertex asymmetry calculated in g4hrs with acceptance cuts applied
     T->Draw("ev.Th>>hth_origin",ACC*"rate");   // vertex theta calculated in g4hrs with acceptance cuts applied
     T->Draw("ev.Q2>>hq2_origin",ACC*"rate");   // vertex Q2 calculated in g4hrs with acceptance cuts applied

     double thisBeamE,thisTh,thisRate,p_ztarg,xfp;
     int thisA;
     T->SetBranchAddress("ev.Th",&thisTh);
     T->SetBranchAddress("ev.beamp",&thisBeamE);
     T->SetBranchAddress("ev.nuclA",&thisA);
     T->SetBranchAddress("rate",&thisRate);
     T->SetBranchAddress("p_ztarg",&p_ztarg);
     T->SetBranchAddress("x_fp_tr",&xfp);

     Double_t nentries = T->GetEntries();
     for(int ii=0; ii<nentries; ii++){
        T->GetEntry(ii);

	double thisAsym = Interpolate(thisBeamE,thisTh,0,1);  // vertex asymmetry (should be the same as ev.A)
	double thisACC = 0;
	if(acc_shape == "original" ) thisACC = FindACC(thisTh,accp_angle,accp,accp_err,100);  // original acceptance function
	if(acc_shape == "box" ) thisACC = f_box->Eval(thisTh);	
	if(acc_shape == "gaus" ) thisACC = f_gaus->Eval(thisTh);	
	if(acc_shape == "smear" ) thisACC = FindACC(thisTh,accp_angle_p,accp,accp_err,100); 
	if(acc_shape == "shift" ) thisACC = FindACC(thisTh,accp_angle_sf,accp,accp_err,100); 

	if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

	if( thisA==208 && (p_peak-p_ztarg)<delta_p) hasym_acc1->Fill(thisAsym*ppb,thisRate*thisACC);
     }

     double setBeamE = 0.9534; 
     int nbin_th = 10000;
     double dtheta = (8.0-3.0)/(nbin_th*1.0);  // delta theta in radius

     for(int ii=0; ii<nbin_th; ii++){
	double thisAngle = 3.0 + dtheta*ii;
	double thisAsym = Interpolate(setBeamE,thisAngle,0,1);  // elemental asymmetry
	double thisXS = Interpolate(setBeamE,thisAngle,0,0);  // elemental cross section

	double deg_to_rad = TMath::Pi()/180.;
	double thisEp = setBeamE/(1.+setBeamE/(208.*amu_c2)*(1.-TMath::Cos(thisAngle*deg_to_rad)));
	double thisQ2 = 4.*setBeamE*thisEp*pow(TMath::Sin(thisAngle*deg_to_rad/2.),2);

	double thisACC=0;
	if(acc_shape == "original" ) thisACC = FindACC(thisAngle,accp_angle,accp,accp_err,100);  // original acceptance function
	if(acc_shape == "box" ) thisACC = f_box->Eval(thisAngle);	
	if(acc_shape == "gaus" ) thisACC = f_gaus->Eval(thisAngle);	
	if(acc_shape == "smear" ) thisACC = FindACC(thisAngle,accp_angle_p,accp,accp_err,100);  
	if(acc_shape == "shift" ) thisACC = FindACC(thisAngle,accp_angle_sf,accp,accp_err,100);  

	if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

	double omega = TMath::Sin(thisAngle*deg_to_rad)*dtheta*deg_to_rad;
	hasym_e_acc->Fill(thisAsym*ppb,thisXS*omega*thisACC);   // elemetal asymmetry with the acceptance function applied
	hth_e_acc->Fill(thisAngle,thisXS*omega*thisACC);   // elemetal angle with the acceptance function applied
	hq2_e_acc->Fill(thisQ2,thisXS*omega*thisACC);   // elemetal Q2 with the acceptance function applied
     } 

     cout<<"================================="<<endl;
     cout<<"orignal vertex asymmetry mean:     "<<setw(8)<<hasym_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_origin->GetRMS()<<endl;
     cout<<"reproduce vertex asymmetry mean:   "<<setw(8)<<hasym_acc1->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_acc1->GetRMS()<<endl;
     cout<<"elemental asymmetry mean:          "<<setw(8)<<hasym_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;
     cout<<"orignal vertex angle mean:   "<<setw(8)<<hth_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hth_origin->GetRMS()<<endl;
     cout<<"elemental angle mean:        "<<setw(8)<<hth_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hth_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;
     cout<<"orignal vertex Q2 mean:   "<<setw(8)<<hq2_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hq2_origin->GetRMS()<<endl;
     cout<<"elemental Q2 mean:        "<<setw(8)<<hq2_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hq2_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;

     Double_t inte1 = hasym_origin->GetMaximum();
     Double_t inte2 = hasym_e_acc->GetMaximum();
     Double_t inte3 = hasym_acc1->GetMaximum();
     hasym_e_acc->Scale(inte1/inte2);
     hasym_acc1->Scale(inte1/inte3);

     TCanvas *c1 = new TCanvas("c1","c1",800,600);
     hasym_origin->SetLineColor(1);
     hasym_acc1->SetLineColor(6);
     hasym_e_acc->SetLineColor(2);

     hasym_origin->SetLineWidth(2);
     hasym_acc1->SetLineWidth(2);
     hasym_e_acc->SetLineWidth(2);

     if(inte1>inte3){
       hasym_origin->Draw("HIST");    
       hasym_acc1->Draw("HIST sames");    
       hasym_e_acc->Draw("HIST sames");
       hasym_origin->SetTitle(Form("Asymmetry with acceptance (dp_inc/p=%.3f);asym(ppb)",delta_p_percent));
     }
     else{
       hasym_acc1->Draw("HIST");    
       hasym_origin->Draw("HIST sames");    
       hasym_e_acc->Draw("HIST sames");
       hasym_acc1->SetTitle("Asymmetry with acceptance; asym(ppb)");
     }

     TLegend *leg = new TLegend(0.14,0.65,0.34,0.80);
     leg->SetTextSize(0.06);
     leg->SetBorderSize(0);
     leg->SetFillStyle(0);
     leg->SetFillColor(0);
     leg->AddEntry(hasym_origin,"origin","L");
     leg->AddEntry(hasym_acc1,"reproduce","L");
     leg->AddEntry(hasym_e_acc,"elemental","L");
     leg->Draw();

     gPad->Update();
     TPaveStats* statA1 = (TPaveStats*)hasym_origin->FindObject("stats");
     TPaveStats* statA2 = (TPaveStats*)hasym_acc1->FindObject("stats");
     TPaveStats* statA3 = (TPaveStats*)hasym_e_acc->FindObject("stats");
     statA1->SetY1NDC(0.90);
     statA1->SetY2NDC(0.75);
     statA2->SetY1NDC(0.75);
     statA2->SetY2NDC(0.60);
     statA3->SetY1NDC(0.60);
     statA3->SetY2NDC(0.45);
     statA1->SetTextColor(1);
     statA2->SetTextColor(6);
     statA3->SetTextColor(2);
     gPad->Modified();

     TLatex tex1;
     tex1.SetTextAlign(11);
     tex1.SetTextSize(0.05);
     if(acc_shape=="smear") tex1.DrawLatexNDC(0.65,0.3,Form("smear %.2f %s",rms_percent*100.,"%"));
     else {
	if(acc_shape=="shift") tex1.DrawLatexNDC(0.65,0.3,Form("shift %.2f deg",shift_th));
        else tex1.DrawLatexNDC(0.65,0.3,Form("#splitline{%s acceptance}{%.1f deg shift}",acc_shape.Data(),value));
     }
    
     inte1 = hth_origin->GetMaximum();
     inte2 = hth_e_acc->GetMaximum();
     hth_e_acc->Scale(inte1/inte2);
     c1->SaveAs(Form("./temp1/asymmetry_%s_%.1f.pdf",acc_shape.Data(),value));

     TCanvas *c2 = new TCanvas("c2","c2",800,600);
     c2->Divide(2,1);
     c2->cd(1);
     hth_origin->SetLineColor(1);
     hth_e_acc->SetLineColor(2);

     hth_origin->SetLineWidth(2);
     hth_e_acc->SetLineWidth(2);

     hth_origin->Draw("HIST");    
     hth_e_acc->Draw("HIST sames");
     hth_origin->SetTitle("angle with acceptance; angle");

     TLegend *leg1 = new TLegend(0.48,0.65,0.68,0.75);
     leg1->SetTextSize(0.06);
     leg1->SetBorderSize(0);
     leg1->SetFillStyle(0);
     leg1->SetFillColor(0);
     leg1->AddEntry(hth_origin,"origin","L");
     leg1->AddEntry(hth_e_acc,"elemental","L");
     leg1->Draw();

     gPad->Update();
     TPaveStats* statB1 = (TPaveStats*)hth_origin->FindObject("stats");
     TPaveStats* statB2 = (TPaveStats*)hth_e_acc->FindObject("stats");
     statB1->SetY1NDC(0.90);
     statB1->SetY2NDC(0.75);
     statB2->SetY1NDC(0.75);
     statB2->SetY2NDC(0.60);
     statB1->SetTextColor(1);
     statB2->SetTextColor(2);
     gPad->Modified();
    
     c2->cd(2);
     inte1 = hq2_origin->GetMaximum();
     inte2 = hq2_e_acc->GetMaximum();
     hq2_e_acc->Scale(inte1/inte2);

     hq2_origin->SetLineColor(1);
     hq2_e_acc->SetLineColor(2);

     hq2_origin->SetLineWidth(2);
     hq2_e_acc->SetLineWidth(2);

     hq2_origin->Draw("HIST");    
     hq2_e_acc->Draw("HIST sames");
     hq2_origin->SetTitle("Q2 with acceptance; Q2");

     TLegend *leg2 = new TLegend(0.48,0.65,0.68,0.75);
     leg2->SetTextSize(0.06);
     leg2->SetBorderSize(0);
     leg2->SetFillStyle(0);
     leg2->SetFillColor(0);
     leg2->AddEntry(hq2_origin,"origin","L");
     leg2->AddEntry(hq2_e_acc,"elemental","L");
     leg2->Draw();

     gPad->Update();
     TPaveStats* statC1 = (TPaveStats*)hq2_origin->FindObject("stats");
     TPaveStats* statC2 = (TPaveStats*)hq2_e_acc->FindObject("stats");
     statC1->SetY1NDC(0.90);
     statC1->SetY2NDC(0.75);
     statC2->SetY1NDC(0.75);
     statC2->SetY2NDC(0.60);
     statC1->SetTextColor(1);
     statC2->SetTextColor(2);
     gPad->Modified();
    
     TLatex tex2;
     tex2.SetTextAlign(11);
     tex2.SetTextSize(0.05);
     if(acc_shape=="smear") tex2.DrawLatexNDC(0.60,0.3,Form("smear %.2f %s",rms_percent*100.,"%"));
     else {
	if(acc_shape=="shift") tex2.DrawLatexNDC(0.60,0.3,Form("shift %.2f deg",shift_th));
        else tex2.DrawLatexNDC(0.5,0.5,Form("#splitline{%s acceptance}{%.1f deg shift}",acc_shape.Data(),value));
     }
  
     c2->SaveAs(Form("./temp1/angle_qsq_%s_%.1f.pdf",acc_shape.Data(),value));

     if(check_ACC){
       TCanvas *c3 = new TCanvas("c3","c3",800,600);
       hACC_check->Draw();
     c3->SaveAs(Form("./temp1/acc_shape_%s_%.1f.pdf",acc_shape.Data(),value));       
     }
     ofstream outfile(Form("./TextFiles/acceptance_%s_shiftTH.csv",acc_shape.Data()),ios_base::app);
     outfile<<acc_shape<<"\t"<<value<<"\t"<<hasym_origin->GetMean()<<"\t"<<hasym_acc1->GetMean()<<"\t"<<hasym_e_acc->GetMean()<<"\t"<<hth_origin->GetMean()<<"\t"<<hth_e_acc->GetMean()<<"\t"<<hq2_origin->GetMean()<<"\t"<<hq2_e_acc->GetMean()<<"\t"<<hasym_origin->GetRMS()<<"\t"<<hasym_acc1->GetRMS()<<"\t"<<hasym_e_acc->GetRMS()<<"\t"<<hth_origin->GetRMS()<<"\t"<<hth_e_acc->GetRMS()<<"\t"<<hq2_origin->GetRMS()<<"\t"<<hq2_e_acc->GetRMS()<<endl;
     outfile.close();

}
