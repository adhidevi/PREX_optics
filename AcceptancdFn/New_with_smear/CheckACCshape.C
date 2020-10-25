#include "LoadACC.h"

void CheckACCshape(){

     // load acceptance table
     double accp_angle[100]={0}, accp[100]={0}, accp_err[100]={0},accp_angle_err[100]={0};
     double accp_angle_1[100]={0};
     double accp_angle_2[100]={0};
     int status = LoadACC("accfunction.csv",accp_angle, accp, accp_err);
     if(status==0) exit(0);  // asymmetry table doesn't exist

     TH1F *hACC = new TH1F("hACC","acceptance function",100,3,8);
     TH1F *hACC_1 = new TH1F("hACC_1","acceptance function rms -3%",100,3,8);
     TH1F *hACC_2 = new TH1F("hACC_2","acceptance function rms +3%",100,3,8);

     int nbin_th = 100000;
     double dtheta = (8.0-3.0)/(nbin_th*1.0);  // delta theta in radius
     for(int ii=0; ii<nbin_th; ii++){
        double thisAngle = 3.0 + dtheta*ii;
        double thisACC = FindACC(thisAngle,accp_angle,accp,accp_err,100);  // acceptance to the corresponding angle
        if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

	hACC->Fill(thisAngle,thisACC/1000.);
     }
   
     Double_t mean = hACC->GetMean();
     Double_t RMS = hACC->GetRMS();
     Double_t mean_bin = hACC->FindBin(mean);
     Double_t peak = hACC->GetBinContent(mean_bin);///hACC_N->GetBinContent(mean_bin);;
cout<<"orignal:  "<<peak<<"  "<<mean<<"  "<<RMS<<endl;

     // gaus shape
     TF1 *f_gaus = new TF1("f_gaus","gaus",3,8);
     f_gaus->SetParameters(peak, mean, RMS);
    
     // box shape
     double th_1 = mean - sqrt(3.)*RMS;
     double th_2 = mean + sqrt(3.)*RMS;
     TF1 *f_box = new TF1("f_box",Form("%f*((x>=%f && x<=%f)? 1 : 0)",peak,th_1,th_2),3,8);
     f_box->SetLineColor(6);
     f_gaus->SetLineColor(7);

     // smearing
     for(int ii=0; ii<100; ii++){
	accp_angle_1[ii] = mean + (accp_angle[ii]-mean)*0.97;
	accp_angle_2[ii] = mean + (accp_angle[ii]-mean)*1.03;
     }

     for(int ii=0; ii<nbin_th; ii++){
        double thisAngle = 3.0 + dtheta*ii;
        double thisACC1 = FindACC(thisAngle,accp_angle_1,accp,accp_err,100);  // acceptance to the corresponding angle
        if(thisACC1<0) {printf("Something wrong here: ACC= %f\n",thisACC1); exit(0);}

        double thisACC2 = FindACC(thisAngle,accp_angle_2,accp,accp_err,100);  // acceptance to the corresponding angle
        if(thisACC2<0) {printf("Something wrong here: ACC= %f\n",thisACC2); exit(0);}

	hACC_1->Fill(thisAngle,thisACC1/1000.);
	hACC_2->Fill(thisAngle,thisACC2/1000.);
     }
     mean = hACC_1->GetMean();
     RMS = hACC_1->GetRMS();
     mean_bin = hACC_1->FindBin(mean);
     peak = hACC_1->GetBinContent(mean_bin);///hACC_1_N->GetBinContent(mean_bin);;
cout<<"-3%:  "<<peak<<"  "<<mean<<"  "<<RMS<<endl;

     mean = hACC_2->GetMean();
     RMS = hACC_2->GetRMS();
     mean_bin = hACC_2->FindBin(mean);
     peak = hACC_2->GetBinContent(mean_bin);///hACC_2_N->GetBinContent(mean_bin);;
cout<<"+3%   "<<peak<<"  "<<mean<<"  "<<RMS<<endl;
    
     TCanvas *c1 = new TCanvas("c1","c1",1000,500);
     hACC->SetLineColor(1);
     hACC_1->SetLineColor(2);
     hACC_2->SetLineColor(4);

     hACC->SetLineWidth(2);
     hACC_1->SetLineWidth(2);
     hACC_2->SetLineWidth(2);

     hACC->Draw();
     hACC_1->Draw("same");
     hACC_2->Draw("same");
     f_gaus->Draw("same");
     f_box->Draw("same");

     TLegend *leg = new TLegend(0.65,0.50,0.85,0.75);
     leg->AddEntry(hACC,"original","L");
     leg->AddEntry(hACC_1,"rms -3%","L");
     leg->AddEntry(hACC_2,"rms +3%","L");
     leg->AddEntry(f_gaus,"gaus","L");
     leg->AddEntry(f_box,"box","L");
     leg->Draw();

     ofstream outfile_gaus;
     outfile_gaus.open("accfunction_gaus.csv");
     ofstream outfile_box;
     outfile_box.open("accfunction_box.csv");

     outfile_gaus<<"vertex angle,acceptance,stat_err"<<endl;
     outfile_box<<"vertex angle,acceptance,stat_err"<<endl;

     for(int ii=0; ii<100; ii++){
	double tmp_th = accp_angle[ii];
	double acc_gaus = f_gaus->Eval(tmp_th);	
	double acc_box = f_box->Eval(tmp_th);	

 	outfile_gaus<<tmp_th<<","<<acc_gaus<<","<<0<<endl;
 	outfile_box<<tmp_th<<","<<acc_box<<","<<0<<endl;
     }
     outfile_gaus.close();
     outfile_box.close();
}
