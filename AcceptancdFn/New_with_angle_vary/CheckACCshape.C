#include "LoadACC.h"
//acc_shape=="smear" for rms change or "shift" for angle shift or "box" or "gaus"
//value==3 for +/,value,"%")-3% rms change, 0.5 for +/-0.5deg angle change
void CheckACCshape(TString acc_shape, double value){

     // load acceptance table
     double accp_angle[100]={0}, accp[100]={0}, accp_err[100]={0},accp_angle_err[100]={0};
     double accp_angle_1[100]={0};
     double accp_angle_2[100]={0};
     double accp_angle_3[100]={0};
     double accp_angle_4[100]={0};
     int status = LoadACC("accfunction.csv",accp_angle, accp, accp_err);
     if(status==0) exit(0);  // asymmetry table doesn't exist

     TH1F *hACC = new TH1F("hACC","acceptance function",100,3,8);
     TH1F *hACC_1 = new TH1F("hACC_1",Form("acceptance function rms -%.1f %s",value,"%"),100,3,8);
     TH1F *hACC_2 = new TH1F("hACC_2",Form("acceptance function rms +%.1f %s",value,"%"),100,3,8);
     TH1F *hACC_3 = new TH1F("hACC_3",Form("acceptance function shift +%.1f deg",value),100,3,8);
     TH1F *hACC_4 = new TH1F("hACC_4",Form("acceptance function shift -%.1f deg",value),100,3,8);

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
     f_gaus->SetLineColor(2);
     TF1* f_gausM = new TF1("f_gausM","gaus",3,8);
     f_gausM->SetParameters(peak, mean-value,RMS);
     f_gausM->SetLineColor(4);
     TF1* f_gausP = new TF1("f_gausP","gaus",3,8);
     f_gausP->SetParameters(peak, mean+value,RMS);
     f_gausP->SetLineColor(6);
    
     // box shape
     double th_1 = mean - sqrt(3.)*RMS;
     double th_2 = mean + sqrt(3.)*RMS;
     TF1 *f_box = new TF1("f_box",Form("%f*((x>=%f && x<=%f)? 1 : 0)",peak,th_1,th_2),3,8);
     f_box->SetLineColor(2);
     TF1 *f_boxM = new TF1("f_boxM",Form("%f*((x>=%f && x<=%f)? 1 : 0)",peak,th_1-value,th_2-value),3,8);
     f_boxM->SetLineColor(4);
     TF1 *f_boxP = new TF1("f_boxP",Form("%f*((x>=%f && x<=%f)? 1 : 0)",peak,th_1+value,th_2+value),3,8);
     f_boxP->SetLineColor(6);

     // smearing
     for(int ii=0; ii<100; ii++){
	accp_angle_1[ii] = mean + (accp_angle[ii]-mean)*(1-value/100.);
	accp_angle_2[ii] = mean + (accp_angle[ii]-mean)*(1+value/100.);
	accp_angle_3[ii] = accp_angle[ii]+value;
	accp_angle_4[ii] = accp_angle[ii]-value;
     }

     for(int ii=0; ii<nbin_th; ii++){
        double thisAngle = 3.0 + dtheta*ii;
        double thisACC1 = FindACC(thisAngle,accp_angle_1,accp,accp_err,100);  // acceptance to the corresponding angle
        if(thisACC1<0) {printf("Something wrong here 1: ACC= %f\n",thisACC1); exit(0);}

        double thisACC2 = FindACC(thisAngle,accp_angle_2,accp,accp_err,100);  // acceptance to the corresponding angle
        if(thisACC2<0) {printf("Something wrong here 2: ACC= %f\n",thisACC2); exit(0);}

        double thisACC3 = FindACC(thisAngle,accp_angle_3,accp,accp_err,100);  // acceptance to the corresponding angle
        if(thisACC3<0) {printf("Something wrong here 3: ACC= %f\n",thisACC3); exit(0);}

        double thisACC4 = FindACC(thisAngle,accp_angle_4,accp,accp_err,100);  // acceptance to the corresponding angle
        if(thisACC4<0) {printf("Something wrong here: ACC= %f\n",thisACC4); exit(0);}

	hACC_1->Fill(thisAngle,thisACC1/1000.);
	hACC_2->Fill(thisAngle,thisACC2/1000.);
	hACC_3->Fill(thisAngle,thisACC3/1000.);
	hACC_4->Fill(thisAngle,thisACC4/1000.);
     }
     mean = hACC_1->GetMean();
     RMS = hACC_1->GetRMS();
     mean_bin = hACC_1->FindBin(mean);
     peak = hACC_1->GetBinContent(mean_bin);///hACC_1_N->GetBinContent(mean_bin);
    if(acc_shape=="smear")
cout<<Form("-%.1f %s:  ",value,"%")<<peak<<"  "<<mean<<"  "<<RMS<<endl;

     mean = hACC_2->GetMean();
     RMS = hACC_2->GetRMS();
     mean_bin = hACC_2->FindBin(mean);
     peak = hACC_2->GetBinContent(mean_bin);///hACC_2_N->GetBinContent(mean_bin);
    if(acc_shape=="smear")
cout<<Form("+%.1f %s   ",value,"%")<<peak<<"  "<<mean<<"  "<<RMS<<endl;

     mean = hACC_3->GetMean();
     RMS = hACC_3->GetRMS();
     mean_bin = hACC_3->FindBin(mean);
     peak = hACC_3->GetBinContent(mean_bin);///hACC_3_N->GetBinContent(mean_bin);
    if(acc_shape=="shift")
cout<<Form("+%.1f   ",value)<<peak<<"  "<<mean<<"  "<<RMS<<endl;

     mean = hACC_4->GetMean();
     RMS = hACC_4->GetRMS();
     mean_bin = hACC_4->FindBin(mean);
     peak = hACC_4->GetBinContent(mean_bin);///hACC_4_N->GetBinContent(mean_bin);
    if(acc_shape=="shift")
cout<<Form("-%.1f   ",value)<<peak<<"  "<<mean<<"  "<<RMS<<endl;
    
     TCanvas *c1 = new TCanvas("c1","c1",1000,500);
     hACC->SetLineColor(1);
     hACC_1->SetLineColor(2);
     hACC_2->SetLineColor(4);
     hACC_3->SetLineColor(2);
     hACC_4->SetLineColor(4);

     hACC->SetLineWidth(2);
     hACC_1->SetLineWidth(2);
     hACC_2->SetLineWidth(2);
     hACC_3->SetLineWidth(2);
     hACC_4->SetLineWidth(2);

     hACC->Draw();
     if(acc_shape=="smear"){
     hACC_1->Draw("same");
     hACC_2->Draw("same");
     }else if(acc_shape=="shift"){
     hACC_3->Draw("same");
     hACC_4->Draw("same");
     }else if(acc_shape=="gaus"){
     f_gaus->Draw("same");
     f_gausM->Draw("same");
     f_gausP->Draw("same");
     }else if(acc_shape=="box"){
     f_box->Draw("same");
     f_boxM->Draw("same");
     f_boxP->Draw("same");
     }
     TLegend *leg = new TLegend(0.70,0.50,0.80,0.70);
     leg->AddEntry(hACC,"original","L");
     leg->SetTextSize(0.06);
     leg->SetBorderSize(0);
     leg->SetFillStyle(0);
     leg->SetFillColor(0);
     if(acc_shape=="smear"){
     leg->AddEntry(hACC_1,Form("rms -%.1f %s",value,"%"),"L");
     leg->AddEntry(hACC_2,Form("rms +%.1f %s",value,"%"),"L");
     }else if(acc_shape=="shift"){
     leg->AddEntry(hACC_3,Form("shift +%.1f deg",value),"L");
     leg->AddEntry(hACC_4,Form("shift -%.1f deg",value),"L");
     }else if(acc_shape=="box"){
     leg->AddEntry(f_box,"Original Box","L");
     leg->AddEntry(f_boxP,Form("shift +%.1f deg",value),"L");
     leg->AddEntry(f_boxM,Form("shift -%.1f deg",value),"L");
     }else if(acc_shape=="gaus"){
     leg->AddEntry(f_gaus,"Original Gaus","L");
     leg->AddEntry(f_gausP,Form("shift +%.1f deg",value),"L");
     leg->AddEntry(f_gausM,Form("shift -%.1f deg",value),"L");
     }
     leg->Draw();
     c1->SaveAs(Form("./temp1/acceptanceFn_%s_%.1f.pdf",acc_shape.Data(),value));
/*
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
*/
}
