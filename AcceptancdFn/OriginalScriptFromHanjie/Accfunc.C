#include "SetCut.h"

void Accfunc(){
     TChain *T = new TChain("T");
     T->Add("/lustre19/expphy/volatile/halla/parity/ryanrich/Ebeam953/TargNominal/Zero1_SandwichLHRS_PREX_0.0_1.root");
     T->Add("/lustre19/expphy/volatile/halla/parity/ryanrich/Ebeam953/TargNominal/Zero1_SandwichLHRS_PREX_0.0_2.root");
     T->Add("/lustre19/expphy/volatile/halla/parity/ryanrich/Ebeam953/TargNominal/Zero1_SandwichLHRS_PREX_0.0_3.root");
     T->Add("/lustre19/expphy/volatile/halla/parity/ryanrich/Ebeam953/TargNominal/Zero1_SandwichLHRS_PREX_0.0_4.root");
     T->Add("/lustre19/expphy/volatile/halla/parity/ryanrich/Ebeam953/TargNominal/Zero1_SandwichLHRS_PREX_0.0_5.root");

     // elastic peak position
     TH1F *hpztarg = new TH1F("hpztarg","p_ztarg distribution with cut",200,945,954);
     T->Draw("p_ztarg>>hpztarg",(colCut+isPb)*"rate");
     hpztarg->Draw("HIST");

     Int_t pbin = hpztarg->GetMaximumBin();
     Double_t p_peak = hpztarg->GetBinCenter(pbin);
     cout<<"p peak:  "<<p_peak<<endl; 

     // dp cut for the accepted events
     TString dpcut = Form("(%f-p_ztarg)<2.2",p_peak);
     TCut DP = Form("%s",dpcut.Data());

     // dp cut for the incident events
     double delta_p_percent = 2./100.;
     double delta_p = p_peak*delta_p_percent;
     TString dpcut_inc = Form("(%f-p_ztarg)<%f",p_peak,delta_p);
     TCut DP_INC = Form("%s",dpcut_inc.Data());

     TCut ACC = DP+isPb+colCut+XCUT;
     TCut INC = isPb+DP_INC;

     int nbin = 100;

     TH1F *hinc_v = new TH1F("hinc_v","incident events angle distribution",nbin,3,8); 
     TH1F *hacc_v = new TH1F("hacc_v","accepted events angle distribution",nbin,3,8); 

     hinc_v->Sumw2();
     hacc_v->Sumw2();

     T->Draw("ev.Th>>hinc_v",INC*"rate");
     T->Draw("ev.Th>>hacc_v",ACC*"rate");

     TGraphErrors *gACC_v = new TGraphErrors();

     ofstream outfile;
     outfile.open("accfunction.csv");
     outfile<<"vertex angle,acceptance,stat_err"<<endl;

     int nn = 0;
     for(int ii=1; ii<nbin+1; ii++){
	Double_t ninc = hinc_v->GetBinContent(ii);
	Double_t th_bin = hinc_v->GetBinCenter(ii);
	if(ninc==0) {
	   outfile<<th_bin<<","<<0<<","<<0<<endl;
	   continue;
	}

	Double_t nacc = hacc_v->GetBinContent(ii);
	Double_t tmp_ratio = nacc/ninc;

        Double_t inc_sumw2_sqrt = hinc_v->GetBinError(ii); // sqrt(sumw2)
	Double_t tmp_err2 = tmp_ratio * (1.-tmp_ratio) * pow(inc_sumw2_sqrt,2)/pow(ninc,2);
	Double_t tmp_err = sqrt(tmp_err2);

	gACC_v->SetPoint(nn,th_bin,tmp_ratio);	
	gACC_v->SetPointError(nn,0,tmp_err);	

	outfile<<th_bin<<","<<tmp_ratio<<","<<tmp_err<<endl;
 	nn++;
     } 
     outfile.close();


     TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
     c1->Draw();
     TPad *p1 = new TPad("p1","p1",0.0,0.0,0.8,1.);
     p1->SetRightMargin(0.1);
     p1->Draw();
     p1->cd();

     gACC_v->SetMarkerStyle(8);
     gACC_v->SetMarkerColor(4);
     gACC_v->Draw("AP");
     gACC_v->SetTitle("acceptance function;ev.Th (deg)");

     c1->cd(0);
     TPad *p2 = new TPad("p2","p2",0.8,0.0,1.,1.);
     p2->SetLeftMargin(0.);
     p2->Draw();
     p2->cd();
     TLatex tex;
     tex.SetTextAlign(11);
     tex.SetTextSize(0.05);
     tex.DrawLatexNDC(0.0,0.8,"thrown-in cut:");
     tex.DrawLatexNDC(0.1,0.75,ispb);
     tex.DrawLatexNDC(0.1,0.7,dpcut_inc);
     tex.DrawLatexNDC(0.0,0.65,"accepted cut:");
     tex.DrawLatexNDC(0.1,0.6,ispb);
     tex.DrawLatexNDC(0.1,0.55,colcut);
     tex.DrawLatexNDC(0.1,0.5,dpcut);
     tex.DrawLatexNDC(0.1,0.45,xcut);

     TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
     c2->Divide(2,1);
     c2->cd(1);
     hinc_v->Draw();
     hinc_v->SetLineColor(4);
     hinc_v->SetTitle("incident angle distribution;ev.Th (deg);");

     c2->cd(2);
     hacc_v->Draw();
     hacc_v->SetLineColor(4);
     hacc_v->SetTitle("accepted angle distribution;ev.Th (deg);");


}
