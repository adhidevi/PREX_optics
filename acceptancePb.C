//run as analyzer -l 'acceptancePb.C(2120,"208Pb")' for example
void acceptancePb(int run,TString target){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetTitleFontSize(0.08);
	gROOT->ForceStyle();
	TGaxis::SetMaxDigits(3);

//LHRS upstream adc cut
	double upadc_cutL = 480;
//RHRS upstream adc cut
	double upadc_cutR = 502;
	TChain* T = new TChain("T");
	if(run<10000){
//	T->Add(Form("/lustre19/expphy/volatile/halla/parity/ryanrich/prexRootFiles/prexLHRS_%d_-*.root",run));
	T->Add(Form("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexLHRS_%d_-*.root",run));
//trigger cut
//	TCut trig_cut = "";
	TCut trig_cut = "(P.evtypebits&2)==2";
//VDC track cut
	TCut vdc_cut = "L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1";
//track cut on theta and phi
	TCut tr_cut = "(L.tr.th[0]<0.05&&L.tr.th[0]>-0.2&&L.tr.ph[0]<0.1&&L.tr.ph[0]>-0.1)";
//track cut on target theta and target phi
	TCut tg_cut = "(L.tr.tg_th[0]<0.055&&L.tr.tg_th[0]>-0.055&&L.tr.tg_ph[0]>-0.018&&L.tr.tg_ph[0]<0.026)";
//FASTBUS adc cut on upstream cut
	TCut adc_cut_up = Form("P.upQadcL>%f",upadc_cutL);
	TCut adc_cut_up_ped = Form("P.upQadcL<%f",upadc_cutL);
	TCut cut_basic = trig_cut+vdc_cut;
	TCut cut = cut_basic+tr_cut+tg_cut;
	TCut cut_wadc = cut+adc_cut_up;
	TCut cut_wped = cut+adc_cut_up_ped;
	TCut x_cut = "(L.tr.x[0]+0.9*L.tr.th[0]) > -0.0665";
	
	double x_low = 940;
	double x_hi = 952;
	double nbins = 200;
	double bin_size = (x_hi-x_low)/nbins;
	TH1F* hp = new TH1F("hp",Form("LHRS %s spectrum no adc cut (run%d);L.gold.p (MeV/c)",target.Data(),run),nbins,x_low,x_hi);
	TH1F* hpc = new TH1F("hpc",Form("LHRS %s spectrum with adc cut (run%d);L.gold.p (MeV/c)",target.Data(),run),nbins,x_low,x_hi);
	TH1F* hpp = new TH1F("hpp",Form("LHRS %s spectrum with < adc cut (run%d);L.gold.p (MeV/c)",target.Data(),run),nbins,x_low,x_hi);
	TH1F* hxcut = new TH1F("hxcut",Form("LHRS %s spectrum with adc and x cut (run%d);L.gold.p (MeV/c)",target.Data(),run),nbins,x_low,x_hi);
	T->Draw("1000*L.gold.p>>hp",cut);
	T->Draw("1000*L.gold.p>>hpc",cut_wadc);
	T->Draw("1000*L.gold.p>>hpp",cut_wped);
	T->Draw("1000*L.gold.p>>hxcut",cut_wadc+x_cut);
	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	c1->Divide(1,3,0.0001,0.00001);
	c1->cd(2);
	gPad->SetGrid();
	gPad->SetLogy();
	gPad->SetBottomMargin(0);
	gPad->SetTopMargin(0);
	hpc->SetLineColor(kBlack);
	hpp->SetLineColor(kOrange);
	hpc->Draw("HISTO");
	hpp->Draw("HISTO SAME");
	hxcut->SetLineColor(kOrange);
//	hxcut->Draw("same");
	hpc->SetLabelSize(0.05,"y");
	double e_peak = hpc->GetBinCenter(hpc->GetMaximumBin());
	const int nline = 5;
	TLine* line[nline];
	TPaveLabel* lbl[nline];
	TString level[] = {"Elastic","1^{st} excited 3^{-}","2^{nd} excited 5^{-}","3^{rd} excited 5^{-}","4^{th} excited 5^{-}"};
	double firstE = 2.615;
	double secondE = 3.198;
	double thirdE = 3.709;
	double fourthE = 3.961;
	line[0] = new TLine(e_peak,0,e_peak,hp->GetMaximum());
	line[1] = new TLine(e_peak-firstE,0,e_peak-firstE,hp->GetMaximum());
	line[2] = new TLine(e_peak-secondE,0,e_peak-secondE,hp->GetMaximum());
	line[3] = new TLine(e_peak-thirdE,0,e_peak-thirdE,hp->GetMaximum());
	line[4] = new TLine(e_peak-fourthE,0,e_peak-fourthE,hp->GetMaximum());
	lbl[0] = new TPaveLabel(0.65,0.80,0.95,0.90,level[0].Data(),"NDC");
	lbl[1] = new TPaveLabel(0.65,0.70,0.95,0.80,level[1].Data(),"NDC");
	lbl[2] = new TPaveLabel(0.65,0.60,0.95,0.70,level[2].Data(),"NDC");
	lbl[3] = new TPaveLabel(0.65,0.50,0.95,0.60,level[3].Data(),"NDC");
	lbl[4] = new TPaveLabel(0.65,0.40,0.95,0.50,level[4].Data(),"NDC");
	for(int iline=0;iline<nline;iline++){
	if(iline==3){
	line[iline]->SetLineColor(30);
	lbl[iline]->SetTextColor(30);
	}else{
	line[iline]->SetLineColor(iline+2);
	lbl[iline]->SetTextColor(iline+2);
	}
	lbl[iline]->SetBorderSize(0);
	lbl[iline]->SetFillStyle(0);
	lbl[iline]->SetFillColor(0);
	line[iline]->SetLineWidth(2);
	line[iline]->Draw();
	}
	TPaveLabel* accepted = new TPaveLabel(0.65,0.75,0.90,0.85,"Accepted by Quartz","NDC");
	accepted->SetBorderSize(0);
	accepted->SetFillStyle(0);
	accepted->SetFillColor(0);
	accepted->SetTextColor(kBlack);
	accepted->Draw();
	TPaveLabel* missed = new TPaveLabel(0.65,0.65,0.90,0.75,"Missed by Quartz","NDC");
	missed->SetBorderSize(0);
	missed->SetFillStyle(0);
	missed->SetFillColor(0);
	missed->SetTextColor(kOrange);
	missed->Draw();
	TPaveLabel* x_accept = new TPaveLabel(0.60,0.75,0.95,0.85,"adc+quartz edge cut","NDC");
	x_accept->SetBorderSize(0);
	x_accept->SetFillStyle(0);
	x_accept->SetFillColor(0);
	x_accept->SetTextColor(kOrange);
//	x_accept->Draw();
	cout<<Form("Elastic: %1.4f",e_peak)<<endl;
	c1->cd(1);
	gPad->SetGrid();
	gPad->SetLogy();
	gPad->SetBottomMargin(0);
	hp->SetLineColor(kBlack);
	hp->Draw("HISTO");
	hp->SetLabelSize(0.05,"y");
	for(int iline=0;iline<nline;iline++){
	line[iline]->SetLineWidth(2);
	line[iline]->Draw();
	lbl[iline]->Draw();
	}
	c1->cd(3);
	gPad->SetGrid();
	gPad->SetLogy();
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0.15);
//	hpc->Scale(1/hpc->Integral("width"));
//	hp->Scale(1/hp->Integral("width"));
//	TH1F* hACC = (TH1F*)hpc->Clone();
//	hACC->Divide(hp);
//	hACC->Draw("HISTO");
	TH1F* hACC = (TH1F*)hpc->Clone();
	TH1F* hACC_xcut = (TH1F*)hxcut->Clone();
	hACC->SetTitle("Fraction (Accepted/Total);L.gold.p (MeV/c)");
	hACC->SetTitleSize(0.08,"x");
	hACC->SetLabelSize(0.05,"x;y");
	hACC->Divide(hp);
	hACC_xcut->Divide(hp);
	hACC->SetLineColor(kBlack);
	hACC->SetMarkerStyle(20);
	hACC->SetMarkerSize(0.4);
	hACC->SetTitle(";R.gold.p (MeV/c);");
	hACC->Draw("prof");
	TPaveLabel* title3 = new TPaveLabel(0.55,0.60,0.90,0.70,"Fraction (Accepted/Total)","NDC");
	title3->SetBorderSize(0);
	title3->SetFillColor(0);
	title3->SetFillStyle(0);
	title3->SetTextColor(kBlack);
	title3->Draw();
//	hACC_xcut->Draw("same");
	hACC->SetTitleOffset(0.7,"x");
	for(int iline=0;iline<nline;iline++){
	line[iline]->SetLineWidth(2);
	line[iline]->Draw();
	}
	c1->SaveAs(Form("./plots/acceptance_function_%s_run%d.pdf",target.Data(),run));

	TCanvas* c0 = new TCanvas("c0","c0",1000,700);
	TH1F* hpCpy = (TH1F*)hp->Clone("hpCpy");
	hpCpy->SetTitle(Form("LHRS %s spectrum no adc cut (run%d);L.gold.p (MeV/c);",target.Data(),run));
	hpCpy->Draw("HISTO");
	TF1* fitG = new TF1("fitG","gaus",945.9,947.1);
	hpCpy->Fit(fitG,"R0");
	fitG->SetLineColor(kRed);
	fitG->SetLineWidth(2);
	fitG->Draw("same");
	
	hpCpy->Scale(1./hpCpy->Integral(),"width");
	double firstEfr = hpCpy->Interpolate(e_peak-firstE)/hpCpy->Interpolate(e_peak);
	double secondEfr = hpCpy->Interpolate(e_peak-secondE)/hpCpy->Interpolate(e_peak);
	double thirdEfr = hpCpy->Interpolate(e_peak-thirdE)/hpCpy->Interpolate(e_peak);
	double fourthEfr = hpCpy->Interpolate(e_peak-fourthE)/hpCpy->Interpolate(e_peak);
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.06);
	latex.SetTextColor(3);
	latex.DrawLatex(0.7,0.85,Form("1st inel = %2.4f",firstEfr));
	latex.SetTextColor(4);
	latex.DrawLatex(0.7,0.8,Form("2nd inel = %2.4f",secondEfr));
	latex.SetTextColor(30);
	latex.DrawLatex(0.7,0.75,Form("3rd inel = %2.4f",thirdEfr));
	latex.SetTextColor(6);
	latex.DrawLatex(0.7,0.7,Form("4th inel = %2.4f",fourthEfr));
	for(int iline=0;iline<nline;iline++){
	line[iline]->Draw();
	}
	}else{
	T->Add(Form("./prex_counting/prexRHRS_%d_-1*.root",run));
//trigger cut
//	TCut trig_cut = "";
	TCut trig_cut = "(P.evtypebits&2)==2";
//VDC track cut
	TCut vdc_cut = "R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1";
//track cut on theta and phi
	TCut tr_cut = "(R.tr.th[0]<0.05&&R.tr.th[0]>-0.2&&R.tr.ph[0]<0.1&&R.tr.ph[0]>-0.1)";
//track cut on target theta and target phi
	TCut tg_cut = "(R.tr.tg_th[0]<0.055&&R.tr.tg_th[0]>-0.055&&R.tr.tg_ph[0]>-0.018&&R.tr.tg_ph[0]<0.026)";
//FASTBUS adc cut on upstream cut
	TCut adc_cut_up = Form("P.upQadcR>%f",upadc_cutR);
	TCut adc_cut_up_ped = Form("P.upQadcR<%f",upadc_cutR);
	TCut cut_basic = trig_cut+vdc_cut;
	TCut cut = cut_basic+tr_cut+tg_cut;
	TCut cut_wadc = cut+adc_cut_up;
	TCut cut_wped = cut+adc_cut_up_ped;
	TCut x_cut = "(R.tr.x[0]+0.9*R.tr.th[0]) > -0.0725";

	TH1F* hp = new TH1F("hp",Form("RHRS %s spectrum no adc cut (run%d);R.gold.p (MeV/c)",target.Data(),run),200,940,952);
	TH1F* hpc = new TH1F("hpc",Form("RHRS %s spectrum with adc cut (run%d);R.gold.p (MeV/c)",target.Data(),run),200,940,952);
	TH1F* hpp = new TH1F("hpp",Form("RHRS %s spectrum with < adc cut (run%d);R.gold.p (MeV/c)",target.Data(),run),200,940,952);
	TH1F* hxcut = new TH1F("hxcut",Form("RHRS %s spectrum with adc and x cut (run%d);R.gold.p (MeV/c)",target.Data(),run),200,940,952);
	T->Draw("1000*R.gold.p>>hp",cut);
	T->Draw("1000*R.gold.p>>hpc",cut_wadc);
	T->Draw("1000*R.gold.p>>hpp",cut_wped);
	T->Draw("1000*R.gold.p>>hxcut",cut_wadc+x_cut);
	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	c1->Divide(1,3,0.0001,0.00001);
	c1->cd(2);
	gPad->SetGrid();
	gPad->SetLogy();
	gPad->SetBottomMargin(0);
	gPad->SetTopMargin(0);
	hpc->SetLineColor(kBlack);
	hpp->SetLineColor(kOrange);
	hpc->Draw("HISTO");
	hpp->Draw("HISTO SAME");
	hxcut->SetLineColor(kOrange);
//	hxcut->Draw("same");
	hpc->SetLabelSize(0.05,"y");
	double e_peak = hpc->GetBinCenter(hpc->GetMaximumBin());
	const int nline = 5;
	TLine* line[nline];
	TPaveLabel* lbl[nline];
	TString level[] = {"Elastic","1^{st} excited 3^{-}","2^{nd} excited 5^{-}","3^{rd} excited 5^{-}","4^{th} excited 5^{-}"};
	line[0] = new TLine(e_peak,0,e_peak,1e3);
	line[1] = new TLine(e_peak-2.615,0,e_peak-2.615,1e3);
	line[2] = new TLine(e_peak-3.198,0,e_peak-3.198,1e3);
	line[3] = new TLine(e_peak-3.709,0,e_peak-3.709,1e3);
	line[4] = new TLine(e_peak-3.961,0,e_peak-3.961,1e3);
	lbl[0] = new TPaveLabel(0.65,0.80,0.95,0.90,level[0].Data(),"NDC");
	lbl[1] = new TPaveLabel(0.65,0.70,0.95,0.80,level[1].Data(),"NDC");
	lbl[2] = new TPaveLabel(0.65,0.60,0.95,0.70,level[2].Data(),"NDC");
	lbl[3] = new TPaveLabel(0.65,0.50,0.95,0.60,level[3].Data(),"NDC");
	lbl[4] = new TPaveLabel(0.65,0.40,0.95,0.50,level[4].Data(),"NDC");
	for(int iline=0;iline<nline;iline++){
	if(iline==3){
	line[iline]->SetLineColor(30);
	lbl[iline]->SetTextColor(30);
	}else{
	line[iline]->SetLineColor(iline+2);
	lbl[iline]->SetTextColor(iline+2);
	}
	lbl[iline]->SetBorderSize(0);
	lbl[iline]->SetFillStyle(0);
	lbl[iline]->SetFillColor(0);
	line[iline]->SetLineWidth(2);
	line[iline]->Draw();
	}
        TPaveLabel* accepted = new TPaveLabel(0.65,0.75,0.90,0.85,"Accepted by Quartz","NDC");
        accepted->SetBorderSize(0);
        accepted->SetFillStyle(0);
        accepted->SetFillColor(0);
        accepted->SetTextColor(kBlack);
        accepted->Draw();
        TPaveLabel* missed = new TPaveLabel(0.65,0.65,0.90,0.75,"Missed by Quartz","NDC");
        missed->SetBorderSize(0);
        missed->SetFillStyle(0);
        missed->SetFillColor(0);
        missed->SetTextColor(kOrange);
        missed->Draw();
	TPaveLabel* x_accept = new TPaveLabel(0.60,0.75,0.95,0.85,"adc+quartz edge cut","NDC");
	x_accept->SetBorderSize(0);
	x_accept->SetFillStyle(0);
	x_accept->SetFillColor(0);
	x_accept->SetTextColor(kOrange);
//	x_accept->Draw();
	cout<<Form("Elastic: %1.4f",e_peak)<<endl;
	c1->cd(1);
	gPad->SetGrid();
	gPad->SetLogy();
	gPad->SetBottomMargin(0);
	hp->SetLineColor(kBlack);
	hp->Draw("HISTO");
	hp->SetLabelSize(0.05,"y");
	for(int iline=0;iline<nline;iline++){
	line[iline]->SetLineWidth(2);
	line[iline]->Draw();
	lbl[iline]->Draw();
	}
	c1->cd(3);
	gPad->SetGrid();
	gPad->SetLogy();
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0.15);
//	hpc->Scale(1/hpc->Integral("width"));
//	hp->Scale(1/hp->Integral("width"));
//	TH1F* hACC = (TH1F*)hpc->Clone();
//	hACC->Divide(hp);
//	hACC->Draw("HISTO");
	TH1F* hACC = (TH1F*)hpc->Clone();
	TH1F* hACC_xcut = (TH1F*)hxcut->Clone();
	hACC->SetTitle("Fraction (Accepted/Total);R.gold.p (MeV/c)");
	hACC->SetTitleSize(0.08,"x");
	hACC->SetLabelSize(0.05,"x;y");
	hACC->Divide(hp);
	hACC_xcut->Divide(hp);
	hACC->SetLineColor(kBlack);
	hACC->Draw("HISTO");
//	hACC_xcut->Draw("same");
	hACC->SetTitleOffset(0.7,"x");
	for(int iline=0;iline<nline;iline++){
	line[iline]->SetLineWidth(2);
	line[iline]->Draw();
	}
	c1->SaveAs(Form("./plots/acceptance_function_%s_run%d.pdf",target.Data(),run));
	}
}
