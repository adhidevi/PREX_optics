//run as analyzer -l 'acceptancePbNew.C(2052,"208Pb")' for example
void acceptancePbNew(int run,TString target){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetTitleFontSize(0.08);
	gROOT->ForceStyle();
	TGaxis::SetMaxDigits(3);

//LHRS upstream adc cut
	double upadc_cutL_approx = 480;
//RHRS upstream adc cut
	double upadc_cutR = 502;
	TChain* T = new TChain("T");
	if(run<10000){
	T->Add(Form("/chafs1/work1/prex_counting/chandan/prexLHRS_%d_-1*.root",run));
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
	
	TCanvas* cADC = new TCanvas("cADC","cADC",500,400);
	TH1F* huq = new TH1F("huq",Form("ADC raw USL (run%d)",run),300,400,700);
	huq->SetLineColor(1);
	T->Draw("P.upQadcL>>huq");
	TH1F* huqCpy = (TH1F*)huq->Clone("huqCpy");
	huqCpy->GetXaxis()->SetRangeUser(upadc_cutL_approx-10,upadc_cutL_approx+10);
	double upadc_cutL = huqCpy->GetXaxis()->GetBinCenter(huqCpy->GetMinimumBin());
	TLine* adcL = new TLine(upadc_cutL,0.0,upadc_cutL,huq->GetMaximum());
	adcL->SetLineColor(2);
	adcL->SetLineWidth(2);
	adcL->Draw();
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.05);
	lat.SetTextColor(2);
	lat.DrawLatex(0.50,0.70,Form("ADC_cut = %3.1f",upadc_cutL));
	TCut adc_cut_up = Form("P.upQadcL>%f",upadc_cutL);
	TCut adc_cut_up_ped = Form("P.upQadcL<%f",upadc_cutL);
	TCut cut_basic = trig_cut+vdc_cut;
	TCut cut = cut_basic+tr_cut+tg_cut;
	TCut cut_wadc = cut+adc_cut_up;
	TCut cut_wped = cut+adc_cut_up_ped;
	TCut x_cut = "(L.tr.x[0]+1.3*L.tr.th[0]) > -0.0665";
	
	double x_low = 940;
	double x_hi = 952;
	double nbins = 200;
	double bin_size = (x_hi-x_low)/nbins;
	TCanvas* cp = new TCanvas("cp","cp",700,500);
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
	hpp->SetLineColor(kBlue+2);
	hpc->Draw("HISTO");
	hpp->Draw("HISTO SAME");
//	hxcut->Draw("same");
	hpc->SetLabelSize(0.05,"y");
	double e_peak = hp->GetBinCenter(hp->GetMaximumBin());
	const int nline = 5;
	TLine* line[nline];
	TPaveLabel* lbl[nline];
	TString level[] = {"Elastic","1^{st} excited 3^{-}","2^{nd} excited 5^{-}","3^{rd} excited 4^{-}","4^{th} excited 5^{-}"};
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
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.07);
	latex.SetTextColor(kMagenta+2);
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
	missed->SetTextColor(kBlue+2);
	missed->Draw();
	TPaveLabel* x_accept = new TPaveLabel(0.60,0.75,0.95,0.85,"adc+quartz edge cut","NDC");
	x_accept->SetBorderSize(0);
	x_accept->SetFillStyle(0);
	x_accept->SetFillColor(0);
//	x_accept->Draw();
	cout<<Form("Elastic: %1.4f",e_peak)<<endl;
	TLine* quartz = new TLine(947.3,0.0,947.3,hp->GetMaximum());
	quartz->SetLineColor(kMagenta+2);
	quartz->SetLineWidth(2);
	quartz->SetLineStyle(10);
	quartz->Draw();
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
	latex.SetTextAngle(90);
	latex.DrawLatex(0.63,0.35,"Quartz Edge");
	quartz->Draw();
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
	hACC->SetTitle(";L.gold.p (MeV/c);");
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
	quartz->Draw();

	c1->SaveAs(Form("./plots/acceptance_function_%s_run%d.png",target.Data(),run));

	TCanvas* c0 = new TCanvas("c0","c0",1000,700);
	gPad->SetLogy();
	TH1F* hpCpy = (TH1F*)hp->Clone("hpCpy");
	hpCpy->SetTitle(Form("LHRS %s spectrum no adc cut (run%d);L.gold.p (MeV/c);",target.Data(),run));
	hpCpy->Draw("HISTO");

/*	TF1* fitP = new TF1("fitP","gaus",hpCpy->GetBinCenter(hpCpy->GetMaximumBin())-0.3,hpCpy->GetBinCenter(hpCpy->GetMaximumBin())+0.3);
	fitP->SetLineColor(kRed);
	fitP->SetLineWidth(2);
	hpCpy->Fit(fitP,"R0");
//	fitP->Draw("same");
	double gP0 = fitP->GetParameter(0);
	double gP1 = fitP->GetParameter(1);
	double gP2 = fitP->GetParameter(2);
	TF1* fitPcrystal = new TF1("fitPcrystal","crystalball",gP1-9*gP2,gP1+3*gP2);
	fitPcrystal->SetParameters(gP0,gP1,gP2,1.64,1.1615);
	fitPcrystal->SetLineColor(kRed);
	fitPcrystal->SetLineWidth(2);
	hpCpy->Fit(fitPcrystal,"R0");
//	fitPcrystal->Draw("same");
	TF1* polfit = new TF1("polfit","pol2",940,942);
	polfit->SetLineColor(3);
	polfit->SetLineWidth(2);
	hpCpy->Fit(polfit,"R0");
//	polfit->Draw("same");

	double crystalpolPar[8];
	TF1* crystalpol = new TF1("crystalpol","crystalball(0)+pol2(5)",hpCpy->GetMinimum(),fitPcrystal->GetXmax());
	fitPcrystal->GetParameters(&crystalpolPar[0]);
	polfit->GetParameters(&crystalpolPar[5]);
	crystalpol->SetParameters(crystalpolPar);
	crystalpol->SetLineColor(6);
	crystalpol->SetLineWidth(2);
	hpCpy->Fit(crystalpol,"R0");
	crystalpol->Draw("same");
*/		
//	hpCpy->Scale(1./hpCpy->Integral(),"width");
	double elasticCount = hpCpy->Interpolate(e_peak);
	double firstCount = hpCpy->Interpolate(e_peak-firstE);
	double secondCount = hpCpy->Interpolate(e_peak-secondE);
	double thirdCount = hpCpy->Interpolate(e_peak-thirdE);
	double fourthCount = hpCpy->Interpolate(e_peak-fourthE);
	cout<<"3 sigma on 1st excited state: "<<3*sqrt(firstCount)<<endl;
	cout<<"3 sigma on 2nd excited state: "<<3*sqrt(secondCount)<<endl;
	cout<<"3 sigma on 3rd excited state: "<<3*sqrt(thirdCount)<<endl;
	cout<<"3 sigma on 4th excited state: "<<3*sqrt(fourthCount)<<endl;

	cout<<"3 sigma bound on 1st excited state: "<<3*sqrt(firstCount)/elasticCount<<endl;
	cout<<"3 sigma bound on 2nd excited state: "<<3*sqrt(secondCount)/elasticCount<<endl;
	cout<<"3 sigma bound on 3rd excited state: "<<3*sqrt(thirdCount)/elasticCount<<endl;
	cout<<"3 sigma bound on 4th excited state: "<<3*sqrt(fourthCount)/elasticCount<<endl;

	double firstEfr = hpCpy->Interpolate(e_peak-firstE)/hpCpy->Interpolate(e_peak);
	double secondEfr = hpCpy->Interpolate(e_peak-secondE)/hpCpy->Interpolate(e_peak);
	double thirdEfr = hpCpy->Interpolate(e_peak-thirdE)/hpCpy->Interpolate(e_peak);
	double fourthEfr = hpCpy->Interpolate(e_peak-fourthE)/hpCpy->Interpolate(e_peak);

	latex.SetTextSize(0.04);
	latex.SetTextAngle(0);
	latex.DrawLatex(0.15,0.85,Form("Relative strength"));
	latex.SetTextColor(3);
	latex.DrawLatex(0.15,0.8,Form("1^{st} excited = %2.4f",firstEfr));
	latex.SetTextColor(4);
	latex.DrawLatex(0.15,0.75,Form("2^{nd} excited = %2.4f",secondEfr));
	latex.SetTextColor(30);
	latex.DrawLatex(0.15,0.7,Form("3^{rd} excited = %2.4f",thirdEfr));
	latex.SetTextColor(6);
	latex.DrawLatex(0.15,0.65,Form("4^{th} excited = %2.4f",fourthEfr));
	for(int iline=0;iline<nline;iline++){
	line[iline]->Draw();
	}
	quartz->Draw();
	latex.SetTextColor(kMagenta+2);
	latex.SetTextAngle(90);
	latex.DrawLatex(0.63,0.35,"Quartz Edge");

	c0->SaveAs(Form("./plots/relative_strength_%s_run%d.png",target.Data(),run));

	TCanvas* cX = new TCanvas("cX","cX",1000,700);
	gPad->SetLogy();
	TH1F* hx1 = new TH1F("hx1",Form("Projected x on quartz (run%d);L.tr.x[0]+1.3*L.tr.th[0];",run),200,-0.10,0.07);
	TH1F* hx2 = new TH1F("hx2",Form("Projected x on quartz (run%d);L.tr.x[0]+1.3*L.tr.th[0];",run),200,-0.10,0.07);
	TH1F* hx3 = new TH1F("hx3",Form("Projected x on quartz (run%d);L.tr.x[0]+1.3*L.tr.th[0];",run),200,-0.10,0.07);
	TH1F* hx4 = new TH1F("hx4",Form("Projected x on quartz (run%d);L.tr.x[0]+1.3*L.tr.th[0];",run),200,-0.10,0.07);
	hx1->SetLineColor(kBlack);
	hx2->SetLineColor(kRed);
	hx3->SetLineColor(kBlue+2);
	hx4->SetLineColor(kOrange);
	double slope = 16.18;
	T->Draw("L.tr.x[0]+1.3*L.tr.th[0]+0.0364>>hx1",cut);
	T->Draw("L.tr.x[0]+1.3*L.tr.th[0]+0.0364>>hx2",cut_wadc,"same");
	T->Draw("L.tr.x[0]+1.3*L.tr.th[0]+0.0364>>hx3",cut_wped,"same");
	T->Draw("L.tr.x[0]+1.3*L.tr.th[0]+0.0364>>hx4",cut_wadc+"L.gold.p<0.9473","same");
	gStyle->SetOptFit(0);
	cout<<Form("peak location in x: %2.4f",hx1->GetBinCenter(hx1->GetMaximumBin()))<<endl;;
	TF1* fit = new TF1("fit","gaus",hx1->GetBinCenter(hx1->GetMaximumBin())-0.01,hx1->GetBinCenter(hx1->GetMaximumBin())+0.01);
	hx1->Fit(fit,"R0");
	double xE0 = -slope*0.0/9.53374/100;
	double xE1 = -slope*firstE/9.53374/100;
	double xE2 = -slope*secondE/9.53374/100;
	double xE3 = -slope*thirdE/9.53374/100;
	double xE4 = -slope*fourthE/9.53374/100;
        line[0] = new TLine(xE0,0,xE0,hx1->GetMaximum());
        line[1] = new TLine(xE1,0,xE1,hx1->GetMaximum()/2.);
        line[2] = new TLine(xE2,0,xE2,hx1->GetMaximum()/2.);
        line[3] = new TLine(xE3,0,xE3,hx1->GetMaximum()/2.);
        line[4] = new TLine(xE4,0,xE4,hx1->GetMaximum()/2.);

	latex.SetTextAngle(0);
        latex.SetTextColor(3);
        latex.DrawLatex(0.2,0.85,Form("1^{st} excited (3^{-})"));
        latex.SetTextColor(4);
        latex.DrawLatex(0.2,0.8,Form("2^{nd} excited (5^{-})"));
        latex.SetTextColor(30);
        latex.DrawLatex(0.2,0.75,Form("3rd excited (4^{-})"));
        latex.SetTextColor(6);
        latex.DrawLatex(0.2,0.7,Form("4th excited (5^{-})"));

        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.6,0.85,Form("All Events"));
        latex.SetTextColor(kRed);
        latex.DrawLatex(0.6,0.8,Form("Accepted By Quartz"));
        latex.SetTextColor(kBlue+2);
        latex.DrawLatex(0.6,0.75,Form("Missed By Quartz"));
        latex.SetTextColor(kOrange);
        latex.DrawLatex(0.6,0.7,Form("Accepted but p < 947.3 MeV"));

	for(int iline=0;iline<nline;iline++){
	if(iline==3){
	line[iline]->SetLineColor(30);
	}else{
	line[iline]->SetLineColor(iline+2);
	}
	line[iline]->SetLineWidth(2);
	line[iline]->Draw();
	}
	quartz = new TLine(-0.035,0.0,-0.035,hx1->GetMaximum());
	quartz->SetLineColor(kMagenta+2);
	quartz->SetLineWidth(2);
	quartz->SetLineStyle(10);
	quartz->Draw();
	latex.SetTextColor(kMagenta+2);
	latex.SetTextAngle(90);
	latex.DrawLatex(0.45,0.40,"Quartz Edge");

	cX->SaveAs(Form("./plots/projected_X_%s_run%d.png",target.Data(),run));

	TCanvas* cXfp = new TCanvas("cXfp","cXfp",1000,700);
	gPad->SetLogy();
	TH1F* hxfp1 = new TH1F("hxfp1",Form("Transport x on VDC (run%d);L.tr.x[0];",run),200,-0.10,0.07);
	TH1F* hxfp2 = new TH1F("hxfp2",Form("Transport x on VDC (run%d);L.tr.x[0];",run),200,-0.10,0.07);
	TH1F* hxfp3 = new TH1F("hxfp3",Form("Transport x on VDC (run%d);L.tr.x[0];",run),200,-0.10,0.07);
	TH1F* hxfp4 = new TH1F("hxfp4",Form("Transport x on VDC (run%d);L.tr.x[0];",run),200,-0.10,0.07);
	hxfp1->SetLineColor(kBlack);
	hxfp2->SetLineColor(kRed);
	hxfp3->SetLineColor(kBlue+2);
	hxfp4->SetLineColor(kOrange);
	double slopefp = 13.98;
	T->Draw("L.tr.x[0]+0.02799>>hxfp1",cut);
	T->Draw("L.tr.x[0]+0.02799>>hxfp2",cut_wadc,"same");
	T->Draw("L.tr.x[0]+0.02799>>hxfp3",cut_wped,"same");
	T->Draw("L.tr.x[0]+0.02799>>hxfp4",cut_wadc+"L.gold.p<0.9473","same");
	gStyle->SetOptFit(0);
	cout<<Form("peak location in x: %2.4f",hxfp1->GetBinCenter(hxfp1->GetMaximumBin()))<<endl;
	TF1* fitfp = new TF1("fitfp","gaus",hxfp1->GetBinCenter(hxfp1->GetMaximumBin())-0.01,hxfp1->GetBinCenter(hxfp1->GetMaximumBin())+0.01);
	hxfp1->Fit(fitfp,"R0");
	double xfpE0 = -slopefp*0.0/9.53374/100;
	double xfpE1 = -slopefp*firstE/9.53374/100;
	double xfpE2 = -slopefp*secondE/9.53374/100;
	double xfpE3 = -slopefp*thirdE/9.53374/100;
	double xfpE4 = -slopefp*fourthE/9.53374/100;
        line[0] = new TLine(xfpE0,0,xfpE0,hxfp1->GetMaximum());
        line[1] = new TLine(xfpE1,0,xfpE1,hxfp1->GetMaximum()/2.);
        line[2] = new TLine(xfpE2,0,xfpE2,hxfp1->GetMaximum()/2.);
        line[3] = new TLine(xfpE3,0,xfpE3,hxfp1->GetMaximum()/2.);
        line[4] = new TLine(xfpE4,0,xfpE4,hxfp1->GetMaximum()/2.);

	latex.SetTextAngle(0);
        latex.SetTextColor(3);
        latex.DrawLatex(0.2,0.85,Form("1^{st} excited (3^{-})"));
        latex.SetTextColor(4);
        latex.DrawLatex(0.2,0.8,Form("2^{nd} excited (5^{-})"));
        latex.SetTextColor(30);
        latex.DrawLatex(0.2,0.75,Form("3rd excited (4^{-})"));
        latex.SetTextColor(6);
        latex.DrawLatex(0.2,0.7,Form("4th excited (5^{-})"));

        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.6,0.85,Form("All Events"));
        latex.SetTextColor(kRed);
        latex.DrawLatex(0.6,0.8,Form("Accepted By Quartz"));
        latex.SetTextColor(kBlue+2);
        latex.DrawLatex(0.6,0.75,Form("Missed By Quartz"));
        latex.SetTextColor(kOrange);
        latex.DrawLatex(0.6,0.7,Form("Accepted but p < 947.3 MeV"));

	for(int iline=0;iline<nline;iline++){
	if(iline==3){
	line[iline]->SetLineColor(30);
	}else{
	line[iline]->SetLineColor(iline+2);
	}
	line[iline]->SetLineWidth(2);
	line[iline]->Draw();
	}
	quartz = new TLine(-0.029,0.0,-0.029,hxfp1->GetMaximum());
	quartz->SetLineColor(kMagenta+2);
	quartz->SetLineWidth(2);
	quartz->SetLineStyle(10);
	quartz->Draw();
	latex.SetTextColor(kMagenta+2);
	latex.SetTextAngle(90);
	latex.DrawLatex(0.45,0.40,"Quartz Edge");

	cXfp->SaveAs(Form("./plots/VDC_X_%s_run%d.png",target.Data(),run));

	}else{
	T->Add(Form("/chafs1/work1/prex_counting/chandan/prexRHRS_%d_-1*.root",run));
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
	TCut x_cut = "(R.tr.x[0]+1.3*R.tr.th[0]) > -0.0725";

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
	hpp->SetLineColor(kBlue+2);
	hpc->Draw("HISTO");
	hpp->Draw("HISTO SAME");
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
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.04);
	latex.SetTextColor(kMagenta+2);
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
        missed->SetTextColor(kBlue+2);
        missed->Draw();
	TPaveLabel* x_accept = new TPaveLabel(0.60,0.75,0.95,0.85,"adc+quartz edge cut","NDC");
	x_accept->SetBorderSize(0);
	x_accept->SetFillStyle(0);
	x_accept->SetFillColor(0);
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
	c1->SaveAs(Form("./plots/acceptance_function_%s_run%d.png",target.Data(),run));
	}
}
