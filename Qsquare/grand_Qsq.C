void grand_Qsq(){
	TGaxis::SetMaxDigits(3);
	gStyle->SetOptFit(0110);
	gStyle->SetLabelSize(0.05,"x");
	ifstream infileL("../alignment/TextFiles/qsquare_diagnostic_quantitiesL.csv");
	ifstream infileR("../alignment/TextFiles/qsquare_diagnostic_quantitiesR.csv");
	if(!infileL||!infileR){
	cout<<"No input file! Missing file"<<endl;
	exit(0);
	}
	string dateL, dateR, tgtL, tgtR;
	double runL, adc_cL, x_peakL, x_cutL, p_peakL, p_cutL, qsq_adcL, rms_adcL, ent_adcL, qsq_xL, rms_xL, ent_xL, qsq_pL, rms_pL, ent_pL, qsq_xyloL, rms_xyloL, ent_xyloL, qsq_xyhiL, rms_xyhiL, ent_xyhiL;
	double runR, adc_cR, x_peakR, x_cutR, p_peakR, p_cutR, qsq_adcR, rms_adcR, ent_adcR, qsq_xR, rms_xR, ent_xR, qsq_pR, rms_pR, ent_pR, qsq_xyloR, rms_xyloR, ent_xyloR, qsq_xyhiR, rms_xyhiR, ent_xyhiR;

	vector<string>DateL, DateR, TargetL, TargetR;
	vector<Float_t> RunL, RunR, QsqL, QsqR, ErrR, ErrL;
	infileL.ignore(10000,'\n');
	while(infileL>>runL>>tgtL>>dateL>>adc_cL>>x_peakL>>x_cutL>>p_peakL>>p_cutL>>qsq_adcL>>rms_adcL>>ent_adcL>>qsq_xL>>rms_xL>>ent_xL>>qsq_pL>>rms_pL>>ent_pL>>qsq_xyloL>>rms_xyloL>>ent_xyloL>>qsq_xyhiL>>rms_xyhiL>>ent_xyhiL){
	DateL.push_back(dateL);
	TargetL.push_back(tgtL);
	RunL.push_back(runL);
	QsqL.push_back(qsq_adcL);
	ErrL.push_back(rms_adcL/TMath::Sqrt(ent_adcL));
	}
	infileL.close();

	infileR.ignore(10000,'\n');
	while(infileR>>runR>>tgtR>>dateR>>adc_cR>>x_peakR>>x_cutR>>p_peakR>>p_cutR>>qsq_adcR>>rms_adcR>>ent_adcR>>qsq_xR>>rms_xR>>ent_xR>>qsq_pR>>rms_pR>>ent_pR>>qsq_xyloR>>rms_xyloR>>ent_xyloR>>qsq_xyhiR>>rms_xyhiR>>ent_xyhiR){
	DateR.push_back(dateR);
	TargetR.push_back(tgtR);
	RunR.push_back(runR);
	QsqR.push_back(qsq_adcR);
	ErrR.push_back(rms_adcR/TMath::Sqrt(ent_adcR));
	}
	infileR.close();

	int nptL = RunL.size();
	int nptR = RunR.size();
	cout<<"Counted: "<<nptL<<" runs in LHRS and "<<nptR<<" runs in RHRS"<<endl;

	Float_t EntryL[nptL];
	Float_t EntryR[nptR];
	for(int ipt=0;ipt<nptL;ipt++)
	EntryL[ipt] = ipt;
	for(int ipt=0;ipt<nptR;ipt++)
	EntryR[ipt] = ipt;
	
	TCanvas* c1 = new TCanvas("c1","c1",1000,600);
//	c1->SetTopMargin(0.05);
//	c1->SetBottomMargin(0.15);
	gPad->SetGridx();
	TGraphErrors* grL = new TGraphErrors(nptL,EntryL,&QsqL[0],0,&ErrL[0]);
	grL->SetMarkerStyle(20);
	grL->Draw("AP");
	grL->Fit("pol0","W");
	grL->SetLineColor(4);
	grL->SetMarkerColor(4);
	TF1* fitL = grL->GetFunction("pol0");
	fitL->SetLineColor(4);

	TGraphErrors* grR = new TGraphErrors(nptR,EntryR,&QsqR[0],0,&ErrR[0]);
	grR->SetMarkerStyle(29);
	grR->SetMarkerSize(1.5);
	grR->Draw("AP");
	grR->Fit("pol0","W");
	grR->SetLineColor(2);
	grR->SetMarkerColor(2);
	TF1* fitR = grR->GetFunction("pol0");
	fitR->SetLineColor(2);
	
	TMultiGraph* mg = new TMultiGraph("mg","PREX-2 Q^{2} over time;;Q^{2} (GeV/c)^{2}");
	mg->Add(grL,"P");
	mg->Add(grR,"P");
	mg->Draw("AP");
	mg->GetYaxis()->CenterTitle();
	mg->GetYaxis()->SetRangeUser(0.00635,0.0066);
	gPad->Update();

	mg->GetXaxis()->Set(nptL,-0.5,nptL-0.5);
	TString x_label[] = {"#splitline{#color[4]{1983}}{#color[2]{21108}}","#splitline{#color[4]{1996}}{#color[2]{21121}}","#splitline{#color[4]{2052}}{#color[2]{21185}}","#splitline{#color[4]{2199}}{#color[2]{21344}}","#color[4]{2219}","#color[4]{2292}","#color[4]{2293}","#color[4]{2294}"};
	for(int ipt=0;ipt<nptL;ipt++){
	mg->GetXaxis()->SetBinLabel(ipt+1,Form("%s",x_label[ipt].Data()));
	}
	mg->GetXaxis()->SetLabelSize(0.055);
	gPad->Modified();
	TLegend* lag = new TLegend(0.15,0.50,0.25,0.60);
	lag->SetFillStyle(0);
	lag->AddEntry(grL,"LHRS","ep");
	lag->AddEntry(grR,"RHRS","ep");
	lag->Draw();
	gPad->Update();
	TPaveStats* statL = (TPaveStats*)grL->FindObject("stats");
	TPaveStats* statR = (TPaveStats*)grR->FindObject("stats");
	statL->SetY2NDC(0.60);
	statL->SetY1NDC(0.55);
	statL->SetX2NDC(0.90);
	statL->SetX1NDC(0.65);
	statL->SetTextColor(4);
	
	statR->SetY2NDC(0.55);
	statR->SetY1NDC(0.50);
	statR->SetX2NDC(0.90);
	statR->SetX1NDC(0.65);
	statR->SetTextColor(2);
	gPad->Modified();
	c1->SaveAs(Form("./plots/grand_Qsq.pdf"));

}
