void grand_Qsq(){
	TGaxis::SetMaxDigits(3);
	gStyle->SetOptFit(0110);
	gStyle->SetLabelSize(0.05,"x");
	ifstream infileL("./TextFiles/qsq_calculatedL.csv");
	ifstream infileR("./TextFiles/qsq_calculatedR.csv");
	if(!infileL||!infileR){
	cout<<"No input file! Missing file"<<endl;
	exit(0);
	}
	string dateL, dateR, targetL, targetR;
	Float_t runL, runR, qsqL, qsqR, errL, errR;
	vector<string>DateL, DateR, TargetL, TargetR;
	vector<Float_t> RunL, RunR, QsqL, QsqR, ErrR, ErrL;
	while(infileL>>dateL>>targetL>>runL>>qsqL>>errL){
	DateL.push_back(dateL);
	TargetL.push_back(targetL);
	RunL.push_back(runL);
	QsqL.push_back(qsqL);
	ErrL.push_back(errL);
	}
	infileL.close();

	while(infileR>>dateR>>targetR>>runR>>qsqR>>errR){
	DateR.push_back(dateR);
	TargetR.push_back(targetR);
	RunR.push_back(runR);
	QsqR.push_back(qsqR);
	ErrR.push_back(errR);
	}
	infileR.close();

	int nptL = RunL.size();
	int nptR = RunR.size();

	Float_t EntryL[nptL];
	Float_t EntryR[nptR];
	for(int ipt=0;ipt<nptL;ipt++)
	EntryL[ipt] = ipt;
	for(int ipt=0;ipt<nptR;ipt++)
	EntryR[ipt] = ipt;
	
	TCanvas* c1 = new TCanvas("c1","c1",1000,600);
	c1->SetTopMargin(0.05);
	c1->SetBottomMargin(0.15);
	gPad->SetGridx();
	TGraphErrors* grL = new TGraphErrors(nptL,EntryL,&QsqL[0],0,&ErrL[0]);
	grL->SetMarkerStyle(20);
	grL->Draw("AP");
	grL->Fit("pol0");
	grL->SetLineColor(4);
	grL->SetMarkerColor(4);
	TF1* fitL = grL->GetFunction("pol0");
	fitL->SetLineColor(4);

	TGraphErrors* grR = new TGraphErrors(nptR,EntryR,&QsqR[0],0,&ErrR[0]);
	grR->SetMarkerStyle(29);
	grR->SetMarkerSize(1.5);
	grR->Draw("AP");
	grR->Fit("pol0");
	grR->SetLineColor(2);
	grR->SetMarkerColor(2);
	TF1* fitR = grR->GetFunction("pol0");
	fitR->SetLineColor(2);
	
	TMultiGraph* mg = new TMultiGraph("mg","mg");
	mg->Add(grL,"P");
	mg->Add(grR,"P");
	mg->Draw("AP");
	mg->SetTitle(";;Q^{2} (GeV/c)^{2}");
	mg->GetYaxis()->CenterTitle();
	mg->GetYaxis()->SetRangeUser(0.00610,0.00645);
	gPad->Update();

	mg->GetXaxis()->Set(nptR,-0.5,nptR-0.5);
	for(int ipt=0;ipt<nptR;ipt++)
	mg->GetXaxis()->SetBinLabel(ipt+1,Form("#splitline{#color[4]{%5.0f}}{#splitline{#color[2]{%5.0f}}{#splitline{%s}{%s}}}",RunL[ipt],RunR[ipt],TargetL[ipt].c_str(),DateL[ipt].c_str()));
//	mg->GetXaxis()->SetLabelSize(0.031);
	mg->GetXaxis()->SetLabelSize(0.036);
	gPad->Modified();
	TLegend* lag = new TLegend(0.12,0.60,0.20,0.70);
	lag->SetFillStyle(0);
	lag->AddEntry(grL,"LHRS","ep");
	lag->AddEntry(grR,"RHRS","ep");
	lag->Draw();
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.05);
	latex.SetTextColor(1);
	latex.DrawLatex(0.25,0.70,"PREX-2 Q^{2} for various targets and dates");
	latex.DrawLatex(0.25,0.65,Form("Error bare are 0.3%s of the Q^{2}","%"));
	latex.SetTextColor(4);
	latex.SetTextSize(0.025);
	latex.DrawLatex(0.9,0.125,"L run");
	latex.SetTextColor(2);
	latex.DrawLatex(0.9,0.10,"R run");
	latex.SetTextColor(1);
	latex.DrawLatex(0.9,0.075,"Target");
	latex.SetTextColor(1);
	latex.DrawLatex(0.9,0.05,"Date");
	gPad->Update();
	TPaveStats* statL = (TPaveStats*)grL->FindObject("stats");
	TPaveStats* statR = (TPaveStats*)grR->FindObject("stats");
	statL->SetY2NDC(0.60);
	statL->SetY1NDC(0.55);
	statL->SetX2NDC(0.90);
	statL->SetX1NDC(0.65);
	statL->SetTextColor(4);
	
	statR->SetY2NDC(0.65);
	statR->SetY1NDC(0.60);
	statR->SetX2NDC(0.90);
	statR->SetX1NDC(0.65);
	statR->SetTextColor(2);
	gPad->Modified();
	c1->SaveAs(Form("./plots/grand_Qsq.pdf"));

	TCanvas* c2 = new TCanvas("c2","c2",600,600);
	TH2F* qsqLR = new TH2F("qsqLR","Qsq LHRS vs Qsq RHRS;RHRS;LHRS",100,0.0063,0.0064,100,0.0061,0.0063);
	qsqLR->GetYaxis()->SetTitleOffset(1.3);
	qsqLR->GetXaxis()->SetLabelSize(0.03);
	qsqLR->SetMarkerStyle(20);
	Float_t* LQ2 = QsqL.data();
	Float_t* RQ2 = QsqR.data();
	for(int iptL=0;iptL<nptL;iptL++){
	qsqLR->Fill(RQ2[iptL],LQ2[iptL]);
	}
	qsqLR->Draw();
	c2->SaveAs("./plots/grand_Qsq_LvsR.pdf");
	gSystem->Exec(Form("pdfunite ./plots/grand_Qsq.pdf ./plots/grand_Qsq_LvsR.pdf ./plots/good_runs_Qsq.pdf"));
	gSystem->Exec(Form("rm -rf ./plots/grand_Qsq.pdf ./plots/grand_Qsq_LvsR.pdf"));
}
