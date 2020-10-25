void grand_Asym(){
	TGaxis::SetMaxDigits(3);
	gStyle->SetOptFit(0110);
	gStyle->SetLabelSize(0.05,"x");
	ifstream infile("./TextFiles/output_thCut_allRuns.csv");
	if(!infile){
	cout<<"No input file! Missing file"<<endl;
	exit(0);
	}
	Float_t run, theta, thetaEr, asym, asymEr, sasym, sasymEr, q2, q2Er, sens, sensEr;
	string date, target;
	vector<string>Date, Target;
	vector<Float_t>Run, Theta, ThetaEr, Asym, AsymEr, SAsym, SAsymEr, Q2, Q2Er, Sens, SensEr;
	while(infile>>run>>date>>target>>theta>>thetaEr>>asym>>asymEr>>sasym>>sasymEr>>q2>>q2Er>>sens>>sensEr){
	Run.push_back(run);
	Date.push_back(date);
	Target.push_back(target);
	Theta.push_back(theta);
	ThetaEr.push_back(thetaEr);
	Q2.push_back(q2);
	Q2Er.push_back(q2Er);
	Asym.push_back(asym);
	AsymEr.push_back(asymEr);
	}
	infile.close();

	int npt = Run.size();
	Float_t Entry[npt];
	for(int ipt=0;ipt<npt;ipt++)
	Entry[ipt] = ipt;
	
	TCanvas* c1 = new TCanvas("c1","c1",800,500);
	c1->SetTopMargin(0.05);
	c1->SetBottomMargin(0.15);
	TGraphErrors* gr = new TGraphErrors(npt,Entry,&Q2[0],0,&Q2Er[0]);
	gr->SetMarkerStyle(20);
	gr->Draw("AP");
	gr->Fit("pol0");
	gr->SetLineColor(4);
	gr->SetMarkerColor(4);
	TF1* fit = gr->GetFunction("pol0");
	fit->SetLineColor(4);

	gPad->Update();

	gr->GetXaxis()->Set(npt,-0.5,npt-0.5);
	for(int ipt=0;ipt<npt;ipt++)
	gr->GetXaxis()->SetBinLabel(ipt+1,Form("#splitline{#color[2]{%5.0f}}{#splitline{%s}{%s}}",Run[ipt],Target[ipt].c_str(),Date[ipt].c_str()));
	gPad->Modified();

	TLegend* lag = new TLegend(0.2,0.20,0.3,0.30);
	lag->SetFillStyle(0);
	lag->AddEntry(gr,"LHRS","ep");
	lag->Draw();
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.05);
	latex.SetTextColor(1);
	latex.DrawLatex(0.25,0.85,"PREX-2 Q^{2} for various targets and dates");
	latex.DrawLatex(0.25,0.80,"Error bare are errors in mean");
	latex.SetTextColor(4);
	latex.SetTextSize(0.04);
	latex.DrawLatex(0.9,0.11,"L run");
	latex.SetTextColor(2);
	latex.DrawLatex(0.9,0.08,"R run");
	latex.SetTextColor(1);
	latex.DrawLatex(0.9,0.045,"Target");
	latex.SetTextColor(1);
	latex.DrawLatex(0.9,0.015,"Date");
	gPad->Update();
/*	TPaveStats* statL = (TPaveStats*)gr->FindObject("stats");
	TPaveStats* statR = (TPaveStats*)grR->FindObject("stats");
	statL->SetY2NDC(0.30);
	statL->SetY1NDC(0.25);
	statL->SetX2NDC(0.90);
	statL->SetX1NDC(0.65);
	statL->SetTextColor(4);
	
	statR->SetY2NDC(0.25);
	statR->SetY1NDC(0.20);
	statR->SetX2NDC(0.90);
	statR->SetX1NDC(0.65);
	statR->SetTextColor(2);
	gPad->Modified();
*/	c1->SaveAs(Form("./plots/grand_Asym.pdf"));
}
