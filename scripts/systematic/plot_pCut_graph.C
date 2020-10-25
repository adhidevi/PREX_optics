void plot_pCut_graph(int run){
	gStyle->SetTitleYOffset(1.3);
	ifstream infile(Form("./TextFiles/output_pCut_run%d.csv",run));
	if(!infile){
	cout<<"Input file doesn't exits! Quiting...."<<endl;
	exit(0);
	}
	
	double qedge, pcut, thM, thR, asM, asR, sasM, sasR, q2M, q2R, senM, senR;
	vector<double> Qedge, Pcut, ThM, ThR, AsM, AsR, SasM, SasR, Q2M, Q2R, SenM, SenR;
	while(infile>>qedge>>pcut>>thM>>thR>>asM>>asR>>sasM>>sasR>>q2M>>q2R>>senM>>senR){
	Qedge.push_back(qedge);
	Pcut.push_back(pcut);
	ThM.push_back(thM);
	ThR.push_back(thR);
	AsM.push_back(asM);
	AsR.push_back(asR);
	SasM.push_back(sasM);
	SasR.push_back(sasR);
	Q2M.push_back(q2M);
	Q2R.push_back(q2R);
	SenM.push_back(senM);
	SenR.push_back(senR);
	}
	infile.close();
	int npt = Pcut.size();
	
	for(int ipt=0;ipt<npt;ipt++)
	cout<<Pcut[ipt]<<endl;

	TCanvas* c1 = new TCanvas("c1","c1",700,500);	
	TGraph* grL = new TGraph(npt,&Pcut[0],&ThM[0]);
	grL->SetMarkerStyle(20);
	grL->SetTitle(Form("#theta_{lab} vs pCut (run%d);pCut (GeV/c);#theta_{lab} (deg)",run));
	grL->Draw("AP");
	double max = grL->GetHistogram()->GetMaximum();
	double min = grL->GetHistogram()->GetMinimum();
	TLine* qlineL = new TLine(Qedge[0],min,Qedge[0],max);
	qlineL->SetLineColor(kOrange);
	qlineL->SetLineWidth(2);
	qlineL->Draw();

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.05);
	latex.SetTextColor(kOrange);
	latex.DrawLatex(0.6,0.8,"q_edge");
}
