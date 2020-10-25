void plot_grand_yloCut(){
//	TGaxis::SetMaxDigits(3);
	gStyle->SetTitleYOffset(1.3);
	ifstream runfile("./TextFiles/runlistR.list");
	int run;
	string date, target;
	vector<int>Run;
	vector<string>Date,Target;
	while(runfile>>date>>target>>run){
	     Date.push_back(date);
	     Target.push_back(target);
	     Run.push_back(run);
	}
	runfile.close();
	int nrun = Run.size();
	int* myRun = Run.data();	
	double qedge, ylocut, theta, asym, sasym, qsq, sens;
	vector<double>Qedge, Ylocut, Theta, Asym, Sasym, Qsq, Sens;
	TGraphErrors* gr1[nrun];
	TGraphErrors* gr2[nrun];
	TGraphErrors* gr3[nrun];
	TGraphErrors* gr4[nrun];
	TGraphErrors* gr5[nrun];
	TString filename;
	
	for(int irun=0;irun<nrun;irun++){
	   filename = Form("./TextFiles/output_yloCut_run%d.csv",myRun[irun]);
 	   ifstream infile(Form("%s",filename.Data()));
	   if(infile.is_open())
	   cout<<Form("reading file: %s",filename.Data())<<endl;
	   if(!infile){
	        cout<<Form("File '''%s'''does not exist! Quiting...",filename.Data())<<endl;
		exit(0);
	   }
	   while(infile>>qedge>>ylocut>>theta>>asym>>sasym>>qsq>>sens){
		Qedge.push_back(qedge);
		Ylocut.push_back(ylocut);
		Theta.push_back(theta);
		Asym.push_back(asym);
		Sasym.push_back(sasym);
		Qsq.push_back(qsq);
		Sens.push_back(sens);
	}
	infile.close();
	int sze = Qedge.size();
	for(int j=0;j<sze;j++)
	cout<<Ylocut[j]<<"\t"<<Theta[j]<<"\t"<<Asym[j]<<"\t"<<Qsq[j]<<endl;

	gr1[irun] = new TGraphErrors(sze,&Ylocut[0],&Theta[0]);
	gr2[irun] = new TGraphErrors(sze,&Ylocut[0],&Asym[0]);
	gr3[irun] = new TGraphErrors(sze,&Ylocut[0],&Qsq[0]);
	gr4[irun] = new TGraphErrors(sze,&Ylocut[0],&Sasym[0]);
	gr5[irun] = new TGraphErrors(sze,&Ylocut[0],&Sens[0]);
	gr1[irun]->SetMarkerStyle(20);
	gr2[irun]->SetMarkerStyle(20);
	gr3[irun]->SetMarkerStyle(20);
	gr4[irun]->SetMarkerStyle(20);
	gr5[irun]->SetMarkerStyle(20);
	if(irun==4){
	gr1[irun]->SetMarkerColor(49);
	gr2[irun]->SetMarkerColor(49);
	gr3[irun]->SetMarkerColor(49);
	gr4[irun]->SetMarkerColor(49);
	gr5[irun]->SetMarkerColor(49);
	}else if(irun==9){
	gr1[irun]->SetMarkerColor(41);
	gr2[irun]->SetMarkerColor(41);
	gr3[irun]->SetMarkerColor(41);
	gr4[irun]->SetMarkerColor(41);
	gr5[irun]->SetMarkerColor(41);
	}else{
	gr1[irun]->SetMarkerColor(irun+1);
	gr2[irun]->SetMarkerColor(irun+1);
	gr3[irun]->SetMarkerColor(irun+1);
	gr4[irun]->SetMarkerColor(irun+1);
	gr5[irun]->SetMarkerColor(irun+1);
	}
	gr1[irun]->Draw("APL");
	gr2[irun]->Draw("APL");
	gr3[irun]->Draw("APL");
	gr4[irun]->Draw("APL");
	gr5[irun]->Draw("APL");
	Qedge.clear();
	Ylocut.clear();
	Theta.clear();
	Asym.clear();
	Sasym.clear();
	Qsq.clear();
	Sens.clear();
	}

	TMultiGraph* mg1 = new TMultiGraph("mg1","mg1");
	TMultiGraph* mg2 = new TMultiGraph("mg2","mg2");
	TMultiGraph* mg3 = new TMultiGraph("mg3","mg3");
	TMultiGraph* mg4 = new TMultiGraph("mg4","mg4");
	TMultiGraph* mg5 = new TMultiGraph("mg5","mg5");
	mg1->SetTitle("Variation of scattering angle (#theta) with yloCut;yloCut (m);Theta (degrees)");
	mg2->SetTitle("Variation of Asymmetry with yloCut;yloCut (m);Asymmetry (ppm)");
	mg3->SetTitle("Variation of Q-square with yloCut;yloCut (m);Qsq ((GeV/c)^{2})");
	mg4->SetTitle("Variation of Stretched Asymmetry with yloCut;yloCut (m);Stretched Asymmetry (ppm)");
	mg5->SetTitle("Variation of Sensitivity with yloCut;yloCut (m);Sensitivity");
	for(int irun=0;irun<nrun;irun++){
	   mg1->Add(gr1[irun],"lp");
	   mg2->Add(gr2[irun],"lp");
	   mg3->Add(gr3[irun],"lp");
	   mg4->Add(gr4[irun],"lp");
	   mg5->Add(gr5[irun],"lp");
	}
	TCanvas* c1 = new TCanvas("c1","c1",900,500);
	c1->Divide(2,1,0.00001,0.00001);
	c1->cd(1);
	gPad->SetPad(0.0,0.0,0.75,1.0);
	gPad->SetRightMargin(0.02);
	gPad->Update();
	mg1->GetXaxis()->CenterTitle();
	mg1->GetYaxis()->CenterTitle();
	TLegend* leg1 = new TLegend(0.00,0.50,1.1,0.90);
	leg1->SetTextSize(0.07);
	mg1->Draw("APL");
        gPad->Update();
        TLine* q_edge = new TLine(-0.013,-10000,-0.013,10000);
        q_edge->SetLineColor(kOrange);
        q_edge->SetLineWidth(2);
        q_edge->Draw();
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.05);
        latex.SetTextColor(kOrange);
        latex.DrawLatex(0.70,0.30,"q_edge");
	c1->cd(2);
	gPad->SetPad(0.76,0.0,1.0,1.0);
	gPad->SetLeftMargin(0.0);
	for(int irun=0;irun<nrun;irun++){
	   leg1->AddEntry(gr1[irun],Form("run%d, %s, %s",myRun[irun],Date[irun].c_str(),Target[irun].c_str()),"p");
	   leg1->SetFillStyle(0);
	   leg1->SetBorderSize(0);
	   leg1->Draw();
	}
	c1->SaveAs("./temp1/gr1.pdf");
	TCanvas* c2 = new TCanvas("c2","c2",900,500);
	c2->Divide(2,1,0.00001,0.00001);
	c2->cd(1);
	gPad->SetPad(0.0,0.0,0.75,1.0);
	gPad->SetRightMargin(0.05);
	gPad->Update();
	mg2->GetXaxis()->CenterTitle();
	mg2->GetYaxis()->CenterTitle();
	mg2->Draw("APL");
        q_edge->Draw();
        latex.DrawLatex(0.70,0.30,"q_edge");
	c2->cd(2);
	gPad->SetPad(0.76,0.0,1.0,1.0);
	gPad->SetLeftMargin(0.0);
	TLegend* leg2 = new TLegend(0.00,0.50,1.1,0.90);
	leg2->SetTextSize(0.07);
	for(int irun=0;irun<nrun;irun++){
	   leg2->AddEntry(gr2[irun],Form("run%d, %s, %s",myRun[irun],Date[irun].c_str(),Target[irun].c_str()),"p");
	   leg2->SetFillStyle(0);
	   leg2->SetBorderSize(0);
	   leg2->Draw();
	}
	c2->SaveAs("./temp1/gr2.pdf");
	TCanvas* c3 = new TCanvas("c3","c3",900,500);
	c3->Divide(2,1,0.00001,0.00001);
	c3->cd(1);
	gPad->SetPad(0.0,0.0,0.75,1.0);
	gPad->SetRightMargin(0.02);
	gPad->Update();
	mg3->GetXaxis()->CenterTitle();
	mg3->GetYaxis()->CenterTitle();
	mg3->Draw("APL");
        q_edge->Draw();
        latex.DrawLatex(0.70,0.30,"q_edge");
	c3->cd(2);
	gPad->SetPad(0.76,0.0,1.0,1.0);
	gPad->SetLeftMargin(0.0);
	TLegend* leg3 = new TLegend(0.00,0.50,1.1,0.90);
	leg3->SetTextSize(0.07);
	for(int irun=0;irun<nrun;irun++){
	   leg3->AddEntry(gr3[irun],Form("run%d, %s, %s",myRun[irun],Date[irun].c_str(),Target[irun].c_str()),"p");
	   leg3->SetFillStyle(0);
	   leg3->SetBorderSize(0);
	   leg3->Draw();
	}
	c3->SaveAs("./temp1/gr3.pdf");
	TCanvas* c4 = new TCanvas("c4","c4",900,500);
	c4->Divide(2,1,0.00001,0.00001);
	c4->cd(1);
	gPad->SetPad(0.0,0.0,0.75,1.0);
	gPad->SetRightMargin(0.02);
	gPad->Update();
	mg4->GetXaxis()->CenterTitle();
	mg4->GetYaxis()->CenterTitle();
	mg4->Draw("APL");
        q_edge->Draw();
        latex.DrawLatex(0.70,0.30,"q_edge");
	c4->cd(2);
	gPad->SetPad(0.76,0.0,1.0,1.0);
	gPad->SetLeftMargin(0.0);
	gPad->SetFillStyle(0);
	TLegend* leg4 = new TLegend(0.00,0.50,1.1,0.90);
	leg4->SetTextSize(0.07);
	for(int irun=0;irun<nrun;irun++){
	   leg4->AddEntry(gr4[irun],Form("run%d, %s, %s",myRun[irun],Date[irun].c_str(),Target[irun].c_str()),"p");
	   leg4->SetFillStyle(0);
	   leg4->SetBorderSize(0);
	   leg4->Draw();
	}
	c4->SaveAs("./temp1/gr4.pdf");
	TCanvas* c5 = new TCanvas("c5","c5",900,500);
	c5->Divide(2,1,0.00001,0.00001);
	c5->cd(1);
	gPad->SetPad(0.0,0.0,0.75,1.0);
	gPad->SetRightMargin(0.02);
	gPad->Update();
	mg5->GetXaxis()->CenterTitle();
	mg5->GetYaxis()->CenterTitle();
	mg5->Draw("APL");
        q_edge->Draw();
        latex.DrawLatex(0.70,0.30,"q_edge");
	c5->cd(2);
	gPad->SetPad(0.76,0.0,1.0,1.0);
	gPad->SetLeftMargin(0.0);
	TLegend* leg5 = new TLegend(0.00,0.50,1.1,0.90);
	leg5->SetTextSize(0.07);
	for(int irun=0;irun<nrun;irun++){
	   leg5->AddEntry(gr5[irun],Form("run%d, %s, %s",myRun[irun],Date[irun].c_str(),Target[irun].c_str()),"p");
	   leg5->SetFillStyle(0);
	   leg5->SetBorderSize(0);
	   leg5->Draw();
	}
	c5->SaveAs("./temp1/gr5.pdf");
	gSystem->Exec(Form("pdfunite ./temp1/gr*.pdf ./plots/yloCut_allRunsR.pdf"));
	gSystem->Exec(Form("rm -rf ./temp1/gr*.pdf"));
}
