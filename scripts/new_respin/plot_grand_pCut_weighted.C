void plot_grand_pCut_weighted(TString arm){
	TGaxis::SetMaxDigits(3);
	gStyle->SetTitleYOffset(1.3);
	TString runlist;
	runlist = "./TextFiles/good_run.list";
	ifstream runfile(Form("%s",runlist.Data()));
	int runL, runR;
	double ebeam;
	string date, target;
	vector<int>RunL, RunR;
	vector<double>Ebeam;
	vector<string>Date,Target;
	while(runfile>>date>>target>>runL>>runR>>ebeam){
	     Date.push_back(date);
	     Target.push_back(target);
	     RunL.push_back(runL);
	     RunR.push_back(runR);
	     Ebeam.push_back(ebeam);
	}
	runfile.close();
	int nrun = RunL.size();
	int* myRunL = RunL.data();
	int* myRunR = RunR.data();
	double qedge, pcut, qsq, qsq_rms, edge_rel, acc_fr;
	vector<double>Qedge, Pcut, Qsq, Qsq_rms, Edge_rel, Acc_fr, Qsq_central;
	TString filename;
	TString edge_file_name;
	   if(arm=="R"){
	   edge_file_name = "./TextFiles/Qsq_edgeCut_allRunsR.csv";
	   }else if(arm=="L"){
	   edge_file_name = "./TextFiles/Qsq_edgeCut_allRunsL.csv";
	   }
	   ifstream edge_file(Form("%s",edge_file_name.Data()));

	   if(edge_file.is_open())
	   cout<<Form("reading file: %s",edge_file_name.Data())<<endl;
	   if(!edge_file){
	        cout<<Form("File '''%s'''does not exist! Quiting...",edge_file_name.Data())<<endl;
		exit(0);
	   }

	   while(edge_file>>qedge>>pcut>>qsq>>qsq_rms>>edge_rel>>acc_fr){
		Qsq_central.push_back(qsq);
		}

	for(int irun=0;irun<nrun;irun++){
	   if(arm=="R"){
	   filename = Form("./TextFiles/Qsq_pCut_run%d.csv",myRunR[irun]);
	   }else if(arm=="L"){
	   filename = Form("./TextFiles/Qsq_pCut_run%d.csv",myRunL[irun]);
	   }
 	   ifstream infile(Form("%s",filename.Data()));

	   if(infile.is_open())
	   cout<<Form("reading file: %s",filename.Data())<<endl;
	   if(!infile){
	        cout<<Form("File '''%s'''does not exist! Quiting...",filename.Data())<<endl;
		exit(0);
	   }
	   while(infile>>qedge>>pcut>>qsq>>qsq_rms>>edge_rel>>acc_fr){
		Qedge.push_back(qedge);
		Pcut.push_back(pcut);
		Qsq.push_back(qsq);
		Qsq_rms.push_back(qsq_rms);
		Edge_rel.push_back(edge_rel);
		Acc_fr.push_back(acc_fr);
		}

	infile.close();
	int sze = Qedge.size();
	for(int j=0;j<sze;j++)
	cout<<Pcut[j]<<"\t"<<Qsq[j]<<endl;

	gr[irun] = new TGraphErrors(sze,&Qsq_central[0],&Qsq[0]);
	gr[irun]->SetMarkerStyle(20);
	if(irun==4){
	gr[irun]->SetMarkerColor(49);
	}else if(irun==9){
	gr[irun]->SetMarkerColor(41);
	}else{
	gr[irun]->SetMarkerColor(irun+1);
	}
	gr[irun]->Draw("APL");
	Qedge.clear();
	Pcut.clear();
	Qsq.clear();
	Qsq_rms.clear();
	Edge_rel.clear();
	Acc_fr.clear();
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
	TLine* q_edge;
	if(arm=="RHRS")
	q_edge = new TLine(0.9454,-10000,0.9454,10000);
	if(arm=="LHRS")
	q_edge = new TLine(0.9474,-10000,0.9474,10000);
	q_edge->SetLineColor(kOrange);
	q_edge->SetLineWidth(2);
	q_edge->Draw();
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.05);
	latex.SetTextColor(kOrange);
	if(arm=="RHRS")
	latex.DrawLatex(0.75,0.20,"q_edge");
	if(arm=="LHRS")
	latex.DrawLatex(0.75,0.80,"q_edge");
	c1->cd(2);
	gPad->SetPad(0.76,0.0,1.0,1.0);
	gPad->SetLeftMargin(0.0);
	for(int irun=0;irun<nrun;irun++){
	   leg1->AddEntry(gr[irun],Form("run%d, %s, %s",myRun[irun],Date[irun].c_str(),Target[irun].c_str()),"p");
	   leg1->SetFillStyle(0);
	   leg1->SetBorderSize(0);
	   leg1->Draw();
	}
	c1->SaveAs("./temp/gr.pdf");
	TCanvas* c2 = new TCanvas("c2","c2",900,500);
	c2->Divide(2,1,0.00001,0.00001);
	c2->cd(1);
	gPad->SetPad(0.0,0.0,0.75,1.0);
	gPad->SetRightMargin(0.02);
	gPad->Update();
	mg2->GetXaxis()->CenterTitle();
	mg2->GetYaxis()->CenterTitle();
	mg2->Draw("APL");
	gPad->Update();
	q_edge->SetLineColor(kOrange);
	q_edge->SetLineWidth(2);
	q_edge->Draw();
	if(arm=="RHRS")
	latex.DrawLatex(0.75,0.20,"q_edge");
	if(arm=="LHRS")
	latex.DrawLatex(0.75,0.80,"q_edge");
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
	c2->SaveAs("./temp/gr2.pdf");
	TCanvas* c3 = new TCanvas("c3","c3",900,500);
	c3->Divide(2,1,0.00001,0.00001);
	c3->cd(1);
	gPad->SetPad(0.0,0.0,0.75,1.0);
	gPad->SetRightMargin(0.02);
	gPad->Update();
	mg3->GetXaxis()->CenterTitle();
	mg3->GetYaxis()->CenterTitle();
	mg3->Draw("APL");
	gPad->Update();
	q_edge->SetLineColor(kOrange);
	q_edge->SetLineWidth(2);
	q_edge->Draw();

	if(arm=="RHRS")
	latex.DrawLatex(0.75,0.20,"q_edge");
	if(arm=="LHRS")
	latex.DrawLatex(0.75,0.80,"q_edge");
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
	c3->SaveAs("./temp/gr3.pdf");
	TCanvas* c4 = new TCanvas("c4","c4",900,500);
	c4->Divide(2,1,0.00001,0.00001);
	c4->cd(1);
	gPad->SetPad(0.0,0.0,0.75,1.0);
	gPad->SetRightMargin(0.02);
	gPad->Update();
	mg4->GetXaxis()->CenterTitle();
	mg4->GetYaxis()->CenterTitle();
	mg4->Draw("APL");
	gPad->Update();
	q_edge->SetLineColor(kOrange);
	q_edge->SetLineWidth(2);
	q_edge->Draw();
	if(arm=="RHRS")
	latex.DrawLatex(0.75,0.20,"q_edge");
	if(arm=="LHRS")
	latex.DrawLatex(0.75,0.80,"q_edge");
	c4->cd(2);
	gPad->SetPad(0.76,0.0,1.0,1.0);
	gPad->SetLeftMargin(0.0);
	TLegend* leg4 = new TLegend(0.00,0.50,1.1,0.90);
	leg4->SetTextSize(0.07);
	for(int irun=0;irun<nrun;irun++){
	   leg4->AddEntry(gr4[irun],Form("run%d, %s, %s",myRun[irun],Date[irun].c_str(),Target[irun].c_str()),"p");
	   leg4->SetFillStyle(0);
	   leg4->SetBorderSize(0);
	   leg4->Draw();
	}
	c4->SaveAs("./temp/gr4.pdf");
	TCanvas* c5 = new TCanvas("c5","c5",900,500);
	c5->Divide(2,1,0.00001,0.00001);
	c5->cd(1);
	gPad->SetPad(0.0,0.0,0.75,1.0);
	gPad->SetRightMargin(0.02);
	gPad->Update();
	mg5->GetXaxis()->CenterTitle();
	mg5->GetYaxis()->CenterTitle();
	mg5->Draw("APL");
	gPad->Update();
	q_edge->SetLineColor(kOrange);
	q_edge->SetLineWidth(2);
	q_edge->Draw();
	if(arm=="RHRS")
	latex.DrawLatex(0.75,0.20,"q_edge");
	if(arm=="LHRS")
	latex.DrawLatex(0.75,0.80,"q_edge");
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
	c5->SaveAs("./temp/gr5.pdf");
	if(arm=="RHRS")
	gSystem->Exec(Form("pdfunite ./temp/gr*.pdf ./plots/pCut_allRunsR.pdf"));
	if(arm=="LHRS")
	gSystem->Exec(Form("pdfunite ./temp/gr*.pdf ./plots/pCut_allRunsL.pdf"));
	gSystem->Exec(Form("rm -rf ./temp/gr*.pdf"));
}
