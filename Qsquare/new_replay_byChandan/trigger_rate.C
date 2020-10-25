void trigger_rate(){
	ifstream infile("../TextFiles/good_runs.list");
	string date, tgt;
	int runL, runR;
	vector<string>Date,Tgt;
	vector<int>RunL,RunR;
	while(infile>>date>>tgt>>runL>>runR){
	Date.push_back(date);
	Tgt.push_back(tgt);
	RunL.push_back(runL);
	RunR.push_back(runR);
	}
	infile.close();
	int nrun = RunL.size();
	cout<<"Total runs = "<<nrun<<endl;
	
	TCanvas* c1 = new TCanvas("c1","c1",900,600);
	for(int irun=0;irun<nrun;irun++){
	TChain* T = new TChain("TSLeft");
	T->Add(Form("/chafs1/work1/prex_counting/QsqRootFiles/prexLHRS_%d_-1*.root",RunL[irun]));
	T->Draw(Form("LeftT1_r:Entry$>>h_%d",irun),"","goff");
	TH2D* h1 = (TH2D*)gDirectory->FindObject(Form("h_%d",irun));
	h1->SetTitle(Form("LeftT1_r vs Entry$ (run%d);Entry$;LeftT1_r",RunL[irun]));	
	h1->SetMarkerStyle(20);
	h1->Draw();
	c1->SaveAs(Form("./temp/triggerT1_run%d.pdf",RunL[irun]));
	}
	gSystem->Exec(Form("pdfunite ./temp/trigger* ./plots/trigger_T1_all_runs.pdf"));
	gSystem->Exec(Form("rm -rf ./temp/trigger*"));
}
