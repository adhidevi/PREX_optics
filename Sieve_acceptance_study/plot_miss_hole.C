void plot_miss_hole(){
	gStyle->SetOptStat(0);
	TChain* T = new TChain("T");
	for(int i=1;i<6;i++)
	T->Add(Form("/work/halla/parity/disk1/ryanrich/sandwich/SandwichLHRS_PREX_0.0_%d.root",i));
	TH2F* h2 = new TH2F("h2","Theta vs Phi sim Pb;Phi;Theta",100,-0.05,0.05,100,-0.06,0.06);
	T->Draw("th_ztarg_tr:ph_ztarg_tr>>h2","x_vdc_tr!=-333&&!(th_ztarg_tr>-0.0038&&th_ztarg_tr<0.0032&&ph_ztarg_tr>-0.0175&&ph_ztarg_tr<-0.0095)&&!(th_ztarg_tr>0.027&&th_ztarg_tr<0.0375&&ph_ztarg_tr>0.012&&ph_ztarg_tr<0.023)&&!(th_ztarg_tr>-0.0375&&th_ztarg_tr<-0.027&&ph_ztarg_tr>0.012&&ph_ztarg_tr<0.023)");
	TBox* boxIN = new TBox(-0.0175,-0.0038,-0.0095,0.0032);
	boxIN->SetLineColor(kRed);
	boxIN->SetFillColor(0);
	boxIN->SetFillStyle(0);
	boxIN->SetLineWidth(2);
	boxIN->Draw();
	TBox* boxOT = new TBox(0.012,0.027,0.023,0.0375);
	boxOT->SetLineColor(kRed);
	boxOT->SetFillColor(0);
	boxOT->SetFillStyle(0);
	boxOT->SetLineWidth(2);
	boxOT->Draw();
	TBox* boxOB = new TBox(0.012,-0.027,0.023,-0.0375);
	boxOB->SetLineColor(kRed);
	boxOB->SetFillColor(0);
	boxOB->SetFillStyle(0);
	boxOB->SetLineWidth(2);
	boxOB->Draw();

}
