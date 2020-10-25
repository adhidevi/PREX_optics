void drawS0L(int run_num){
	TChain* T = new TChain("T");
	T->Add(Form("/lustre19/expphy/volatile/halla/ryanrich/prexRootFiles/prexLHRS_%d_-*.root",run_num));
        TCanvas* c3n = new TCanvas("c3n","S0A and S0B hits",1000,800);
        c3n->Divide(2,2);
        c3n->cd(1);
        gPad->SetLogy(1);
        T->Draw("Ndata.P.S0AtdcL");
        c3n->cd(2);
        gPad->SetLogy(1);
        T->Draw("Ndata.P.S0BtdcL");
	int nhitAmax = T->GetMaximum("Ndata.P.S0AtdcL");
	int nhitBmax = T->GetMaximum("Ndata.P.S0BtdcL");
	c3n->cd(3);
	gPad->SetLogy(1);
	gStyle->SetOptStat(0);
	cout<<Form("nhitAmax: %d, nhitBmax: %d",nhitAmax, nhitBmax)<<endl;
	const int nhit = nhitAmax;
	TH1F* hA[4];
	for(int i=0;i<4;i++){
	T->Draw(Form("P.S0AtdcL[%d]>>h%d",i,i),"","goff");
	hA[i] = (TH1F*)gDirectory->FindObject(Form("h%d",i));
	}
	T->Draw("P.S0AtdcL>>hAA","","goff");
	TH1F* hAA = (TH1F*)gDirectory->FindObject("hAA");
	hAA->SetLineColor(30);
	hAA->Draw();
	for(int i=0;i<4;i++){
	hA[i]->SetLineColor(i+1);
	hA[i]->Draw("same");
	}
	TLatex latex;
	latex.SetTextSize(0.04);
	latex.SetNDC();
	latex.SetTextColor(30);
	latex.DrawLatex(0.6,0.85,"P.S0AtdcL");
	latex.SetTextColor(1);
	latex.DrawLatex(0.6,0.80,"P.S0AtdcL[0]");
	latex.SetTextColor(2);
	latex.DrawLatex(0.6,0.75,"P.S0AtdcL[1]");
	latex.SetTextColor(3);
	latex.DrawLatex(0.6,0.70,"P.S0AtdcL[2]");
	latex.SetTextColor(4);
	latex.DrawLatex(0.6,0.65,"P.S0AtdcL[3]");
	c3n->cd(4);
	gPad->SetLogy(1);
	TH1F* hB[4];
	for(int i=0;i<4;i++){
	T->Draw(Form("P.S0BtdcL[%d]>>hB%d",i,i),"","goff");
	hB[i] = (TH1F*)gDirectory->FindObject(Form("hB%d",i));
	}
	T->Draw("P.S0BtdcL>>hBB","","goff");
	TH1F* hBB = (TH1F*)gDirectory->FindObject("hBB");
	hBB->SetLineColor(30);
	hBB->Draw();
	for(int i=0;i<4;i++){
	hB[i]->SetLineColor(i+1);
	hB[i]->Draw("same");
	}
	latex.SetTextColor(30);
	latex.DrawLatex(0.6,0.85,"P.S0BtdcL");
	latex.SetTextColor(1);
	latex.DrawLatex(0.6,0.80,"P.S0BtdcL[0]");
	latex.SetTextColor(2);
	latex.DrawLatex(0.6,0.75,"P.S0BtdcL[1]");
	latex.SetTextColor(3);
	latex.DrawLatex(0.6,0.70,"P.S0BtdcL[2]");
	latex.SetTextColor(4);
	latex.DrawLatex(0.6,0.65,"P.S0BtdcL[3]");
	TCanvas* ct = new TCanvas("ct","t_A-t_B vs VDC",600,400);
	T->Draw("P.S0AtdcL-P.S0BtdcL:L.tr.x[0]>>hABX(200,-0.5,0.5,200,0,2500)");

        TCanvas* c3 = new TCanvas("c3","S0A and S0B hits",1200,800);
        c3->Divide(3,2);
        c3->cd(1);
        gPad->SetLogy(1);
        T->Draw("Ndata.P.S0AtdcL");
        c3->cd(2);
	gPad->SetLogy(1);
	gStyle->SetOptStat(0);
	T->SetLineColor(1);
	T->Draw("P.S0AtdcL[0]","L.vdc.u1.rawtime>860 && L.vdc.u1.rawtime<1850");
	c3->cd(3);
	T->Draw("P.S0AtdcL[0]>>hAA0","L.vdc.u1.rawtime>860 && L.vdc.u1.rawtime<1850");
	TH1F* hAA0 = (TH1F*)gDirectory->FindObject("hAA0");
	hAA0->GetYaxis()->SetRangeUser(-100,1000);
	hAA0->Draw();
	c3->cd(4);
        gPad->SetLogy(1);
        T->Draw("Ndata.P.S0BtdcL");
	c3->cd(5);
	gPad->SetLogy(1);
	T->Draw("P.S0BtdcL[0]","L.vdc.u1.rawtime>860 && L.vdc.u1.rawtime<1850");
	gPad->Update();
	c3->cd(6);
	T->Draw("P.S0BtdcL[0]>>hBB0","L.vdc.u1.rawtime>860 && L.vdc.u1.rawtime<1850");
	TH1F* hBB0 = (TH1F*)gDirectory->FindObject("hBB0");
	hBB0->GetYaxis()->SetRangeUser(-100,1000);
	hBB0->Draw();
	TCanvas* ct1 = new TCanvas("ct1","t_A-t_B vs VDC",1000,1000);
	ct1->Divide(2,2);
	ct1->cd(1);
	T->Draw("P.S0AtdcL[0]-P.S0BtdcL[0]:L.tr.x[0]>>hABX1(200,-1,1,200,-100,100)","L.vdc.u1.rawtime>860 && L.vdc.u1.rawtime<1850");
	ct1->cd(2);
	T->Draw("P.S0AtdcL[0]-P.S0BtdcL[0]>>hAB1(200,-100,100)","L.vdc.u1.rawtime>860 && L.vdc.u1.rawtime<1850");
	ct1->cd(3);
	T->Draw("P.S0AtdcL[0]-P.S0BtdcL[0]:L.tr.x[0]>>hABX2(200,-1,1,200,-500,500)","L.vdc.u1.rawtime>860 && L.vdc.u1.rawtime<1850");
	ct1->cd(4);
	T->Draw("P.S0AtdcL[0]-P.S0BtdcL[0]>>hAB2(200,-500,500)","L.vdc.u1.rawtime>860 && L.vdc.u1.rawtime<1850");
}
