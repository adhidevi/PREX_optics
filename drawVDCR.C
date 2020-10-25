void drawVDCR(int run_num){
	TChain* T = new TChain("T");
	T->Add(Form("/lustre19/expphy/volatile/halla/parity/ryanrich/prexRootFiles/prexRHRS_%d_-1*.root",run_num));
	TCanvas* c1 = new TCanvas("c1","VDC ncluster",1000,600);
	c1->Divide(2,2);
	c1->cd(1);
	T->Draw("R.vdc.u1.nclust");
	c1->cd(2);
	T->Draw("R.vdc.v1.nclust");
	c1->cd(3);
	T->Draw("R.vdc.u2.nclust");
	c1->cd(4);
	T->Draw("R.vdc.v2.nclust");
	
	TLine* u1time1 = new TLine(1850,0.0,1850,100000);
	u1time1->SetLineColor(2);
	TLine* u1time2 = new TLine(860,0.0,860,100000);
	u1time2->SetLineColor(2);

	TCanvas* c2 = new TCanvas("c2","VDC time",1000,600);
 	c2->Divide(2,2);
 	c2->cd(1);
	T->SetLineColor(1);
 	T->Draw("R.vdc.u1.rawtime>>hvdctu1(200,700,2500)");
	T->SetLineColor(2);
 	T->Draw("R.vdc.u1.rawtime[0]>>hvdctu1_0(200,700,2500)","","same");
	T->SetLineColor(3);
 	T->Draw("R.vdc.u1.rawtime[1]>>hvdctu1_1(200,700,2500)","","same");
	T->SetLineColor(4);
 	T->Draw("R.vdc.u1.rawtime[2]>>hvdctu1_2(200,700,2500)","","same");
	T->SetLineColor(6);
 	T->Draw("R.vdc.u1.rawtime[3]>>hvdctu1_3(200,700,2500)","","same");
	T->SetLineColor(7);
 	T->Draw("R.vdc.u1.rawtime[4]>>hvdctu1_4(200,700,2500)","","same");
	
	u1time1->Draw();
	u1time2->Draw();

 	c2->cd(2);
	T->SetLineColor(1);
 	T->Draw("R.vdc.v1.rawtime>>hvdctv1(200,700,2500)");
	T->SetLineColor(2);
 	T->Draw("R.vdc.v1.rawtime[0]>>hvdctv1_0(200,700,2500)","","same");
	T->SetLineColor(3);
 	T->Draw("R.vdc.v1.rawtime[1]>>hvdctv1_1(200,700,2500)","","same");
	T->SetLineColor(4);
 	T->Draw("R.vdc.v1.rawtime[2]>>hvdctv1_2(200,700,2500)","","same");
	T->SetLineColor(6);
 	T->Draw("R.vdc.v1.rawtime[3]>>hvdctv1_3(200,700,2500)","","same");
	T->SetLineColor(7);
 	T->Draw("R.vdc.v1.rawtime[4]>>hvdctv1_4(200,700,2500)","","same");
	u1time1->Draw();
	u1time2->Draw();

 	c2->cd(3);
	T->SetLineColor(1);
 	T->Draw("R.vdc.u2.rawtime>>hvdctu2(200,700,2500)");
	T->SetLineColor(2);
 	T->Draw("R.vdc.u2.rawtime[0]>>hvdctu2_0(200,700,2500)","","same");
	T->SetLineColor(3);
 	T->Draw("R.vdc.u2.rawtime[1]>>hvdctu2_1(200,700,2500)","","same");
	T->SetLineColor(4);
 	T->Draw("R.vdc.u2.rawtime[2]>>hvdctu2_2(200,700,2500)","","same");
	T->SetLineColor(6);
 	T->Draw("R.vdc.u2.rawtime[3]>>hvdctu2_3(200,700,2500)","","same");
	T->SetLineColor(7);
 	T->Draw("R.vdc.u2.rawtime[4]>>hvdctu2_4(200,700,2500)","","same");
	u1time1->Draw();
	u1time2->Draw();

 	c2->cd(4);
	T->SetLineColor(1);
	T->Draw("R.vdc.v2.rawtime>>hvdctv2(200,700,2500)");
	T->SetLineColor(2);
	T->Draw("R.vdc.v2.rawtime[0]>>hvdctv2_0(200,700,2500)","","same");
	T->SetLineColor(3);
	T->Draw("R.vdc.v2.rawtime[1]>>hvdctv2_1(200,700,2500)","","same");
	T->SetLineColor(4);
	T->Draw("R.vdc.v2.rawtime[2]>>hvdctv2_2(200,700,2500)","","same");
	T->SetLineColor(6);
	T->Draw("R.vdc.v2.rawtime[3]>>hvdctv2_3(200,700,2500)","","same");
	T->SetLineColor(7);
	T->Draw("R.vdc.v2.rawtime[4]>>hvdctv2_4(200,700,2500)","","same");
	u1time1->Draw();
	u1time2->Draw();

	TCanvas* c10 = new TCanvas("c10","VDC wire",1000,600);
	c10->Divide(2,2);
	c10->cd(1);
	T->SetLineColor(kBlue);
	T->Draw("R.vdc.u1.wire");
	c10->cd(2);
	T->Draw("R.vdc.v1.wire");
	c10->cd(3);
	T->Draw("R.vdc.u2.wire");
	c10->cd(4);
	T->Draw("R.vdc.v2.wire");

	TCanvas* c3 = new TCanvas("c3","S0A vs S0B",1000,600);
	c3->Divide(2,2);
	c3->cd(1);
	T->Draw("P.S0AadcR>>hS0Aadc(200,0,2500)");
	c3->cd(2);
	T->Draw("P.S0BadcR>>hS0Badc(200,0,2500)");
	c3->cd(3);
	T->Draw("P.S0AadcR:P.S0BadcR>>hS0ABadc(200,0,2500,200,0,2500)","","colz");
	
	TCanvas* c3n = new TCanvas("c3n","S0A and S0B hits",1000,600);
	c3n->Divide(1,2);
	gPad->SetLogy(1);
	c3n->cd(1);
	T->Draw("Ndata.P.S0AtdcR");
	c3n->cd(2);
	gPad->SetLogy(1);
	T->Draw("Ndata.P.S0BtdcR");

	TCanvas* c3t = new TCanvas("c3t","S0A vs S0B",1000,900);
	c3t->Divide(2,3);
	c3t->cd(1);
	gPad->SetLogy(1);
	T->Draw("P.S0AtdcR>>hS0Atdc(200,0,2500)");
	c3t->cd(2);
	gPad->SetLogy(1);
	T->Draw("P.S0BtdcR>>hS0Btdc(200,0,2500)");
	c3t->cd(3);
	gPad->SetLogy(1);
	T->Draw("P.S0AtdcR>>hS0Atdc1(200,0,2500)");
	c3t->cd(4);
	gPad->SetLogy(1);
	T->Draw("P.S0BtdcR>>hS0Btdc1(200,0,2500)");
	c3t->cd(5);
	T->Draw("P.S0AtdcR:P.S0BtdcR>>hS0ABtdc(200,0,2500,200,0,2500)","","colz");

	TCanvas* c4 = new TCanvas("c4","S0 vs VDC position",1000,600);
	c4->Divide(2,2);
	c4->cd(1);
	T->Draw("P.S0AadcR:R.tr.x[0]>>hAx(200,-0.5,0.5,200,0,2500)");
	c4->cd(2);
	T->Draw("P.S0AadcR:R.tr.y[0]>>hAy(200,-0.5,0.5,200,0,2500)");
	c4->cd(3);
	T->Draw("P.S0BadcR:R.tr.x[0]>>hBx(200,-0.5,0.5,200,0,2500)");
	c4->cd(4);
	T->Draw("P.S0BadcR:R.tr.y[0]>>hBy(200,-0.5,0.5,200,0,2500)");

//use detL5.C for S3 adc distribution	

	TCanvas* c5 = new TCanvas("c5","VDC nhit",1000,600);
	c5->Divide(2,2);
	c5->cd(1);
	gPad->SetLogy(1);
	T->Draw("R.vdc.u1.nhit>>hu1hit(30,0,30)");
	c5->cd(2);
	gPad->SetLogy(1);
	T->Draw("R.vdc.v1.nhit>>hv1hit(30,0,30)");
	c5->cd(3);
	gPad->SetLogy(1);
	T->Draw("R.vdc.u2.nhit>>hu2hit(30,0,30)");
	c5->cd(4);
	gPad->SetLogy(1);
	T->Draw("R.vdc.v2.nhit>>hv2hit(30,0,30)");

	TCanvas* c6 = new TCanvas("c6","VDC vdc time vs wire",1000,600);
	c6->Divide(2,2);
	c6->cd(1);
	T->Draw("R.vdc.u1.rawtime:R.vdc.u1.wire>>h1(200,0,400,200,0,3000)");
	c6->cd(2);
	T->Draw("R.vdc.v1.rawtime:R.vdc.v1.wire>>h2(200,0,400,200,0,3000)");
	c6->cd(3);
	T->Draw("R.vdc.u2.rawtime:R.vdc.u2.wire>>h3(200,0,400,200,0,3000)");
	c6->cd(4);
	T->Draw("R.vdc.v2.rawtime:R.vdc.v2.wire>>h4(200,0,400,200,0,3000)");
}
