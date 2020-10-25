void hamc_data(){
	gStyle->SetOptStat(0);
	TFile* hamc_file = new TFile("../../hamc/hamc.root");
	TFile* data_file = new TFile("/lustre19/expphy/volatile/halla/parity/devi/rootfiles/prexRHRS_21188_-1.root");
	TTree* hamc = (TTree*)hamc_file->Get("hamc");
	TTree* data = (TTree*)data_file->Get("T");
	TH1F* hs = new TH1F("hs","PREX-2 p;p;",200,940,954);
	TH1F* hd = new TH1F("hd","PREX-2 p;p;",200,940,954);
	hs->SetLineColor(kBlue);
	hd->SetLineColor(kRed);

	hamc->Draw("950*(1+dpp)>>hs","","goff");
	data->Draw("1000*R.gold.p>>hd","(R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1)&&(R.tr.th[0]<0.05&&R.tr.th[0]>-0.2&&R.tr.ph[0]<0.1&&R.tr.ph[0]>-0.1)&&(R.tr.tg_th[0]<0.055&&R.tr.tg_th[0]>-0.055&&R.tr.tg_ph[0]>-0.018&&R.tr.tg_ph[0]<0.026)&&(P.evtypebits&2)==2","goff");
	hs->Draw("hist");
	hd->Draw("hist same");
	hs->Scale(1./hs->Integral());
	hd->Scale(1./hd->Integral());
	double firstIP = 2.615;
	double secondIP = 3.198;
	double thirdIP = 3.475;
	double fourthIP = 3.708;
	double elasticPeak = hs->GetBinCenter(hs->GetMaximumBin());
	double firstIPPeak = hs->GetBinCenter(hs->GetMaximumBin())-firstIP;
	double secondIPPeak = hs->GetBinCenter(hs->GetMaximumBin())-secondIP;
	double thirdIPPeak = hs->GetBinCenter(hs->GetMaximumBin())-thirdIP;
	double fourthIPPeak = hs->GetBinCenter(hs->GetMaximumBin())-fourthIP;
	cout<<hs->Interpolate(firstIPPeak)/hs->GetMaximum()<<endl;
	cout<<hs->Interpolate(secondIPPeak)/hs->GetMaximum()<<endl;
	cout<<hs->Interpolate(thirdIPPeak)/hs->GetMaximum()<<endl;
	cout<<hs->Interpolate(fourthIPPeak)/hs->GetMaximum()<<endl;
	int nline = 5;
	TLine* line[nline];
	line[0] = new TLine(elasticPeak,0.0,elasticPeak,hs->GetMaximum());
	line[1] = new TLine(firstIPPeak,0.0,firstIPPeak,hs->GetMaximum());
	line[2] = new TLine(secondIPPeak,0.0,secondIPPeak,hs->GetMaximum());
	line[3] = new TLine(thirdIPPeak,0.0,thirdIPPeak,hs->GetMaximum());
	line[4] = new TLine(fourthIPPeak,0.0,fourthIPPeak,hs->GetMaximum());
	for(int iline = 0;iline<nline;iline++){
	line[iline]->SetLineColor(iline+2);
	line[iline]->Draw();
	}
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.04);
	latex.SetTextColor(kBlue);
	latex.DrawLatex(0.2,0.80,"HAMC");
	latex.SetTextColor(kRed);
	latex.DrawLatex(0.2,0.75,"DATA");
}
