void momRcut(){
 TChain* T = new TChain("T");
 T->Add("/w/halla-scifs17exp/parity/disk1/bob/podd2016/analyzer/Afile_6908.root");
 TH1F *hmomR = new TH1F("hmomR","R-HRS momentum, thin Pb",200,1055,1063);
 TH2F *hthph = new TH2F("hthph","theta-phi acceptance",200,-0.06,0.06,100,-0.06,0.06);
 TH2F *hxy = new TH2F("hxy","X-Y in focal plane",100,-0.8,0.3,100,-0.08,0.08);
 
// TCut cut = "R.tr.n==1&&R.tr.tg_th>-0.02&&R.tr.tg_th<0.02&&R.tr.tg_ph>-0.015&&R.tr.tg_ph<0.015";
 TCut cut = "R.tr.n==1&&R.tr.th<0.005&&R.tr.th>-0.02&&R.tr.ph<0.015&&R.tr.ph>-0.015";
 T->Draw("1063*(1.+R.tr.tg_dp)>>hmomR",cut);
} 
