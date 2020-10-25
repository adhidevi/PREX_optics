{
// Run 27825 is a thin Pb LHRS run at 1 nA taken on June 12, 2010

 TChain* T = new TChain("T");
 T->Add("/w/halla-scifs17exp/parity/disk1/bob/podd2016/analyzer/Afile_27825.root");
 TH1F *hmomL = new TH1F("hmomL","L-HRS momentum, thin Pb",200,1055,1063);
 TH2F *hthph = new TH2F("hthph","theta-phi acceptance",200,-0.06,0.06,100,-0.06,0.06);
 TH2F *hxy = new TH2F("hxy","X-Y in focal plane",100,-0.8,0.3,100,-0.08,0.08);

 //T->Draw("L.tr.tg_th:L.tr.tg_ph>>hthph","L.tr.tg_th>-0.02&&L.tr.tg_th<0.02&&L.tr.tg_ph>-0.015&&L.tr.tg_ph<0.015");
// T->Draw("L.tr.y:L.tr.x>>hxy","L.tr.y>-0.01&&L.tr.y<0.01");
// TCut cut = "L.tr.n==1&&L.tr.tg_th>-0.02&&L.tr.tg_th<0.02&&L.tr.tg_ph>-0.015&&L.tr.tg_ph<0.015";
 TCut cut = "L.tr.n==1&&L.tr.tg_th>-0.013&&L.tr.tg_ph<0.02&&L.tr.tg_ph>-0.01&&L.tr.tg_ph<0.004";
 T->Draw("1063*(1.+L.tr.tg_dp)>>hmomL",cut);
 
}
