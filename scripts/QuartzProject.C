

//Author - Ryan Richards 
//Macro to look at 1D Distributions at Detector Plane to evaluate ADC cuts and x cuts

void QuartzProject(int run, double val){

TChain *T = new TChain("T");
T->Add(Form("prexLHRS_%d_-1.root",run));


TCut vdc_cut = "L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1";
TCut adc_cut = Form("P.upQadcL>%f",val);
TCut adc_miss = Form("P.upQadcL<%f",val);

TCut totalCut = vdc_cut&&adc_cut;
TCut totalMiss = vdc_cut&&adc_miss;

TLine *adc = new TLine(val,0,val,10000);
adc->SetLineColor(kRed);


T->SetAlias("xq","L.tr.x+0.9*L.tr.th");
T->SetAlias("yq","L.tr.y+0.9*L.tr.ph");

//0 is all the projected tracks, 1 is above ADC cut, 2 below ADC
TH1F *hx[3], *hy[3];
int color[3] = {1,2,4};

//Better to draw ADC spectra as well -- zooming in I think
TH1F *spec = new TH1F("spec","Upstream Quartz ADC Spectrum",200,300,1000);
T->Project(spec->GetName(),"P.upQadcL","");


for(int i = 0; i < 3; i++){ 

 hx[i] = new TH1F(Form("hx[%i]",i),"X at Detector Plane",600,-0.5,0.1);
 hy[i] = new TH1F(Form("hy[%i]",i),"Y at Detector Plane",600,-0.1,0.1);

 hx[i]->SetLineColor(color[i]); hy[i]->SetLineColor(color[i]);
 hx[i]->GetXaxis()->SetTitle("L.tr.x+0.9*L.tr.th (m)");
 hy[i]->GetXaxis()->SetTitle("L.tr.y+0.9*L.tr.ph (m)");


 if( i == 0 ) { T->Project(hx[i]->GetName(),"xq",vdc_cut); T->Project(hy[i]->GetName(),"yq",vdc_cut); }
 else if( i == 1 ) { T->Project(hx[i]->GetName(),"xq",totalCut); T->Project(hy[i]->GetName(),"yq",totalCut); }
 else { T->Project(hx[i]->GetName(),"xq",totalMiss); T->Project(hy[i]->GetName(),"yq",totalMiss); }


}


TCanvas *c = new TCanvas();
c->SetLogy();
spec->Draw();
adc->Draw("same");

TCanvas *c1 = new TCanvas();
hx[0]->Draw();
hx[1]->Draw("same");
hx[2]->Draw("same");

/*
TCanvas *c2 = new TCanvas();
hy[0]->Draw();
hy[1]->Draw("same");
hy[2]->Draw("same");
*/

TCanvas *c3 = new TCanvas();
hx[0]->Draw();
hx[1]->Draw("same");
hx[2]->Draw("same");
c3->SetLogy();

}


