int LoadACC(TString filename, double angle[], double acc[], double acc_err[]){
     ifstream infile;
     infile.open(filename);

     if(!infile.is_open()){
	printf("Can't open file accfunction file\n");
	return 0;
     }

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;

    while(tmp.ReadLine(infile)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,",");
          angle[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,",");
          acc[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          acc_err[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
    infile.close();
    return 1;
}

double FindACC(double theta, double angle[], double acc[], double acc_err[],int size){
    if(angle[0]<0) return -1;   // acceptance table is wrong
    if(theta<angle[0] || theta>=angle[size-1]) return 0;  // exceed the acceptance table
	
    double acceptance=-2;
    for(int ii=0; ii<size-1; ii++){
	if( angle[ii]<=theta && theta<angle[ii+1] ){
	   double acc_low, acc_high;
	   if(acc[ii]<acc_err[ii]) acc_low=0;  // get rid of statistical fluctuation
	   else acc_low = acc[ii];

	   if(acc[ii+1]<acc_err[ii+1]) acc_high=0;
	   else acc_high = acc[ii+1];

	   double a = (acc_high-acc_low)/(angle[ii+1]-angle[ii]);
	   double b = acc_low - a*angle[ii];

	   acceptance = a*theta+b;
	   break;
	}
    }
	
    if(acceptance==-2) printf("FindACC: Something is wrong \n ");
    return acceptance;

}

void GetACCshape(TString filename, Double_t &mean, Double_t &RMS, Double_t &peak){

     // load acceptance table
     double accp_angle[100]={0}, accp[100]={0}, accp_err[100]={0},accp_angle_err[100]={0};
     int status = LoadACC(filename,accp_angle, accp, accp_err);
     if(status==0) exit(0);  // asymmetry table doesn't exist

     TH1F *hACC = new TH1F("hACC","acceptance function",100,3,8);

     int nbin_th = 100000;
     double dtheta = (8.0-3.0)/(nbin_th*1.0);  // delta theta in radius
     for(int ii=0; ii<nbin_th; ii++){
        double thisAngle = 3.0 + dtheta*ii;
        double thisACC = FindACC(thisAngle,accp_angle,accp,accp_err,100);  // acceptance to the corresponding angle
        if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}
        
        hACC->Fill(thisAngle,thisACC/1000.);
     }
   
     mean = hACC->GetMean();
     RMS = hACC->GetRMS();
     Double_t mean_bin = hACC->FindBin(mean);
     peak = hACC->GetBinContent(mean_bin);

     delete hACC;
     return;
}

