void combine_weighted_csv(TString arm){
	TString runlistname = "./TextFiles/good_run.list";
	ifstream runlist(Form("%s",runlistname.Data()));
	if(!runlist){
	cout<<"The file "<<runlistname<<" does not exits! Quiting..."<<endl;
	exit(0);
	}
	string date, tgt;
	int runL, runR;
	double beamE;
	vector<string>Date, Tgt;
	vector<int>RunL, RunR;
	vector<double>BeamE;
	while(runlist>>date>>tgt>>runL>>runR>>beamE){
	Date.push_back(date);
	Tgt.push_back(tgt);
	RunL.push_back(runL);
	RunR.push_back(runR);
	BeamE.push_back(beamE);
	}
	runlist.close();
	int nrun = RunL.size();
	cout<<"I counted "<<nrun<<" runs in each HRS."<<endl;
        int* myRunL = RunL.data();
        int* myRunR = RunR.data();

	TString filename;
        double qedge, pcut, qsq, qsq_rms, edge_rel, acc_fr;
        vector<double>Qedge, Pcut, Qsq, Qsq_rms, Edge_rel, Acc_fr, Qsq_central;
	
	ofstream outfile(Form("./TextFiles/acc_weighted_qsq_all_%sHRS_runs.csv",arm.Data()));

        for(int irun=0;irun<nrun;irun++){
           if(arm=="R"){
           filename = Form("./TextFiles/Qsq_pCut_W_run%d.csv",myRunR[irun]);
           }else if(arm=="L"){
           filename = Form("./TextFiles/Qsq_pCut_W_run%d.csv",myRunL[irun]);
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
	if(arm=="R"){
	cout<<"Run = "<<myRunR[irun]<<endl;
	outfile<<"Run = "<<myRunR[irun]<<endl;
	}else if(arm=="L"){
	cout<<"Run = "<<myRunL[irun]<<endl;
	outfile<<"Run = "<<myRunL[irun]<<endl;
	}
        outfile<<"q_edge"<<"\t"<<"pCut"<<"\t"<<"QsqW"<<"\t"<<"QsqW_rms"<<"\t"<<"q_dge_rel"<<"\t"<<"Acc_fr"<<endl;
        for(int j=0;j<sze;j++){
        cout<<Qedge[j]<<"\t"<<Pcut[j]<<"\t"<<Qsq[j]<<"\t"<<Qsq_rms[j]<<"\t"<<Edge_rel[j]<<"\t"<<Acc_fr[j]<<endl;
        outfile<<Qedge[j]<<"\t"<<Pcut[j]<<"\t"<<Qsq[j]<<"\t"<<Qsq_rms[j]<<"\t"<<Edge_rel[j]<<"\t"<<Acc_fr[j]<<endl;
	}
	Qedge.clear();
	Pcut.clear();
	Qsq.clear();
	Qsq_rms.clear();
	Edge_rel.clear();
	Edge_rel.clear();
	Acc_fr.clear();
	}

        ifstream qfile(Form("./TextFiles/Qsq_edgeCut_W_allRuns%s.csv",arm.Data()));

        if(!qfile){
        cout<<"File does not exist! Quiting..."<<endl;
        exit(0);
        }
        while(qfile>>qedge>>pcut>>qsq>>qsq_rms>>edge_rel>>acc_fr){
		Qedge.push_back(qedge);
               	Pcut.push_back(pcut);
                Qsq.push_back(qsq);
                Qsq_rms.push_back(qsq_rms);
                Edge_rel.push_back(edge_rel);
                Acc_fr.push_back(acc_fr);
                }
        qfile.close();
	
	if(arm=="R")
        outfile<<"RunR"<<"\t"<<"q_edge"<<"\t"<<"QsqW"<<"\t"<<"QsqW_rms"<<"\t"<<"q_dge_rel"<<"\t"<<"Acc_fr"<<endl;
	else if(arm=="L")
        outfile<<"RunL"<<"\t"<<"q_edge"<<"\t"<<"QsqW"<<"\t"<<"QsqW_rms"<<"\t"<<"q_dge_rel"<<"\t"<<"Acc_fr"<<endl;
        for(int j=0;j<nrun;j++){
	if(arm=="R"){
        cout<<RunR[j]<<"\t"<<Qedge[j]<<"\t"<<Qsq[j]<<"\t"<<Qsq_rms[j]<<"\t"<<Edge_rel[j]<<"\t"<<Acc_fr[j]<<endl;
        outfile<<RunR[j]<<"\t"<<Qedge[j]<<"\t"<<Qsq[j]<<"\t"<<Qsq_rms[j]<<"\t"<<Edge_rel[j]<<"\t"<<Acc_fr[j]<<endl;
	}else if(arm=="L"){
        cout<<RunL[j]<<"\t"<<Qedge[j]<<"\t"<<Qsq[j]<<"\t"<<Qsq_rms[j]<<"\t"<<Edge_rel[j]<<"\t"<<Acc_fr[j]<<endl;
        outfile<<RunL[j]<<"\t"<<Qedge[j]<<"\t"<<Qsq[j]<<"\t"<<Qsq_rms[j]<<"\t"<<Edge_rel[j]<<"\t"<<Acc_fr[j]<<endl;
	}
	}
	Qedge.clear();
	Pcut.clear();
	Qsq.clear();
	Qsq_rms.clear();
	Edge_rel.clear();
	Edge_rel.clear();
	Acc_fr.clear();

}
