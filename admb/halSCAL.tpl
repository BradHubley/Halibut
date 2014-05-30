//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Brad Hubley													 
//Date:	May 2014														 
//Purpose:SCAL for Atlantic Halibut.											 
//Notes: 			 
//							 
//																 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	int sim;
	int rseed;
	
	LOCAL_CALCS
		sim=0;
		rseed=0;
		int on,opt;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			sim=1;
			rseed=atoi(ad_comm::argv[on+1]);
		}
	END_CALCS

	// dimensions
	init_int syr;
	init_int eyr;
	init_int nage;
	init_int nsexes;
	init_int nlbin;
	//!! int nlbin2=nlbin;
	init_int nlbin2;
	init_number minbin;
	init_number stepbin;
	init_vector like_weight(1,12);				// weighting vector for likelihood (ala KT)	
	
	// Survey data
	init_vector RVindex(syr,eyr);				// summer RV survey index in stratified total numbers	
	init_matrix RVcatlen(syr,eyr,1,nlbin);		// summer RV survey catch-at-length in stratified total numbers	
	init_int HSsyr;
	init_vector HSindex(HSsyr,eyr);				// halibut survey index glm results
	init_matrix HSMcatlen(HSsyr,eyr,1,nlbin);	// halibut survey catch-at-length males
	init_matrix HSFcatlen(HSsyr,eyr,1,nlbin);	// halibut survey catch-at-length females
	
	// Fishery data
	init_vector LLctM(syr,eyr);					// longline catch, males (tonnes)
	init_vector LLctF(syr,eyr);					// longline catch, females (tonnes)
	init_vector OTct(syr,eyr);					// otter trawl catch (tonnes)
	init_int LLsyr;
	init_matrix LLMcatlen1(LLsyr,eyr,1,nlbin2);	// longline catch-at-length, males
	init_matrix LLFcatlen1(LLsyr,eyr,1,nlbin2);	// longline catch-at-length, females
	init_int OTsyr;
	init_matrix OTcatlen1(OTsyr,eyr,1,nlbin2);	// otter trawl catch-at-length
	
	// strips off last 4 bins
	matrix LLMcatlen(LLsyr,eyr,1,nlbin);
	!! LLMcatlen=trans(trans(LLMcatlen1).sub(1,nlbin));
	matrix LLFcatlen(LLsyr,eyr,1,nlbin);
	!! LLFcatlen=trans(trans(LLFcatlen1).sub(1,nlbin));
	matrix OTcatlen(OTsyr,eyr,1,nlbin);
	!! OTcatlen=trans(trans(OTcatlen1).sub(1,nlbin));


	// LH parameters
	init_vector lwa(1,nsexes);
	init_vector lwb(1,nsexes);
	init_vector linf(1,nsexes);
	init_vector vbk(1,nsexes);
	init_vector t0(1,nsexes);
	init_vector laaM(1,nage);					// length-at-age, males
	init_vector laaF(1,nage);					// length-at-age, females
	init_vector laaM_sigma(1,nage);
	init_vector laaF_sigma(1,nage);
	init_vector mata(1,nsexes);
	init_vector matb(1,nsexes);
	
	// Initilaization parameters
	init_number iro;
	init_number icr;
	init_number irbar;
	init_number ifbar;
	init_number iahat;
	init_number ighat;
	init_number ifull;
	init_number iSDR;
	init_number iSDL;
	init_int eof;
	int iter;
	!!iter=0;
	vector age(1,nage);
	!!age.fill_seqadd(1,1);
	vector iyr(syr,eyr);
	!!iyr.fill_seqadd(syr,1);
	number mM;
	!!mM=1.2*vbk(1);
	number mF;
	!!mF=1.2*vbk(2);
	LOCAL_CALCS
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}else{
			cout<<"Reading data successful."<<endl;
		}
	END_CALCS
	
	// LH parameters
	vector ma(1,nage);					//maturity-at-age
	vector waM(1,nage);					//weight-at-age male
	vector waF(1,nage);					//weight-at-age female
	vector fa(1,nage);					//fecundity-at-age
	LOCAL_CALCS
		for (int a=1;a<=nage;a++){
			ma(a)=1/(1+exp(-matb(2)*((linf(2)-(linf(2)-t0(2))*exp(-vbk(2)*a))-mata(2))));
		}
		waM=lwa(1)*pow(laaF,lwb(1))/1000; 	//convert to kg
		waF=lwa(2)*pow(laaM,lwb(2))/1000; 
		fa=elem_prod(waF,ma);
	END_CALCS
	
	// convert numbers-at-lengths to proportions
	matrix pRVcatlen(syr,eyr,1,nlbin);		
	matrix pHSMcatlen(HSsyr,eyr,1,nlbin);
	matrix pHSFcatlen(HSsyr,eyr,1,nlbin);
	matrix pOTcatlen(OTsyr,eyr,1,nlbin);		
	matrix pLLMcatlen(LLsyr,eyr,1,nlbin);
	matrix pLLFcatlen(LLsyr,eyr,1,nlbin);
	LOCAL_CALCS
		for(int t=syr;t<=eyr;t++){
			pRVcatlen(t) = RVcatlen(t)/sum(RVcatlen(t));
		}
		for(int t=HSsyr;t<=eyr;t++){
			pHSMcatlen(t) = HSMcatlen(t)/sum(HSMcatlen(t));
			pHSFcatlen(t) = HSFcatlen(t)/sum(HSFcatlen(t));
		}
		for(int t=OTsyr;t<=eyr;t++){
			pOTcatlen(t) = OTcatlen(t)/sum(OTcatlen(t));
		}
		for(int t=LLsyr;t<=eyr;t++){
			pLLMcatlen(t) = LLMcatlen(t)/sum(LLMcatlen(t));
			pLLFcatlen(t) = LLFcatlen(t)/sum(LLFcatlen(t));
		}
	END_CALCS
	
	
	// ****Simulations****

PARAMETER_SECTION
	
	// Recruitment parameters
	init_number log_ro(4);
	init_number log_cr(5);
	init_number log_rbar;
	init_bounded_dev_vector wt(syr-nage,eyr,-10.,10.,3);
	
	// Fishing mortalities
	init_number log_LLF_fbar;
	init_number log_LLM_fbar;
	init_number log_OT_fbar;
	init_bounded_dev_vector log_LLF_ft_dev(syr,eyr,-10.,10.,2);
	init_bounded_dev_vector log_LLM_ft_dev(syr,eyr,-10.,10.,2);
	init_bounded_dev_vector log_OT_ft_dev(syr,eyr,-10.,10.,2);
	
	// Selectivities
	init_bounded_number HSsel_half(0,nage,4);
	init_bounded_number HSsel_SD(0,5,4);
	init_bounded_number RVsel_full(0,nage,4);
	init_bounded_number RVsel_SDR(0,5,4);
	init_bounded_number RVsel_SDL(0,5,4);
	
	init_bounded_number LLsel_half(0,nage,4);
	init_bounded_number LLsel_SD(0,5,4);
	init_bounded_number OTsel_full(0,nage,4);
	init_bounded_number OTsel_SDR(0,5,4);
	init_bounded_number OTsel_SDL(0,5,4);
	
	
	init_bounded_number rho(0,1,-1);
	init_bounded_number cvgrow(0,1);
	
	init_number ptdev;
	objective_function_value ofv;
	sdreport_number tdev;
	
	// initialization
	!!log_ro=log(iro);
	!!log_cr=log(icr);
	!!log_rbar=log(irbar);
	!!log_LLF_fbar=log(ifbar);
	!!log_LLM_fbar=log(ifbar);
	!!log_OT_fbar=log(ifbar);
	!!rho=0.5;
	!!cvgrow=0.1;
	!!ptdev=log(1./0.08);
	!!HSsel_half=iahat;
	!!HSsel_SD=ighat;
	!!RVsel_full=ifull;
	!!RVsel_SDR=iSDR;
	!!RVsel_SDL=iSDL;
	!!LLsel_half=iahat;
	!!LLsel_SD=ighat;
	!!OTsel_full=ifull;
	!!OTsel_SDR=iSDR;
	!!OTsel_SDL=iSDL;
	
	number ro;
	number cr;
	number rbar;
	number fbar;
	number so;
	number beta;
	number RVq;
	number HSq;
	number fmsy;
	number msy;
	number bmsy;
	number sig;
	number tau;

	vector HSM_va(1,nage);
	vector HSF_va(1,nage);
	vector RV_va(1,nage);
	vector LLM_va(1,nage);
	vector LLF_va(1,nage);
	vector OT_va(1,nage);
	vector LLM_ft(syr,eyr);
	vector LLF_ft(syr,eyr);
	vector OT_ft(syr,eyr);
	vector M_bt(syr,eyr+1);
	vector F_bt(syr,eyr+1);
	vector LLF_ct_hat(syr,eyr);
	vector LLM_ct_hat(syr,eyr);
	vector OT_ct_hat(syr,eyr);
	vector LLF_ct_resid(syr,eyr);
	vector LLM_ct_resid(syr,eyr);
	vector OT_ct_resid(syr,eyr);
	vector RV_resid(syr,eyr);
	vector HS_resid(HSsyr,eyr);
	vector rt(syr+1,eyr);
	vector rt_resid(syr+1,eyr);
	vector HS_pred(syr,eyr+1);
	vector RV_pred(syr,eyr+1);
	vector lbins(1,nlbin);

	matrix M_Nat(syr,eyr+1,1,nage);
	matrix F_Nat(syr,eyr+1,1,nage);
	matrix M_Fat(syr,eyr,1,nage);
	matrix F_Fat(syr,eyr,1,nage);
	matrix OT_Fat(syr,eyr,1,nage);
	matrix M_Zat(syr,eyr,1,nage);
	matrix F_Zat(syr,eyr,1,nage);
	matrix F_alk(1,nlbin,1,nage);
	matrix M_alk(1,nlbin,1,nage);
	matrix RV_plhat(syr,eyr,1,nlbin);
	matrix HSM_plhat(HSsyr,eyr,1,nlbin);
	matrix HSF_plhat(HSsyr,eyr,1,nlbin);
	matrix OT_plhat(OTsyr,eyr,1,nlbin);
	matrix LLM_plhat(LLsyr,eyr,1,nlbin);
	matrix LLF_plhat(LLsyr,eyr,1,nlbin);

PRELIMINARY_CALCS_SECTION
//  if(sim)
//  {
//  	run_data_simulation();
//  }

PROCEDURE_SECTION
	initialization();
	statedynamics();
	observation_model();
	stock_recruit_model();
	objective_function();

//	if(mceval_phase())
//	{ 
//		forecast();
//		get_CF(value(fmsy),value(msy),value(bmsy));
//		mcmc_output();
//	}
//	if(last_phase())
//	{
//		forecast();
//		get_CF(value(fmsy),value(msy),value(bmsy));
//	} 

FUNCTION initialization
	
	// vulnerabilities-at-age (survey)
	RV_va=sel_dhn(1,nage,RVsel_full,RVsel_SDR,RVsel_SDL);
	HSM_va=plogis(age,HSsel_half,HSsel_SD);
	HSF_va=plogis(age,HSsel_half,HSsel_SD);
	
	// vulnerabilities-at-age (fishery)
	LLM_va=plogis(age,LLsel_half,LLsel_SD);
	LLF_va=plogis(age,LLsel_half,LLsel_SD);
	OT_va=sel_dhn(1,nage,OTsel_full,OTsel_SDR,OTsel_SDL);
	
	// fishing mortality 
	LLM_ft=mfexp(log_LLM_fbar+log_LLM_ft_dev);
	LLF_ft=mfexp(log_LLF_fbar+log_LLF_ft_dev);
	OT_ft=mfexp(log_OT_fbar+log_OT_ft_dev);
	M_Fat=outer_prod(LLM_ft,LLM_va);
	F_Fat=outer_prod(LLF_ft,LLF_va);
	OT_Fat=outer_prod(OT_ft,OT_va);

	// total mortality
	M_Zat=mM+M_Fat+OT_Fat;
	F_Zat=mF+F_Fat+OT_Fat;
	


FUNCTION statedynamics

	// Males
	dvar_vector lxoM=pow(mfexp(-mM),age-1.);
	lxoM(nage)/=(1.-mfexp(-mM));
	M_Nat(syr,1)=mfexp(log_rbar+wt(syr-1))/2;
	for(int j=2;j<=nage;j++) M_Nat(syr,j)=mfexp(log_rbar+wt(syr-j))*lxoM(j)/2;
	for(int i=syr;i<=eyr;i++)
	{
		M_Nat(i+1,1)=mfexp(log_rbar+wt(i))/2;
		M_Nat(i+1)(2,nage)=++elem_prod(M_Nat(i)(1,nage-1),mfexp(-M_Zat(i)(1,nage-1)));
		M_Nat(i+1,nage)+=M_Nat(i,nage)*mfexp(-M_Zat(i,nage));
	}
	M_bt=M_Nat*elem_prod(waM,LLM_va);

	// Females
	dvar_vector lxoF=pow(mfexp(-mF),age-1.);
	lxoF(nage)/=(1.-mfexp(-mF));
	F_Nat(syr,1)=mfexp(log_rbar+wt(syr-1))/2;
	for(int j=2;j<=nage;j++) F_Nat(syr,j)=mfexp(log_rbar+wt(syr-j))*lxoF(j)/2;
	for(int i=syr;i<=eyr;i++)
	{
		F_Nat(i+1,1)=mfexp(log_rbar+wt(i))/2;
		F_Nat(i+1)(2,nage)=++elem_prod(F_Nat(i)(1,nage-1),mfexp(-F_Zat(i)(1,nage-1)));
		F_Nat(i+1,nage)+=F_Nat(i,nage)*mfexp(-F_Zat(i,nage));
	}
	F_bt=F_Nat*elem_prod(waF,LLF_va);


FUNCTION observation_model
	
	// Fishery
	// longline females
	dvar_matrix LLMC(syr,eyr,1,nage);
	LLMC=elem_prod(elem_div(M_Fat,M_Zat),elem_prod(1.-mfexp(-M_Zat),M_Nat));
	LLM_ct_hat=LLMC*waM/1000;
	LLF_ct_resid=log(LLctF)-log(LLF_ct_hat);
	// longline females
	dvar_matrix LLFC(syr,eyr,1,nage);
	LLFC=elem_prod(elem_div(F_Fat,F_Zat),elem_prod(1.-mfexp(-F_Zat),F_Nat));
	LLF_ct_hat=LLFC*waF/1000;
	LLM_ct_resid=log(LLctM)-log(LLM_ct_hat);
	// ottertrawl
	dvar_matrix OTMC(syr,eyr,1,nage);
	dvar_matrix OTFC(syr,eyr,1,nage);
	OTMC=elem_prod(elem_div(OT_Fat,M_Zat),elem_prod(1.-mfexp(-M_Zat),M_Nat));
	OTFC=elem_prod(elem_div(OT_Fat,F_Zat),elem_prod(1.-mfexp(-F_Zat),F_Nat));
	OT_ct_hat=OTFC*waF/1000+OTMC*waM/1000;
	OT_ct_resid=log(OTct)-log(OT_ct_hat);
	
	//Survey
	//RV residuals (walters and ludwig 1994)
	RV_pred=(F_Nat+M_Nat)*RV_va;
	RV_resid=log(RVindex)-log(RV_pred(syr,eyr));
	RVq=mfexp(mean(RV_resid));
	RV_resid-=mean(RV_resid);
	//HS residuals (walters and ludwig 1994)
	HS_pred=(F_Nat*HSF_va+M_Nat*HSM_va);
	HS_resid=log(HSindex)-log(HS_pred(HSsyr,eyr));
	HSq=mfexp(mean(HS_resid));
	HS_resid-=mean(HS_resid);
	
	//Age Length transition matrix
	lbins.fill_seqadd(minbin,stepbin*2.);
	dvar_vector Fsdl=laaF*cvgrow;
	dvar_vector Msdl=laaM*cvgrow;
	dvariable z1;
	dvariable z2;
	F_alk.initialize();
	M_alk.initialize();
	for(int i=1;i<=nage;i++) //loop over ages
	{
		 for(int j=1;j<=nlbin;j++) //loop over length bins
		{
			// Males
			z1=((lbins(j)-stepbin)-laaM(i))/Msdl(i);
			z2=((lbins(j)+stepbin)-laaM(i))/Msdl(i);
			M_alk(j,i)=cumd_norm(z2)-cumd_norm(z1);
			// Females
			z1=((lbins(j)-stepbin)-laaF(i))/Fsdl(i);
			z2=((lbins(j)+stepbin)-laaF(i))/Fsdl(i);
			F_alk(j,i)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
	}//end nage
	//F_alk/=sum(F_alk);
	//ofstream ofs("testFALK.txt");
	//ofs<<setprecision(2)<<F_alk<<endl;
	//ad_exit(1);
	
	
	// Predicted proportions at length
	dvar_vector paM;
	dvar_vector paF;
	dvar_vector pl;
	dvar_vector plF;
	dvar_vector plM;
	
	for(int i=syr;i<=eyr;i++)
	{
		// RV survey
		paM=elem_prod(M_Nat(i),RV_va);
		plM=M_alk*paM;
		paF=elem_prod(F_Nat(i),RV_va);
		plF=F_alk*paF;
		pl=plF+plM;
		RV_plhat(i)=pl/sum(pl);
	}
	
	for(int i=HSsyr;i<=eyr;i++)
	{
		// halibut survey females
		paF=elem_prod(F_Nat(i),HSF_va);
		plF=F_alk*paF;
		HSF_plhat(i)=plF/sum(plF);
		// halibut survey males
		paM=elem_prod(M_Nat(i),HSM_va);
		plM=M_alk*paM;
		HSM_plhat(i)=plM/sum(plM);
		}
	
	for(int i=LLsyr;i<=eyr;i++)
	{
		// longline females
		paF=LLFC(i);
		plF=F_alk*paF;
		LLF_plhat(i)=plF/sum(plF);
		// longline males
		paM=LLMC(i);
		plM=M_alk*paM;
		LLM_plhat(i)=plM/sum(plM);
	}
	
	for(int i=OTsyr;i<=eyr;i++)
	{
		// otter trawl
		paM=OTMC(i);
		plM=M_alk*paM;
		paF=OTFC(i);
		plF=F_alk*paF;
		pl=plF+plM;
		OT_plhat(i)=pl/sum(pl);
	}

FUNCTION stock_recruit_model
	ro=mfexp(log_ro);
	cr=mfexp(log_cr)+1.;
	dvar_vector lxo=pow(mfexp(-mF),age-1.);
	lxo(nage)/=(1.-mfexp(-mF));
	dvariable phieo=lxo*fa;
	so=cr/phieo;
	beta=(cr-1.)/(ro*phieo);
	dvar_vector sbt=F_Nat*fa;
	dvar_vector nmr=so*sbt(syr,eyr-1);
	dvar_vector den=(1.+beta*sbt(syr,eyr-1));
	dvar_vector tmp_rt=++elem_div(nmr,den);
	rt=column(F_Nat.sub(syr+1,eyr),1)+column(M_Nat.sub(syr+1,eyr),1);
	rt_resid=log(tmp_rt)-log(rt);

FUNCTION objective_function 
	
	// Negative Log Likelihoods
	dvar_vector nll_vec(1,12);
	tdev=sqrt(1./mfexp(ptdev));
	sig=sqrt(rho*1./mfexp(ptdev));//process error sd.dev
	tau=sqrt((1.-rho)*1./mfexp(ptdev));
	double tau2;
	dvar_matrix nu(syr,eyr,1,nlbin);
	double pmin=0.0025;
	
	nll_vec.initialize();
	nll_vec(1)=dnorm(LLF_ct_resid,0.1);
	nll_vec(2)=dnorm(LLF_ct_resid,0.1);
	nll_vec(3)=dnorm(OT_ct_resid,0.1);
	nll_vec(4)=dnorm(RV_resid,tau);
	nll_vec(5)=dnorm(HS_resid,tau);
	nll_vec(6)=dnorm(rt_resid,sig);
	nll_vec(7)=dmvlogistic(pRVcatlen,RV_plhat,nu,tau2,pmin);
	nll_vec(8)=dmvlogistic(pHSMcatlen,HSM_plhat,nu,tau2,pmin);
	nll_vec(9)=dmvlogistic(pHSFcatlen,HSF_plhat,nu,tau2,pmin);
	nll_vec(10)=dmvlogistic(pOTcatlen,OT_plhat,nu,tau2,pmin);
	nll_vec(11)=dmvlogistic(pLLMcatlen,LLM_plhat,nu,tau2,pmin);
	nll_vec(12)=dmvlogistic(pLLFcatlen,LLF_plhat,nu,tau2,pmin);
	
	//cout<<tau2<<endl;
	//ad_exit(1);
	
	// Priors
	dvar_vector p_vec(1,5);
	dvariable h=cr/(4.+cr);
	p_vec.initialize();
	p_vec(1)=dbeta((h-0.2)/0.8,2,2);
	if(last_phase())
	{
		p_vec(3)=dnorm(wt,2);
		p_vec(4)=dnorm(log_LLF_ft_dev,2);
	}
	else
	{
		p_vec(3)=100.*norm2(wt);
		p_vec(4)=100.*norm2(log_LLF_ft_dev);
	}
	p_vec(5)=dbeta(rho,50,50);
	
	

    for(int i=1;i<=12;i++)
     {
     ofv+= nll_vec(i)*like_weight(i);
     }

	//ofv=sum(nll_vec);	//+sum(p_vec);

FUNCTION dvar_vector sel_dhn(const int& minage, const int& maxage, dvariable& full, dvariable& varR, dvariable& varL)
	
	dvar_vector selection(minage,maxage);
	selection.initialize();
	
	for (int a=minage;a<=maxage;a++)
	{
	if (full>a)
	{
		selection(a)=exp(-1.0*(square(a-full))/square(varL));             
	}
	else
	{
		selection(a)=exp(-1.0*(square(a-full))/square(varR));             
	}
	}	
	selection = selection,max(selection);

	return(selection);

FUNCTION mcmc_output
//	if(iter==0)
//	{
//		ofstream ofs("refpar.mcmc");
//		ofs<<"fmsy\t bmsy\t msy\t b/bmsy\t f/fmsy"<<endl;
//		//ofs<<"r\t k\t q\t sig\t"<<endl;
//	}
//	iter++;
//	double fratio=value(ft(eyr)/fmsy);
//	double bratio=value(Nt(eyr)*wa/bmsy);
//	ofstream ofs("refpar.mcmc",ios::app);
//	ofs<<fmsy<<"\t"<<bmsy<<"\t"<<msy<<"\t"<<bratio<<"\t"<<fratio<<endl;

FUNCTION forecast

FUNCTION run_data_simulation


//=============================================================================
REPORT_SECTION
	
	// Data
	REPORT(syr);
	REPORT(eyr);
	REPORT(nage);
	REPORT(nsexes);
	REPORT(nlbin);
	REPORT(nlbin2);
	REPORT(minbin);
	REPORT(stepbin);
	REPORT(RVindex);
	REPORT(RVcatlen);
	REPORT(HSsyr);
	REPORT(HSindex);
	REPORT(HSMcatlen);
	REPORT(HSFcatlen);
	REPORT(LLctM);
	REPORT(LLctF);
	REPORT(OTct);
	REPORT(LLsyr);
	REPORT(LLMcatlen);
	REPORT(LLFcatlen);
	REPORT(OTsyr);
	REPORT(OTcatlen);
	REPORT(lwa);
	REPORT(lwb);
	REPORT(linf);
	REPORT(vbk);
	REPORT(t0);
	REPORT(laaM);
	REPORT(laaF);
	REPORT(laaM_sigma);
	REPORT(laaF_sigma);
	REPORT(mata);
	REPORT(matb);

	// Output
	REPORT(F_Nat);
	REPORT(M_Nat);
	
	// selectivity
	REPORT(RV_va);
	REPORT(HSM_va);
	REPORT(HSF_va);
	REPORT(OT_va);
	REPORT(LLM_va);
	REPORT(LLF_va);

	// fishing mortalities
	REPORT(LLM_ft);
	REPORT(LLF_ft);
	REPORT(OT_ft);
	REPORT(F_Zat);
	REPORT(M_Zat);
	REPORT(F_Fat);
	REPORT(M_Fat);
	REPORT(OT_Fat);

	// recruitment
	REPORT(rbar);
	REPORT(rt);
	REPORT(wt);

	// error terms
	REPORT(tau);
	REPORT(sig);
	REPORT(rho);

	// predicted indices
	REPORT(RVq);
	REPORT(RV_pred);
	REPORT(HSq);
	REPORT(HS_pred);
	REPORT(LLF_ct_hat);
	REPORT(LLM_ct_hat);
	REPORT(OT_ct_hat);

	// proportions at length
	REPORT(pRVcatlen);
	REPORT(pHSMcatlen);
	REPORT(pHSFcatlen);
	REPORT(pOTcatlen);
	REPORT(pLLMcatlen);
	REPORT(pLLFcatlen);
	REPORT(RV_plhat);
	REPORT(HSM_plhat);
	REPORT(HSF_plhat);
	REPORT(OT_plhat);
	REPORT(LLM_plhat);
	REPORT(LLF_plhat);

	//double fmsy;
	//double msy;
	//double bmsy;
	//if(last_phase())
	//{
	//	get_CF(fmsy,msy,bmsy);
	//}
//	REPORT(ro);
//	REPORT(cr);
//	REPORT(mfexp(log_fbar));
//	REPORT(mfexp(log_fbar+log_ft_dev));
//	REPORT(plt);
//	REPORT(RV_plhat);
//	REPORT(ahat);
//	REPORT(ghat);
//	REPORT(rho);
//	REPORT(cvgrow);
//	REPORT(ct);
//	REPORT(ct_hat);
//	REPORT(fmsy);
//	REPORT(msy);
//	REPORT(bmsy);

TOP_OF_MAIN_SECTION
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;

