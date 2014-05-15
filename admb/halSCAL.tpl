//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Brad Hubley													 
//Date:	May 2014														 
//Purpose:SCAL for Atlantic Halibut.											 
//Notes: 	males=1, females=2			 
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
//	!!ad_comm::change_datafile_name("halSCAL.ctl");
//	init_number simsig;
//	init_number simtau;
//	init_number simqo;
//	init_number simao;
//	init_number simro;
//	init_number simcr;
//	init_number simrbar;
//	init_number simahat;
//	init_number simghat;
//	init_number simcvgrow;
//	init_vector simF(syr,eyr);
//	init_int eofc;
//	LOCAL_CALCS
//		if(eofc!=999)
//		{
//			cout<<"Error reading control file.\n Fix it."<<endl;
//			ad_exit(1);
//		}
//	END_CALCS

PARAMETER_SECTION
	
	// estimated parameters
	init_number log_ro(4);
	init_number log_cr(5);
	init_number log_rbar;
	
	init_number log_LLF_fbar;
	init_number log_LLM_fbar;
	init_number log_OT_fbar;
	init_bounded_dev_vector log_LLF_ft_dev(syr,eyr,-10.,10.,3);
	init_bounded_dev_vector log_LLM_ft_dev(syr,eyr,-10.,10.,3);
	init_bounded_dev_vector log_OT_ft_dev(syr,eyr,-10.,10.,3);
	
	init_number ptdev;
	init_bounded_number ahat(0,nage);
	init_bounded_number ghat(0,5);
	init_bounded_number rho(0,1);
	init_bounded_number cvgrow(0,1);
	init_bounded_dev_vector wt(syr-nage,eyr,-10.,10.,2);
	
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
	!!ahat=iahat;
	!!ghat=ighat;
	
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
	matrix HSM_plhat(syr,eyr,1,nlbin);
	matrix HSF_plhat(syr,eyr,1,nlbin);
	matrix OT_plhat(syr,eyr,1,nlbin);
	matrix LLM_plhat(syr,eyr,1,nlbin);
	matrix LLF_plhat(syr,eyr,1,nlbin);

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
	RV_va=plogis(age,ahat,ghat);
	HSM_va=plogis(age,ahat,ghat);
	HSF_va=plogis(age,ahat,ghat);
	
	// vulnerabilities-at-age (fishery)
	LLM_va=plogis(age,ahat,ghat);
	LLF_va=plogis(age,ahat,ghat);
	OT_va=plogis(age,ahat,ghat);
	
	// fishing mortality 
	LLM_ft=mfexp(log_LLM_fbar+log_LLM_ft_dev);
	LLF_ft=mfexp(log_LLF_fbar+log_LLF_ft_dev);
	OT_ft=mfexp(log_OT_fbar+log_OT_ft_dev);
	M_Fat=outer_prod(LLM_ft,LLM_va);
	F_Fat=outer_prod(LLF_ft,LLF_va);
	OT_Fat=outer_prod(OT_ft,OT_va)*0.5;

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
	LLM_ct_hat=LLMC*waM;
	LLF_ct_resid=log(LLctF)-log(LLF_ct_hat);
	// longline females
	dvar_matrix LLFC(syr,eyr,1,nage);
	LLFC=elem_prod(elem_div(F_Fat,F_Zat),elem_prod(1.-mfexp(-F_Zat),F_Nat));
	LLF_ct_hat=LLFC*waF;
	LLM_ct_resid=log(LLctM)-log(LLM_ct_hat);
	// ottertrawl
	dvar_matrix OTMC(syr,eyr,1,nage);
	dvar_matrix OTFC(syr,eyr,1,nage);
	OTMC=elem_prod(elem_div(OT_Fat,M_Zat),elem_prod(1.-mfexp(-M_Zat),M_Nat));
	OTFC=elem_prod(elem_div(OT_Fat,F_Zat),elem_prod(1.-mfexp(-F_Zat),F_Nat));
	OT_ct_hat=OTFC*waF+OTMC*waM;
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
		// longline females
		paF=LLFC(i);
		plF=F_alk*paF;
		LLF_plhat(i)=plF/sum(plF);
		// longline males
		paM=LLMC(i);
		plM=M_alk*paM;
		LLM_plhat(i)=plM/sum(plM);
		// otter trawl
		paM=OTMC(i);
		plM=M_alk*paM;
		paF=OTFC(i);
		plF=F_alk*paF;
		pl=plF+plM;
		OT_plhat(i)=pl/sum(pl);
		
		// halibut survey females
		paF=elem_prod(F_Nat(i),HSF_va);
		plF=F_alk*paF;
		HSF_plhat(i)=plF/sum(plF);
		// halibut survey males
		paM=elem_prod(M_Nat(i),HSM_va);
		plM=M_alk*paM;
		HSM_plhat(i)=plM/sum(plM);
		// RV survey
		paM=elem_prod(M_Nat(i),RV_va);
		plM=M_alk*paM;
		paF=elem_prod(F_Nat(i),RV_va);
		plF=F_alk*paF;
		pl=plF+plM;
		RV_plhat(i)=pl/sum(pl);
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
	
	nll_vec.initialize();
	nll_vec(1)=dnorm(LLF_ct_resid,0.05);
	nll_vec(2)=dnorm(LLF_ct_resid,0.05);
	nll_vec(3)=dnorm(OT_ct_resid,0.05);
	nll_vec(4)=dnorm(RV_resid,tau);
	nll_vec(5)=dnorm(HS_resid,tau);
	nll_vec(6)=dnorm(rt_resid,sig);
	nll_vec(7)=dmvlogistic(RVcatlen,RV_plhat,tau2);
	nll_vec(8)=dmvlogistic(HSMcatlen,HSM_plhat,tau2);
	nll_vec(9)=dmvlogistic(HSFcatlen,HSF_plhat,tau2);
	nll_vec(10)=dmvlogistic(OTcatlen,OT_plhat,tau2);
	nll_vec(11)=dmvlogistic(LLMcatlen,LLM_plhat,tau2);
	nll_vec(12)=dmvlogistic(LLFcatlen,LLF_plhat,tau2);
	
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
	
	
	ofv=sum(nll_vec);	//+sum(p_vec);

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
//	random_number_generator rng(rseed);
//	dmatrix C(syr,eyr,1,nage);
//	dvector tmp(syr-nage,eyr);
//	dvector eps(syr,eyr);
//	dvector simqt(syr,eyr);
//	tmp.fill_randn(rng);
//	eps.fill_randn(rng);
//	wt=tmp*simsig;
//	eps*=simtau;
//	log_ro=log(simro);
//	log_cr=log(simcr);
//	log_rbar=log(simrbar);
//	ahat=simahat;
//	ghat=simghat;
//	ro=mfexp(log_ro);
//	cr=mfexp(log_cr);
//	rbar=mfexp(log_rbar);
//	dvector lxo=pow(mfexp(-m),age-1.);
//	lxo(nage)/=(1.-mfexp(-m));
//	double phieo=lxo*fa;
//	double phibo=lxo*wa;
//	double Bo=value(ro)*phibo;
//	so=cr/phieo;
//	beta=(cr-1.)/(ro*phieo);
//	//selectivity
//	va=plogis(age,ahat,ghat);
//	//Make some fish
//	//initial numbers
//	Nt(syr,1)=mfexp(log_rbar+wt(syr-1));
//	for(int j=2;j<=nage;j++)Nt(syr,j)=mfexp(log_rbar+wt(syr-j))*lxo(j);
//	ft=simF;
//
//	//p(l|a)
//	dvector sdl=la*simcvgrow;
//	//cout<<sdl<<endl;
//	//ad_exit(1);
//	dmatrix pla(1,nlbin,1,nage);
//	dvector lbins(1,nlbin);
//	lbins.fill_seqadd(minbin,stepbin*2.);
//	double z1;
//	double z2;
//	pla.initialize();
//	for(int i=1;i<=nage;i++) //loop over ages
//	{
//		 for(int j=1;j<=nlbin;j++) //loop over length bins
//		{
//			z1=((lbins(j)-stepbin)-la(i))/sdl(i);
//			z2=((lbins(j)+stepbin)-la(i))/sdl(i);
//			pla(j,i)=cumd_norm(z2)-cumd_norm(z1);
//		}//end nbins
//	}//end nage
//	//cout<<pla<<endl;
//	//ad_exit(1);
//	for(int i=syr;i<=eyr;i++)
//	{
//		dvector ba=value(elem_prod(Nt(i),wa));
//		Zt(i)=m+ft(i)*va;
//		//update numbers
//		double sbt=value(Nt(i)*fa);
//		simqt(i)=simqo/(simao+(1.-simao)*sum(ba)/Bo);
//		Nt(i+1,1)=so*sbt/(1.+beta*sbt)*mfexp(wt(i));
//		Nt(i+1)(2,nage)=++elem_prod(Nt(i)(1,nage-1),mfexp(-Zt(i)(1,nage-1)));
//		Nt(i+1,nage)+=Nt(i,nage)*mfexp(-Zt(i,nage));
//		dvector zttmp=value(Zt(i));
//		C(i)=elem_prod(elem_div(value(ft(i)*va),zttmp),elem_prod(1.-mfexp(-zttmp),value(Nt(i))));
//		// get proportions at age in the catch
//		dvector pa=C(i)/sum(C(i));
//		// probability of being a length given any age
//		dvector pl=pla*pa;
//		//cout<<pl<<endl;
//		//ad_exit(1);
//		pl/=sum(pl);
//		plt(i)=rmvlogistic(pl,0.3,rseed+i);
//	}
//	ct=C*wa;
//	bt=Nt*elem_prod(wa,va);
//	yt=elem_prod(simqt,value(elem_prod(bt(syr,eyr),mfexp(eps))));
//	wt=0;
//	//cout<<yt<<endl;
//	//cout<<""<<endl;
//	//cout<<ct<<endl;
//	//cout<<""<<endl;
//	//cout<<plt<<endl;
//	//ad_exit(1);
//
//
//
FUNCTION void calc_partials(const double& fe, double& phie, double& phif, double& phiq, double& dphif_df, double& dphiq_df, double& dRe_df)
//	//Use this function to calculate the partial derivatives (Table 2 in Martell et al. 2008)
//	//Arguments: fe=fishing rate
//	int i;
//	dvector lx=(pow(exp(-m),age-1.));
//	lx(nage)/=(1.-exp(-(m)));
//	dvector lz(1,nage);
//	dvector za=value(m+fe*va);
//	dvector sa=1.-exp(-za);
//	dvector qa=elem_prod(elem_div(value(va),za),sa);
//	double dlz_df=0;
//
//	lz[1]=1.0; 
//	dphiq_df=0; dphif_df=0;
//	phie=(sum(elem_prod(lx,fa)));
//	for(i=1; i<=nage; i++)
//	{
//		if(i>1) lz[i]=lz[i-1]*exp(-za[i-1]);
//		if(i>1) dlz_df=dlz_df*exp(-za[i-1]) - lz[i-1]*value(va[i-1])*exp(-za[i-1]);
//		if(i==nage){ //6/11/2007 added plus group.
//			lz[i]/=(1.-mfexp(-za[i]));
//			//dlz_df=dlz_df*mfexp(-za[i-1]) - lz[i-1]*va[i-1]*mfexp(-za[i-1])/(1.-mfexp(-za[i]))
//			dlz_df=value(dlz_df/(1.-mfexp(-za[i]))
//					-lz[i-1]*mfexp(-za[i-1])*va[i]*mfexp(-za[i])
//			/((1.-mfexp(-za[i]))*(1.-mfexp(-za[i]))));
//		}	
//		dphif_df=dphif_df+(fa[i])*dlz_df;
//		dphiq_df=dphiq_df+(wa[i]*qa[i]*dlz_df+(lz[i]*wa[i]*value(va[i]*va[i]))/za[i]*(exp(-za[i])-sa[i]/za[i]));
//	}
//	phif=sum(elem_prod(lz,(fa)));
//	phiq=sum(elem_prod(elem_prod(lz,(wa)),qa));
//	dRe_df=value(ro/(cr-1.))*phie/square(phif)*dphif_df;
//	//dphif_df=sum(elem_prod(fa,dlz_df));
//	//dvector t2=elem_div(elem_prod(elem_prod(lz,value(va)),wa),za);
//	//dvector t3=exp(-za)-elem_div(sa,za);
//	//dphiq_df=sum(elem_prod(elem_prod(wa,qa),dlz_df)+elem_prod(t2,t3));
//
//
//
FUNCTION void get_CF(double& fe, double& msy,double& bmsy)
//	//This function uses Newton-Raphson method to iteratively solve for F*
//	//Then calculates C* given F* (See eq 1.3 in Martell 2008 CJFAS)
//	int iter;
//	double dy,ddy,re;
//	double phie,phif,phiq,dphif_df,dphiq_df,dRe_df;
//	fe=(m);  //initial guess for Fmsy
//	for(iter= 1; iter<=50; iter++)
//	{
//		calc_partials(fe,phie,phif,phiq,dphif_df,dphiq_df,dRe_df);
//		re=value(ro*(cr-phie/phif)/(cr-1.));
//		dy=re*phiq+fe*phiq*dRe_df+fe*re*dphiq_df;
//		ddy=phiq*dRe_df+re*dphiq_df;
//		//Newton update
//		fe=fe-dy/ddy;
//		if(sfabs(dy)<1.e-10)break;
//		//cout<<"Fe dy\t"<<fe<<" "<<dy<<" "<<fe-dy/ddy<<endl;
//	}
//	msy=fe*re*phiq;
//	bmsy=re*phif;
//	//cout<<"Fe "<<fe<<endl;
//
//
//
REPORT_SECTION
	
	
	REPORT(RVq)
	REPORT(RV_pred)
	REPORT(RVindex)
	REPORT(HSq)
	REPORT(HS_pred)
	REPORT(HSindex)

	REPORT(LLF_ct_hat)
	REPORT(LLctF)
	REPORT(LLM_ct_hat)
	REPORT(LLctM)
	REPORT(OT_ct_hat)
	REPORT(OTct)
	
	REPORT(rbar)
	REPORT(rt)
	


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

