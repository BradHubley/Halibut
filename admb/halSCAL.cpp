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
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <halscal.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
		sim=0;
		rseed=0;
		int on,opt;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			sim=1;
			rseed=atoi(ad_comm::argv[on+1]);
		}
  syr.allocate("syr");
  eyr.allocate("eyr");
  nage.allocate("nage");
  nsexes.allocate("nsexes");
  nlbin.allocate("nlbin");
  nlbin2.allocate("nlbin2");
  minbin.allocate("minbin");
  stepbin.allocate("stepbin");
  RVindex.allocate(syr,eyr,"RVindex");
  RVcatlen.allocate(syr,eyr,1,nlbin,"RVcatlen");
  HSsyr.allocate("HSsyr");
  HSindex.allocate(HSsyr,eyr,"HSindex");
  HSMcatlen.allocate(HSsyr,eyr,1,nlbin,"HSMcatlen");
  HSFcatlen.allocate(HSsyr,eyr,1,nlbin,"HSFcatlen");
  LLctM.allocate(syr,eyr,"LLctM");
  LLctF.allocate(syr,eyr,"LLctF");
  OTct.allocate(syr,eyr,"OTct");
  LLsyr.allocate("LLsyr");
  LLMcatlen1.allocate(LLsyr,eyr,1,nlbin2,"LLMcatlen1");
  LLFcatlen1.allocate(LLsyr,eyr,1,nlbin2,"LLFcatlen1");
  OTsyr.allocate("OTsyr");
  OTcatlen1.allocate(OTsyr,eyr,1,nlbin2,"OTcatlen1");
  LLMcatlen.allocate(LLsyr,eyr,1,nlbin);
 LLMcatlen=trans(trans(LLMcatlen1).sub(1,nlbin));
  LLFcatlen.allocate(LLsyr,eyr,1,nlbin);
 LLFcatlen=trans(trans(LLFcatlen1).sub(1,nlbin));
  OTcatlen.allocate(OTsyr,eyr,1,nlbin);
 OTcatlen=trans(trans(OTcatlen1).sub(1,nlbin));
  lwa.allocate(1,nsexes,"lwa");
  lwb.allocate(1,nsexes,"lwb");
  linf.allocate(1,nsexes,"linf");
  vbk.allocate(1,nsexes,"vbk");
  t0.allocate(1,nsexes,"t0");
  laaM.allocate(1,nage,"laaM");
  laaF.allocate(1,nage,"laaF");
  laaM_sigma.allocate(1,nage,"laaM_sigma");
  laaF_sigma.allocate(1,nage,"laaF_sigma");
  mata.allocate(1,nsexes,"mata");
  matb.allocate(1,nsexes,"matb");
  iro.allocate("iro");
  icr.allocate("icr");
  irbar.allocate("irbar");
  ifbar.allocate("ifbar");
  iahat.allocate("iahat");
  ighat.allocate("ighat");
  ifull.allocate("ifull");
  iSDR.allocate("iSDR");
  iSDL.allocate("iSDL");
  eof.allocate("eof");
iter=0;
  age.allocate(1,nage);
age.fill_seqadd(1,1);
  iyr.allocate(syr,eyr);
iyr.fill_seqadd(syr,1);
mM=1.2*vbk(1);
mF=1.2*vbk(2);
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}else{
			cout<<"Reading data successful."<<endl;
		}
  ma.allocate(1,nage);
  waM.allocate(1,nage);
  waF.allocate(1,nage);
  fa.allocate(1,nage);
		for (int a=1;a<=nage;a++){
			ma(a)=1/(1+exp(-matb(2)*((linf(2)-(linf(2)-t0(2))*exp(-vbk(2)*a))-mata(2))));
		}
		waM=lwa(1)*pow(laaF,lwb(1))/1000; 	//convert to kg
		waF=lwa(2)*pow(laaM,lwb(2))/1000; 
		fa=elem_prod(waF,ma);
  pRVcatlen.allocate(syr,eyr,1,nlbin);
  pHSMcatlen.allocate(HSsyr,eyr,1,nlbin);
  pHSFcatlen.allocate(HSsyr,eyr,1,nlbin);
  pOTcatlen.allocate(OTsyr,eyr,1,nlbin);
  pLLMcatlen.allocate(LLsyr,eyr,1,nlbin);
  pLLFcatlen.allocate(LLsyr,eyr,1,nlbin);
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
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_ro.allocate(4,"log_ro");
  log_cr.allocate(5,"log_cr");
  log_rbar.allocate("log_rbar");
  wt.allocate(syr-nage,eyr,-10.,10.,2,"wt");
  log_LLF_fbar.allocate("log_LLF_fbar");
  log_LLM_fbar.allocate("log_LLM_fbar");
  log_OT_fbar.allocate("log_OT_fbar");
  log_LLF_ft_dev.allocate(syr,eyr,-10.,10.,3,"log_LLF_ft_dev");
  log_LLM_ft_dev.allocate(syr,eyr,-10.,10.,3,"log_LLM_ft_dev");
  log_OT_ft_dev.allocate(syr,eyr,-10.,10.,3,"log_OT_ft_dev");
  HSsel_half.allocate(0,nage,4,"HSsel_half");
  HSsel_SD.allocate(0,5,4,"HSsel_SD");
  RVsel_full.allocate(0,nage,4,"RVsel_full");
  RVsel_SDR.allocate(0,5,4,"RVsel_SDR");
  RVsel_SDL.allocate(0,5,4,"RVsel_SDL");
  LLsel_half.allocate(0,nage,4,"LLsel_half");
  LLsel_SD.allocate(0,5,4,"LLsel_SD");
  OTsel_full.allocate(0,nage,4,"OTsel_full");
  OTsel_SDR.allocate(0,5,4,"OTsel_SDR");
  OTsel_SDL.allocate(0,5,4,"OTsel_SDL");
  rho.allocate(0,1,"rho");
  cvgrow.allocate(0,1,"cvgrow");
  ptdev.allocate("ptdev");
  ofv.allocate("ofv");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  tdev.allocate("tdev");
log_ro=log(iro);
log_cr=log(icr);
log_rbar=log(irbar);
log_LLF_fbar=log(ifbar);
log_LLM_fbar=log(ifbar);
log_OT_fbar=log(ifbar);
rho=0.5;
cvgrow=0.1;
ptdev=log(1./0.08);
HSsel_half=iahat;
HSsel_SD=ighat;
RVsel_full=ifull;
RVsel_SDR=iSDR;
RVsel_SDL=iSDL;
LLsel_half=iahat;
LLsel_SD=ighat;
OTsel_full=ifull;
OTsel_SDR=iSDR;
OTsel_SDL=iSDL;
  ro.allocate("ro");
  #ifndef NO_AD_INITIALIZE
  ro.initialize();
  #endif
  cr.allocate("cr");
  #ifndef NO_AD_INITIALIZE
  cr.initialize();
  #endif
  rbar.allocate("rbar");
  #ifndef NO_AD_INITIALIZE
  rbar.initialize();
  #endif
  fbar.allocate("fbar");
  #ifndef NO_AD_INITIALIZE
  fbar.initialize();
  #endif
  so.allocate("so");
  #ifndef NO_AD_INITIALIZE
  so.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  RVq.allocate("RVq");
  #ifndef NO_AD_INITIALIZE
  RVq.initialize();
  #endif
  HSq.allocate("HSq");
  #ifndef NO_AD_INITIALIZE
  HSq.initialize();
  #endif
  fmsy.allocate("fmsy");
  #ifndef NO_AD_INITIALIZE
  fmsy.initialize();
  #endif
  msy.allocate("msy");
  #ifndef NO_AD_INITIALIZE
  msy.initialize();
  #endif
  bmsy.allocate("bmsy");
  #ifndef NO_AD_INITIALIZE
  bmsy.initialize();
  #endif
  sig.allocate("sig");
  #ifndef NO_AD_INITIALIZE
  sig.initialize();
  #endif
  tau.allocate("tau");
  #ifndef NO_AD_INITIALIZE
  tau.initialize();
  #endif
  HSM_va.allocate(1,nage,"HSM_va");
  #ifndef NO_AD_INITIALIZE
    HSM_va.initialize();
  #endif
  HSF_va.allocate(1,nage,"HSF_va");
  #ifndef NO_AD_INITIALIZE
    HSF_va.initialize();
  #endif
  RV_va.allocate(1,nage,"RV_va");
  #ifndef NO_AD_INITIALIZE
    RV_va.initialize();
  #endif
  LLM_va.allocate(1,nage,"LLM_va");
  #ifndef NO_AD_INITIALIZE
    LLM_va.initialize();
  #endif
  LLF_va.allocate(1,nage,"LLF_va");
  #ifndef NO_AD_INITIALIZE
    LLF_va.initialize();
  #endif
  OT_va.allocate(1,nage,"OT_va");
  #ifndef NO_AD_INITIALIZE
    OT_va.initialize();
  #endif
  LLM_ft.allocate(syr,eyr,"LLM_ft");
  #ifndef NO_AD_INITIALIZE
    LLM_ft.initialize();
  #endif
  LLF_ft.allocate(syr,eyr,"LLF_ft");
  #ifndef NO_AD_INITIALIZE
    LLF_ft.initialize();
  #endif
  OT_ft.allocate(syr,eyr,"OT_ft");
  #ifndef NO_AD_INITIALIZE
    OT_ft.initialize();
  #endif
  M_bt.allocate(syr,eyr+1,"M_bt");
  #ifndef NO_AD_INITIALIZE
    M_bt.initialize();
  #endif
  F_bt.allocate(syr,eyr+1,"F_bt");
  #ifndef NO_AD_INITIALIZE
    F_bt.initialize();
  #endif
  LLF_ct_hat.allocate(syr,eyr,"LLF_ct_hat");
  #ifndef NO_AD_INITIALIZE
    LLF_ct_hat.initialize();
  #endif
  LLM_ct_hat.allocate(syr,eyr,"LLM_ct_hat");
  #ifndef NO_AD_INITIALIZE
    LLM_ct_hat.initialize();
  #endif
  OT_ct_hat.allocate(syr,eyr,"OT_ct_hat");
  #ifndef NO_AD_INITIALIZE
    OT_ct_hat.initialize();
  #endif
  LLF_ct_resid.allocate(syr,eyr,"LLF_ct_resid");
  #ifndef NO_AD_INITIALIZE
    LLF_ct_resid.initialize();
  #endif
  LLM_ct_resid.allocate(syr,eyr,"LLM_ct_resid");
  #ifndef NO_AD_INITIALIZE
    LLM_ct_resid.initialize();
  #endif
  OT_ct_resid.allocate(syr,eyr,"OT_ct_resid");
  #ifndef NO_AD_INITIALIZE
    OT_ct_resid.initialize();
  #endif
  RV_resid.allocate(syr,eyr,"RV_resid");
  #ifndef NO_AD_INITIALIZE
    RV_resid.initialize();
  #endif
  HS_resid.allocate(HSsyr,eyr,"HS_resid");
  #ifndef NO_AD_INITIALIZE
    HS_resid.initialize();
  #endif
  rt.allocate(syr+1,eyr,"rt");
  #ifndef NO_AD_INITIALIZE
    rt.initialize();
  #endif
  rt_resid.allocate(syr+1,eyr,"rt_resid");
  #ifndef NO_AD_INITIALIZE
    rt_resid.initialize();
  #endif
  HS_pred.allocate(syr,eyr+1,"HS_pred");
  #ifndef NO_AD_INITIALIZE
    HS_pred.initialize();
  #endif
  RV_pred.allocate(syr,eyr+1,"RV_pred");
  #ifndef NO_AD_INITIALIZE
    RV_pred.initialize();
  #endif
  lbins.allocate(1,nlbin,"lbins");
  #ifndef NO_AD_INITIALIZE
    lbins.initialize();
  #endif
  M_Nat.allocate(syr,eyr+1,1,nage,"M_Nat");
  #ifndef NO_AD_INITIALIZE
    M_Nat.initialize();
  #endif
  F_Nat.allocate(syr,eyr+1,1,nage,"F_Nat");
  #ifndef NO_AD_INITIALIZE
    F_Nat.initialize();
  #endif
  M_Fat.allocate(syr,eyr,1,nage,"M_Fat");
  #ifndef NO_AD_INITIALIZE
    M_Fat.initialize();
  #endif
  F_Fat.allocate(syr,eyr,1,nage,"F_Fat");
  #ifndef NO_AD_INITIALIZE
    F_Fat.initialize();
  #endif
  OT_Fat.allocate(syr,eyr,1,nage,"OT_Fat");
  #ifndef NO_AD_INITIALIZE
    OT_Fat.initialize();
  #endif
  M_Zat.allocate(syr,eyr,1,nage,"M_Zat");
  #ifndef NO_AD_INITIALIZE
    M_Zat.initialize();
  #endif
  F_Zat.allocate(syr,eyr,1,nage,"F_Zat");
  #ifndef NO_AD_INITIALIZE
    F_Zat.initialize();
  #endif
  F_alk.allocate(1,nlbin,1,nage,"F_alk");
  #ifndef NO_AD_INITIALIZE
    F_alk.initialize();
  #endif
  M_alk.allocate(1,nlbin,1,nage,"M_alk");
  #ifndef NO_AD_INITIALIZE
    M_alk.initialize();
  #endif
  RV_plhat.allocate(syr,eyr,1,nlbin,"RV_plhat");
  #ifndef NO_AD_INITIALIZE
    RV_plhat.initialize();
  #endif
  HSM_plhat.allocate(syr,eyr,1,nlbin,"HSM_plhat");
  #ifndef NO_AD_INITIALIZE
    HSM_plhat.initialize();
  #endif
  HSF_plhat.allocate(syr,eyr,1,nlbin,"HSF_plhat");
  #ifndef NO_AD_INITIALIZE
    HSF_plhat.initialize();
  #endif
  OT_plhat.allocate(syr,eyr,1,nlbin,"OT_plhat");
  #ifndef NO_AD_INITIALIZE
    OT_plhat.initialize();
  #endif
  LLM_plhat.allocate(syr,eyr,1,nlbin,"LLM_plhat");
  #ifndef NO_AD_INITIALIZE
    LLM_plhat.initialize();
  #endif
  LLF_plhat.allocate(syr,eyr,1,nlbin,"LLF_plhat");
  #ifndef NO_AD_INITIALIZE
    LLF_plhat.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
}

void model_parameters::userfunction(void)
{
  ofv =0.0;
	initialization();
	statedynamics();
	observation_model();
	stock_recruit_model();
	objective_function();
}

void model_parameters::initialization(void)
{
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
	OT_Fat=outer_prod(OT_ft,OT_va)*0.5;
	// total mortality
	M_Zat=mM+M_Fat+OT_Fat;
	F_Zat=mF+F_Fat+OT_Fat;
}

void model_parameters::statedynamics(void)
{
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
}

void model_parameters::observation_model(void)
{
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
}

void model_parameters::stock_recruit_model(void)
{
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
}

void model_parameters::objective_function(void)
{
	// Negative Log Likelihoods
	dvar_vector nll_vec(1,12);
	tdev=sqrt(1./mfexp(ptdev));
	sig=sqrt(rho*1./mfexp(ptdev));//process error sd.dev
	tau=sqrt((1.-rho)*1./mfexp(ptdev));
	double tau2;
	dvar_matrix nu(syr,eyr,1,nlbin);
	double pmin=0.0025;
	nll_vec.initialize();
	nll_vec(1)=dnorm(LLF_ct_resid,0.05);
	nll_vec(2)=dnorm(LLF_ct_resid,0.05);
	nll_vec(3)=dnorm(OT_ct_resid,0.05);
	nll_vec(4)=dnorm(RV_resid,tau);
	nll_vec(5)=dnorm(HS_resid,tau);
	nll_vec(6)=dnorm(rt_resid,sig);
	nll_vec(7)=dmvlogistic(RVcatlen,RV_plhat,nu,tau2,pmin);
	nll_vec(8)=dmvlogistic(HSMcatlen,HSM_plhat,nu,tau2,pmin);
	nll_vec(9)=dmvlogistic(HSFcatlen,HSF_plhat,nu,tau2,pmin);
	nll_vec(10)=dmvlogistic(OTcatlen,OT_plhat,nu,tau2,pmin);
	nll_vec(11)=dmvlogistic(LLMcatlen,LLM_plhat,nu,tau2,pmin);
	nll_vec(12)=dmvlogistic(LLFcatlen,LLF_plhat,nu,tau2,pmin);
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
}

dvar_vector model_parameters::sel_dhn(const int& minage, const int& maxage, dvariable& full, dvariable& varR, dvariable& varL)
{
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
}

void model_parameters::mcmc_output(void)
{
}

void model_parameters::forecast(void)
{
}

void model_parameters::run_data_simulation(void)
{
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
	REPORT(syr)
	REPORT(HSsyr)
	REPORT(eyr)
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
	REPORT(LLM_ft)
	REPORT(LLF_ft)
	REPORT(OT_ft)
	REPORT(rbar)
	REPORT(rt)
	REPORT(tau)
	REPORT(sig)
	//REPORT(tau2)
	REPORT(RV_va)
	REPORT(HSM_va)
	REPORT(HSF_va)
	REPORT(OT_va)
	REPORT(LLM_va)
	REPORT(LLF_va)
	REPORT(pRVcatlen)
	REPORT(pHSMcatlen)
	REPORT(pHSFcatlen)
	REPORT(pOTcatlen)
	REPORT(pLLMcatlen)
	REPORT(pLLFcatlen)
	//double fmsy;
	//double msy;
	//double bmsy;
	//if(last_phase())
	//{
	//	get_CF(fmsy,msy,bmsy);
	//}
}

void model_parameters::final_calcs()
{
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
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
