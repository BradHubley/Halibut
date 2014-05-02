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
#include <halSCAL.htp>

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
  minbin.allocate("minbin");
  stepbin.allocate("stepbin");
  wa.allocate(1,nsexes,"wa");
  wb.allocate(1,nsexes,"wb");
  linf.allocate(1,nsexes,"linf");
  vbk.allocate(1,nsexes,"vbk");
  t0.allocate(1,nsexes,"t0");
  laaF.allocate(1,nage,"laaF");
  laaM.allocate(1,nage,"laaM");
  laaF_sigma.allocate(1,nage,"laaF_sigma");
  laaM_sigma.allocate(1,nage,"laaM_sigma");
  RVindex.allocate(syr,eyr,"RVindex");
  RVcatlen.allocate(syr,eyr,1,nlbin,"RVcatlen");
  HSsyr.allocate("HSsyr");
  HSindex.allocate(HSsyr,eyr,"HSindex");
  HScatlenF.allocate(HSsyr,eyr,1,nlbin,"HScatlenF");
  HScatlenM.allocate(HSsyr,eyr,1,nlbin,"HScatlenM");
  LLctF.allocate(syr,eyr,"LLctF");
  LLctM.allocate(syr,eyr,"LLctM");
  OTct.allocate(syr,eyr,"OTct");
  LLsyr.allocate("LLsyr");
  LLcatlenF.allocate(LLsyr,eyr,1,nlbin,"LLcatlenF");
  LLcatlenM.allocate(LLsyr,eyr,1,nlbin,"LLcatlenM");
  OTsyr.allocate("OTsyr");
  OTcatlen.allocate(OTsyr,eyr,1,nlbin,"OTcatlen");
  iro.allocate("iro");
  icr.allocate("icr");
  ct.allocate(syr,eyr,"ct");
  irbar.allocate("irbar");
  ifbar.allocate("ifbar");
  iahat.allocate("iahat");
  ighat.allocate("ighat");
  eof.allocate("eof");
iter=0;
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}else{
			cout<<"Reading data successful."<<endl;
			ad_exit(1);
		}
  age.allocate(1,nage);
  la.allocate(1,nage);
  wa.allocate(1,nage);
  fa.allocate(1,nage);
		m=1.2*vbk;
		age.fill_seqadd(1,1);
		la=linf*(1.-mfexp(-vbk*age));
		wa=wla*pow(la,wlb);
		//weight at age * maturity at age
		fa=elem_prod(wa,plogis(age,log(3)/vbk,0.1*log(3)/vbk));
ad_comm::change_datafile_name("SCAL.ctl");
  simsig.allocate("simsig");
  simtau.allocate("simtau");
  simqo.allocate("simqo");
  simao.allocate("simao");
  simro.allocate("simro");
  simcr.allocate("simcr");
  simrbar.allocate("simrbar");
  simahat.allocate("simahat");
  simghat.allocate("simghat");
  simcvgrow.allocate("simcvgrow");
  simF.allocate(syr,eyr,"simF");
  eofc.allocate("eofc");
		if(eofc!=999)
		{
			cout<<"Error reading control file.\n Fix it."<<endl;
			ad_exit(1);
		}
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_ro.allocate(4,"log_ro");
  log_cr.allocate(5,"log_cr");
  log_rbar.allocate("log_rbar");
  log_fbar.allocate("log_fbar");
  ahat.allocate(0,nage,"ahat");
  ghat.allocate(0,5,"ghat");
  ptdev.allocate("ptdev");
  rho.allocate(0,1,"rho");
  cvgrow.allocate(0,1,"cvgrow");
  wt.allocate(syr-nage,eyr,-10.,10.,2,"wt");
  log_ft_dev.allocate(syr,eyr,-10.,10.,3,"log_ft_dev");
  nll.allocate("nll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
log_ro=log(iro);
log_cr=log(icr);
log_rbar=log(irbar);
log_fbar=log(ifbar);
rho=0.5;
cvgrow=0.2;
ptdev=log(1./0.08);
ahat=iahat;
ghat=ighat;
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
  q.allocate("q");
  #ifndef NO_AD_INITIALIZE
  q.initialize();
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
  tdev.allocate("tdev");
  va.allocate(1,nage,"va");
  #ifndef NO_AD_INITIALIZE
    va.initialize();
  #endif
  ft.allocate(syr,eyr,"ft");
  #ifndef NO_AD_INITIALIZE
    ft.initialize();
  #endif
  bt.allocate(syr,eyr+1,"bt");
  #ifndef NO_AD_INITIALIZE
    bt.initialize();
  #endif
  ct_hat.allocate(syr,eyr,"ct_hat");
  #ifndef NO_AD_INITIALIZE
    ct_hat.initialize();
  #endif
  ct_resid.allocate(syr,eyr,"ct_resid");
  #ifndef NO_AD_INITIALIZE
    ct_resid.initialize();
  #endif
  yt_resid.allocate(syr,eyr,"yt_resid");
  #ifndef NO_AD_INITIALIZE
    yt_resid.initialize();
  #endif
  rt.allocate(syr+1,eyr,"rt");
  #ifndef NO_AD_INITIALIZE
    rt.initialize();
  #endif
  rt_resid.allocate(syr+1,eyr,"rt_resid");
  #ifndef NO_AD_INITIALIZE
    rt_resid.initialize();
  #endif
  Nt.allocate(syr,eyr+1,1,nage,"Nt");
  #ifndef NO_AD_INITIALIZE
    Nt.initialize();
  #endif
  Zt.allocate(syr,eyr,1,nage,"Zt");
  #ifndef NO_AD_INITIALIZE
    Zt.initialize();
  #endif
  plhat.allocate(syr,eyr,1,nlbin,"plhat");
  #ifndef NO_AD_INITIALIZE
    plhat.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
  if(sim)
  {
  	run_data_simulation();
  }
}

void model_parameters::userfunction(void)
{
  nll =0.0;
	initialization();
	statedynamics();
	observation_model();
	stock_recruit_model();
	objective_function();
	if(mceval_phase())
	{ 
		forecast();
		get_CF(value(fmsy),value(msy),value(bmsy));
		mcmc_output();
	}
	if(last_phase())
	{
		forecast();
		get_CF(value(fmsy),value(msy),value(bmsy));
	} 
}

void model_parameters::initialization(void)
{
	va=plogis(age,ahat,ghat);
	ft=mfexp(log_fbar+log_ft_dev);
	Zt=m+outer_prod(ft,va);
}

void model_parameters::statedynamics(void)
{
	dvar_vector lxo=pow(mfexp(-m),age-1.);
	lxo(nage)/=(1.-mfexp(-m));
	Nt(syr,1)=mfexp(log_rbar+wt(syr-1));
	for(int j=2;j<=nage;j++) Nt(syr,j)=mfexp(log_rbar+wt(syr-j))*lxo(j);
	for(int i=syr;i<=eyr;i++)
	{
		Nt(i+1,1)=mfexp(log_rbar+wt(i));
		Nt(i+1)(2,nage)=++elem_prod(Nt(i)(1,nage-1),mfexp(-Zt(i)(1,nage-1)));
		Nt(i+1,nage)+=Nt(i,nage)*mfexp(-Zt(i,nage));
	}
	bt=Nt*elem_prod(wa,va);
}

void model_parameters::observation_model(void)
{
	dvar_matrix C(syr,eyr,1,nage);
	dvar_matrix F(syr,eyr,1,nage);
	F=outer_prod(ft,va);
	C=elem_prod(elem_div(F,Zt),elem_prod(1.-mfexp(-Zt),Nt));
	ct_hat=C*wa;
	//catch residuals
	ct_resid=log(ct)-log(ct_hat);
	//cpue residuals ala waters and ludwig 1994
	yt_resid=log(yt)-log(bt(syr,eyr));
	q=mfexp(mean(yt_resid));
	yt_resid-=mean(yt_resid);
	//predicted proportions in the catch
		//p(l|a)
	dvar_vector sdl=la*cvgrow;
	//cout<<sdl<<endl;
	//ad_exit(1);
	dvar_matrix pla(1,nlbin,1,nage);
	dvar_vector lbins(1,nlbin);
	lbins.fill_seqadd(minbin,stepbin*2.);
	dvariable z1;
	dvariable z2;
	pla.initialize();
	for(int i=1;i<=nage;i++) //loop over ages
	{
		 for(int j=1;j<=nlbin;j++) //loop over length bins
		{
			z1=((lbins(j)-stepbin)-la(i))/sdl(i);
			z2=((lbins(j)+stepbin)-la(i))/sdl(i);
			pla(j,i)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
	}//end nage
	for(int i=syr;i<=eyr;i++)
	{
		dvar_vector pa=C(i)/sum(C(i));
		dvar_vector pl=pla*pa;
		plhat(i)=pl/sum(pl);
	}
	//cout<<plhat(syr)<<endl;
	//ad_exit(1);
}

void model_parameters::stock_recruit_model(void)
{
	ro=mfexp(log_ro);
	cr=mfexp(log_cr)+1.;
	dvar_vector lxo=pow(mfexp(-m),age-1.);
	lxo(nage)/=(1.-mfexp(-m));
	dvariable phieo=lxo*fa;
	so=cr/phieo;
	beta=(cr-1.)/(ro*phieo);
	dvar_vector sbt=Nt*fa;
	dvar_vector nmr=so*sbt(syr,eyr-1);
	dvar_vector den=(1.+beta*sbt(syr,eyr-1));
	dvar_vector tmp_rt=++elem_div(nmr,den);
	rt=column(Nt.sub(syr+1,eyr),1);
	rt_resid=log(tmp_rt)-log(rt);
}

void model_parameters::objective_function(void)
{
	dvar_vector nll_vec(1,4);
	tdev=sqrt(1./mfexp(ptdev));
	sig=sqrt(rho*1./mfexp(ptdev));//process error sd.dev
	tau=sqrt((1.-rho)*1./mfexp(ptdev));
	nll_vec.initialize();
	nll_vec(1)=dnorm(ct_resid,0.05);
	nll_vec(2)=dnorm(yt_resid,tau);
	nll_vec(3)=dnorm(rt_resid,sig);
	double tau2;
	nll_vec(4)=dmvlogistic(plt,plhat,tau2);
	dvar_vector p_vec(1,5);
	p_vec.initialize();
	dvariable h=cr/(4.+cr);
	p_vec(1)=dbeta((h-0.2)/0.8,2,2);
	if(last_phase())
	{
		p_vec(3)=dnorm(wt,2);
		p_vec(4)=dnorm(log_ft_dev,2);
	}
	else
	{
		p_vec(3)=100.*norm2(wt);
		p_vec(4)=100.*norm2(log_ft_dev);
	}
	p_vec(5)=dbeta(rho,50,50);
	nll=sum(nll_vec)+sum(p_vec);
}

void model_parameters::mcmc_output(void)
{
	if(iter==0)
	{
		ofstream ofs("refpar.mcmc");
		ofs<<"fmsy\t bmsy\t msy\t b/bmsy\t f/fmsy"<<endl;
		//ofs<<"r\t k\t q\t sig\t"<<endl;
	}
	iter++;
	double fratio=value(ft(eyr)/fmsy);
	double bratio=value(Nt(eyr)*wa/bmsy);
	ofstream ofs("refpar.mcmc",ios::app);
	ofs<<fmsy<<"\t"<<bmsy<<"\t"<<msy<<"\t"<<bratio<<"\t"<<fratio<<endl;
}

void model_parameters::forecast(void)
{
}

void model_parameters::run_data_simulation(void)
{
	random_number_generator rng(rseed);
	dmatrix C(syr,eyr,1,nage);
	dvector tmp(syr-nage,eyr);
	dvector eps(syr,eyr);
	dvector simqt(syr,eyr);
	tmp.fill_randn(rng);
	eps.fill_randn(rng);
	wt=tmp*simsig;
	eps*=simtau;
	log_ro=log(simro);
	log_cr=log(simcr);
	log_rbar=log(simrbar);
	ahat=simahat;
	ghat=simghat;
	ro=mfexp(log_ro);
	cr=mfexp(log_cr);
	rbar=mfexp(log_rbar);
	dvector lxo=pow(mfexp(-m),age-1.);
	lxo(nage)/=(1.-mfexp(-m));
	double phieo=lxo*fa;
	double phibo=lxo*wa;
	double Bo=value(ro)*phibo;
	so=cr/phieo;
	beta=(cr-1.)/(ro*phieo);
	//selectivity
	va=plogis(age,ahat,ghat);
	//Make some fish
	//initial numbers
	Nt(syr,1)=mfexp(log_rbar+wt(syr-1));
	for(int j=2;j<=nage;j++)Nt(syr,j)=mfexp(log_rbar+wt(syr-j))*lxo(j);
	ft=simF;
	//p(l|a)
	dvector sdl=la*simcvgrow;
	//cout<<sdl<<endl;
	//ad_exit(1);
	dmatrix pla(1,nlbin,1,nage);
	dvector lbins(1,nlbin);
	lbins.fill_seqadd(minbin,stepbin*2.);
	double z1;
	double z2;
	pla.initialize();
	for(int i=1;i<=nage;i++) //loop over ages
	{
		 for(int j=1;j<=nlbin;j++) //loop over length bins
		{
			z1=((lbins(j)-stepbin)-la(i))/sdl(i);
			z2=((lbins(j)+stepbin)-la(i))/sdl(i);
			pla(j,i)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
	}//end nage
	//cout<<pla<<endl;
	//ad_exit(1);
	for(int i=syr;i<=eyr;i++)
	{
		dvector ba=value(elem_prod(Nt(i),wa));
		Zt(i)=m+ft(i)*va;
		//update numbers
		double sbt=value(Nt(i)*fa);
		simqt(i)=simqo/(simao+(1.-simao)*sum(ba)/Bo);
		Nt(i+1,1)=so*sbt/(1.+beta*sbt)*mfexp(wt(i));
		Nt(i+1)(2,nage)=++elem_prod(Nt(i)(1,nage-1),mfexp(-Zt(i)(1,nage-1)));
		Nt(i+1,nage)+=Nt(i,nage)*mfexp(-Zt(i,nage));
		dvector zttmp=value(Zt(i));
		C(i)=elem_prod(elem_div(value(ft(i)*va),zttmp),elem_prod(1.-mfexp(-zttmp),value(Nt(i))));
		// get proportions at age in the catch
		dvector pa=C(i)/sum(C(i));
		// probability of being a length given any age
		dvector pl=pla*pa;
		//cout<<pl<<endl;
		//ad_exit(1);
		pl/=sum(pl);
		plt(i)=rmvlogistic(pl,0.3,rseed+i);
	}
	ct=C*wa;
	bt=Nt*elem_prod(wa,va);
	yt=elem_prod(simqt,value(elem_prod(bt(syr,eyr),mfexp(eps))));
	wt=0;
	//cout<<yt<<endl;
	//cout<<""<<endl;
	//cout<<ct<<endl;
	//cout<<""<<endl;
	//cout<<plt<<endl;
	//ad_exit(1);
}

void model_parameters::calc_partials(const double& fe, double& phie, double& phif, double& phiq, double& dphif_df, double& dphiq_df, double& dRe_df)
{
	//Use this function to calculate the partial derivatives (Table 2 in Martell et al. 2008)
	//Arguments: fe=fishing rate
	int i;
	dvector lx=(pow(exp(-m),age-1.));
	lx(nage)/=(1.-exp(-(m)));
	dvector lz(1,nage);
	dvector za=value(m+fe*va);
	dvector sa=1.-exp(-za);
	dvector qa=elem_prod(elem_div(value(va),za),sa);
	double dlz_df=0;
	lz[1]=1.0; 
	dphiq_df=0; dphif_df=0;
	phie=(sum(elem_prod(lx,fa)));
	for(i=1; i<=nage; i++)
	{
		if(i>1) lz[i]=lz[i-1]*exp(-za[i-1]);
		if(i>1) dlz_df=dlz_df*exp(-za[i-1]) - lz[i-1]*value(va[i-1])*exp(-za[i-1]);
		if(i==nage){ //6/11/2007 added plus group.
			lz[i]/=(1.-mfexp(-za[i]));
			//dlz_df=dlz_df*mfexp(-za[i-1]) - lz[i-1]*va[i-1]*mfexp(-za[i-1])/(1.-mfexp(-za[i]))
			dlz_df=value(dlz_df/(1.-mfexp(-za[i]))
					-lz[i-1]*mfexp(-za[i-1])*va[i]*mfexp(-za[i])
			/((1.-mfexp(-za[i]))*(1.-mfexp(-za[i]))));
		}	
		dphif_df=dphif_df+(fa[i])*dlz_df;
		dphiq_df=dphiq_df+(wa[i]*qa[i]*dlz_df+(lz[i]*wa[i]*value(va[i]*va[i]))/za[i]*(exp(-za[i])-sa[i]/za[i]));
	}
	phif=sum(elem_prod(lz,(fa)));
	phiq=sum(elem_prod(elem_prod(lz,(wa)),qa));
	dRe_df=value(ro/(cr-1.))*phie/square(phif)*dphif_df;
	//dphif_df=sum(elem_prod(fa,dlz_df));
	//dvector t2=elem_div(elem_prod(elem_prod(lz,value(va)),wa),za);
	//dvector t3=exp(-za)-elem_div(sa,za);
	//dphiq_df=sum(elem_prod(elem_prod(wa,qa),dlz_df)+elem_prod(t2,t3));
}

void model_parameters::get_CF(double& fe, double& msy,double& bmsy)
{
	//This function uses Newton-Raphson method to iteratively solve for F*
	//Then calculates C* given F* (See eq 1.3 in Martell 2008 CJFAS)
	int iter;
	double dy,ddy,re;
	double phie,phif,phiq,dphif_df,dphiq_df,dRe_df;
	fe=(m);  //initial guess for Fmsy
	for(iter= 1; iter<=50; iter++)
	{
		calc_partials(fe,phie,phif,phiq,dphif_df,dphiq_df,dRe_df);
		re=value(ro*(cr-phie/phif)/(cr-1.));
		dy=re*phiq+fe*phiq*dRe_df+fe*re*dphiq_df;
		ddy=phiq*dRe_df+re*dphiq_df;
		//Newton update
		fe=fe-dy/ddy;
		if(sfabs(dy)<1.e-10)break;
		//cout<<"Fe dy\t"<<fe<<" "<<dy<<" "<<fe-dy/ddy<<endl;
	}
	msy=fe*re*phiq;
	bmsy=re*phif;
	//cout<<"Fe "<<fe<<endl;
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
	//double fmsy;
	//double msy;
	//double bmsy;
	//if(last_phase())
	//{
	//	get_CF(fmsy,msy,bmsy);
	//}
	REPORT(ro);
	REPORT(cr);
	REPORT(mfexp(log_rbar));
	REPORT(mfexp(log_rbar+wt));
	REPORT(mfexp(log_fbar));
	REPORT(mfexp(log_fbar+log_ft_dev));
	REPORT(plt);
	REPORT(plhat);
	REPORT(ahat);
	REPORT(ghat);
	REPORT(rho);
	REPORT(cvgrow);
	REPORT(ct);
	REPORT(ct_hat);
	REPORT(fmsy);
	REPORT(msy);
	REPORT(bmsy);
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
