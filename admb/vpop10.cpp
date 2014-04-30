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
#include <vpop10.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_report1 = new ofstream("mcout.dat",ios::app);  // you must delete this file before a new mceval run;
  pad_tmp = new ofstream("tmp.dat",ios::app);  // you must delete this file before a new mceval run;
  debug.allocate("debug");
if (debug==1){char mark; cout<<"starting debug /n"; cin>>mark; }
  first_year.allocate("first_year");
  last_year.allocate("last_year");
  first_age.allocate("first_age");
  last_age.allocate("last_age");
  first_length.allocate("first_length");
  length_inc.allocate("length_inc");
  N_lens.allocate("N_lens");
  N_sexes.allocate("N_sexes");
  Flength.allocate("Flength");
  fminbins.allocate("fminbins");
  fmaxbins.allocate("fmaxbins");
  mminbins.allocate("mminbins");
  mmaxbins.allocate("mmaxbins");
  fminbins2.allocate("fminbins2");
  fmaxbins2.allocate("fmaxbins2");
  mminbins2.allocate("mminbins2");
  mmaxbins2.allocate("mmaxbins2");
  fminbins3.allocate("fminbins3");
  fmaxbins3.allocate("fmaxbins3");
  mminbins3.allocate("mminbins3");
  mmaxbins3.allocate("mmaxbins3");
  Frecyr.allocate("Frecyr");
  Minit.allocate("Minit");
  Mopt.allocate("Mopt");
  N1opt.allocate("N1opt");
  N1_prior.allocate(1,7,"N1_prior");
  year1.allocate(first_age,last_age,"year1");
  start_Z_prior.allocate(1,7,"start_Z_prior");
  age1_prior.allocate(1,7,"age1_prior");
  R0_prior.allocate(1,7,"R0_prior");
  beta_prior.allocate(1,7,"beta_prior");
  alpha_prior.allocate(1,7,"alpha_prior");
  p1_prior.allocate(1,7,"p1_prior");
  p2_prior.allocate(1,7,"p2_prior");
  p3_prior.allocate(1,7,"p3_prior");
  log_RecDev_prior.allocate(1,7,"log_RecDev_prior");
  M_prior.allocate(1,7,"M_prior");
  ll_select_switch.allocate("ll_select_switch");
  old_summer_rv_Sfull_prior.allocate(1,7,"old_summer_rv_Sfull_prior");
  old_summer_rv_varLest_prior.allocate(1,7,"old_summer_rv_varLest_prior");
  old_summer_rv_varRest_prior.allocate(1,7,"old_summer_rv_varRest_prior");
  halibut_fixed_Sfull_prior.allocate(1,N_sexes,1,7,"halibut_fixed_Sfull_prior");
  halibut_fixed_varLest_prior.allocate(1,N_sexes,1,7,"halibut_fixed_varLest_prior");
  ll_Sfull_prior.allocate(1,N_sexes,1,7,"ll_Sfull_prior");
  ll_varLest_prior.allocate(1,N_sexes,1,7,"ll_varLest_prior");
  ot_Sfull_prior.allocate(1,7,"ot_Sfull_prior");
  ot_varLest_prior.allocate(1,7,"ot_varLest_prior");
  ot_varRest_prior.allocate(1,7,"ot_varRest_prior");
  log_q_old_summer_rv_prior.allocate(1,7,"log_q_old_summer_rv_prior");
  log_q_halibut_fixed_prior.allocate(1,7,"log_q_halibut_fixed_prior");
  u_ll_prior.allocate(1,7,"u_ll_prior");
  like_weight.allocate(1,12,"like_weight");
  N_catlen_commercial.allocate("N_catlen_commercial");
  N_catlen_surveys.allocate("N_catlen_surveys");
  catlen_LikeType.allocate(1,N_catlen_surveys,"catlen_LikeType");
  maxssizes.allocate("maxssizes");
  bi.allocate(1,N_sexes,"bi");
  bii.allocate(1,N_sexes,"bii");
  linf.allocate(1,N_sexes,"linf");
  k.allocate(1,N_sexes,"k");
  t0.allocate(1,N_sexes,"t0");
  sigma_a.allocate(1,N_sexes,"sigma_a");
  sigma_b.allocate(1,N_sexes,"sigma_b");
  matage.allocate(1,N_sexes,"matage");
  mat_a.allocate(1,N_sexes,"mat_a");
  mat_b.allocate(1,N_sexes,"mat_b");
  summer_rv_n_years.allocate("summer_rv_n_years");
if (debug==1){char mark; cout<<"summer_rv_n_years "<<summer_rv_n_years<<endl; cin>>mark; }
  summer_rv_years.allocate(1,summer_rv_n_years,"summer_rv_years");
  summer_rv_total.allocate(1,summer_rv_n_years,"summer_rv_total");
  halibut_fixed_n_years.allocate("halibut_fixed_n_years");
  halibut_fixed_years.allocate(1,halibut_fixed_n_years,"halibut_fixed_years");
  halibut_fixed_total.allocate(1,halibut_fixed_n_years,"halibut_fixed_total");
  halibut_fixed_matureF.allocate(1,halibut_fixed_n_years,"halibut_fixed_matureF");
if (debug==1){char mark; cout<<"halibut_fixed_total "<<halibut_fixed_total<<endl; cin>>mark; }
  N_catlen_rv.allocate("N_catlen_rv");
  obs_catlen_rv.allocate(1,N_catlen_rv,-3,N_lens,"obs_catlen_rv");
if (debug==1){char mark; cout<<"obs_catlen_rv "<<obs_catlen_rv(106)<<endl; cin>>mark; }
  Canada.allocate(1960,last_year,"Canada");
  Foreign.allocate(1960,last_year,"Foreign");
  Total.allocate(1960,last_year,"Total");
  ll_F.allocate(1960,last_year,"ll_F");
  ll_M.allocate(1960,last_year,"ll_M");
  ot.allocate(1960,last_year,"ot");
  comm_total_n_years.allocate("comm_total_n_years");
  comm_total_years.allocate(1,comm_total_n_years,"comm_total_years");
  comm_total.allocate(1,comm_total_n_years,"comm_total");
  N_lens_comm.allocate("N_lens_comm");
  N_catlen_comm.allocate("N_catlen_comm");
if (debug==1){char mark; cout<<"N_catlen_comm "<<N_catlen_comm<<endl; cin>>mark; }
  obs_catlen_comm.allocate(1,N_catlen_comm,-3,N_lens_comm,"obs_catlen_comm");
if (debug==1){char mark; cout<<"obs_catlen_comm "<<obs_catlen_comm(66)<<endl; cin>>mark; }
  pva_horizon.allocate("pva_horizon");
  length_mort_vec.allocate("length_mort_vec");
  mort_vec.allocate(1,length_mort_vec,"mort_vec");
  lenatage.allocate(1,N_sexes,first_age,last_age,"lenatage");
  F_sigma.allocate(first_age,last_age,"F_sigma");
if (debug==1){char mark; cout<<"F_sigma "<<F_sigma<<endl; cin>>mark; }
  M_sigma.allocate(first_age,last_age,"M_sigma");
  dud.allocate("dud");
if (debug==1){char mark; cout<<"What number are we thinking of? "<<dud<<endl; cin>>mark; }
  w.allocate(1,N_sexes,first_age,last_age);
  obs_catlenP_rv.allocate(1,N_catlen_rv,1,N_lens);
  out_obs_catlen_rv.allocate(1,N_catlen_rv,-3,N_lens);
  obs_catlenP_comm.allocate(1,N_catlen_comm,1,N_lens_comm);
  out_obs_catlen_comm.allocate(1,N_catlen_comm,-3,N_lens_comm);
  mat.allocate(1,N_sexes,first_age,last_age);
annualsigma=log_RecDev_prior(6);  // cv for recruitment deviates
biasRannual=mfexp(-0.5*square(annualsigma));
srows=last_year+2-first_year;
  F_agelen.allocate(first_age,last_age,1,N_lens);
  M_agelen.allocate(first_age,last_age,1,N_lens);
  F_agelen_comm.allocate(first_age,last_age,1,N_lens_comm);
  M_agelen_comm.allocate(first_age,last_age,1,N_lens_comm);
  F_lenage.allocate(first_age,last_age,1,N_lens);
  M_lenage.allocate(first_age,last_age,1,N_lens);
  F_lenage_comm.allocate(first_age,last_age,1,N_lens_comm);
  M_lenage_comm.allocate(first_age,last_age,1,N_lens_comm);
  lens.allocate(1,N_lens);
lens.fill_seqadd(first_length,length_inc);
  lens_comm.allocate(1,N_lens_comm);
lens_comm.fill_seqadd(first_length,length_inc);
  ages.allocate(first_age,last_age);
ages.fill_seqadd(first_age,1);
  years.allocate(first_year,last_year+1);
years.fill_seqadd(first_year,1);
    int nopt=0;
    int on;
    char * end;
    if( (on=option_match(ad_comm::argc,ad_comm::argv,"-mcont",nopt)) >-1) Minit = strtod(ad_comm::argv[on+1], &end);
    if( (on=option_match(ad_comm::argc,ad_comm::argv,"-boss",nopt)) >-1) 
     {
     fminbins=strtod(ad_comm::argv[on+=1], &end);     
     fmaxbins=strtod(ad_comm::argv[on+=1], &end);     
     mminbins=strtod(ad_comm::argv[on+=1], &end);     
     mmaxbins=strtod(ad_comm::argv[on+=1], &end);     
     fminbins2=strtod(ad_comm::argv[on+=1], &end);     
     fmaxbins2=strtod(ad_comm::argv[on+=1], &end);     
     mminbins2=strtod(ad_comm::argv[on+=1], &end);     
     mmaxbins2=strtod(ad_comm::argv[on+=1], &end);     
     fminbins3=strtod(ad_comm::argv[on+=1], &end);     
     fmaxbins3=strtod(ad_comm::argv[on+=1], &end);     
     mminbins3=strtod(ad_comm::argv[on+=1], &end);     
     mmaxbins3=strtod(ad_comm::argv[on+=1], &end);     
     }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  dummy.allocate("dummy");
  #ifndef NO_AD_INITIALIZE
  dummy.initialize();
  #endif
phz=N1_prior(1);lb=N1_prior(2);ub=N1_prior(3);
  N1.allocate(lb,ub,phz,"N1");
phz=start_Z_prior(1);lb=start_Z_prior(2);ub=start_Z_prior(3);
  start_Z.allocate(lb,ub,phz,"start_Z");
phz=age1_prior(1);lb=age1_prior(2);ub=age1_prior(3);
  age1.allocate(first_year,last_year-2,lb,ub,phz,"age1");
phz=R0_prior(1);lb=R0_prior(2);ub=R0_prior(3);
  R0.allocate(lb,ub,phz,"R0");
phz=beta_prior(1);lb=beta_prior(2);ub=beta_prior(3);
  beta.allocate(lb,ub,phz,"beta");
phz=alpha_prior(1);lb=alpha_prior(2);ub=alpha_prior(3);
  alpha.allocate(lb,ub,phz,"alpha");
phz=p1_prior(1);lb=p1_prior(2);ub=p1_prior(3);
  p1.allocate(lb,ub,phz,"p1");
phz=p2_prior(1);lb=p2_prior(2);ub=p2_prior(3);
  p2.allocate(lb,ub,phz,"p2");
phz=p3_prior(1);lb=p3_prior(2);ub=p3_prior(3);
  p3.allocate(lb,ub,phz,"p3");
phz=M_prior(1);lb=M_prior(2);ub=M_prior(3);
  M.allocate(lb,ub,phz,"M");
phz=log_RecDev_prior(1);lb=log_RecDev_prior(2);ub=log_RecDev_prior(3);
  log_RecDev.allocate(Frecyr,last_year,lb,ub,phz,"log_RecDev");
  all_log_RecDev.allocate(first_year,last_year+1,"all_log_RecDev");
  #ifndef NO_AD_INITIALIZE
    all_log_RecDev.initialize();
  #endif
phz=old_summer_rv_Sfull_prior(1);lb=old_summer_rv_Sfull_prior(2);ub=old_summer_rv_Sfull_prior(3);
  old_summer_rv_Sfull.allocate(lb,ub,phz,"old_summer_rv_Sfull");
phz=old_summer_rv_varLest_prior(1);lb=old_summer_rv_varLest_prior(2);ub=old_summer_rv_varLest_prior(3);
  old_summer_rv_varLest.allocate(lb,ub,phz,"old_summer_rv_varLest");
phz=old_summer_rv_varRest_prior(1);lb=old_summer_rv_varRest_prior(2);ub=old_summer_rv_varRest_prior(3);
  old_summer_rv_varRest.allocate(lb,ub,phz,"old_summer_rv_varRest");
phz=halibut_fixed_Sfull_prior(1,1);lb=halibut_fixed_Sfull_prior(1,2);ub=halibut_fixed_Sfull_prior(1,3);
  halibut_fixed_Sfull_F.allocate(lb,ub,phz,"halibut_fixed_Sfull_F");
phz=halibut_fixed_varLest_prior(1,1);lb=halibut_fixed_varLest_prior(1,2);ub=halibut_fixed_varLest_prior(1,3);
  halibut_fixed_varLest_F.allocate(lb,ub,phz,"halibut_fixed_varLest_F");
phz=halibut_fixed_Sfull_prior(2,1);lb=halibut_fixed_Sfull_prior(2,2);ub=halibut_fixed_Sfull_prior(2,3);
  halibut_fixed_Sfull_M.allocate(lb,ub,phz,"halibut_fixed_Sfull_M");
phz=halibut_fixed_varLest_prior(2,1);lb=halibut_fixed_varLest_prior(2,2);ub=halibut_fixed_varLest_prior(2,3);
  halibut_fixed_varLest_M.allocate(lb,ub,phz,"halibut_fixed_varLest_M");
phz=ll_Sfull_prior(1,1);lb=ll_Sfull_prior(1,2);ub=ll_Sfull_prior(1,3);
  ll_Sfull_F.allocate(lb,ub,phz,"ll_Sfull_F");
phz=ll_varLest_prior(1,1);lb=ll_varLest_prior(1,2);ub=ll_varLest_prior(1,3);
  ll_varLest_F.allocate(lb,ub,phz,"ll_varLest_F");
phz=ll_Sfull_prior(2,1);lb=ll_Sfull_prior(2,2);ub=ll_Sfull_prior(2,3);
  ll_Sfull_M.allocate(lb,ub,phz,"ll_Sfull_M");
phz=ll_varLest_prior(2,1);lb=ll_varLest_prior(2,2);ub=ll_varLest_prior(2,3);
  ll_varLest_M.allocate(lb,ub,phz,"ll_varLest_M");
phz=ot_Sfull_prior(1);lb=ot_Sfull_prior(2);ub=ot_Sfull_prior(3);
  ot_Sfull.allocate(lb,ub,phz,"ot_Sfull");
phz=ot_varLest_prior(1);lb=ot_varLest_prior(2);ub=ot_varLest_prior(3);
  ot_varLest.allocate(lb,ub,phz,"ot_varLest");
phz=ot_varRest_prior(1);lb=ot_varRest_prior(2);ub=ot_varRest_prior(3);
  ot_varRest.allocate(lb,ub,phz,"ot_varRest");
phz=log_q_old_summer_rv_prior(1);lb=log_q_old_summer_rv_prior(2);ub=log_q_old_summer_rv_prior(3);
  log_q_old_summer_rv.allocate(lb,ub,phz,"log_q_old_summer_rv");
phz=log_q_halibut_fixed_prior(1);lb=log_q_halibut_fixed_prior(2);ub=log_q_halibut_fixed_prior(3);
  log_q_halibut_fixed.allocate(lb,ub,phz,"log_q_halibut_fixed");
phz=u_ll_prior(1);lb=u_ll_prior(2);ub=u_ll_prior(3);
  u_ll_2003.allocate(lb,ub,phz,"u_ll_2003");
  summer_rv_resid.allocate(1,summer_rv_n_years,"summer_rv_resid");
  #ifndef NO_AD_INITIALIZE
    summer_rv_resid.initialize();
  #endif
  halibut_fixed_resid.allocate(1,halibut_fixed_n_years,"halibut_fixed_resid");
  #ifndef NO_AD_INITIALIZE
    halibut_fixed_resid.initialize();
  #endif
  sex_rat_rv_resid.allocate(1,N_lens,"sex_rat_rv_resid");
  #ifndef NO_AD_INITIALIZE
    sex_rat_rv_resid.initialize();
  #endif
  sex_rat_comm_resid.allocate(1,N_lens_comm,"sex_rat_comm_resid");
  #ifndef NO_AD_INITIALIZE
    sex_rat_comm_resid.initialize();
  #endif
  pred_catlen_rv.allocate(1,N_catlen_rv,1,N_lens,"pred_catlen_rv");
  #ifndef NO_AD_INITIALIZE
    pred_catlen_rv.initialize();
  #endif
  pred_catlenP_rv.allocate(1,N_catlen_rv,1,N_lens,"pred_catlenP_rv");
  #ifndef NO_AD_INITIALIZE
    pred_catlenP_rv.initialize();
  #endif
  pred_catlen_comm.allocate(1,N_catlen_comm,1,N_lens_comm,"pred_catlen_comm");
  #ifndef NO_AD_INITIALIZE
    pred_catlen_comm.initialize();
  #endif
  pred_catlenP_comm.allocate(1,N_catlen_comm,1,N_lens_comm,"pred_catlenP_comm");
  #ifndef NO_AD_INITIALIZE
    pred_catlenP_comm.initialize();
  #endif
  catlen_survey_resid.allocate(1,N_catlen_rv,1,N_lens,"catlen_survey_resid");
  #ifndef NO_AD_INITIALIZE
    catlen_survey_resid.initialize();
  #endif
  catlen_comm_resid.allocate(1,N_catlen_comm,1,N_lens_comm,"catlen_comm_resid");
  #ifndef NO_AD_INITIALIZE
    catlen_comm_resid.initialize();
  #endif
  catlen_comm_kernal_resid.allocate(1,N_catlen_comm,1,N_lens_comm,"catlen_comm_kernal_resid");
  #ifndef NO_AD_INITIALIZE
    catlen_comm_kernal_resid.initialize();
  #endif
  catlen_survey_kernal_resid.allocate(1,N_catlen_rv,1,N_lens,"catlen_survey_kernal_resid");
  #ifndef NO_AD_INITIALIZE
    catlen_survey_kernal_resid.initialize();
  #endif
  Females.allocate(first_year,last_year+1,first_age,last_age,"Females");
  #ifndef NO_AD_INITIALIZE
    Females.initialize();
  #endif
  Males.allocate(first_year,last_year+1,first_age,last_age,"Males");
  #ifndef NO_AD_INITIALIZE
    Males.initialize();
  #endif
  naa2009.allocate(first_age,last_age,"naa2009");
  Total_N.allocate(first_year,last_year+1,"Total_N");
  #ifndef NO_AD_INITIALIZE
    Total_N.initialize();
  #endif
  Total_W.allocate(first_year,last_year+1,"Total_W");
  #ifndef NO_AD_INITIALIZE
    Total_W.initialize();
  #endif
  Total_F.allocate(first_year,last_year+1,"Total_F");
  #ifndef NO_AD_INITIALIZE
    Total_F.initialize();
  #endif
  Total_M.allocate(first_year,last_year+1,"Total_M");
  #ifndef NO_AD_INITIALIZE
    Total_M.initialize();
  #endif
  SSN.allocate(first_year,last_year+1,"SSN");
  #ifndef NO_AD_INITIALIZE
    SSN.initialize();
  #endif
  SSN2.allocate(first_year,last_year+1,"SSN2");
  #ifndef NO_AD_INITIALIZE
    SSN2.initialize();
  #endif
  SSB.allocate(first_year,last_year+1,"SSB");
  #ifndef NO_AD_INITIALIZE
    SSB.initialize();
  #endif
  SSB2.allocate(first_year,last_year+1,"SSB2");
  SSB2009.allocate("SSB2009");
  FF2009.allocate("FF2009");
  FM2009.allocate("FM2009");
  rec.allocate(first_year,last_year+1,"rec");
  #ifndef NO_AD_INITIALIZE
    rec.initialize();
  #endif
  dudfem.allocate(first_age,last_age,"dudfem");
  #ifndef NO_AD_INITIALIZE
    dudfem.initialize();
  #endif
  dudmale.allocate(first_age,last_age,"dudmale");
  #ifndef NO_AD_INITIALIZE
    dudmale.initialize();
  #endif
  surv.allocate(1,N_sexes,first_age,last_age,"surv");
  #ifndef NO_AD_INITIALIZE
    surv.initialize();
  #endif
  surv50.allocate(1,N_sexes,first_age,last_age,"surv50");
  #ifndef NO_AD_INITIALIZE
    surv50.initialize();
  #endif
  surv25.allocate(1,N_sexes,first_age,last_age,"surv25");
  #ifndef NO_AD_INITIALIZE
    surv25.initialize();
  #endif
  surv75.allocate(1,N_sexes,first_age,last_age,"surv75");
  #ifndef NO_AD_INITIALIZE
    surv75.initialize();
  #endif
  w_F.allocate(first_age,last_age,"w_F");
  #ifndef NO_AD_INITIALIZE
    w_F.initialize();
  #endif
  pred_summer_rv_total.allocate(1970,last_year,"pred_summer_rv_total");
  #ifndef NO_AD_INITIALIZE
    pred_summer_rv_total.initialize();
  #endif
  pred_halibut_fixed_total.allocate(1998,last_year,"pred_halibut_fixed_total");
  #ifndef NO_AD_INITIALIZE
    pred_halibut_fixed_total.initialize();
  #endif
  summer_RVage.allocate(1970,last_year,first_age,last_age,"summer_RVage");
  #ifndef NO_AD_INITIALIZE
    summer_RVage.initialize();
  #endif
  halibut_fixed_age.allocate(1998,last_year,first_age,last_age,"halibut_fixed_age");
  #ifndef NO_AD_INITIALIZE
    halibut_fixed_age.initialize();
  #endif
  VB_Canada.allocate(1960,last_year+1,"VB_Canada");
  #ifndef NO_AD_INITIALIZE
    VB_Canada.initialize();
  #endif
  VB_Foreign.allocate(1960,last_year+1,"VB_Foreign");
  #ifndef NO_AD_INITIALIZE
    VB_Foreign.initialize();
  #endif
  VB_ll.allocate(1960,last_year+1,"VB_ll");
  #ifndef NO_AD_INITIALIZE
    VB_ll.initialize();
  #endif
  VB_ot.allocate(1960,last_year+1,"VB_ot");
  #ifndef NO_AD_INITIALIZE
    VB_ot.initialize();
  #endif
  TB_Canada.allocate(1960,last_year+1,"TB_Canada");
  #ifndef NO_AD_INITIALIZE
    TB_Canada.initialize();
  #endif
  TB_Foreign.allocate(1960,last_year+1,"TB_Foreign");
  #ifndef NO_AD_INITIALIZE
    TB_Foreign.initialize();
  #endif
  TB_ll.allocate(1960,last_year+1,"TB_ll");
  #ifndef NO_AD_INITIALIZE
    TB_ll.initialize();
  #endif
  TB_ot.allocate(1960,last_year+1,"TB_ot");
  #ifndef NO_AD_INITIALIZE
    TB_ot.initialize();
  #endif
  C_ll_F.allocate(1960,last_year,"C_ll_F");
  #ifndef NO_AD_INITIALIZE
    C_ll_F.initialize();
  #endif
  C_ll_M.allocate(1960,last_year,"C_ll_M");
  #ifndef NO_AD_INITIALIZE
    C_ll_M.initialize();
  #endif
  C_ot.allocate(1960,last_year,"C_ot");
  #ifndef NO_AD_INITIALIZE
    C_ot.initialize();
  #endif
  u_Canada.allocate(1960,last_year,"u_Canada");
  #ifndef NO_AD_INITIALIZE
    u_Canada.initialize();
  #endif
  u_Foreign.allocate(1960,last_year,"u_Foreign");
  #ifndef NO_AD_INITIALIZE
    u_Foreign.initialize();
  #endif
  u_ll.allocate(1960,last_year,"u_ll");
  #ifndef NO_AD_INITIALIZE
    u_ll.initialize();
  #endif
  u_ll_pdf.allocate("u_ll_pdf");
  #ifndef NO_AD_INITIALIZE
  u_ll_pdf.initialize();
  #endif
  u_ot.allocate(1960,last_year,"u_ot");
  #ifndef NO_AD_INITIALIZE
    u_ot.initialize();
  #endif
  u_ll_F.allocate(first_year,last_year,"u_ll_F");
  u_ll_M.allocate(first_year,last_year,"u_ll_M");
  #ifndef NO_AD_INITIALIZE
    u_ll_M.initialize();
  #endif
  old_summer_rv_F_sel.allocate(first_age,last_age,"old_summer_rv_F_sel");
  #ifndef NO_AD_INITIALIZE
    old_summer_rv_F_sel.initialize();
  #endif
  old_summer_rv_M_sel.allocate(first_age,last_age,"old_summer_rv_M_sel");
  #ifndef NO_AD_INITIALIZE
    old_summer_rv_M_sel.initialize();
  #endif
  halibut_fixed_F_sel.allocate(first_age,last_age,"halibut_fixed_F_sel");
  #ifndef NO_AD_INITIALIZE
    halibut_fixed_F_sel.initialize();
  #endif
  halibut_fixed_M_sel.allocate(first_age,last_age,"halibut_fixed_M_sel");
  #ifndef NO_AD_INITIALIZE
    halibut_fixed_M_sel.initialize();
  #endif
  ll_F_sel.allocate(first_age,last_age,"ll_F_sel");
  #ifndef NO_AD_INITIALIZE
    ll_F_sel.initialize();
  #endif
  ll_M_sel.allocate(first_age,last_age,"ll_M_sel");
  #ifndef NO_AD_INITIALIZE
    ll_M_sel.initialize();
  #endif
  ot_F_sel.allocate(first_age,last_age,"ot_F_sel");
  #ifndef NO_AD_INITIALIZE
    ot_F_sel.initialize();
  #endif
  ot_M_sel.allocate(first_age,last_age,"ot_M_sel");
  #ifndef NO_AD_INITIALIZE
    ot_M_sel.initialize();
  #endif
  varL.allocate("varL");
  #ifndef NO_AD_INITIALIZE
  varL.initialize();
  #endif
  varR.allocate("varR");
  #ifndef NO_AD_INITIALIZE
  varR.initialize();
  #endif
  Sfull.allocate("Sfull");
  #ifndef NO_AD_INITIALIZE
  Sfull.initialize();
  #endif
  F_utotal.allocate(1960,last_year,first_age,last_age,"F_utotal");
  #ifndef NO_AD_INITIALIZE
    F_utotal.initialize();
  #endif
  M_utotal.allocate(1960,last_year,first_age,last_age,"M_utotal");
  #ifndef NO_AD_INITIALIZE
    M_utotal.initialize();
  #endif
  F_ucatch.allocate(1960,last_year,first_age,last_age,"F_ucatch");
  #ifndef NO_AD_INITIALIZE
    F_ucatch.initialize();
  #endif
  M_ucatch.allocate(1960,last_year,first_age,last_age,"M_ucatch");
  #ifndef NO_AD_INITIALIZE
    M_ucatch.initialize();
  #endif
  prop_landed.allocate(first_year,last_year+1,1,2,first_age,last_age,"prop_landed");
  #ifndef NO_AD_INITIALIZE
    prop_landed.initialize();
  #endif
  prop_fmort.allocate(first_year,last_year+1,1,2,first_age,last_age,"prop_fmort");
  #ifndef NO_AD_INITIALIZE
    prop_fmort.initialize();
  #endif
  S0.allocate("S0");
  #ifndef NO_AD_INITIALIZE
  S0.initialize();
  #endif
  summer_rv_like.allocate("summer_rv_like");
  #ifndef NO_AD_INITIALIZE
  summer_rv_like.initialize();
  #endif
  halibut_fixed_like.allocate("halibut_fixed_like");
  #ifndef NO_AD_INITIALIZE
  halibut_fixed_like.initialize();
  #endif
  catlen_survey_like.allocate(1,N_catlen_surveys,"catlen_survey_like");
  #ifndef NO_AD_INITIALIZE
    catlen_survey_like.initialize();
  #endif
  catlen_comm_like.allocate(1,N_catlen_commercial,"catlen_comm_like");
  #ifndef NO_AD_INITIALIZE
    catlen_comm_like.initialize();
  #endif
  catlen_survey_kernal_like.allocate(1,N_catlen_surveys,"catlen_survey_kernal_like");
  #ifndef NO_AD_INITIALIZE
    catlen_survey_kernal_like.initialize();
  #endif
  catlen_comm_kernal_like.allocate(1,N_catlen_commercial,"catlen_comm_kernal_like");
  #ifndef NO_AD_INITIALIZE
    catlen_comm_kernal_like.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  ffpen.allocate("ffpen");
  #ifndef NO_AD_INITIALIZE
  ffpen.initialize();
  #endif
  upen.allocate("upen");
  #ifndef NO_AD_INITIALIZE
  upen.initialize();
  #endif
  reccheck.allocate("reccheck");
  #ifndef NO_AD_INITIALIZE
  reccheck.initialize();
  #endif
  recpen.allocate("recpen");
  #ifndef NO_AD_INITIALIZE
  recpen.initialize();
  #endif
  recpen2.allocate("recpen2");
  #ifndef NO_AD_INITIALIZE
  recpen2.initialize();
  #endif
  r_penF.allocate("r_penF");
  #ifndef NO_AD_INITIALIZE
  r_penF.initialize();
  #endif
  alphapen.allocate("alphapen");
  #ifndef NO_AD_INITIALIZE
  alphapen.initialize();
  #endif
  age1_pen.allocate("age1_pen");
  #ifndef NO_AD_INITIALIZE
  age1_pen.initialize();
  #endif
  age1_last_year_pdf.allocate("age1_last_year_pdf");
  #ifndef NO_AD_INITIALIZE
  age1_last_year_pdf.initialize();
  #endif
  recdiff_pdf.allocate("recdiff_pdf");
  #ifndef NO_AD_INITIALIZE
  recdiff_pdf.initialize();
  #endif
  recdiff_sigma.allocate("recdiff_sigma");
  #ifndef NO_AD_INITIALIZE
  recdiff_sigma.initialize();
  #endif
  temp1.allocate("temp1");
  #ifndef NO_AD_INITIALIZE
  temp1.initialize();
  #endif
  gm.allocate("gm");
  #ifndef NO_AD_INITIALIZE
  gm.initialize();
  #endif
  pred_f.allocate(1,N_catlen_rv,1,N_lens,"pred_f");
  #ifndef NO_AD_INITIALIZE
    pred_f.initialize();
  #endif
  pred_m.allocate(1,N_catlen_rv,1,N_lens,"pred_m");
  #ifndef NO_AD_INITIALIZE
    pred_m.initialize();
  #endif
  pred_f_comm.allocate(1,N_catlen_comm,1,N_lens_comm,"pred_f_comm");
  #ifndef NO_AD_INITIALIZE
    pred_f_comm.initialize();
  #endif
  pred_m_comm.allocate(1,N_catlen_comm,1,N_lens_comm,"pred_m_comm");
  #ifndef NO_AD_INITIALIZE
    pred_m_comm.initialize();
  #endif
  pred_sex_rat.allocate(1,N_catlen_rv,1,N_lens,"pred_sex_rat");
  #ifndef NO_AD_INITIALIZE
    pred_sex_rat.initialize();
  #endif
  pred_sex_rat_comm.allocate(1,N_catlen_comm,1,N_lens_comm,"pred_sex_rat_comm");
  #ifndef NO_AD_INITIALIZE
    pred_sex_rat_comm.initialize();
  #endif
  fcomp.allocate(1,12,"fcomp");
  #ifndef NO_AD_INITIALIZE
    fcomp.initialize();
  #endif
  ofv.allocate("ofv");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  SPRF0.allocate("SPRF0");
  #ifndef NO_AD_INITIALIZE
  SPRF0.initialize();
  #endif
  ssbeq.allocate("ssbeq");
  #ifndef NO_AD_INITIALIZE
  ssbeq.initialize();
  #endif
  receq.allocate("receq");
  #ifndef NO_AD_INITIALIZE
  receq.initialize();
  #endif
  obs_catlen_rv_short.allocate(1,N_catlen_rv,fminbins2,fmaxbins2,"obs_catlen_rv_short");
  #ifndef NO_AD_INITIALIZE
    obs_catlen_rv_short.initialize();
  #endif
  pred_catlen_rv_short.allocate(1,N_catlen_rv,fminbins2,fmaxbins2,"pred_catlen_rv_short");
  #ifndef NO_AD_INITIALIZE
    pred_catlen_rv_short.initialize();
  #endif
  obs_catlen_ll_short_female.allocate(1,N_catlen_comm,fminbins,fmaxbins,"obs_catlen_ll_short_female");
  #ifndef NO_AD_INITIALIZE
    obs_catlen_ll_short_female.initialize();
  #endif
  pred_catlen_ll_short_female.allocate(1,N_catlen_comm,fminbins,fmaxbins,"pred_catlen_ll_short_female");
  #ifndef NO_AD_INITIALIZE
    pred_catlen_ll_short_female.initialize();
  #endif
  obs_catlen_ll_short_male.allocate(1,N_catlen_comm,mminbins,mmaxbins,"obs_catlen_ll_short_male");
  #ifndef NO_AD_INITIALIZE
    obs_catlen_ll_short_male.initialize();
  #endif
  pred_catlen_ll_short_male.allocate(1,N_catlen_comm,mminbins,mmaxbins,"pred_catlen_ll_short_male");
  #ifndef NO_AD_INITIALIZE
    pred_catlen_ll_short_male.initialize();
  #endif
  obs_catlen_ot_short.allocate(1,N_catlen_comm,fminbins3,fmaxbins3,"obs_catlen_ot_short");
  #ifndef NO_AD_INITIALIZE
    obs_catlen_ot_short.initialize();
  #endif
  pred_catlen_ot_short.allocate(1,N_catlen_comm,fminbins3,fmaxbins3,"pred_catlen_ot_short");
  #ifndef NO_AD_INITIALIZE
    pred_catlen_ot_short.initialize();
  #endif
  obs_catlen_hsurvey_short_female.allocate(1,N_catlen_rv,fminbins,fmaxbins,"obs_catlen_hsurvey_short_female");
  #ifndef NO_AD_INITIALIZE
    obs_catlen_hsurvey_short_female.initialize();
  #endif
  pred_catlen_hsurvey_short_female.allocate(1,N_catlen_rv,fminbins,fmaxbins,"pred_catlen_hsurvey_short_female");
  #ifndef NO_AD_INITIALIZE
    pred_catlen_hsurvey_short_female.initialize();
  #endif
  obs_catlen_hsurvey_short_male.allocate(1,N_catlen_rv,mminbins,mmaxbins,"obs_catlen_hsurvey_short_male");
  #ifndef NO_AD_INITIALIZE
    obs_catlen_hsurvey_short_male.initialize();
  #endif
  pred_catlen_hsurvey_short_male.allocate(1,N_catlen_rv,mminbins,mmaxbins,"pred_catlen_hsurvey_short_male");
  #ifndef NO_AD_INITIALIZE
    pred_catlen_hsurvey_short_male.initialize();
  #endif
  obs_catlen_hCI_short_female.allocate(1,N_catlen_rv,fminbins,fmaxbins,"obs_catlen_hCI_short_female");
  #ifndef NO_AD_INITIALIZE
    obs_catlen_hCI_short_female.initialize();
  #endif
  pred_catlen_hCI_short_female.allocate(1,N_catlen_rv,fminbins,fmaxbins,"pred_catlen_hCI_short_female");
  #ifndef NO_AD_INITIALIZE
    pred_catlen_hCI_short_female.initialize();
  #endif
  obs_catlen_hCI_short_male.allocate(1,N_catlen_rv,mminbins,mmaxbins,"obs_catlen_hCI_short_male");
  #ifndef NO_AD_INITIALIZE
    obs_catlen_hCI_short_male.initialize();
  #endif
  pred_catlen_hCI_short_male.allocate(1,N_catlen_rv,mminbins,mmaxbins,"pred_catlen_hCI_short_male");
  #ifndef NO_AD_INITIALIZE
    pred_catlen_hCI_short_male.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
  //VH
  int i1=0;
  int i2=0;
  i1 = (option_match(ad_comm::argc,ad_comm::argv,"-phase",i2)); 
  if(i1<=0) initialize_parameters(); 
  for(obs=1;obs<=N_catlen_comm;obs++)
    {
    if(obs_catlen_comm(obs,-2)==1)	// select longline data
     {
     for(i=1;i<=8;i++)		// make all values < 34 cm zero
       {
       obs_catlen_comm(obs,i)=0;
       }
     }
    }
  setup_fixed_parameters();
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1e-4}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{200,200,400,500,500,50000,50000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::userfunction(void)
{
  ofv =0.0;
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  dummy=-5;
  fpen=0.;
  ffpen=0.;
  recpen=0;
  recpen2=0;
  r_penF=0;
  reccheck=0;
  upen=0;
  alphapen=0.;
  age1_pen=0.;
  ofv=0;
  setup_variable_parameters();
  initial_conditions();
  dynamics();
 // evaluating all the likelihoods
   Summer_RV_Likelihood();
   Halibut_Fixed_Likelihood();
   Catlen_Comm_Likelihood();
   Catlen_Survey_Likelihood();
   Age1_pen();
   Age1_last_year_pdf();
   Recdiff_pdf();
  if(dummy>1)
    {
    ofv=square(dummy);
    }
  else
    {
  // adding up all the log-likelihoods
  ofv+=ffpen*10; // pen on estimating F iteratively
  ofv+=recdiff_pdf;      //random walk on recruitment	
    for(i=1;i<=12;i++)
     {
     ofv+= fcomp(i)*like_weight(i);
     }
   if(mceval_phase()) 
   {
   MCWrite(); 
   }
  } //close else
}

void model_parameters::initialize_parameters(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  N1=N1_prior(7);
  start_Z=start_Z_prior(7);
  R0=R0_prior(7);
  beta=beta_prior(7);
  alpha=alpha_prior(7);
  p1=p1_prior(7);
  p2=p2_prior(7);
  p3=p3_prior(7);
  M=M_prior(7);
  log_RecDev=log_RecDev_prior(7);
  all_log_RecDev = 0.;
  age1(1970)=580662.6;
  age1(1971)=580662.6;
  age1(1972)=580662.6;
  age1(1973)=580662.6;
  age1(1974)=580662.6;
  age1(1975)=580662.6;
  age1(1976)=580662.6;
  age1(1977)=580662.6;
  age1(1978)=580662.6;
  age1(1979)=580662.6;
  age1(1980)=580662.6;
  age1(1981)=580662.6;
  age1(1982)=580662.6;
  age1(1983)=580662.6;
  age1(1984)=580662.6;
  age1(1985)=580662.6;
  age1(1986)=580662.6;
  age1(1987)=580662.6;
  age1(1988)=958432;
  age1(1989)=864676;
  age1(1990)=843690;
  age1(1991)=818227;
  age1(1992)=757501;
  age1(1993)=685773;
  age1(1994)=644483;
  age1(1995)=607888;
  age1(1996)=574659;
  age1(1997)=559200;
  age1(1998)=459269;
  age1(1999)=343876;
  age1(2000)=263401;
  age1(2001)=190435;
  age1(2002)=138422;
  age1(2003)=580662.6;
  age1(2004)=580662.6;
  age1(2005)=580662.6;
  age1(2006)=580662.6;
  age1(2007)=580662.6;
  old_summer_rv_Sfull=old_summer_rv_Sfull_prior(7);
  old_summer_rv_varLest=old_summer_rv_varLest_prior(7);
  old_summer_rv_varRest=old_summer_rv_varRest_prior(7);
  halibut_fixed_Sfull_F=halibut_fixed_Sfull_prior(1,7);
  halibut_fixed_varLest_F=halibut_fixed_varLest_prior(1,7);
  halibut_fixed_Sfull_M=halibut_fixed_Sfull_prior(2,7);
  halibut_fixed_varLest_M=halibut_fixed_varLest_prior(2,7);
  ll_Sfull_F=ll_Sfull_prior(1,7);
  ll_varLest_F=ll_varLest_prior(1,7);
  ll_Sfull_M=ll_Sfull_prior(2,7);
  ll_varLest_M=ll_varLest_prior(2,7);
  ot_Sfull=ot_Sfull_prior(7);
  ot_varLest=ot_varLest_prior(7);
  ot_varRest=ot_varRest_prior(7);
  ot_Sfull=ot_Sfull_prior(7);
  ot_varLest=ot_varLest_prior(7);
  ot_varRest=ot_varRest_prior(7);
  log_q_old_summer_rv=log_q_old_summer_rv_prior(7);
  log_q_halibut_fixed=log_q_halibut_fixed_prior(7);
  u_ll_2003=u_ll_prior(7);
}

void model_parameters::setup_fixed_parameters(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  F_agelen=0.;
  M_agelen=0.;
  F_agelen_comm=0.;
  M_agelen_comm=0.;
  F_lenage=0.;
  M_lenage=0.;
  F_lenage_comm=0.;
  M_lenage_comm=0.;
  mat=0.;
  out_obs_catlen_rv = 0.;
  obs_catlenP_rv = 0.;
  out_obs_catlen_comm = 0.;
  obs_catlenP_comm = 0.;
  double pi=3.14159265358;
  for (sex=1;sex<=N_sexes;sex++)
  {
  for (age=first_age;age<=last_age;age++)
  {
  mat(sex,age)=1/(1+exp(-mat_b(sex)*((linf(sex)-(linf(sex)-t0(sex))*exp(-k(sex)*age))-mat_a(sex))));
  if(age>=5){F_sigma(age)=10;}
  if(age>=4){M_sigma(age)=10;}
  }
  }
  for (sex=1;sex<=N_sexes;sex++)
  {
  w(sex)=bi(sex)*pow(lenatage(sex),bii(sex))/1000;			// convert to kg
  }
  w_F=w(1);
  int mls_bin=25;
  dvariable discard_mort=0.23;
  dvector temp(first_age,last_age);
  dvector temp2_F(1,N_lens);
  dvector temp2_M(1,N_lens);
  dvector temp3_F(1,N_lens_comm);
  dvector temp3_M(1,N_lens_comm);
  temp=0;
  temp2_F=0;
  temp2_M=0;
  temp3_F=0;
  temp3_M=0;
  for (age=first_age;age<=last_age; age++)
  {
  F_agelen(age) = (length_inc/(sqrt(2*pi)*F_sigma(age)))*mfexp(-square(lens-lenatage(1,age))/(2*square(F_sigma(age))));
  M_agelen(age) = (length_inc/(sqrt(2*pi)*M_sigma(age)))*mfexp(-square(lens-lenatage(2,age))/(2*square(M_sigma(age))));
  temp=rowsum(F_agelen);
  F_agelen(age)/=temp(age); 
  temp=rowsum(M_agelen);
  M_agelen(age)/=temp(age);
  F_agelen_comm(age) = (length_inc/(sqrt(2*pi)*F_sigma(age)))*mfexp(-square(lens_comm-lenatage(1,age))/(2*square(F_sigma(age))));
  M_agelen_comm(age) = (length_inc/(sqrt(2*pi)*M_sigma(age)))*mfexp(-square(lens_comm-lenatage(2,age))/(2*square(M_sigma(age))));
  temp=rowsum(F_agelen_comm);
  F_agelen_comm(age)/=temp(age); 
  temp=rowsum(M_agelen_comm);
  M_agelen_comm(age)/=temp(age);
  F_lenage(age) = (length_inc/(sqrt(2*pi)*F_sigma(age)))*mfexp(-square(lens-lenatage(1,age))/(2*square(F_sigma(age))));
  M_lenage(age) = (length_inc/(sqrt(2*pi)*M_sigma(age)))*mfexp(-square(lens-lenatage(2,age))/(2*square(M_sigma(age))));
  F_lenage_comm(age) = (length_inc/(sqrt(2*pi)*F_sigma(age)))*mfexp(-square(lens_comm-lenatage(1,age))/(2*square(F_sigma(age))));
  M_lenage_comm(age) = (length_inc/(sqrt(2*pi)*M_sigma(age)))*mfexp(-square(lens_comm-lenatage(2,age))/(2*square(M_sigma(age))));
    for (sex=1;sex<=2;sex++)
    {
       for (year=first_year; year<=last_year+1; year++)
       {
          if (year<=1993)
          {
            prop_landed(year,1,age)= 1;
            prop_landed(year,2,age)= 1;
            prop_fmort(year,1,age)= 1;
            prop_fmort(year,2,age)= 1;
          }
          else {
            prop_landed(year,1,age)= sum(F_agelen(age)(mls_bin,N_lens));
            prop_landed(year,2,age)= sum(M_agelen(age)(mls_bin,N_lens));
            prop_fmort(year,1,age)= sum(F_agelen(age)(1,mls_bin)*discard_mort) + sum(F_agelen(age)(mls_bin,N_lens));
            prop_fmort(year,2,age)= sum(M_agelen(age)(1,mls_bin)*discard_mort) + sum(M_agelen(age)(mls_bin,N_lens));
          }
       }
    }
  }   // close age
  for (i=1;i<=N_lens;i++)
  {
  for (age=first_age;age<=last_age; age++)
  {
  temp2_F(i)+=F_lenage(age,i);
  temp2_M(i)+=M_lenage(age,i);
  }
  }
  for (i=1;i<=N_lens;i++)
  {
  for (age=first_age;age<=last_age; age++)
  {
  F_lenage(age,i)/=temp2_F(i); 
  M_lenage(age,i)/=temp2_M(i); 
  }
  }
  for (i=1;i<=N_lens_comm;i++)
  {
  for (age=first_age;age<=last_age; age++)
  {
  temp3_F(i)+=F_lenage_comm(age,i);
  temp3_M(i)+=M_lenage_comm(age,i);
  }
  }
  for (i=1;i<=N_lens_comm;i++)
  {
  for (age=first_age;age<=last_age; age++)
  {
  F_lenage_comm(age,i)/=temp3_F(i); 
  M_lenage_comm(age,i)/=temp3_M(i); 
  }
  }
  for (obs=1;obs<=N_catlen_rv;obs++)
  {
   obs_catlenP_rv(obs)(Flength,N_lens) = obs_catlen_rv(obs)(Flength,N_lens) / sum(obs_catlen_rv(obs)(Flength,N_lens));
   out_obs_catlen_rv(obs)(-3,0) = obs_catlen_rv(obs)(-3,0);
   out_obs_catlen_rv(obs)(-2,0) = obs_catlen_rv(obs)(-2,0);
   out_obs_catlen_rv(obs)(1,N_lens) = obs_catlenP_rv(obs)(1,N_lens);
  }
 for (obs=1;obs<=N_catlen_comm;obs++)
  {
   obs_catlenP_comm(obs)(Flength,N_lens_comm) = obs_catlen_comm(obs)(Flength,N_lens_comm) / sum(obs_catlen_comm(obs)(Flength,N_lens_comm));
   out_obs_catlen_comm(obs)(-3,0) = obs_catlen_comm(obs)(-3,0);
   out_obs_catlen_comm(obs)(-2,0) = obs_catlen_comm(obs)(-2,0);
   out_obs_catlen_comm(obs)(1,N_lens_comm) = obs_catlenP_comm(obs)(1,N_lens_comm);
  }
  //not working. need to include discard mortality
  for (year=1960;year<=2009;year++) 
  {
   C_ll_M(year) = ll_M(year);
   C_ll_F(year) = ll_F(year);
   C_ot(year) = ot(year);
  }
}

void model_parameters::setup_variable_parameters(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  old_summer_rv_F_sel = Make_select_dhn(first_age, last_age,old_summer_rv_Sfull,old_summer_rv_varLest,old_summer_rv_varRest);
  ot_F_sel = Make_select_dhn(first_age, last_age,ot_Sfull,ot_varLest,ot_varRest);
  halibut_fixed_F_sel = Make_select_logistic(first_age, last_age,halibut_fixed_Sfull_F,halibut_fixed_varLest_F);
  halibut_fixed_M_sel = Make_select_logistic(first_age, last_age,halibut_fixed_Sfull_M,halibut_fixed_varLest_M);
  ll_F_sel = Make_select_logistic(first_age, last_age,ll_Sfull_F,ll_varLest_F);
  ll_M_sel = Make_select_logistic(first_age, last_age,ll_Sfull_M,ll_varLest_M);
  surv = mfexp(-M);
  surv50 = mfexp(-M*0.5);
  surv25 = mfexp(-M*0.25);
  surv75 = mfexp(-M*0.75);
  all_log_RecDev(Frecyr,last_year)=log_RecDev;  //Frecyr
}

void model_parameters::initial_conditions(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  Females.initialize();
  Males.initialize();
  Total_N.initialize();
  Total_W.initialize();
  Total_F.initialize();
  Total_M.initialize();
  SSN.initialize();
  SSN2.initialize();
  SSB.initialize();
  SSB2.initialize();
  Females=0;
  Males=0;
  Total_N=0;
  Total_W=0;
  Total_F=0;
  Total_M=0;
  SSN=0;
  SSN2=0;
  SSB=0;
  SSB2=0;
  SSB2009=0;
  FF2009=0;
  FM2009=0;
  //SPR in biomass of female biomass (not numbers)   
  dudfem(first_age)=0.5;   //age-1 1:1 sex ratio
  for (age=first_age+1;age<=last_age-1;age++)
     {
     dudfem(age)=dudfem(age-1)*surv(1,age-1);
     } 
   dudfem(last_age) = dudfem(last_age-1)*surv(1,last_age-1)/(1-surv(1,last_age));
   SPRF0=sum(elem_prod(elem_prod(dudfem,mat(1)),w(1)));
  // B-H  see Gibsion + Myers 2003 alewife Can. Tec. Report 2468
  ssbeq=p1*(-log(p2/(2.71828*p1*SPRF0)));   //  Ricker reparameterized
   receq=2.71828*p2/p1*ssbeq*exp(-ssbeq/p1);	// Ricker reparamertized   
 if(N1opt==0)  
  {
  Females(first_year,first_age)=receq*0.5;	//	assume 50:50
  Males(first_year,first_age)=receq*0.5;
  }
 if(N1opt==1)
  {
  reccheck=receq-N1;
  reccheck=posfun(reccheck,1,recpen);
  Females(first_year,first_age)= N1*0.5;
  Males(first_year,first_age)= N1*0.5;
  }
 if(N1opt==2)
  {
  for (year=first_year;year<=last_year;year++)
  {
  Females(year,first_age)=age1(year)*0.5;
  Males(year,first_age)=age1(year)*0.5;
  }
  for (age=first_age+1;age<=last_age;age++)
  {
  Females(first_year,age)=year1(age)*0.5;
  Males(first_year,age)=year1(age)*0.5;
  }
  }
 if(N1opt==3)
  {
  for (year=first_year;year<=last_year-2;year++)
  {
  Females(year,first_age)=age1(year)*0.5;
  Males(year,first_age)=age1(year)*0.5;
  }
  for (age=first_age+1;age<=last_age;age++)
  {
  Females(first_year,age)=Females(first_year,age-1)*exp(-start_Z);
  Males(first_year,age)=Males(first_year,age-1)*exp(-start_Z);
  }
  }
  Total_N(first_year)=sum(Females(first_year))+sum(Males(first_year));
  Total_F(first_year)=(sum(elem_prod(Females(first_year),w(1))))/1000; 	// in t
  Total_M(first_year)=(sum(elem_prod(Males(first_year),w(2))))/1000; 	// in t
  Total_W(first_year)=Total_F(first_year)+Total_M(first_year);
  SSN(first_year) = sum(elem_prod(Females(first_year),mat(1)));
  SSN2(first_year) = sum(elem_prod(Females(first_year),mat(1))) +sum(elem_prod(Males(first_year),mat(2)));
  SSB(first_year) = sum(elem_prod(elem_prod(Females(first_year),mat(1)),w(1)))/1000;  // convert to tons
  SSB2(first_year) = (sum(elem_prod(elem_prod(Females(first_year),mat(1)),w(1))) + sum(elem_prod(elem_prod(Males(first_year),mat(2)),w(2))))/1000;  // convert to tons
  VB_ll(first_year)=((elem_prod(elem_prod(elem_prod(prop_landed(first_year,1),Females(first_year)),surv50(1)),ll_F_sel)*w(1))+(elem_prod(elem_prod(Males(first_year),surv50(2)),elem_prod(ll_M_sel,prop_landed(first_year,2)))*w(2)))/1000; 
  TB_Canada(first_year)=(elem_prod(Females(first_year),surv50(1))*w(1)+
    					 elem_prod(  Males(first_year),surv50(2))*w(2))/1000; 
}

void model_parameters::dynamics(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  F_utotal.initialize();
  M_utotal.initialize();
  dvariable WF; 
  dvariable CF;
  dvariable WM; 
  dvariable CM;
  dvariable COT;
  dvar_vector test_sel(1,20);
  test_sel=1;
  test_sel(1)=0;
  test_sel(2)=0;
  test_sel(3)=0;
  test_sel(4)=.5;
  for (year=first_year;year<=last_year;year++)
  {
  WF=Total_F(year);
  WM=Total_M(year);
  if(year<1988)
   {
   CF=C_ll_F(year)/(1-0.055);		//kluge. bump up landings to total catch
   CM=C_ll_M(year)/(1-0.055);
   }else{  
   CF=C_ll_F(year);
   CM=C_ll_M(year);
   }
  if(year<1984)
   {
   COT=C_ot(year)/(1-0.055);
   }else{
   COT=C_ot(year);
   }   
  u_ll_F(year)=GetF_forward(WF,CF,M);
  u_ll_M(year)=GetF_forward(WM,CM,M);
  u_ot(year)=GetF_forward(WF+WM,COT,M);
      for (age=first_age;age<=last_age-1;age++)
        {
         F_utotal(year,age)=prop_fmort(year,1,age)*(u_ll_F(year)*ll_F_sel(age)+ u_ot(year)*ot_F_sel(age)*.5);
         M_utotal(year,age)=prop_fmort(year,2,age)*(u_ll_M(year)*ll_M_sel(age)+ u_ot(year)*ot_F_sel(age)*.5);
         F_ucatch(year,age)=(u_ll_F(year)*ll_F_sel(age));
         M_ucatch(year,age)=(u_ll_M(year)*ll_M_sel(age));
         Females(year+1,age+1)=Females(year,age)*exp(-F_utotal(year,age))*surv(1,age);
         Males(year+1,age+1)=Males(year,age)*exp(-M_utotal(year,age))*surv(2,age);
        }
      //Max Age
         F_utotal(year,last_age)=prop_fmort(year,1,last_age)*(u_ll_F(year)*ll_F_sel(last_age)+ u_ot(year)*ot_F_sel(age)*.5);
         M_utotal(year,last_age)=prop_fmort(year,2,last_age)*(u_ll_M(year)*ll_M_sel(last_age)+ u_ot(year)*ot_F_sel(age)*.5);
         F_ucatch(year,last_age)=(u_ll_F(year)*ll_F_sel(last_age));
         M_ucatch(year,last_age)=(u_ll_M(year)*ll_M_sel(last_age));
         Females(year+1,last_age)=Females(year,last_age)*exp(-F_utotal(year,last_age))*surv(1,last_age);
         Males(year+1,last_age)=Males(year,last_age)*exp(-M_utotal(year,last_age))*surv(2,last_age);
  SSB(year) = sum(elem_prod(elem_prod(Females(year),mat(1)),w(1)))/1000;	 // convert to tons
  SSB2(year) = (sum(elem_prod(elem_prod(Females(year),mat(1)),w(1))) + sum(elem_prod(elem_prod(Males(year),mat(2)),w(2))))/1000;  // convert to tons
  SSB2009=SSB2(last_year);
  naa2009=Females(last_year)+Males(last_year);
  FF2009=F_utotal(last_year,10);
  FM2009=M_utotal(last_year,10);
  Females(year+1,first_age)=age1(year+1)*.5;
  Males(year+1,first_age)=age1(year+1)*.5;
        if(log_RecDev_prior(1)>0)
        {
        }
       else
        {
        }
  int index;
  index=3;
  temp1=1.0;
  double temp2;
  temp2=1/(double)index;
    for(i=2;i<=index+1;i++)
    {
    temp1*=Females(last_year-i,first_age);
    }
    gm=pow(temp1,temp2);
    for(i=-1;i<=index-2;i++)
    {
    Females(last_year-i,first_age)=gm;
    Males(last_year-i,first_age)=gm;
    }
    Females(last_year,first_age+1)=gm*exp(-F_utotal(last_year-1,first_age))*surv(1,age);
    Males(last_year,first_age+1)=gm*exp(-M_utotal(last_year-1,first_age))*surv(2,age);
  rec=column(Females,1)*2;  
   VB_ll(year+1)=((elem_prod(elem_prod(Females(year+1),surv50(1)),elem_prod(prop_landed(year+1,1),ll_F_sel))*w(1))+
                   (elem_prod(elem_prod(Males(year+1),surv50(2)),elem_prod(prop_landed(year+1,2),ll_M_sel))*w(2)))/1000;
   TB_Canada(year+1)=(elem_prod(Females(year+1),surv50(1))*w(1)+
   						  elem_prod(  Males(year+1),surv50(2))*w(2))/1000; 
  dvar_vector temp_F(first_year,last_year);
  dvar_vector temp_M(first_year,last_year);
  temp_F=0;
  temp_M=0;
      for (age=first_age;age<=last_age;age++)
       {
       temp_F(year)+=Females(year+1,age)*w(1,age)/1000;
       temp_M(year)+=Males(year+1,age)*w(2,age)/1000;
       }  
   Total_N(year+1)=sum(Females(year+1))+sum(Males(year+1));
   Total_F(year+1)=temp_F(year);
   Total_M(year+1)=temp_M(year);
   Total_W(year+1)=Total_F(year+1)+Total_M(year+1);
   SSN(year+1) = sum(elem_prod(Females(year+1),mat(1)));             
   SSN2(year+1) = sum(elem_prod(Females(year+1),mat(1))) + sum(elem_prod(Males(year+1),mat(2)));             
   }
  for (year=first_year;year<=last_year;year++)
      {
      for (age=first_age;age<=last_age;age++)
       {
       summer_RVage(year,age)=Females(year,age)*exp(-(F_utotal(year,age)+M)*.5)*exp(log_q_old_summer_rv)*old_summer_rv_F_sel(age)  +  
                               Males(year,age)*exp(-(M_utotal(year,age)+M)*.5)*exp(log_q_old_summer_rv)*old_summer_rv_F_sel(age);
       }
      }
  for (year=1998;year<=last_year;year++)
      {
      for (age=first_age;age<=last_age;age++)
       {
        halibut_fixed_age(year,age)=Females(year,age)*w(1,age)*exp(-(F_utotal(year,age)+M)*.5)*exp(log_q_halibut_fixed)*halibut_fixed_F_sel(age) +  
                               Males(year,age)*w(2,age)*exp(-(M_utotal(year,age)+M)*.5)*exp(log_q_halibut_fixed)*halibut_fixed_M_sel(age);
       }
      }
  pred_summer_rv_total=rowsum(summer_RVage);
  pred_halibut_fixed_total=rowsum(halibut_fixed_age);
  recpen2=pow(log(Females(first_year,first_age))-log(Females(first_year+1,first_age)),2);
}

void model_parameters::Summer_RV_Likelihood(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
   dvar_vector junk(1,summer_rv_n_years);
    for(i=1;i<=summer_rv_n_years;i++)
     {
     junk(i)=pred_summer_rv_total(summer_rv_years(i));   
     }
    summer_rv_like=LogNormal_LogLike(summer_rv_total+.01,junk+.01,summer_rv_n_years,.1,2,summer_rv_resid);
    //  cout<<summer_rv_like<<endl;
	// ad_exit(1);
    fcomp(1)=summer_rv_like;
}

void model_parameters::Halibut_Fixed_Likelihood(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  dvar_vector junk(1,halibut_fixed_n_years);
    for(i=1;i<=halibut_fixed_n_years;i++)
     {
     junk(i)=pred_halibut_fixed_total(halibut_fixed_years(i));   
     }
   halibut_fixed_like=LogNormal_LogLike(halibut_fixed_total+.01,junk+.01,halibut_fixed_n_years,0.1,2,halibut_fixed_resid);
   fcomp(5)=halibut_fixed_like;
}

void model_parameters::Catlen_Comm_Likelihood(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  pred_catlen_comm.initialize();
  catlen_comm_like=0.;
  int ssex;
  for(obs=1;obs<=N_catlen_comm;obs++)
  {
    series = obs_catlen_comm(obs,-2);
    ssex = obs_catlen_comm(obs,-0);
    if(obs_catlen_comm(obs,-3)>maxssizes)		
    {
    maxs = maxssizes;
    }
    else
    {
    maxs = obs_catlen_comm(obs,-3);
    }
    year=obs_catlen_comm(obs,-1);
    //1=LL
    //2=OT
    switch (series)
    {
    case 1:
	if(ssex==1)
	{
	pred_catlen_comm(obs) =elem_prod(elem_prod(exp(-(F_utotal(year)+M)),Females(year)),ll_F_sel)*F_agelen_comm;
	}
	else
	{
	pred_catlen_comm(obs) =elem_prod(elem_prod(exp(-(M_utotal(year)+M)),Males(year)),ll_M_sel)*M_agelen_comm;
	}
    break;
    case 2:
	pred_catlen_comm(obs)+=elem_prod(elem_prod(exp(-(F_utotal(year)+M)),Females(year)),ot_F_sel)*F_agelen_comm;
	pred_catlen_comm(obs)+=elem_prod(elem_prod(exp(-(M_utotal(year)+M)),Males(year)),ot_F_sel)*M_agelen_comm;
    break;
   }  //case loop
   switch (series)
     {
     case 1:
     if(ssex==1)
      {
       pred_catlen_ll_short_female(obs)=pred_catlen_comm(obs)(fminbins,fmaxbins);
       obs_catlen_ll_short_female(obs)=obs_catlen_comm(obs)(fminbins,fmaxbins);
       // Turn into proportions
       pred_catlen_ll_short_female(obs)=pred_catlen_ll_short_female(obs)/sum(pred_catlen_ll_short_female(obs));
       obs_catlen_ll_short_female(obs)=obs_catlen_ll_short_female(obs)/sum(obs_catlen_ll_short_female(obs));
       pred_catlenP_comm(obs)(fminbins,fmaxbins)=pred_catlen_ll_short_female(obs)(fminbins,fmaxbins);
      }
     if(ssex==2)
      {
       pred_catlen_ll_short_male(obs)=pred_catlen_comm(obs)(mminbins,mmaxbins);
       obs_catlen_ll_short_male(obs)=obs_catlen_comm(obs)(mminbins,mmaxbins);
       // Turn into proportions
       pred_catlen_ll_short_male(obs)=pred_catlen_ll_short_male(obs)/sum(pred_catlen_ll_short_male(obs));
       obs_catlen_ll_short_male(obs)=obs_catlen_ll_short_male(obs)/sum(obs_catlen_ll_short_male(obs));
       pred_catlenP_comm(obs)(mminbins,mmaxbins)=pred_catlen_ll_short_male(obs)(mminbins,mmaxbins);
      }
  break;
    case 2:
       pred_catlen_ot_short(obs)=pred_catlen_comm(obs)(fminbins3,fmaxbins3);
       obs_catlen_ot_short(obs)=obs_catlen_comm(obs)(fminbins3,fmaxbins3);
       // Turn into proportions
       pred_catlen_ot_short(obs)=pred_catlen_ot_short(obs)/sum(pred_catlen_ot_short(obs));
       obs_catlen_ot_short(obs)=obs_catlen_ot_short(obs)/sum(obs_catlen_ot_short(obs));
       pred_catlenP_comm(obs)(fminbins3,fmaxbins3)=pred_catlen_ot_short(obs)(fminbins3,fmaxbins3);
    break;
    }  //case loop
  data_int n_short;
  n_short=0;
  switch (series)
  {
  case 1:
  if(ssex==1)
   {
   n_short=fmaxbins-fminbins+1;
   catlen_comm_like(series)+=Multinomial(obs_catlen_ll_short_female(obs)(fminbins,fmaxbins), pred_catlen_ll_short_female(obs)(fminbins,fmaxbins),n_short,maxs,catlen_LikeType(series),catlen_comm_resid(obs)(fminbins,fmaxbins));
   }
  if(ssex==2)
   {
   n_short=mmaxbins-mminbins+1;
   catlen_comm_like(series)+=Multinomial(obs_catlen_ll_short_male(obs)(mminbins,mmaxbins),pred_catlen_ll_short_male(obs)(mminbins,mmaxbins),n_short,maxs,catlen_LikeType(series),catlen_comm_resid(obs)(mminbins,mmaxbins));
   }
  break;
  case 2:
   n_short=fmaxbins3-fminbins3+1;
   catlen_comm_like(series)+=Multinomial(obs_catlen_ot_short(obs)(fminbins3,fmaxbins3), pred_catlen_ot_short(obs)(fminbins3,fmaxbins3),n_short,maxs,catlen_LikeType(series),catlen_comm_resid(obs)(fminbins3,fmaxbins3));  
  break;
  }  // close case loop
  }  // close obs loop
  fcomp(9)=sum(catlen_comm_like);
}

void model_parameters::Catlen_Survey_Likelihood(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  pred_catlen_rv.initialize();
  catlen_survey_like=0.;
  int ssex;
  for(obs=1;obs<=N_catlen_rv;obs++)
  {
  series = obs_catlen_rv(obs,-2);
  ssex = obs_catlen_rv(obs,-0);
  if(obs_catlen_rv(obs,-3)>maxssizes)		
  {
  maxs = maxssizes;
  }
  else
  {
  maxs = obs_catlen_rv(obs,-3);
  }
  year=obs_catlen_rv(obs,-1);
  //1=summer_rv
  //7=halibut_fixed
  switch (series)
  {
  case 1:
	if(ssex==1)
	{
    pred_catlen_rv(obs) =elem_prod(elem_prod(exp(-(F_utotal(year)+M)*.5),Females(year)),old_summer_rv_F_sel)*F_agelen;
	}
	else
	{
    pred_catlen_rv(obs) =elem_prod(elem_prod(exp(-(M_utotal(year)+M)*.5),Males(year)),old_summer_rv_F_sel)*M_agelen;
	}
  break;
   case 7:
	if(ssex==1)
	{
	pred_catlen_rv(obs) =elem_prod(elem_prod(exp(-(F_utotal(year)+M)*.5),Females(year)),halibut_fixed_F_sel)*F_agelen;
	}
	else
	{
	pred_catlen_rv(obs) =elem_prod(elem_prod(exp(-(M_utotal(year)+M)*.5),Males(year)),halibut_fixed_M_sel)*M_agelen;
	}
    break;
  }
  switch (series)
  {
  case 1:
     if(ssex==1)
      {
       pred_catlen_rv_short(obs)=pred_catlen_rv(obs)(fminbins2,fmaxbins2);
       obs_catlen_rv_short(obs)=obs_catlen_rv(obs)(fminbins2,fmaxbins2);
       // Turn into proportions
       pred_catlen_rv_short(obs)=pred_catlen_rv_short(obs)/sum(pred_catlen_rv_short(obs));
       obs_catlen_rv_short(obs)=obs_catlen_rv_short(obs)/sum(obs_catlen_rv_short(obs));
       pred_catlenP_rv(obs)(fminbins2,fmaxbins2)=pred_catlen_rv_short(obs)(fminbins2,fmaxbins2);
      }
     if(ssex==2)
      {
       pred_catlen_rv_short(obs)=pred_catlen_rv(obs)(mminbins2,mmaxbins2);
       obs_catlen_rv_short(obs)=obs_catlen_rv(obs)(mminbins2,mmaxbins2);
       // Turn into proportions
       pred_catlen_rv_short(obs)=pred_catlen_rv_short(obs)/sum(pred_catlen_rv_short(obs));
       obs_catlen_rv_short(obs)=obs_catlen_rv_short(obs)/sum(obs_catlen_rv_short(obs));
       pred_catlenP_rv(obs)(mminbins2,mmaxbins2)=pred_catlen_rv_short(obs)(mminbins2,mmaxbins2);
      }
  break;
  case 7:
     if(ssex==1)
      {
       newmax=fmaxbins;
       if(fmaxbins>72){newmax=72;}
       pred_catlen_hsurvey_short_female(obs)(fminbins,newmax)=pred_catlen_rv(obs)(fminbins,newmax);
       obs_catlen_hsurvey_short_female(obs)(fminbins,newmax)=obs_catlen_rv(obs)(fminbins,newmax);
       // Turn into proportions
       pred_catlen_hsurvey_short_female(obs)=pred_catlen_hsurvey_short_female(obs)/sum(pred_catlen_hsurvey_short_female(obs));
       obs_catlen_hsurvey_short_female(obs)=obs_catlen_hsurvey_short_female(obs)/sum(obs_catlen_hsurvey_short_female(obs));
       pred_catlenP_rv(obs)(fminbins,fmaxbins)=pred_catlen_hsurvey_short_female(obs)(fminbins,fmaxbins);
      }
     if(ssex==2)
      {
       newmax=mmaxbins;
       if(mmaxbins>72){newmax=72;}
       pred_catlen_hsurvey_short_male(obs)(mminbins,newmax)=pred_catlen_rv(obs)(mminbins,newmax);
       obs_catlen_hsurvey_short_male(obs)(mminbins,newmax)=obs_catlen_rv(obs)(mminbins,newmax);
       // Turn into proportions
       pred_catlen_hsurvey_short_male(obs)=pred_catlen_hsurvey_short_male(obs)/sum(pred_catlen_hsurvey_short_male(obs));
       obs_catlen_hsurvey_short_male(obs)=obs_catlen_hsurvey_short_male(obs)/sum(obs_catlen_hsurvey_short_male(obs));
       pred_catlenP_rv(obs)(mminbins,mmaxbins)=pred_catlen_hsurvey_short_male(obs)(mminbins,mmaxbins);
      }
  break;
  }
  int ix;
  ix=series;
  if(series>3) ix=series-2;
  if(like_weight(ix)>=1)
  {
  data_int n_short;
  n_short=0;
  switch (series)
  {
  case 1:
  if(ssex==1)
   {
   n_short=fmaxbins2-fminbins2+1;
   catlen_survey_like(series)+=Multinomial(obs_catlen_rv_short(obs)(fminbins2,fmaxbins2), pred_catlen_rv_short(obs)(fminbins2,fmaxbins2),n_short,maxs,catlen_LikeType(series),catlen_survey_resid(obs)(fminbins2,fmaxbins2));
   }
  if(ssex==2)
   {
   n_short=mmaxbins2-mminbins2+1;
   catlen_survey_like(series)+=Multinomial(obs_catlen_rv_short(obs)(mminbins2,mmaxbins2),pred_catlen_rv_short(obs)(mminbins2,mmaxbins2),n_short,maxs,catlen_LikeType(series),catlen_survey_resid(obs)(mminbins2,mmaxbins2));
   }
  break;
  case 7:
  if(ssex==1)
   {
   newmax=fmaxbins;
   if(fmaxbins>72){newmax=72;}else{newmax=fmaxbins;}
   n_short=newmax-fminbins+1;
   catlen_survey_like(series)+=Multinomial(obs_catlen_hsurvey_short_female(obs)(fminbins,newmax),
   		pred_catlen_hsurvey_short_female(obs)(fminbins,newmax),n_short,maxs,catlen_LikeType(series),catlen_survey_resid(obs)(fminbins,newmax));
   }
  if(ssex==2)
   {
   newmax=mmaxbins;
   if(mmaxbins>72){newmax=72;}else{newmax=mmaxbins;}
   n_short=newmax-mminbins+1;
   catlen_survey_like(series)+=Multinomial(obs_catlen_hsurvey_short_male(obs)(mminbins,newmax),
	   pred_catlen_hsurvey_short_male(obs)(mminbins,newmax),n_short,maxs,catlen_LikeType(series),catlen_survey_resid(obs)(mminbins,newmax));
   }
  break;
  }
  }//close series
  }//close obs
  fcomp(10)=sum(catlen_survey_like);
}

void model_parameters::MCWrite(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  report1 
    <<" old_summer_rv_Sfull "<<old_summer_rv_Sfull
    <<" old_summer_rv_varRest "<<old_summer_rv_varRest
    <<" halibut_fixed_Sfull_F "<<halibut_fixed_Sfull_F
    <<" halibut_fixed_Sfull_M "<<halibut_fixed_Sfull_M
    <<" ll_Sfull_F "<<ll_Sfull_F
    <<" ll_Sfull_M "<<ll_Sfull_M
    <<" ot_Sfull "<<ot_Sfull
    <<" ot_varRest "<<ot_varRest
    <<" log_q_old_summer_rv "<<log_q_old_summer_rv
    <<" log_q_halibut_fixed "<<log_q_halibut_fixed
    <<" SSB2009 "<<SSB2009                        
    <<" FF2009 "<<FF2009                        
    <<" FM2009 "<<FM2009                        
    <<" naa2009 "<<naa2009                        
  <<" ofv "<<ofv
  <<endl;
}

dvariable model_parameters::GetF_forward(const dvariable& N, const dvariable& C, const dvariable& M)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  dvariable Ctmp=C;
  dvariable fest=Ctmp/N;
  double minsurv=0.01;
  if(fest> 0.99)
  {
    Ctmp =N*(1-posfun(1-fest,minsurv,ffpen));
    cout<<"posfun for "<<fest<<" "<<C<<" "<<N<<" "<<ffpen<<endl;
    fest=Ctmp/N;
  }
  dvariable ztmp;
  dvariable Pcat;
  dvariable df;
  for (int iter=1;iter<=10;iter++)
  {
    ztmp=fest+M;
    Pcat=fest/ztmp*(1-mfexp(-ztmp))*N;
    df=N/ztmp*(mfexp(-ztmp)*(-1+fest/ztmp+fest)+1-fest/ztmp);
    fest+= -1*(Pcat-Ctmp)/df;
  }
  return(fest);
}

dvariable model_parameters::ffromcn(const dvariable& N, const dvariable& C, const dvariable& M)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  dvariable Ctmp=C;
  dvariable fest=M;
  dvariable ztmp;
  for (int iter=1;iter<=10;iter++)
  {
    ztmp=fest+M;
    fest=Ctmp*ztmp/N*(1-exp(-1*ztmp));
  }
  return(fest);
}

dvar_vector model_parameters::Make_select_dhn(int minage, int maxage, dvariable full, dvariable varL, dvariable varR)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  dvar_vector selection(minage,maxage);
  selection.initialize();
  for (age=minage;age<=maxage;age++)
  {
	if (full>age)
	{
		selection(age)=exp(-1.0*(square(age-full))/varL);             
	}
	else
	{
		selection(age)=exp(-1.0*(square(age-full))/varR);             
	}
  }	
  selection = selection/max(selection);
  return(selection);
}

dvar_vector model_parameters::Make_select_logistic(int minage, int maxage, dvariable full, dvariable varL)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  dvar_vector selection(minage,maxage);
  selection.initialize();
  for (age=minage;age<=maxage;age++)
  {
  selection(age)=1/(1+exp(-varL*(age-full))); 
  }	
  selection = selection/max(selection);
  return(selection);
}

void model_parameters::Age1_pen(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  dvar_vector junk(1,summer_rv_n_years-2);
  dvar_vector R(1,summer_rv_n_years-2);
    for(i=1;i<=summer_rv_n_years-2;i++)
     {
     R(i)=posfun(2.71828*p2/p1*(SSB(summer_rv_years(i))-p3)*exp(-(SSB(summer_rv_years(i))-p3)/p1)*biasRannual,10,r_penF);
     // both sexes
     }
    for(i=1;i<=summer_rv_n_years-2;i++)
     {
     junk(i)=rec(summer_rv_years(i+1));   
     }
      age1_pen=LogNormal_LogLike(junk,R,summer_rv_n_years-2,.6,2);  
}

void model_parameters::Age1_last_year_pdf(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  dvariable age1_last_year;
  double age1_last_year_std;
  age1_last_year=0;
  age1_last_year_pdf=0;
  for(year=2005;year<=2007;year++)
    {
    age1_last_year+=(rec(year));
    }
    age1_last_year=log(age1_last_year/3);
  age1_last_year_std=.1;
  age1_last_year_pdf=Normal(log(rec(last_year)+1),age1_last_year,age1_last_year_std*1);
}

void model_parameters::Recdiff_pdf(void)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  recdiff_pdf=0;
  recdiff_sigma=.5;
  dvar_vector recdiff(first_year,last_year-1);
  for(year=first_year;year<=last_year-1;year++)
    {
    recdiff(year)=log(rec(year)/rec(year+1));
    }
  recdiff_pdf=norm2(recdiff/pow(recdiff_sigma,2));
}

dvariable model_parameters::LogNormal_LogLike(dvar_vector observed, dvar_vector predicted, double len,double sigma, double type)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
   dvariable neglogLike;
   dvariable sumlogr;
   dvariable len2;
   dvariable sumsqrlogratio;
   double pi=3.1415926535897;
   sumlogr=0;
   sumsqrlogratio=0;
   len2=0;    
   for(i=1;i<=len;i++)
     {
       if(observed(i)>-0.000001)
        {
         sumlogr+=log(observed(i));
         sumsqrlogratio+=square(log(observed(i))-log(predicted(i)));
         len2+=1;
        }
     }
     if(type==2)
     {
     neglogLike = 1/(2*square(sigma))*sumsqrlogratio   + sumlogr   + len2*log(sigma)*sqrt(2*pi);  //specified sigma used
     }
    if(type==1 &last_phase())
    {
    neglogLike = 0.5*len2*log( (2*pi/len2)*sumsqrlogratio) + sumlogr + len2/2;  //est. sigma used
    }
  return(neglogLike);
}

dvariable model_parameters::LogNormal_LogLike(dvar_vector observed, dvar_vector predicted, double len,double sigma, double type,dvar_vector resid)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
   dvariable neglogLike;
   dvariable sumlogr;
   dvariable len2;
   dvariable sumsqrlogratio;
   double pi=3.1415926535897;
   if(type==1)
   {
     cout<<" code not designed for type 1 and calculate resids "<<endl;
   }
   else
   {
     sumlogr=0;
     sumsqrlogratio=0;
     len2=0;    
     for(i=1;i<=len;i++)
     {
       if(observed(i)>-0.000001)
        {
         sumlogr+=log(observed(i));
         sumsqrlogratio+=square(log(observed(i))-log(predicted(i)));
         len2+=1;
         resid(i)=(log(observed(i))-log(predicted(i)))/sigma;
        }
     }
     neglogLike = 1/(2*square(sigma))*sumsqrlogratio   + sumlogr   + len2*log(sigma)*sqrt(2*pi);  //specified sigma used
    }
  return(neglogLike);
}

dvariable model_parameters::calc_sigma(dvar_vector observed, dvar_vector predicted, double len)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
   dvariable sigma;
   dvariable sumsqrlogratio;
   dvariable len2;
   sumsqrlogratio=0;
   len2=0; 
   for(i=1;i<=len;i++)
     {
       if(observed(i)>-0.000001)
        {
         sumsqrlogratio+=square(log(observed(i))-log(predicted(i)));
         len2+=1;
        }
     }
   sigma = sqrt(sumsqrlogratio/len2);
   return(sigma);
}

dvariable model_parameters::Normal(dvariable observed,dvariable mean,double nsigma)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
   dvariable normal;
   double pi=3.1415926535897;
   dvariable neglogpdf;
   //pdf=(1/(sqrt(2*pi)*nsigma)*exp(-(pow(observed-mean,2)/(2*pow(nsigma,2)))));
   //neglogpdf=-log(pdf);
    neglogpdf=log(sqrt(2*pi)*nsigma)+pow(observed-mean,2)/(2*pow(nsigma,2));
   return(neglogpdf);
}

dvariable model_parameters::Binomial(dvariable observed,double x,double n)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
   dvariable p;
   dvariable pdf;
   dvariable neglogpdf;
   p=observed;				
   pdf=pow(p,x)*pow(1-p,n-x);
   neglogpdf=-log(pdf);
   return(neglogpdf);
}

dvariable model_parameters::negfun(const dvariable&x,const double eps,dvariable& pen)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  {
  if (x<=eps) {
  return x;
  } else {
  pen+=.01*square(x-eps);
  cout<<" neg fun "<<x<<" "<<eps/(2-x/eps)<<" "<<pen<<endl;
  return x/(2-eps/x);
  }
  }
}

dvariable model_parameters::Multinomial(dvar_vector obs,dvar_vector pred,data_int& a,double& b,int& Ltype,dvar_vector resid)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  {
   dvariable log_likelihood;
     resid=elem_div((obs-pred),sqrt(elem_prod(pred,(1.-pred))/b));
      log_likelihood=-1.* sum(elem_prod(b*obs,log(pred+1e-6))) + sum(elem_prod(b*obs,log(obs+1e-6)));
    return(log_likelihood);
   }
}

dvariable model_parameters::Multifan(dvector& obs,dvar_vector& pred,int a,double& b,int Ltype)
{
  ofstream& report1= *pad_report1;
  ofstream& tmp= *pad_tmp;
  {
  if (obs.indexmin() != pred.indexmin() || obs.indexmax() != pred.indexmax() )
  {
    cerr << "Index limits on observed vector are not equal to the Index\n"
      "limits on the predicted vector in robust_p function\n";
  }
  /* Routine to do robust likelihood computation for proportions 
    a = number of ages
    b = tau  effective sample size
  */
  RETURN_ARRAYS_INCREMENT(); //Need this statement because the function
              //returns a variable type
   dvariable log_likelihood;
   if(last_phase())
   {
    dvar_vector v(obs.indexmin(),obs.indexmax());
     switch(Ltype)
     {
     case 7: // Original multifan
     v =   (elem_prod(pred ,1.  - pred )+.1/a);
     break;
     case 8: // Using the observed -- From Maunder and Starr
     v =   (elem_prod(obs ,1.  - obs )+.1/a);
     break;
     }
     dvar_vector l  = exp(-0.5*b* elem_div(square(pred - obs), v ));
     log_likelihood = -1.0*sum(log(l + .01));
     log_likelihood  += 0.5 * sum(log(v));
   } 
   else
   {
     log_likelihood=-1.* sum(elem_prod(b*obs,log(pred+1e-6))) + sum(elem_prod(b*obs,log(obs+1e-6)));
   }
   RETURN_ARRAYS_DECREMENT(); // Need this to decrement the stack increment
                              // caused by RETURN_ARRAYB_INCREMENT();
  return(log_likelihood);
  }
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
    REPORT(first_year);
    REPORT(last_year);
    REPORT(first_age);
    REPORT(last_age);
    REPORT(N_lens);  
    REPORT(N_lens_comm);  
    REPORT(lens);  
    REPORT(lens_comm);  
    REPORT(log_q_old_summer_rv);
    REPORT(log_q_halibut_fixed);
    REPORT(N1);
    REPORT(start_Z);
    REPORT(age1);
    REPORT(R0);
    REPORT(receq);
    REPORT(alpha);
    REPORT(beta);
    REPORT(p2);
    REPORT(p1);
    REPORT(p3);
    REPORT(M);
    REPORT(S0);
    REPORT(lenatage);
    REPORT(w);
    REPORT(w_F);
    REPORT(mat);
    REPORT(SSB);
    REPORT(SSB2);
    REPORT(Total_N);
    REPORT(SSN);
    REPORT(SSN2);
    REPORT(Total_F);
    REPORT(Total_M);
    REPORT(Total_W);
    REPORT(rec);
    REPORT(all_log_RecDev);
    REPORT(F_sigma);
    REPORT(M_sigma);
    REPORT(F_agelen);
    REPORT(M_agelen);
    REPORT(F_agelen_comm);
    REPORT(M_agelen_comm);
    REPORT(F_lenage);
    REPORT(M_lenage);
    REPORT(F_lenage_comm);
    REPORT(M_lenage_comm);
    REPORT(out_obs_catlen_rv);	//proportions
    REPORT(obs_catlen_rv);		//raw
    REPORT(pred_catlen_rv_short);	//proportions
    REPORT(pred_catlenP_rv);	//proportions
    REPORT(pred_catlen_hsurvey_short_female);	//proportions
    REPORT(pred_catlen_hsurvey_short_male);	//proportions
    REPORT(pred_catlen_hCI_short_female);	//proportions
    REPORT(pred_catlen_hCI_short_male);	//proportions
    REPORT(out_obs_catlen_comm);
    REPORT(obs_catlen_comm);	//raw
    REPORT(pred_catlen_ll_short_female);
    REPORT(pred_catlen_ll_short_male);
    REPORT(pred_catlenP_comm);	//proportions
    REPORT(pred_catlen_comm);	//raw
    REPORT(pred_catlen_ot_short);
    REPORT(Females);
    REPORT(Males);
    REPORT(summer_RVage);
    REPORT(halibut_fixed_age);
    REPORT(summer_rv_total);
    REPORT(summer_rv_years);
    REPORT(pred_summer_rv_total);
    REPORT(halibut_fixed_years);
    REPORT(halibut_fixed_total);
    REPORT(pred_halibut_fixed_total);
    REPORT(comm_total);
    REPORT(comm_total_years);
    REPORT(F_utotal);
    REPORT(M_utotal);
    REPORT(old_summer_rv_F_sel);
    REPORT(halibut_fixed_F_sel);
    REPORT(halibut_fixed_M_sel);
    REPORT(ll_M_sel);
    REPORT(ll_F_sel);
    REPORT(ot_F_sel);
    REPORT(ages);
    REPORT(u_Canada);
    REPORT(u_Foreign);
    REPORT(C_ll_F);
    REPORT(C_ll_M);
    REPORT(C_ot);
    REPORT(VB_ll);
    REPORT(VB_ot);
    REPORT(VB_Canada);
    REPORT(TB_Canada);
    REPORT(u_ll_F);
    REPORT(u_ll_M);
    REPORT(u_ot);
    REPORT(Canada);
    REPORT(Foreign);
    REPORT(ssbeq);
    REPORT(old_summer_rv_Sfull);
    REPORT(old_summer_rv_varLest);
    REPORT(old_summer_rv_varRest);
    REPORT(	ll_Sfull_F);
    REPORT(ll_varLest_F);
    REPORT(	ll_Sfull_M);
    REPORT(ll_varLest_M);
    REPORT(ot_Sfull);
    REPORT(ot_varLest);
    REPORT(ot_varRest);
    REPORT(summer_rv_like);
    REPORT(halibut_fixed_like);
    REPORT(sum(catlen_survey_like));
    REPORT(sum(catlen_comm_like));
    REPORT(fpen);
    REPORT(r_penF);
    REPORT(reccheck);
    REPORT(recpen);
    REPORT(recpen2);
    REPORT(recdiff_pdf);
    REPORT(u_ll_pdf);
    REPORT(dudfem);
    REPORT(SPRF0);
    REPORT(ssbeq);
    REPORT(receq);
    REPORT(ofv);
    REPORT(elem_prod(fcomp,like_weight));
    REPORT(summer_rv_resid);
    REPORT(halibut_fixed_resid);
   // if (last_phase())
   // {
   // for (ii=1; ii<=N_catlen_rv; ii++)
   //   report<<ii<<" "<<obs_catlen_rv(ii,-3)<<" "<<obs_catlen_rv(ii,-2)<<" "<<obs_catlen_rv(ii,-1)<<" "<<obs_catlen_rv(ii,0)<<" "<<catlen_survey_resid(ii)<<endl;
   // for ( ii=1; ii<=N_catlen_comm; ii++)
   //   report<<ii<<" "<<obs_catlen_comm(ii,-3)<<" "<<obs_catlen_comm(ii,-2)<<" "<<obs_catlen_comm(ii,-1)<<" "<<obs_catlen_comm(ii,0)<<" "<<catlen_comm_resid(ii)<<endl;
   // }
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
{
  delete pad_report1;
  pad_report1 = NULL;
  delete pad_tmp;
  pad_tmp = NULL;
}

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
  arrmblsize=10000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(4000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000000);
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
