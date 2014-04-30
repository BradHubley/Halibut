// Halibut (Sept 2010) Kurtis Trzcinski
// Halibut (Sept 2014) Brad Hubley
 
// ./vpop10 -mcmc 120000 -mcsave 100 // del first 200 rows in mcout.dat for burn in.
// ./vpop10 -mcmc 1200000 -mcsave 1000 // del first 200 rows in mcout.dat for burn in.

TOP_OF_MAIN_SECTION
  time(&start);
  arrmblsize=10000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(4000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000000);
//  gradient_structure::set_MAX_NVAR_OFFSET(999);

//=====================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++

DATA_SECTION
   !!CLASS ofstream report1("mcout.dat",ios::app);  // you must delete this file before a new mceval run
   !!CLASS ofstream tmp("tmp.dat",ios::app);  // you must delete this file before a new mceval run

//  INPUT DATA FILE 
  init_int debug	// for checking the reading in of the data
  !!if (debug==1){char mark; cout<<"starting debug /n"; cin>>mark; }
  // ** model dimentioning parameters **
  init_int first_year
  init_int last_year
  init_int first_age	// First age for model (recruitment is at this age)
  init_int last_age	// 
  init_int first_length
  init_int length_inc	// size of the length interval
  init_int N_lens 	// number of length intervals
  init_int N_sexes  //1=female, 2=male

  init_int Flength  // First length-class to include in cll fitting
  init_int fminbins  // Use the sub-vector operation for creating a structure that is restricted to the + and - length bin
  init_int fmaxbins  // for LL
  init_int mminbins  // 
  init_int mmaxbins  // 
  init_int fminbins2  // for RV
  init_int fmaxbins2  // 
  init_int mminbins2  // 
  init_int mmaxbins2  // 
  init_int fminbins3  // for OT
  init_int fmaxbins3  // 
  init_int mminbins3  // 
  init_int mmaxbins3  // 
  init_int Frecyr
  init_number Minit
  init_int Mopt // 1=constant over ages, 2 = Mjuv (age1) different, 3 = Mmat different, 4 = both different 
  init_int N1opt //1=on,0=off (equilibrium unfished = starting pop)
//  !!if (debug==1){char mark; cout<<"N1opt "<<N1opt<<endl; cin>>mark; }
  // ** Priors **       
  //Phase -1 is not estimate i.e. fixed
  //Priortype 0(or non defined value)=uniform, 1=normal 2=lognormal
  //phase LB UB Priortype Mean cv defaultvalue
  init_vector N1_prior(1,7)
//  init_vector age1(first_year,last_year)
  init_vector year1(first_age,last_age)
  init_vector start_Z_prior(1,7)
  init_vector age1_prior(1,7)
  init_vector R0_prior(1,7)
  init_vector beta_prior(1,7)
  init_vector alpha_prior(1,7) 
  init_vector p1_prior(1,7)
  init_vector p2_prior(1,7)
  init_vector p3_prior(1,7)
  init_vector log_RecDev_prior(1,7)
  init_vector M_prior(1,7) 

  //survey selectivities
  init_int ll_select_switch

  init_vector old_summer_rv_Sfull_prior(1,7) 				// Summer_rv
  init_vector old_summer_rv_varLest_prior(1,7) 
  init_vector old_summer_rv_varRest_prior(1,7) 

  init_matrix halibut_fixed_Sfull_prior(1,N_sexes,1,7) 		// halibut survey fixed station
  init_matrix halibut_fixed_varLest_prior(1,N_sexes,1,7) 

  //commercial selectivities:
  init_matrix ll_Sfull_prior(1,N_sexes,1,7) 				// LL
  init_matrix ll_varLest_prior(1,N_sexes,1,7) 
  init_vector ot_Sfull_prior(1,7) 							// OT
  init_vector ot_varLest_prior(1,7) 
  init_vector ot_varRest_prior(1,7) 


  //survey q's
  init_vector log_q_old_summer_rv_prior(1,7)
  init_vector log_q_halibut_fixed_prior(1,7)
  init_vector u_ll_prior(1,7)

  // ** Likelihood **
  init_vector like_weight(1,12)
  init_int N_catlen_commercial
  init_int N_catlen_surveys
  init_ivector catlen_LikeType(1,N_catlen_surveys) 
  init_int maxssizes 

  // ** fixed parameters **
  // weight@length: 
  init_vector bi(1,N_sexes)
  init_vector bii(1,N_sexes)
  init_vector linf(1,N_sexes) 
  init_vector k(1,N_sexes)
  init_vector t0(1,N_sexes)
  
  //  Relationship between mean length and std
  init_vector sigma_a(1,N_sexes)
  init_vector sigma_b(1,N_sexes)
  init_vector matage(1,N_sexes) // maturity age for splitting for the CPUE
  init_vector mat_a(1,N_sexes)  // a param for logistic maturity ogive
  init_vector mat_b(1,N_sexes) // b param for logistic maturity ogive

//####################################################################################
  // ** DATA **

  //summer rv
  init_int summer_rv_n_years
!!if (debug==1){char mark; cout<<"summer_rv_n_years "<<summer_rv_n_years<<endl; cin>>mark; }
  init_vector summer_rv_years(1,summer_rv_n_years)
  init_vector summer_rv_total(1,summer_rv_n_years)
  //init_vector summer_rv_matureF(1,summer_rv_n_years)

  init_int halibut_fixed_n_years
  init_vector halibut_fixed_years(1,halibut_fixed_n_years)
  init_vector halibut_fixed_total(1,halibut_fixed_n_years)
  init_vector halibut_fixed_matureF(1,halibut_fixed_n_years)
!!if (debug==1){char mark; cout<<"halibut_fixed_total "<<halibut_fixed_total<<endl; cin>>mark; }

//########################################################################3
  //rv CATCH AT LENGTH
 
  init_int N_catlen_rv // catch-at-length observations - total for all gears and years
//!!if (debug==1){char mark; cout<<"N_catlen_rv "<<N_catlen_rv<<endl; cin>>mark; }
  init_matrix obs_catlen_rv(1,N_catlen_rv,-3,N_lens) //year sex 16 19............. 
!!if (debug==1){char mark; cout<<"obs_catlen_rv "<<obs_catlen_rv(106)<<endl; cin>>mark; }

//############################################################################################
  //COMMERCIAL DATA - total landings by country (mt) -999 are missing values, presumably low enough to treat as zeros - see tpl
  init_vector Canada(1960,last_year)
//  !!if (debug==1){char mark; cout<<"Canada "<<Canada<<endl; cin>>mark; }
  init_vector Foreign(1960,last_year)
  init_vector Total(1960,last_year)

  //All nations landings by gear 
  init_vector ll_F(1960,last_year)			// divided up based on sex ratio in the catch from 1988 to present. mean prportion female used from 1960 to 1987
  init_vector ll_M(1960,last_year)
  init_vector ot(1960,last_year)

  init_int comm_total_n_years
  init_vector comm_total_years(1,comm_total_n_years)
  init_vector comm_total(1,comm_total_n_years)

  //##CANADIAN+Foreign CAL 
  init_int N_lens_comm   
  init_int N_catlen_comm // catch-at-length observations - total for all gears and years
  !!if (debug==1){char mark; cout<<"N_catlen_comm "<<N_catlen_comm<<endl; cin>>mark; }
  init_matrix obs_catlen_comm(1,N_catlen_comm,-3,N_lens_comm) 
  !!if (debug==1){char mark; cout<<"obs_catlen_comm "<<obs_catlen_comm(66)<<endl; cin>>mark; }
 
  //PVA inputs
  init_int pva_horizon // length in years
  init_int length_mort_vec
  init_vector mort_vec(1,length_mort_vec) //- incidental mortality inputs
 
  init_matrix lenatage(1,N_sexes,first_age,last_age)		//mean length at age
  init_vector F_sigma(first_age,last_age)		        //sigma len at age
  !!if (debug==1){char mark; cout<<"F_sigma "<<F_sigma<<endl; cin>>mark; }
  init_vector M_sigma(first_age,last_age)

  init_int dud
  !!if (debug==1){char mark; cout<<"What number are we thinking of? "<<dud<<endl; cin>>mark; }

 
  // ** other parameters **
  matrix w(1,N_sexes,first_age,last_age) // weight-at-age

  matrix obs_catlenP_rv(1,N_catlen_rv,1,N_lens) // proportions 
  matrix out_obs_catlen_rv(1,N_catlen_rv,-3,N_lens) 
  matrix obs_catlenP_comm(1,N_catlen_comm,1,N_lens_comm) // proportions 
  matrix out_obs_catlen_comm(1,N_catlen_comm,-3,N_lens_comm) 

  matrix mat(1,N_sexes,first_age,last_age)
  number pi
  number lb
  number ub
  number annualsigma
  number biasRannual
  number maxs

  !!annualsigma=log_RecDev_prior(6);  // cv for recruitment deviates
  !!biasRannual=mfexp(-0.5*square(annualsigma));

  int age
  int sex
  int year
  int i
  int j
  int phz
  int obs
  int series
  int newmax

  int srows
  !!srows=last_year+2-first_year;
  matrix F_agelen(first_age,last_age,1,N_lens)
  matrix M_agelen(first_age,last_age,1,N_lens)
  matrix F_agelen_comm(first_age,last_age,1,N_lens_comm)
  matrix M_agelen_comm(first_age,last_age,1,N_lens_comm)
  matrix F_lenage(first_age,last_age,1,N_lens)
  matrix M_lenage(first_age,last_age,1,N_lens)
  matrix F_lenage_comm(first_age,last_age,1,N_lens_comm)
  matrix M_lenage_comm(first_age,last_age,1,N_lens_comm)
  vector lens(1,N_lens)
  !!lens.fill_seqadd(first_length,length_inc);
  vector lens_comm(1,N_lens_comm)
  !!lens_comm.fill_seqadd(first_length,length_inc);
  vector ages(first_age,last_age)
  !!ages.fill_seqadd(first_age,1);
  vector years(first_year,last_year+1)
  !!years.fill_seqadd(first_year,1);

//============================================================================
 LOC_CALCS
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

 END_CALCS

//============================================================================
PARAMETER_SECTION

  number dummy;

 !!phz=N1_prior(1);lb=N1_prior(2);ub=N1_prior(3);
  init_bounded_number N1(lb,ub,phz)

  !!phz=start_Z_prior(1);lb=start_Z_prior(2);ub=start_Z_prior(3);
  init_bounded_number start_Z(lb,ub,phz)

  !!phz=age1_prior(1);lb=age1_prior(2);ub=age1_prior(3);
  init_bounded_vector age1(first_year,last_year-2,lb,ub,phz)

  !!phz=R0_prior(1);lb=R0_prior(2);ub=R0_prior(3);
  init_bounded_number R0(lb,ub,phz)

  !!phz=beta_prior(1);lb=beta_prior(2);ub=beta_prior(3);
  init_bounded_number beta(lb,ub,phz)

  !!phz=alpha_prior(1);lb=alpha_prior(2);ub=alpha_prior(3);
  init_bounded_number alpha(lb,ub,phz)

  !!phz=p1_prior(1);lb=p1_prior(2);ub=p1_prior(3);
  init_bounded_number p1(lb,ub,phz)

  !!phz=p2_prior(1);lb=p2_prior(2);ub=p2_prior(3);
  init_bounded_number p2(lb,ub,phz)

  !!phz=p3_prior(1);lb=p3_prior(2);ub=p3_prior(3);
  init_bounded_number p3(lb,ub,phz)

  !!phz=M_prior(1);lb=M_prior(2);ub=M_prior(3);
  init_bounded_number M(lb,ub,phz)

  !!phz=log_RecDev_prior(1);lb=log_RecDev_prior(2);ub=log_RecDev_prior(3);
  init_bounded_vector log_RecDev(Frecyr,last_year,lb,ub,phz)  //Frecyr
  vector all_log_RecDev(first_year,last_year+1)

  !!phz=old_summer_rv_Sfull_prior(1);lb=old_summer_rv_Sfull_prior(2);ub=old_summer_rv_Sfull_prior(3);
  init_bounded_number old_summer_rv_Sfull(lb,ub,phz)
  !!phz=old_summer_rv_varLest_prior(1);lb=old_summer_rv_varLest_prior(2);ub=old_summer_rv_varLest_prior(3);
  init_bounded_number old_summer_rv_varLest(lb,ub,phz)
  !!phz=old_summer_rv_varRest_prior(1);lb=old_summer_rv_varRest_prior(2);ub=old_summer_rv_varRest_prior(3);
  init_bounded_number old_summer_rv_varRest(lb,ub,phz)

  !!phz=halibut_fixed_Sfull_prior(1,1);lb=halibut_fixed_Sfull_prior(1,2);ub=halibut_fixed_Sfull_prior(1,3);
  init_bounded_number halibut_fixed_Sfull_F(lb,ub,phz)
  !!phz=halibut_fixed_varLest_prior(1,1);lb=halibut_fixed_varLest_prior(1,2);ub=halibut_fixed_varLest_prior(1,3);
  init_bounded_number halibut_fixed_varLest_F(lb,ub,phz)

  !!phz=halibut_fixed_Sfull_prior(2,1);lb=halibut_fixed_Sfull_prior(2,2);ub=halibut_fixed_Sfull_prior(2,3);
  init_bounded_number halibut_fixed_Sfull_M(lb,ub,phz)
  !!phz=halibut_fixed_varLest_prior(2,1);lb=halibut_fixed_varLest_prior(2,2);ub=halibut_fixed_varLest_prior(2,3);
  init_bounded_number halibut_fixed_varLest_M(lb,ub,phz)

  !!phz=ll_Sfull_prior(1,1);lb=ll_Sfull_prior(1,2);ub=ll_Sfull_prior(1,3);
  init_bounded_number ll_Sfull_F(lb,ub,phz)
  !!phz=ll_varLest_prior(1,1);lb=ll_varLest_prior(1,2);ub=ll_varLest_prior(1,3);
  init_bounded_number ll_varLest_F(lb,ub,phz)

  !!phz=ll_Sfull_prior(2,1);lb=ll_Sfull_prior(2,2);ub=ll_Sfull_prior(2,3);
  init_bounded_number ll_Sfull_M(lb,ub,phz)
  !!phz=ll_varLest_prior(2,1);lb=ll_varLest_prior(2,2);ub=ll_varLest_prior(2,3);
  init_bounded_number ll_varLest_M(lb,ub,phz)

  !!phz=ot_Sfull_prior(1);lb=ot_Sfull_prior(2);ub=ot_Sfull_prior(3);
  init_bounded_number ot_Sfull(lb,ub,phz)
  !!phz=ot_varLest_prior(1);lb=ot_varLest_prior(2);ub=ot_varLest_prior(3);
  init_bounded_number ot_varLest(lb,ub,phz)
  !!phz=ot_varRest_prior(1);lb=ot_varRest_prior(2);ub=ot_varRest_prior(3);
  init_bounded_number ot_varRest(lb,ub,phz)

   !!phz=log_q_old_summer_rv_prior(1);lb=log_q_old_summer_rv_prior(2);ub=log_q_old_summer_rv_prior(3);
   init_bounded_number log_q_old_summer_rv(lb,ub,phz)
 
   !!phz=log_q_halibut_fixed_prior(1);lb=log_q_halibut_fixed_prior(2);ub=log_q_halibut_fixed_prior(3);
   init_bounded_number log_q_halibut_fixed(lb,ub,phz)

  !!phz=u_ll_prior(1);lb=u_ll_prior(2);ub=u_ll_prior(3);
  init_bounded_number u_ll_2003(lb,ub,phz)

  vector summer_rv_resid(1,summer_rv_n_years)
  vector halibut_fixed_resid(1,halibut_fixed_n_years)
  vector sex_rat_rv_resid(1,N_lens)
  vector sex_rat_comm_resid(1,N_lens_comm)

  matrix pred_catlen_rv(1,N_catlen_rv,1,N_lens)
  matrix pred_catlenP_rv(1,N_catlen_rv,1,N_lens) //the proportions
  matrix pred_catlen_comm(1,N_catlen_comm,1,N_lens_comm)
  matrix pred_catlenP_comm(1,N_catlen_comm,1,N_lens_comm) //the proportions

  matrix catlen_survey_resid(1,N_catlen_rv,1,N_lens)
  matrix catlen_comm_resid(1,N_catlen_comm,1,N_lens_comm)
  matrix catlen_comm_kernal_resid(1,N_catlen_comm,1,N_lens_comm)
  matrix catlen_survey_kernal_resid(1,N_catlen_rv,1,N_lens)

  // ALL THE REST OF THE DERIVED THINGS

  matrix Females(first_year,last_year+1,first_age,last_age)
  matrix Males(first_year,last_year+1,first_age,last_age)
  sdreport_vector naa2009(first_age,last_age)
  vector Total_N(first_year,last_year+1)
  vector Total_W(first_year,last_year+1)
  vector Total_F(first_year,last_year+1)
  vector Total_M(first_year,last_year+1)
  vector SSN(first_year,last_year+1) //start of year			FEMALES ONLY		
  vector SSN2(first_year,last_year+1) //start of year			MALES + FEMALES
  vector SSB(first_year,last_year+1) //start of year	FEMALES ONLY		
  sdreport_vector SSB2(first_year,last_year+1) //start of year			MALES + FEMALES
//  vector SSB(first_year,last_year+1) //start of year
  sdreport_number SSB2009
  sdreport_number FF2009
  sdreport_number FM2009
  vector rec(first_year,last_year+1) // both sexes combined
  vector dudfem(first_age,last_age) // virgin spawning abudnance
  vector dudmale(first_age,last_age)
  matrix surv(1,N_sexes,first_age,last_age)
  matrix surv50(1,N_sexes,first_age,last_age)
  matrix surv25(1,N_sexes,first_age,last_age)
  matrix surv75(1,N_sexes,first_age,last_age)
  vector w_F(first_age,last_age)

  //predicted rv totals
  vector pred_summer_rv_total(1970,last_year)
  vector pred_halibut_fixed_total(1998,last_year)
 
  matrix summer_RVage(1970,last_year,first_age,last_age)
  matrix halibut_fixed_age(1998,last_year,first_age,last_age)

   vector VB_Canada(1960,last_year+1) //filled in below
   vector VB_Foreign(1960,last_year+1) 

   vector VB_ll(1960,last_year+1)
   vector VB_ot(1960,last_year+1)

   vector TB_Canada(1960,last_year+1) //filled in below
   vector TB_Foreign(1960,last_year+1) 

   vector TB_ll(1960,last_year+1)
   vector TB_ot(1960,last_year+1)

   vector C_ll_F(1960,last_year)
   vector C_ll_M(1960,last_year)
   vector C_ot(1960,last_year)

   vector u_Canada(1960,last_year) 
   vector u_Foreign(1960,last_year) 
 
   vector u_ll(1960,last_year)
   number u_ll_pdf
   vector u_ot(1960,last_year)
   sdreport_vector u_ll_F(first_year,last_year);
   vector u_ll_M(first_year,last_year);
  
  vector old_summer_rv_F_sel(first_age,last_age)
  vector old_summer_rv_M_sel(first_age,last_age)

  vector halibut_fixed_F_sel(first_age,last_age)
  vector halibut_fixed_M_sel(first_age,last_age)

  vector ll_F_sel(first_age,last_age)
  vector ll_M_sel(first_age,last_age)

  vector ot_F_sel(first_age,last_age)
  vector ot_M_sel(first_age,last_age)

  number varL
  number varR
  number Sfull

  matrix F_utotal(1960,last_year,first_age,last_age)
  matrix M_utotal(1960,last_year,first_age,last_age)
  matrix F_ucatch(1960,last_year,first_age,last_age)
  matrix M_ucatch(1960,last_year,first_age,last_age)
  3darray prop_landed(first_year,last_year+1,1,2,first_age,last_age)
  3darray prop_fmort(first_year,last_year+1,1,2,first_age,last_age)

  number S0
  
  number summer_rv_like
  number halibut_fixed_like

  vector catlen_survey_like(1,N_catlen_surveys)
  vector catlen_comm_like(1,N_catlen_commercial)  //same number of indices
  vector catlen_survey_kernal_like(1,N_catlen_surveys)
  vector catlen_comm_kernal_like(1,N_catlen_commercial)  //same number of indices

  number fpen
  number ffpen
  number upen
  number reccheck
  number recpen
  number recpen2
  number r_penF
  number alphapen
  number age1_pen
  number age1_last_year_pdf
  number recdiff_pdf
  number recdiff_sigma

  number temp1
  number gm
  
  matrix pred_f(1,N_catlen_rv,1,N_lens);
  matrix pred_m(1,N_catlen_rv,1,N_lens);
  matrix pred_f_comm(1,N_catlen_comm,1,N_lens_comm);
  matrix pred_m_comm(1,N_catlen_comm,1,N_lens_comm);
  matrix pred_sex_rat(1,N_catlen_rv,1,N_lens);
  matrix pred_sex_rat_comm(1,N_catlen_comm,1,N_lens_comm);

  vector fcomp(1,12)
  objective_function_value ofv

  number SPRF0
//  sdreport_number ssbeq
  number ssbeq
  number receq
//  number DBL_EPSILON;

  matrix obs_catlen_rv_short(1,N_catlen_rv,fminbins2,fmaxbins2);
  matrix pred_catlen_rv_short(1,N_catlen_rv,fminbins2,fmaxbins2);

  matrix obs_catlen_ll_short_female(1,N_catlen_comm,fminbins,fmaxbins);
  matrix pred_catlen_ll_short_female(1,N_catlen_comm,fminbins,fmaxbins);
  matrix obs_catlen_ll_short_male(1,N_catlen_comm,mminbins,mmaxbins);
  matrix pred_catlen_ll_short_male(1,N_catlen_comm,mminbins,mmaxbins);

  matrix obs_catlen_ot_short(1,N_catlen_comm,fminbins3,fmaxbins3);
  matrix pred_catlen_ot_short(1,N_catlen_comm,fminbins3,fmaxbins3);

  matrix obs_catlen_hsurvey_short_female(1,N_catlen_rv,fminbins,fmaxbins);
  matrix pred_catlen_hsurvey_short_female(1,N_catlen_rv,fminbins,fmaxbins);
  matrix obs_catlen_hsurvey_short_male(1,N_catlen_rv,mminbins,mmaxbins);
  matrix pred_catlen_hsurvey_short_male(1,N_catlen_rv,mminbins,mmaxbins);

  matrix obs_catlen_hCI_short_female(1,N_catlen_rv,fminbins,fmaxbins);
  matrix pred_catlen_hCI_short_female(1,N_catlen_rv,fminbins,fmaxbins);
  matrix obs_catlen_hCI_short_male(1,N_catlen_rv,mminbins,mmaxbins);
  matrix pred_catlen_hCI_short_male(1,N_catlen_rv,mminbins,mmaxbins);

//============================================================================
PRELIMINARY_CALCS_SECTION

//  DBL_EPSILON=2.2204460492503131e-16;

  //VH
  int i1=0;
  int i2=0;
  i1 = (option_match(ad_comm::argc,ad_comm::argv,"-phase",i2)); 
  if(i1<=0) initialize_parameters(); 

//  summer_rv_total*=10;
//  halibut_fixed_total*=10;

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


//  cout << obs_catlen_comm(31,1) << endl;

//  cout << "Okay after initialize_parameters" << endl;
  setup_fixed_parameters();
//  cout << "Okay after setup_fixed_parameters" << endl;

//=====================================================

RUNTIME_SECTION
  convergence_criteria 1e-4
  maximum_function_evaluations 200,200,400,500,500,50000,50000


//============================================================================
PROCEDURE_SECTION
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
//  cout << "Okay after setup" << endl;
  initial_conditions();
//  cout << "Okay after initial_conditions" << endl;
  dynamics();
//  cout << "Okay after dynamics" << endl;
  
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
//  ofv+=fpen*10; //applied to males and females to keep positive * chages results when off
  ofv+=ffpen*10; // pen on estimating F iteratively
//  ofv+=100*upen; //trying to keep U under .5
//  ofv+=recpen; //constrain N1 to be less that R at equlibrium
//  ofv+=r_penF*10; //Keeps R positive when trying to estimate depensation parameters * chages results when off
  
//  ofv+=age1_pen;      //pen on process error on recruitment	
  ofv+=recdiff_pdf;      //random walk on recruitment	
//  ofv+=age1_last_year_pdf;
//  ofv+=sex_rat_like; // proportion female at len in commercial fishery, halibut survey and commercial index

//													like_weight
//  ofv+= summer_rv_like;		// fcomp(1)			10 
//  ofv+= halibut_fixed_like;		// fcomp(5)		1 
//  ofv+= halibut_comm_index_like;	// fcomp(6)		0 
//  ofv+= halibut_3nop_cpue_like;	// fcomp(7)		0 
//  ofv+= halibut_4vwx_cpue_like;	// fcomp(8)		0 
//  ofv+=sum(catlen_comm_like);		// fcomp(9)		1 
//  ofv+=sum(catlen_survey_like);	// fcomp(10)	1 	
//  ofv+= comm_total_like;		// fcomp(11)		0 
//  ofv+= tagging_like;			// fcomp(12)		0
 
   
    for(i=1;i<=12;i++)
     {
     ofv+= fcomp(i)*like_weight(i);
     }

   if(mceval_phase()) 
   {
   MCWrite(); 
   }
  } //close else

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//============================================================================
FUNCTION initialize_parameters
  N1=N1_prior(7);
  start_Z=start_Z_prior(7);
//  start_Z=Minit;
  R0=R0_prior(7);
  beta=beta_prior(7);
  alpha=alpha_prior(7);
  p1=p1_prior(7);
  p2=p2_prior(7);
  p3=p3_prior(7);
  M=M_prior(7);
//  M=Minit;
  log_RecDev=log_RecDev_prior(7);
  all_log_RecDev = 0.;

// starting values back calc. from the catch at length
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
//  age1(2008)=580662.6;
//  age1(2009)=580662.6;

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

//============================================================================
FUNCTION setup_fixed_parameters
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
//  F_sigma(age)=10;
//  M_sigma(age)=10;
//  if(age<5){F_sigma(age)=5;}
//  if(age<5){M_sigma(age)=5;}
  if(age>=5){F_sigma(age)=10;}
  if(age>=4){M_sigma(age)=10;}
  }
  }

  for (sex=1;sex<=N_sexes;sex++)
  {
  w(sex)=bi(sex)*pow(lenatage(sex),bii(sex))/1000;			// convert to kg
  }
  w_F=w(1);

// AGE-LENGTH TRANSITION MATRIX
// liner relationship
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

// make sum to 1 
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


//normalize to go from len to age
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

//-----------------------------------------------------------------------------
FUNCTION setup_variable_parameters
 
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

//*********************************************************
FUNCTION initial_conditions
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
//-----------------------------------------------------------------------------

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

//====================================================================

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
//  VB_ot(first_year)=((elem_prod(elem_prod(Females(first_year),surv50(1)),elem_prod(prop_landed(first_year,1),ot_F_sel))*w(1))+(elem_prod(elem_prod(Males(first_year),surv50(2)),elem_prod(prop_landed(first_year,2),ot_F_sel))*w(2)))/1000;

  TB_Canada(first_year)=(elem_prod(Females(first_year),surv50(1))*w(1)+
    					 elem_prod(  Males(first_year),surv50(2))*w(2))/1000; 

//-----------------------------------------------------------------------------
FUNCTION dynamics
//NOTE in Bob's simulation model the catch and the virtual population are in 1000s in the halibut assessment model they are in raw numbers
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

// Canadian Landings by gear
  u_ll_F(year)=GetF_forward(WF,CF,M);
  u_ll_M(year)=GetF_forward(WM,CM,M);
  u_ot(year)=GetF_forward(WF+WM,COT,M);
//  u_ll_F(year)=ffromcn(WF,CF,M);
//  u_ll_M(year)=ffromcn(WM,CM,M);

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
  
// estimate recruitment as a free parameter
  Females(year+1,first_age)=age1(year+1)*.5;
  Males(year+1,first_age)=age1(year+1)*.5;
  
// recruitment assume 1:1 sex ratio at age 1
// VH: change structure so bias correction is included when log_RecDev's estimated in run (but maybe not phase)         
        if(log_RecDev_prior(1)>0)
        {
//	  Females(year+1,first_age)=2.71828*p2/p1*SSB(year)*exp(-SSB(year)/p1)*exp(all_log_RecDev(year+1))*biasRannual*0.5;
//	  Males(year+1,first_age)=2.71828*p2/p1*SSB(year)*exp(-SSB(year)/p1)*exp(all_log_RecDev(year+1))*biasRannual*0.5;
        }
       else
        {
//	  Females(year+1,first_age)=2.71828*p2/p1*SSB(year)*exp(-SSB(year)/p1)*exp(all_log_RecDev(year+1))*0.5;
//	  Males(year+1,first_age)=2.71828*p2/p1*SSB(year)*exp(-SSB(year)/p1)*exp(all_log_RecDev(year+1))*0.5;
        }

// calculate the geometric mean of recruitment for last 3 years
  int index;
  index=3;
// temp1 declared in paramter section
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
//   VB_ot(year+1)=(elem_prod(elem_prod(Females(year+1),surv50(1)),elem_prod(prop_landed(year+1,1),ot_F_sel))*w(1))+
//                  (elem_prod(elem_prod(Males(year+1),surv50(2)),elem_prod(prop_landed(year+1,2),ot_F_sel))*w(2))/1000;
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

//  cout<<temp(2008)<<endl;
   Total_N(year+1)=sum(Females(year+1))+sum(Males(year+1));
   Total_F(year+1)=temp_F(year);
   Total_M(year+1)=temp_M(year);
//   Total_F(year+1)=(sum(elem_prod(Females(year+1),w(1))))/1000; 	// in t
//   Total_M(year+1)=(sum(elem_prod(Males(year+1),w(2))))/1000; 		// in t
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION Summer_RV_Likelihood

   dvar_vector junk(1,summer_rv_n_years);

    for(i=1;i<=summer_rv_n_years;i++)
     {
     junk(i)=pred_summer_rv_total(summer_rv_years(i));   
     }

    summer_rv_like=LogNormal_LogLike(summer_rv_total+.01,junk+.01,summer_rv_n_years,.1,2,summer_rv_resid);
    //  cout<<summer_rv_like<<endl;
	// ad_exit(1);

    
    fcomp(1)=summer_rv_like;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION Halibut_Fixed_Likelihood
   
  dvar_vector junk(1,halibut_fixed_n_years);

    for(i=1;i<=halibut_fixed_n_years;i++)
     {
     junk(i)=pred_halibut_fixed_total(halibut_fixed_years(i));   
     }

   halibut_fixed_like=LogNormal_LogLike(halibut_fixed_total+.01,junk+.01,halibut_fixed_n_years,0.1,2,halibut_fixed_resid);
   fcomp(5)=halibut_fixed_like;
//-----------------------------------------------------------------------------
FUNCTION Catlen_Comm_Likelihood
 
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
//   catlen_comm_like(series)+=0;  
   catlen_comm_like(series)+=Multinomial(obs_catlen_ot_short(obs)(fminbins3,fmaxbins3), pred_catlen_ot_short(obs)(fminbins3,fmaxbins3),n_short,maxs,catlen_LikeType(series),catlen_comm_resid(obs)(fminbins3,fmaxbins3));  
  break;
  }  // close case loop
  }  // close obs loop
  fcomp(9)=sum(catlen_comm_like);
//-----------------------------------------------------------------------------
FUNCTION Catlen_Survey_Likelihood

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
// F_utotal to account for mid fishery
	if(ssex==1)
	{
	pred_catlen_rv(obs) =elem_prod(elem_prod(exp(-(F_utotal(year)+M)*.5),Females(year)),halibut_fixed_F_sel)*F_agelen;
	}
	else
	{
	pred_catlen_rv(obs) =elem_prod(elem_prod(exp(-(M_utotal(year)+M)*.5),Males(year)),halibut_fixed_M_sel)*M_agelen;
	}
//  cout<<pred_catlen_rv(obs)<<endl;
    break;
  }

// TRUNCATE LENGTH BINS BY CREATING + AND - GROUPS
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION MCWrite
  report1 
//    <<" p1 "<<p1                        
//    <<" p2 "<<p2                        
//    <<" p3 "<<p3                        
//    <<" age1 "<<age1                        
    <<" old_summer_rv_Sfull "<<old_summer_rv_Sfull
//    <<" old_summer_rv_varLest "<<old_summer_rv_varLest
    <<" old_summer_rv_varRest "<<old_summer_rv_varRest
    <<" halibut_fixed_Sfull_F "<<halibut_fixed_Sfull_F
    <<" halibut_fixed_Sfull_M "<<halibut_fixed_Sfull_M
    <<" ll_Sfull_F "<<ll_Sfull_F
    <<" ll_Sfull_M "<<ll_Sfull_M
//    <<" ll_varLest_F "<<ll_varLest_F
//    <<" ll_varLest_M "<<ll_varLest_M
    <<" ot_Sfull "<<ot_Sfull
//    <<" ot_varLest "<<ot_varLest
    <<" ot_varRest "<<ot_varRest
    <<" log_q_old_summer_rv "<<log_q_old_summer_rv
    <<" log_q_halibut_fixed "<<log_q_halibut_fixed
    <<" SSB2009 "<<SSB2009                        
    <<" FF2009 "<<FF2009                        
    <<" FM2009 "<<FM2009                        
    <<" naa2009 "<<naa2009                        
  <<" ofv "<<ofv
  <<endl;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//estimate F iteratively
//code by VH
FUNCTION dvariable GetF_forward(const dvariable& N, const dvariable& C, const dvariable& M)
  dvariable Ctmp=C;
  dvariable fest=Ctmp/N;
//  dvariable fest=M;
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//estimate F iteratively
//code by Bob M.
FUNCTION dvariable ffromcn(const dvariable& N, const dvariable& C, const dvariable& M)
  dvariable Ctmp=C;
  dvariable fest=M;
//  dvariable fest=C/N;
  dvariable ztmp;
  for (int iter=1;iter<=10;iter++)
  {
    ztmp=fest+M;
    fest=Ctmp*ztmp/N*(1-exp(-1*ztmp));
  }
  return(fest);

//=============================================================================
FUNCTION dvar_vector Make_select_dhn(int minage, int maxage, dvariable full, dvariable varL, dvariable varR)

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
//=============================================================================
FUNCTION dvar_vector Make_select_logistic(int minage, int maxage, dvariable full, dvariable varL)

  dvar_vector selection(minage,maxage);
  selection.initialize();

  for (age=minage;age<=maxage;age++)
  {
  selection(age)=1/(1+exp(-varL*(age-full))); 
  }	
  selection = selection/max(selection);
  return(selection);

//================================================================================================
FUNCTION Age1_pen
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

//================================================================================================
FUNCTION Age1_last_year_pdf
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
 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//================================================================================================
FUNCTION Recdiff_pdf

  recdiff_pdf=0;
  recdiff_sigma=.5;
  dvar_vector recdiff(first_year,last_year-1);

  for(year=first_year;year<=last_year-1;year++)
    {
    recdiff(year)=log(rec(year)/rec(year+1));
    }
  recdiff_pdf=norm2(recdiff/pow(recdiff_sigma,2));
 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//****************************************************************************************************************

FUNCTION dvariable LogNormal_LogLike(dvar_vector observed, dvar_vector predicted, double len,double sigma, double type)

   dvariable neglogLike;
   dvariable sumlogr;
   dvariable len2;
   dvariable sumsqrlogratio;
   double pi=3.1415926535897;
  
//THIS ALLOWS YOU TO REM OUT VALUES IN THE DAT FILE BY MAKING THEM NEGATIVE:
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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//****************************************************************************************************************

FUNCTION dvariable LogNormal_LogLike(dvar_vector observed, dvar_vector predicted, double len,double sigma, double type,dvar_vector resid)

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
//THIS ALLOWS YOU TO REM OUT VALUES IN THE DAT FILE BY MAKING THEM NEGATIVE:
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
// calculate residuals here  
  return(neglogLike);

//======================================================================================================================
FUNCTION dvariable calc_sigma(dvar_vector observed, dvar_vector predicted, double len)

   dvariable sigma;
   dvariable sumsqrlogratio;
   dvariable len2;

   sumsqrlogratio=0;
   len2=0; 

//THIS ALLOWS YOU TO REM OUT VALUES IN THE DAT FILE BY MAKING THEM NEGATIVE:

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

//=====================================================================================================================
FUNCTION dvariable Normal(dvariable observed,dvariable mean,double nsigma)

   dvariable normal;
   double pi=3.1415926535897;
   dvariable neglogpdf;
   //pdf=(1/(sqrt(2*pi)*nsigma)*exp(-(pow(observed-mean,2)/(2*pow(nsigma,2)))));
   //neglogpdf=-log(pdf);
    neglogpdf=log(sqrt(2*pi)*nsigma)+pow(observed-mean,2)/(2*pow(nsigma,2));
   return(neglogpdf);

//=====================================================================================================================
FUNCTION dvariable Binomial(dvariable observed,double x,double n)

   dvariable p;
   dvariable pdf;
   dvariable neglogpdf;
   p=observed;				
   pdf=pow(p,x)*pow(1-p,n-x);

   neglogpdf=-log(pdf);
   return(neglogpdf);


//=====================================================================================================================
// Inverse of posfun
FUNCTION dvariable negfun(const dvariable&x,const double eps,dvariable& pen)
  {
  if (x<=eps) {
  return x;
  } else {
  pen+=.01*square(x-eps);

  cout<<" neg fun "<<x<<" "<<eps/(2-x/eps)<<" "<<pen<<endl;
  return x/(2-eps/x);
  }
  }
//=====================================================================================================================
//----------------------------------------------------------------------------
//FUNCTION dvariable Multinomial(dvector obs,dvar_vector pred,data_int& a,double& b,int& Ltype,dvar_vector& resid)
//hal18//FUNCTION dvariable Multinomial(dvector& obs,dvar_vector& pred,int a,double& b,int Ltype,dvar_vector &resid)
              
FUNCTION dvariable Multinomial(dvar_vector obs,dvar_vector pred,data_int& a,double& b,int& Ltype,dvar_vector resid)
  {
//VH
   dvariable log_likelihood;

     resid=elem_div((obs-pred),sqrt(elem_prod(pred,(1.-pred))/b));
//    return sum(-1.* elem_prod(Nsamp*oprop,log(pprop))+ elem_prod(Nsamp*oprop,log(oprop)));
      log_likelihood=-1.* sum(elem_prod(b*obs,log(pred+1e-6))) + sum(elem_prod(b*obs,log(obs+1e-6)));
    return(log_likelihood);
   }
//=====================================================================================================================


FUNCTION dvariable Multifan(dvector& obs,dvar_vector& pred,int a,double& b,int Ltype)
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
//VH
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
 
//-----------------------------------------------------------------------------
REPORT_SECTION

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

