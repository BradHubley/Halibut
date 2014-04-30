#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <caa1.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_report1 = new ofstream("mcout.dat",ios::app);  // you must delete this file before a new run;
  debug.allocate("debug");
  all_first_year.allocate("all_first_year");
if (debug==1){char mark; cout<<"all_first_year= "<<all_first_year<<endl; cin>>mark; }
  all_last_year.allocate("all_last_year");
  all_first_age.allocate("all_first_age");
  all_last_age.allocate("all_last_age");
  select_switch.allocate("select_switch");
  select_data.allocate(1,4,"select_data");
 first_year=select_data(1);last_year=select_data(2);first_age=select_data(3);last_age=select_data(4);
  F_full_prior.allocate(1,4,"F_full_prior");
if (debug==1){char mark; cout<<"F_full_prior= "<<F_full_prior<<endl; cin>>mark; }
  start_Z_prior.allocate(1,4,"start_Z_prior");
  year1_prior.allocate(1,4,"year1_prior");
  age1_prior.allocate(1,4,"age1_prior");
  mort_switch.allocate("mort_switch");
  M_const_prior.allocate(1,4,"M_const_prior");
  C_Sfullold_prior.allocate(1,4,"C_Sfullold_prior");
  C_log_varLold_prior.allocate(1,4,"C_log_varLold_prior");
  C_log_varRold_prior.allocate(1,4,"C_log_varRold_prior");
  RV_Sfullold_prior.allocate(1,4,"RV_Sfullold_prior");
  RV_log_varLold_prior.allocate(1,4,"RV_log_varLold_prior");
  RV_log_varRold_prior.allocate(1,4,"RV_log_varRold_prior");
  RVspring_Sfull_prior.allocate(1,4,"RVspring_Sfull_prior");
  RVspring_log_varL_prior.allocate(1,4,"RVspring_log_varL_prior");
  RVspring_log_varR_prior.allocate(1,4,"RVspring_log_varR_prior");
  RV_log_qold_prior.allocate(1,4,"RV_log_qold_prior");
  RVspring_log_q_prior.allocate(1,4,"RVspring_log_q_prior");
  sigma1.allocate("sigma1");
  sigma2.allocate("sigma2");
  like_weight.allocate(1,4,"like_weight");
  n_landings.allocate("n_landings");
  landings.allocate(1,n_landings,"landings");
  landing_year.allocate(1,n_landings,"landing_year");
if (debug==1){char mark; cout<<"landing_year= "<<landing_year<<endl; cin>>mark; }
  all_mat.allocate(all_first_year,all_last_year,all_first_age,all_last_age,"all_mat");
  n_all_comcaa_years.allocate("n_all_comcaa_years");
  all_comcaa_years.allocate(1,n_all_comcaa_years,"all_comcaa_years");
if (debug==1){char mark; cout<<"all_comcaa_years= "<<all_comcaa_years<<endl; cin>>mark; }
  all_obs_comcaa.allocate(1988,all_last_year,all_first_age,all_last_age,"all_obs_comcaa");
if (debug==1){char mark; cout<<"all_obs_comcaa= "<<all_obs_comcaa(1988)<<endl; cin>>mark; }
  n_all_rvcaa_years.allocate("n_all_rvcaa_years");
  all_rvcaa_years.allocate(1,n_all_rvcaa_years,"all_rvcaa_years");
if (debug==1){char mark; cout<<"all_rvcaa_years= "<<all_rvcaa_years<<endl; cin>>mark; }
  all_obs_rvcaa.allocate(all_first_year,all_last_year,all_first_age,all_last_age,"all_obs_rvcaa");
if (debug==1){char mark; cout<<"all_obs_rvcaa= "<<all_obs_rvcaa(1970)<<endl; cin>>mark; }
  all_obs_rvwaa.allocate(all_first_year,all_last_year,all_first_age,all_last_age,"all_obs_rvwaa");
if (debug==1){char mark; cout<<"all_obs_rvwaa= "<<all_obs_rvwaa(1970)<<endl; cin>>mark; }
  n_all_rvspring_years.allocate("n_all_rvspring_years");
if (debug==1){char mark; cout<<"n_all_rvspring_years= "<<n_all_rvspring_years<<endl; cin>>mark; }
  all_rvspring_years.allocate(1,n_all_rvspring_years,"all_rvspring_years");
  all_obs_rvspring.allocate(1998,all_last_year,all_first_age,all_last_age,"all_obs_rvspring");
  dummy.allocate("dummy");
if (debug==1){char mark; cout<<"dummy= "<<dummy<<endl; cin>>mark; }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  obs_comcaa.allocate(1988,last_year,first_age,last_age,"obs_comcaa");
  #ifndef NO_AD_INITIALIZE
    obs_comcaa.initialize();
  #endif
  obs_rvcaa.allocate(first_year,last_year,first_age,last_age,"obs_rvcaa");
  #ifndef NO_AD_INITIALIZE
    obs_rvcaa.initialize();
  #endif
  obs_rvwaa.allocate(first_year,last_year,first_age,last_age,"obs_rvwaa");
  #ifndef NO_AD_INITIALIZE
    obs_rvwaa.initialize();
  #endif
  obs_rvspring.allocate(1998,2009,first_age,last_age,"obs_rvspring");
  #ifndef NO_AD_INITIALIZE
    obs_rvspring.initialize();
  #endif
  mat.allocate(first_year,last_year,first_age,last_age,"mat");
  #ifndef NO_AD_INITIALIZE
    mat.initialize();
  #endif
  Nage.allocate(first_year,last_year,first_age,last_age,"Nage");
  #ifndef NO_AD_INITIALIZE
    Nage.initialize();
  #endif
  Cage.allocate(first_year,last_year,first_age,last_age,"Cage");
  #ifndef NO_AD_INITIALIZE
    Cage.initialize();
  #endif
  RVage.allocate(first_year,last_year,first_age,last_age,"RVage");
  #ifndef NO_AD_INITIALIZE
    RVage.initialize();
  #endif
  RVspring_age.allocate(1998,2009,first_age,last_age,"RVspring_age");
  #ifndef NO_AD_INITIALIZE
    RVspring_age.initialize();
  #endif
  Wtage.allocate(first_year,last_year,first_age,last_age,"Wtage");
  #ifndef NO_AD_INITIALIZE
    Wtage.initialize();
  #endif
  matbio.allocate(first_year,last_year,first_age,last_age,"matbio");
  #ifndef NO_AD_INITIALIZE
    matbio.initialize();
  #endif
  obs_Cpage.allocate(1988,last_year,first_age,last_age,"obs_Cpage");
  #ifndef NO_AD_INITIALIZE
    obs_Cpage.initialize();
  #endif
  obs_RVpage.allocate(first_year,last_year,first_age,last_age,"obs_RVpage");
  #ifndef NO_AD_INITIALIZE
    obs_RVpage.initialize();
  #endif
  obs_RVspring_page.allocate(1998,last_year,first_age,last_age,"obs_RVspring_page");
  #ifndef NO_AD_INITIALIZE
    obs_RVspring_page.initialize();
  #endif
  pred_Cpage.allocate(1988,last_year,first_age,last_age,"pred_Cpage");
  #ifndef NO_AD_INITIALIZE
    pred_Cpage.initialize();
  #endif
  pred_RVpage.allocate(first_year,last_year,first_age,last_age,"pred_RVpage");
  #ifndef NO_AD_INITIALIZE
    pred_RVpage.initialize();
  #endif
  pred_RVspring_page.allocate(1998,last_year,first_age,last_age,"pred_RVspring_page");
  #ifndef NO_AD_INITIALIZE
    pred_RVspring_page.initialize();
  #endif
  F.allocate(first_year,last_year,first_age,last_age,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  select.allocate(first_year,last_year,first_age,last_age,"select");
  #ifndef NO_AD_INITIALIZE
    select.initialize();
  #endif
  RV_select.allocate(first_year,last_year,first_age,last_age,"RV_select");
  #ifndef NO_AD_INITIALIZE
    RV_select.initialize();
  #endif
  RVspring_select.allocate(1998,2009,first_age,last_age,"RVspring_select");
  #ifndef NO_AD_INITIALIZE
    RVspring_select.initialize();
  #endif
  obs_total_catch.allocate(1988,last_year,"obs_total_catch");
  #ifndef NO_AD_INITIALIZE
    obs_total_catch.initialize();
  #endif
  pred_total_catch.allocate(first_year,last_year,"pred_total_catch");
  #ifndef NO_AD_INITIALIZE
    pred_total_catch.initialize();
  #endif
  obs_rvtotal_catch.allocate(first_year,last_year,"obs_rvtotal_catch");
  #ifndef NO_AD_INITIALIZE
    obs_rvtotal_catch.initialize();
  #endif
  pred_rvtotal_catch.allocate(first_year,last_year,"pred_rvtotal_catch");
  #ifndef NO_AD_INITIALIZE
    pred_rvtotal_catch.initialize();
  #endif
  obs_rvspring_total_catch.allocate(1998,2009,"obs_rvspring_total_catch");
  #ifndef NO_AD_INITIALIZE
    obs_rvspring_total_catch.initialize();
  #endif
  pred_rvspring_total_catch.allocate(1998,2009,"pred_rvspring_total_catch");
  #ifndef NO_AD_INITIALIZE
    pred_rvspring_total_catch.initialize();
  #endif
  ssb.allocate(first_year,last_year,"ssb");
  RVwt_age.allocate(first_age,last_age,"RVwt_age");
  #ifndef NO_AD_INITIALIZE
    RVwt_age.initialize();
  #endif
phz=F_full_prior(1);lb=F_full_prior(2);ub=F_full_prior(3);
  F_full.allocate(first_year,last_year,lb,ub,phz,"F_full");
phz=start_Z_prior(1);lb=start_Z_prior(2);ub=start_Z_prior(3);
  start_Z.allocate(lb,ub,phz,"start_Z");
phz=year1_prior(1);lb=year1_prior(2);ub=year1_prior(3);
  year1.allocate(first_age,last_age,lb,ub,phz,"year1");
phz=age1_prior(1);lb=age1_prior(2);ub=age1_prior(3);
  age1.allocate(first_year,last_year,lb,ub,phz,"age1");
phz=M_const_prior(1);lb=M_const_prior(2);ub=M_const_prior(3);
  M_const.allocate(lb,ub,phz,"M_const");
phz=C_Sfullold_prior(1);lb=C_Sfullold_prior(2);ub=C_Sfullold_prior(3);
  C_Sfullold.allocate(lb,ub,phz,"C_Sfullold");
phz=C_log_varLold_prior(1);lb=C_log_varLold_prior(2);ub=C_log_varLold_prior(3);
  C_log_varLold.allocate(lb,ub,phz,"C_log_varLold");
phz=C_log_varRold_prior(1);lb=C_log_varRold_prior(2);ub=C_log_varRold_prior(3);
  C_log_varRold.allocate(lb,ub,phz,"C_log_varRold");
phz=RV_Sfullold_prior(1);lb=RV_Sfullold_prior(2);ub=RV_Sfullold_prior(3);
  RV_Sfullold.allocate(lb,ub,phz,"RV_Sfullold");
phz=RV_log_varLold_prior(1);lb=RV_log_varLold_prior(2);ub=RV_log_varLold_prior(3);
  RV_log_varLold.allocate(lb,ub,phz,"RV_log_varLold");
phz=RV_log_varRold_prior(1);lb=RV_log_varRold_prior(2);ub=RV_log_varRold_prior(3);
  RV_log_varRold.allocate(lb,ub,phz,"RV_log_varRold");
phz=RVspring_Sfull_prior(1);lb=RVspring_Sfull_prior(2);ub=RVspring_Sfull_prior(3);
  RVspring_Sfull.allocate(lb,ub,phz,"RVspring_Sfull");
phz=RVspring_log_varL_prior(1);lb=RVspring_log_varL_prior(2);ub=RVspring_log_varL_prior(3);
  RVspring_log_varL.allocate(lb,ub,phz,"RVspring_log_varL");
phz=RVspring_log_varR_prior(1);lb=RVspring_log_varR_prior(2);ub=RVspring_log_varR_prior(3);
  RVspring_log_varR.allocate(lb,ub,phz,"RVspring_log_varR");
phz=RV_log_qold_prior(1);lb=RV_log_qold_prior(2);ub=RV_log_qold_prior(3);
  RV_log_qold.allocate(lb,ub,phz,"RV_log_qold");
phz=RVspring_log_q_prior(1);lb=RVspring_log_q_prior(2);ub=RVspring_log_q_prior(3);
  RVspring_log_q.allocate(lb,ub,phz,"RVspring_log_q");
  mean_rvwaa.allocate(first_age,last_age,"mean_rvwaa");
  #ifndef NO_AD_INITIALIZE
    mean_rvwaa.initialize();
  #endif
  mean_rvwaa1970.allocate(first_age,last_age,"mean_rvwaa1970");
  #ifndef NO_AD_INITIALIZE
    mean_rvwaa1970.initialize();
  #endif
  mean_rvwaa2003.allocate(first_age,last_age,"mean_rvwaa2003");
  #ifndef NO_AD_INITIALIZE
    mean_rvwaa2003.initialize();
  #endif
  M.allocate(first_year,last_year,first_age,last_age,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  ssb2003.allocate("ssb2003");
  ctotal_like.allocate("ctotal_like");
  #ifndef NO_AD_INITIALIZE
  ctotal_like.initialize();
  #endif
  Cpage_like.allocate("Cpage_like");
  #ifndef NO_AD_INITIALIZE
  Cpage_like.initialize();
  #endif
  rvtotal_like.allocate("rvtotal_like");
  #ifndef NO_AD_INITIALIZE
  rvtotal_like.initialize();
  #endif
  RVpage_like.allocate("RVpage_like");
  #ifndef NO_AD_INITIALIZE
  RVpage_like.initialize();
  #endif
  rvspring_total_like.allocate("rvspring_total_like");
  #ifndef NO_AD_INITIALIZE
  rvspring_total_like.initialize();
  #endif
  RVspring_page_like.allocate("RVspring_page_like");
  #ifndef NO_AD_INITIALIZE
  RVspring_page_like.initialize();
  #endif
  ctotal_sigma.allocate("ctotal_sigma");
  #ifndef NO_AD_INITIALIZE
  ctotal_sigma.initialize();
  #endif
  rvtotal_sigma.allocate("rvtotal_sigma");
  #ifndef NO_AD_INITIALIZE
  rvtotal_sigma.initialize();
  #endif
  Cpage_di.allocate(1988,last_year,first_age,last_age,"Cpage_di");
  #ifndef NO_AD_INITIALIZE
    Cpage_di.initialize();
  #endif
  RVpage_di.allocate(first_year,last_year,first_age,last_age,"RVpage_di");
  #ifndef NO_AD_INITIALIZE
    RVpage_di.initialize();
  #endif
  RVspring_page_di.allocate(1998,last_year,first_age,last_age,"RVspring_page_di");
  #ifndef NO_AD_INITIALIZE
    RVspring_page_di.initialize();
  #endif
  age1_last_year_pdf.allocate("age1_last_year_pdf");
  #ifndef NO_AD_INITIALIZE
  age1_last_year_pdf.initialize();
  #endif
  ofv.allocate("ofv");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
  M_const=M_const_prior(4);
  C_Sfullold=C_Sfullold_prior(4);
  C_log_varLold=C_log_varLold_prior(4);
  C_log_varRold=C_log_varRold_prior(4);
 
  RV_log_qold=RV_log_qold_prior(4);
  RVspring_log_q=RV_log_qold_prior(4);
 
  RV_Sfullold=RV_Sfullold_prior(4);
  RV_log_varLold=RV_log_varLold_prior(4);
  RV_log_varRold=RV_log_varRold_prior(4);
  RVspring_Sfull=RVspring_Sfull_prior(4);
  RVspring_log_varL=RVspring_log_varL_prior(4);
  RVspring_log_varR=RVspring_log_varR_prior(4);
  F_full=F_full_prior(4); 
  start_Z=start_Z_prior(4);
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
  age1(2008)=580662.6;
  age1(2009)=580662.6;
   for (age=first_age;age<=last_age;age++)
       {
       year1(age)/=1000;			//numbers in thousands
       }
   for (year=first_year+1;year<=last_year;year++)
       {
       age1(year)/=1000;
       }
   for (age=first_age;age<=last_age;age++)
       {
         for (year=1988;year<=last_year;year++)
           {
           obs_comcaa(year,age)=all_obs_comcaa(year,age)/1000;	// Convert to number in 1000's
           }
       }
   for (age=first_age;age<=last_age;age++)
       {
         for (year=first_year;year<=last_year;year++)
           {
           mat(year,age)=all_mat(year,age);
           obs_rvcaa(year,age)=all_obs_rvcaa(year,age)/1000;		//// Convert to number in 1000's
           obs_rvwaa(year,age)= all_obs_rvwaa(year,age);
           }
       }
   for (age=first_age;age<=last_age;age++)
      {
         for (year=1998;year<=2009;year++){obs_rvspring(year,age)=all_obs_rvspring(year,age);}  // raw numbers
       }
  obs_total_catch=rowsum(obs_comcaa);
  obs_rvtotal_catch=rowsum(obs_rvcaa);
  mean_rvwaa=colsum(obs_rvwaa)/(last_year-first_year+1);  // over all years
  mean_rvwaa1970=obs_rvwaa(1970)*1000;  // convert to grams
  mean_rvwaa2003=obs_rvwaa(2003)*1000;  // convert to grams
  obs_rvspring_total_catch=rowsum(obs_rvspring);
   for (age=first_age;age<=last_age;age++)
    {
    mean_rvwaa(age)=mean_rvwaa(age)*1000;  // convert to grams
    }
   for (age=first_age;age<=last_age;age++)
    {
    RVwt_age(age)=(obs_rvwaa(last_year-2,age)+obs_rvwaa(last_year-1,age)+obs_rvwaa(last_year,age))/3; 
    // over 3 most recent years
    }
   // calculating proportions at age for the catch and RV   
          for (year=1988;year<=last_year;year++)
               {
                for (age=first_age;age<=last_age;age++)
                 {
                  obs_Cpage(year,age)=obs_comcaa(year,age)/obs_total_catch(year);
                 }
                }
          for (year=first_year;year<=last_year;year++)
               {
                for (age=first_age;age<=last_age;age++)
                 {
                  obs_RVpage(year,age)=obs_rvcaa(year,age)/obs_rvtotal_catch(year);
                 }
                }
          for (age=first_age;age<=last_age;age++)
              {
              for (year=1998;year<=last_year;year++)
                  {
                  obs_RVspring_page(year,age)=obs_rvspring(year,age)/obs_rvspring_total_catch(year);
                  }
                }
 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

void model_parameters::userfunction(void)
{
  ofv =0.0;
  ofstream& report1= *pad_report1;
  ofv=0;
  age1_last_year_pdf=0;
  initialize_model();
  Population_Dynamics();
  ctotal_Likelihood();
  rvtotal_Likelihood();
  Cpage_Likelihood();
  RVpage_Likelihood();
  rvspring_total_Likelihood();
  RVspring_page_Likelihood();
  Age1_2007_pdf();
  ofv+=age1_last_year_pdf;
  ofv+=ctotal_like;
  ofv+=rvtotal_like;
  ofv+=Cpage_like;
  ofv+=RVpage_like;
  ofv+=rvspring_total_like;
  ofv+=RVspring_page_like;
  ssb2003=ssb(2003);
  if(mceval_phase()) MCWrite();
}

void model_parameters::initialize_model(void)
{
  ofstream& report1= *pad_report1;
 //fill in top of Nage matrix 
         for (year=first_year;year<=last_year;year++)      
          {
          Nage(year,first_age)=age1(year);
          }    
  //fill in left side of Nage matrix 
         for (age=first_age+1;age<=last_age;age++)      
          {
          Nage(first_year,age)=Nage(first_year,age-1)*exp(-start_Z);
          }    
         for (year=first_year;year<=last_year;year++)      
          {
          RV_select(year) = Make_select_logistic(first_age,last_age,RV_Sfullold,RV_log_varLold);
          select(year) = Make_select_logistic(first_age,last_age,C_Sfullold,C_log_varLold);
          }    
         for (year=1998;year<=last_year;year++)      
          {
          RVspring_select(year) = Make_select_logistic(first_age,last_age,RVspring_Sfull,RVspring_log_varL);
          }    
}

void model_parameters::MCWrite(void)
{
  ofstream& report1= *pad_report1;
  report1 
  << " C.Sfullold " <<  C_Sfullold
  << " C.log.varLold  " <<  C_log_varLold 
  << " C.log.varRold  " <<  C_log_varRold 
  << " RV.log.qold  " <<  RV_log_qold 
  << " RV.Sfullold  " <<  RV_Sfullold 
  << " RV.log.varLold  " <<  RV_log_varLold  
  << " RV.log.varRold   " <<  RV_log_varRold 
  <<" year1 " <<year1
  <<" age1 " <<age1
  <<" F.full " <<F_full
  <<" M " <<M
  <<" ofv " <<ofv
  <<endl;
}

void model_parameters::Population_Dynamics(void)
{
  ofstream& report1= *pad_report1;
   for (year=first_year;year<=1987;year++)		// no catch
     {
    for (age=first_age;age<=last_age;age++)
       {
        F_full(year)=0;
        F(year,age)=0;
       }
     }
   for (year=1988;year<=last_year;year++)
     {
    for (age=first_age;age<=last_age;age++)
       {
       F(year,age)=F_full(year)*select(year,age);
       }
     }
   for(year=first_year;year<=last_year;year++)
      {
    for (age=first_age;age<=last_age;age++)
       {
	// M constant
       if(mort_switch==0){M(year,age)=M_const;}
        }
      }
    for (year=first_year;year<=last_year-1;year++)
     {
      for (age=first_age;age<=last_age-1;age++)
       {
       Nage(year+1,age+1)=Nage(year,age)*exp(-(M(year,age)+F(year,age)));  
       }
     }
    for (year=first_year;year<=last_year;year++)
     {
      for (age=first_age;age<=last_age;age++)
        {
        Wtage(year,age)=Nage(year,age)*obs_rvwaa(year,age);
        matbio(year,age)=Wtage(year,age)*mat(year,age);
        RVage(year,age)=Nage(year,age)*exp(-0.5*(M(year,age)+F(year,age)))*exp(RV_log_qold)*RV_select(year,age);
        }
      }
    for (year=1988;year<=last_year;year++)
     {
      for (age=first_age;age<=last_age;age++)
        {
        Cage(year,age)=Nage(year,age)*exp(-0.5*(M(year,age)))*(1-exp(-F(year,age))); 
        }
      }
    for (year=1998;year<=last_year;year++)
      {
      for (age=first_age;age<=last_age;age++)
        {
        RVspring_age(year,age)=Nage(year,age)*exp(-0.5*(M(year,age)+F(year,age)))*exp(RVspring_log_q)*RVspring_select(year,age);
        }
      }
   pred_total_catch=rowsum(Cage);
   pred_rvtotal_catch=rowsum(RVage);
   pred_rvspring_total_catch=rowsum(RVspring_age);
   ssb=rowsum(matbio);
}

dvariable model_parameters::LogNormal_LogLike(dvar_vector observed, dvar_vector predicted, double len,double sigma, double type)
{
  ofstream& report1= *pad_report1;
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

dvariable model_parameters::calc_sigma(dvar_vector observed, dvar_vector predicted, double len)
{
  ofstream& report1= *pad_report1;
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

dvariable model_parameters::Multinomial_LogLike(dvector observed, dvar_vector predicted, double n)
{
  ofstream& report1= *pad_report1;
 //for the age data, this calcs the like for a single year and a loop over years must be done around the function call  	
 //  observed is the count in each cat
 //  predicted is the predicted proportion in each cat
   dvariable neglogLike;
   dvariable term1;
   dvariable term2;
   dvariable term3;
    term1= gammln(n+1);
    term2= sum(gammln(observed+0.001));
    term3=sum(elem_prod(observed,log(predicted+0.001)));
  neglogLike= -(term1-term2+term3);
 // cout<<neglogLike<<endl;
  return(neglogLike);
}

void model_parameters::ctotal_Likelihood(void)
{
  ofstream& report1= *pad_report1;
  dvar_vector temp_comcaa(1,n_all_comcaa_years);
  dvar_vector obs1(1,n_all_comcaa_years);
  for(i=1;i<=n_all_comcaa_years;i++)
     {
     temp_comcaa(i)=pred_total_catch(all_comcaa_years(i));   
     obs1(i)=obs_total_catch(all_comcaa_years(i));   
     }
    ctotal_like=LogNormal_LogLike(obs1+0.1,temp_comcaa+0.1,n_all_comcaa_years,sigma1,2);
 //   ctotal_sigma=calc_sigma(obs1+0.01,temp_comcaa+0.1,n_all_comcaa_years);
}

void model_parameters::Cpage_Likelihood(void)
{
  ofstream& report1= *pad_report1;
    for(year=1988;year<=last_year;year++)  
     {  
      for (age=first_age;age<=last_age;age++)
        {
        pred_Cpage(year,age)=Cage(year,age)/pred_total_catch(year);
        }
     }
    Cpage_like=-sum(elem_prod(obs_Cpage,log(pred_Cpage+.001)));
    Cpage_di=-2*elem_prod(obs_Cpage,log(pred_Cpage+.001));
}

void model_parameters::rvtotal_Likelihood(void)
{
  ofstream& report1= *pad_report1;
  dvar_vector temp_rvcaa(1,n_all_rvcaa_years);
  dvar_vector obs1(1,n_all_rvcaa_years);
  for(i=1;i<=n_all_rvcaa_years;i++)
     {
     temp_rvcaa(i)=pred_rvtotal_catch(all_rvcaa_years(i));   
     obs1(i)=obs_rvtotal_catch(all_rvcaa_years(i));   
     }
    rvtotal_like=LogNormal_LogLike(obs1+0.1,temp_rvcaa+0.1,n_all_rvcaa_years,sigma2,2);
}

void model_parameters::RVpage_Likelihood(void)
{
  ofstream& report1= *pad_report1;
     for(year=first_year;year<=last_year;year++)  
     {  
      for (age=first_age;age<=last_age;age++)
        {
        pred_RVpage(year,age)=RVage(year,age)/pred_rvtotal_catch(year);
        }
     }
       RVpage_like=-sum(elem_prod(obs_RVpage+.001,log(pred_RVpage+.001)));
       RVpage_di=-2*elem_prod(obs_RVpage+.001,log(pred_RVpage+.001));
}

void model_parameters::rvspring_total_Likelihood(void)
{
  ofstream& report1= *pad_report1;
  dvar_vector temp_rvcaa(1,n_all_rvspring_years);
  dvar_vector obs1(1,n_all_rvspring_years);
  for(i=1;i<=n_all_rvspring_years;i++)
     {
     temp_rvcaa(i)=pred_rvspring_total_catch(all_rvspring_years(i));   
     obs1(i)=obs_rvspring_total_catch(all_rvspring_years(i));   
     }
     rvspring_total_like = sum(pow(log(obs1+.1) - log(temp_rvcaa+.01),2));
}

void model_parameters::RVspring_page_Likelihood(void)
{
  ofstream& report1= *pad_report1;
     for(year=1998;year<=last_year;year++)  
     {  
      for (age=first_age;age<=last_age;age++)
        {
        pred_RVspring_page(year,age)=RVspring_age(year,age)/pred_rvspring_total_catch(year);
        }
     }
       RVspring_page_like=-sum(elem_prod(obs_RVspring_page+.001,log(pred_RVspring_page+.001)));
       RVspring_page_di=-2*elem_prod(obs_RVspring_page+.001,log(pred_RVspring_page+.001));
}

dvariable model_parameters::Normal(dvariable observed,dvariable mean,double nsigma)
{
  ofstream& report1= *pad_report1;
   dvariable normal;
   double pi=3.1415926535897;
   dvariable neglogpdf;
   //pdf=(1/(sqrt(2*pi)*nsigma)*exp(-(pow(observed-mean,2)/(2*pow(nsigma,2)))));
   //neglogpdf=-log(pdf);
    neglogpdf=log(sqrt(2*pi)*nsigma)+pow(observed-mean,2)/(2*pow(nsigma,2));
   return(neglogpdf);
}

dvar_vector model_parameters::Make_select_dhn(int minage, int maxage, dvariable full, dvariable varL, dvariable varR)
{
  ofstream& report1= *pad_report1;
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
  dvar_vector selection(minage,maxage);
  selection.initialize();
  for (age=minage;age<=maxage;age++)
  {
  selection(age)=1/(1+exp(-varL*(age-full))); 
  }	
  selection = selection/max(selection);
  return(selection);
}

void model_parameters::Age1_2007_pdf(void)
{
  ofstream& report1= *pad_report1;
  dvariable age1_last_year;
  double age1_last_year_std;
  dvariable temp;
  dvariable temp2;
  temp=0;
  for(year=last_year-4;year<=last_year;year++)  
    {  
    temp+=Nage(year,1);
    }
  temp2=temp/5;
  age1_last_year=temp2;
  age1_last_year_std=.4;
  age1_last_year_pdf=Normal(log(age1(last_year)),log(age1_last_year),age1_last_year_std);
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
    report << "r1" << endl;
    report << "first.year" << endl;
    report << first_year << endl;
    report << "r1" << endl;
    report << "last.year" << endl;
    report << last_year << endl;
    report << "r1" << endl;
    report << "first.age" << endl;
    report << first_age << endl;
    report << "r1" << endl;
    report << "last.age" << endl;
    report << last_age << endl;
    report << "r1" << endl;
    report << "n.landings" << endl;
    report << n_landings << endl;
    report << "r1" << endl;
    report << "landing.year" << endl;
    report << landing_year << endl;
    report << "r1" << endl;
    report << "landings" << endl;
    report << landings << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "mat" << endl;
    report << mat << endl;
    report << "r"<<last_year-1988+1 << endl;
    report << "obs.comcaa" << endl;
    report << obs_comcaa << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "obs.rvcaa" << endl;
    report << obs_rvcaa << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "obs.rvwaa" << endl;
    report << obs_rvwaa << endl;
    report << "r"<<last_year-1998+1 << endl;
    report << "obs.rvspring" << endl;
    report << obs_rvspring << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "Nage" << endl;
    report << Nage << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "M" << endl;
    report << M << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "Cage" << endl;
    report << Cage << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "RVage" << endl;
    report << RVage << endl;
    report << "r"<<last_year-1998+1 << endl;
    report << "RVspring.age" << endl;
    report << RVspring_age << endl;
    report << "r1" << endl;
    report << "obs.rvspring.total.catch" << endl;
    report << obs_rvspring_total_catch << endl;
    report << "r1" << endl;
    report << "pred.rvspring.total.catch" << endl;
    report << pred_rvspring_total_catch << endl;
    report << "r"<<last_year-1988+1 << endl;
    report << "obs.Cpage" << endl;
    report << obs_Cpage << endl;
    report << "r"<<last_year-1988+1 << endl;
    report << "pred.Cpage" << endl;
    report << pred_Cpage << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "obs.RVpage" << endl;
    report << obs_RVpage << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "pred.RVpage" << endl;
    report << pred_RVpage << endl;
    report << "r"<<last_year-1988+1 << endl;
    report << "Cpage.di" << endl;
    report << Cpage_di << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "RVpage.di" << endl;
    report << RVpage_di << endl;
    report << "r"<<last_year-1998+1 << endl;
    report << "RVspring.page.di" << endl;
    report << RVspring_page_di << endl;
    report << "r"<<last_year-1998+1 << endl;
    report << "obs.RVspring.page" << endl;
    report << obs_RVspring_page << endl;
    report << "r"<<last_year-1998+1 << endl;
    report << "pred.RVspring.page" << endl;
    report << pred_RVspring_page << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "F" << endl;
    report << F << endl;
    report << "r1" << endl;
    report << "F.full" << endl;
    report << F_full << endl;
    report << "r1" << endl;
    report << "ssb" << endl;
    report << ssb << endl;
    report << "r1" << endl;
    report << "C.Sfullold" << endl;
    report << C_Sfullold << endl;
    report << "r1" << endl;
    report << "C.log.varLold" << endl;
    report << C_log_varLold << endl;
    report << "r1" << endl;
    report << "C.log.varRold" << endl;
    report << C_log_varRold << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "select" << endl;
    report << select << endl;
    report << "r1" << endl;
    report << "RV.Sfullold" << endl;
    report << RV_Sfullold << endl;
    report << "r1" << endl;
    report << "RV.log.varLold" << endl;
    report << RV_log_varLold << endl;
    report << "r1" << endl;
    report << "RV.log.varRold" << endl;
    report << RV_log_varRold << endl;
    report << "r1" << endl;
    report << "RVspring.Sfull" << endl;
    report << RVspring_Sfull << endl;
    report << "r1" << endl;
    report << "RVspring.log.varL" << endl;
    report << RVspring_log_varL << endl;
    report << "r1" << endl;
    report << "RVspring.log.varR" << endl;
    report << RVspring_log_varR << endl;
    report << "r1" << endl;
    report << "RV.log.qold" << endl;
    report << RV_log_qold << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "RV.select" << endl;
    report << RV_select << endl;
    report << "r"<<last_year-1998+1 << endl;
    report << "RVspring.select" << endl;
    report << RVspring_select << endl;
    report << "r1" << endl;
    report << "obs.total.catch" << endl;
    report << obs_total_catch << endl;
    report << "r1" << endl;
    report << "pred.total.catch" << endl;
    report << pred_total_catch << endl;
    report << "r1" << endl;
    report << "year1" << endl;
    report << year1 << endl;
    report << "r1" << endl;
    report << "mean.RVwt.age" << endl;
    report << RVwt_age << endl;
    report << "r1" << endl;
    report << "mean.rvwaa" << endl;
    report << mean_rvwaa << endl;
    report << "r"<<last_year-first_year+1 << endl;
    report << "Wtage" << endl;
    report << Wtage << endl;
    report << "r1" << endl;
    report << "ctotal.like" << endl;
    report << ctotal_like << endl;
    report << "r1" << endl;
    report << "Cpage.like" << endl;
    report << Cpage_like << endl;
    report << "r1" << endl;
    report << "rvtotal.like" << endl;
    report << rvtotal_like << endl;
    report << "r1" << endl;
    report << "RVpage.like" << endl;
    report << RVpage_like << endl;
    report << "r1" << endl;
    report << "rvspring_total.like" << endl;
    report << rvspring_total_like << endl;
    report << "r1" << endl;
    report << "RVspring_page.like" << endl;
    report << RVspring_page_like << endl;
    report << "r1" << endl;
    report << "ctotal.sigma" << endl;
    report << ctotal_sigma << endl;
    report << "r1" << endl;
    report << "rvtotal.sigma" << endl;
    report << rvtotal_sigma << endl;
    report << "r1" << endl;
    report << "ofv" << endl;
    report << ofv << endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_report1;
  pad_report1 = NULL;
}

void model_parameters::final_calcs(void){}

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
  arrmblsize=50000000;
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
  
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
