//modifed for halibut Oct 2010  sexs and gears combined
//MODEL FOR 4WVS Cod
//JAMIE GIBSON MAY 25, 2002
//modified by K. TRZCINSKI and Jamie Gibson 2004
//This model is uses selectivity by weight or age set by select_switch

//to do mcmc: 
//delete mcout.dat
//caa1 -mcmc 120000 -mcsave 100 // del first 200 rows in mcout.dat for burn in.
//caa1 -mceval


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TOP_OF_MAIN_SECTION
  time(&start);
  arrmblsize=50000000;
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
  
//======================================================================
DATA_SECTION

!!CLASS ofstream report1("mcout.dat",ios::app);  // you must delete this file before a new run
  init_int debug	// 0=don't debug input; 1=debug input 
  init_int all_first_year
  !!if (debug==1){char mark; cout<<"all_first_year= "<<all_first_year<<endl; cin>>mark; }
  init_int all_last_year
  init_int all_first_age	// First age for model (recruitment is at this age)
  init_int all_last_age	// Last age for model
  init_int select_switch   //0 = by weight; 1 = by age 
  init_vector select_data(1,4)    //years and ages to include
  int first_year
  int last_year
  int first_age
  int last_age 
  !! first_year=select_data(1);last_year=select_data(2);first_age=select_data(3);last_age=select_data(4);

  // ** Bounds and Starting values**       
  //phase LB UB defaultvalue
  init_vector F_full_prior(1,4)
  !!if (debug==1){char mark; cout<<"F_full_prior= "<<F_full_prior<<endl; cin>>mark; }
  init_vector start_Z_prior(1,4)
  init_vector year1_prior(1,4)
  init_vector age1_prior(1,4)
  init_int mort_switch
  init_vector M_const_prior(1,4)
  init_vector C_Sfullold_prior(1,4)
  init_vector C_log_varLold_prior(1,4) 
  init_vector C_log_varRold_prior(1,4) 
  init_vector RV_Sfullold_prior(1,4)
  init_vector RV_log_varLold_prior(1,4) 
  init_vector RV_log_varRold_prior(1,4) 
  init_vector RVspring_Sfull_prior(1,4) 
  init_vector RVspring_log_varL_prior(1,4) 
  init_vector RVspring_log_varR_prior(1,4) 
  init_vector RV_log_qold_prior(1,4) 
  init_vector RVspring_log_q_prior(1,4) 
    
  //**FIXED VALUES**

  init_number sigma1
  init_number sigma2
  init_vector like_weight(1,4)

  // ** DATA **
  init_number n_landings
  init_vector landings(1,n_landings) 
  init_vector landing_year(1,n_landings)
  !!if (debug==1){char mark; cout<<"landing_year= "<<landing_year<<endl; cin>>mark; }
  init_matrix all_mat(all_first_year,all_last_year,all_first_age,all_last_age)   // 
  init_int n_all_comcaa_years
  init_vector all_comcaa_years(1,n_all_comcaa_years)
  !!if (debug==1){char mark; cout<<"all_comcaa_years= "<<all_comcaa_years<<endl; cin>>mark; }
  init_matrix all_obs_comcaa(1988,all_last_year,all_first_age,all_last_age) // numbers caught in the commercial fishery
  !!if (debug==1){char mark; cout<<"all_obs_comcaa= "<<all_obs_comcaa(1988)<<endl; cin>>mark; }
  init_int n_all_rvcaa_years
  init_vector all_rvcaa_years(1,n_all_rvcaa_years)
  !!if (debug==1){char mark; cout<<"all_rvcaa_years= "<<all_rvcaa_years<<endl; cin>>mark; }
  init_matrix all_obs_rvcaa(all_first_year,all_last_year,all_first_age,all_last_age) // numbers caught in the rv survey
  !!if (debug==1){char mark; cout<<"all_obs_rvcaa= "<<all_obs_rvcaa(1970)<<endl; cin>>mark; }
  init_matrix all_obs_rvwaa(all_first_year,all_last_year,all_first_age,all_last_age) // mean weight at age in the rv survey
  !!if (debug==1){char mark; cout<<"all_obs_rvwaa= "<<all_obs_rvwaa(1970)<<endl; cin>>mark; }
  init_int n_all_rvspring_years
  !!if (debug==1){char mark; cout<<"n_all_rvspring_years= "<<n_all_rvspring_years<<endl; cin>>mark; }
  init_vector all_rvspring_years(1,n_all_rvspring_years)
  init_matrix all_obs_rvspring(1998,all_last_year,all_first_age,all_last_age)  
  init_int dummy
  !!if (debug==1){char mark; cout<<"dummy= "<<dummy<<endl; cin>>mark; }
 

  int j
  int i
  int age
  int year
  int phz
  number lb
  number ub
  int len
  int qlen
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PARAMETER_SECTION
  matrix obs_comcaa(1988,last_year,first_age,last_age)
  matrix obs_rvcaa(first_year,last_year,first_age,last_age)
  matrix obs_rvwaa(first_year,last_year,first_age,last_age)
  matrix obs_rvspring(1998,2009,first_age,last_age)
  matrix mat(first_year,last_year,first_age,last_age)
//  sdreport_matrix Nage(first_year,last_year,first_age,last_age) //numbers at age at start of year
  matrix Nage(first_year,last_year,first_age,last_age) //numbers at age at start of year
  matrix Cage(first_year,last_year,first_age,last_age)
  matrix RVage(first_year,last_year,first_age,last_age)
  matrix RVspring_age(1998,2009,first_age,last_age)
  matrix Wtage(first_year,last_year,first_age,last_age) //using RV weights
  matrix matbio(first_year,last_year,first_age,last_age) //using RV weights
  matrix obs_Cpage(1988,last_year,first_age,last_age)   // proportions at age
  matrix obs_RVpage(first_year,last_year,first_age,last_age)  // proportions at age 
  matrix obs_RVspring_page(1998,last_year,first_age,last_age)  // proportions at age 
  matrix pred_Cpage(1988,last_year,first_age,last_age)   // proportions at age
  matrix pred_RVpage(first_year,last_year,first_age,last_age)  // proportions at age 
  matrix pred_RVspring_page(1998,last_year,first_age,last_age)  // proportions at age 
  matrix F(first_year,last_year,first_age,last_age) 
  matrix select(first_year,last_year,first_age,last_age)
  matrix RV_select(first_year,last_year,first_age,last_age)
  matrix RVspring_select(1998,2009,first_age,last_age)


  vector obs_total_catch(1988,last_year)
  vector pred_total_catch(first_year,last_year)
  vector obs_rvtotal_catch(first_year,last_year)
  vector pred_rvtotal_catch(first_year,last_year)
  vector obs_rvspring_total_catch(1998,2009)
  vector pred_rvspring_total_catch(1998,2009)

//  vector ssb(first_year,last_year)  //at start of year
  sdreport_vector ssb(first_year,last_year)  //at start of year
  vector RVwt_age(first_age,last_age)
 
  !!phz=F_full_prior(1);lb=F_full_prior(2);ub=F_full_prior(3);
   init_bounded_vector F_full(first_year,last_year,lb,ub,phz)
  !!phz=start_Z_prior(1);lb=start_Z_prior(2);ub=start_Z_prior(3);
  init_bounded_number start_Z(lb,ub,phz)
  !!phz=year1_prior(1);lb=year1_prior(2);ub=year1_prior(3);
   init_bounded_vector year1(first_age,last_age,lb,ub,phz) 
   !!phz=age1_prior(1);lb=age1_prior(2);ub=age1_prior(3);
   init_bounded_vector age1(first_year,last_year,lb,ub,phz) 

  !!phz=M_const_prior(1);lb=M_const_prior(2);ub=M_const_prior(3);
   init_bounded_number M_const(lb,ub,phz)

  !!phz=C_Sfullold_prior(1);lb=C_Sfullold_prior(2);ub=C_Sfullold_prior(3);
   init_bounded_number C_Sfullold(lb,ub,phz)
  !!phz=C_log_varLold_prior(1);lb=C_log_varLold_prior(2);ub=C_log_varLold_prior(3);
   init_bounded_number C_log_varLold(lb,ub,phz)
  !!phz=C_log_varRold_prior(1);lb=C_log_varRold_prior(2);ub=C_log_varRold_prior(3);
   init_bounded_number C_log_varRold(lb,ub,phz)
  !!phz=RV_Sfullold_prior(1);lb=RV_Sfullold_prior(2);ub=RV_Sfullold_prior(3);
   init_bounded_number RV_Sfullold(lb,ub,phz)
  !!phz=RV_log_varLold_prior(1);lb=RV_log_varLold_prior(2);ub=RV_log_varLold_prior(3);
   init_bounded_number RV_log_varLold(lb,ub,phz)
  !!phz=RV_log_varRold_prior(1);lb=RV_log_varRold_prior(2);ub=RV_log_varRold_prior(3);
   init_bounded_number RV_log_varRold(lb,ub,phz)

  !!phz=RVspring_Sfull_prior(1);lb=RVspring_Sfull_prior(2);ub=RVspring_Sfull_prior(3);
   init_bounded_number RVspring_Sfull(lb,ub,phz)
  !!phz=RVspring_log_varL_prior(1);lb=RVspring_log_varL_prior(2);ub=RVspring_log_varL_prior(3);
   init_bounded_number RVspring_log_varL(lb,ub,phz)
  !!phz=RVspring_log_varR_prior(1);lb=RVspring_log_varR_prior(2);ub=RVspring_log_varR_prior(3);
   init_bounded_number RVspring_log_varR(lb,ub,phz)

  !!phz=RV_log_qold_prior(1);lb=RV_log_qold_prior(2);ub=RV_log_qold_prior(3);
   init_bounded_number RV_log_qold(lb,ub,phz)
  !!phz=RVspring_log_q_prior(1);lb=RVspring_log_q_prior(2);ub=RVspring_log_q_prior(3);
   init_bounded_number RVspring_log_q(lb,ub,phz)

  vector mean_rvwaa(first_age,last_age)
  vector mean_rvwaa1970(first_age,last_age)
  vector mean_rvwaa2003(first_age,last_age)

  matrix M(first_year,last_year,first_age,last_age)

//  likeprof_number M1
//  likeprof_number M2
    likeprof_number ssb2003

//used in likelihoods
  number ctotal_like
  number Cpage_like
  number rvtotal_like
  number RVpage_like

  number rvspring_total_like
  number RVspring_page_like

  number ctotal_sigma
  number rvtotal_sigma

  matrix Cpage_di(1988,last_year,first_age,last_age)
  matrix RVpage_di(first_year,last_year,first_age,last_age)
  matrix RVspring_page_di(1998,last_year,first_age,last_age)

  number age1_last_year_pdf

  objective_function_value ofv 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PRELIMINARY_CALCS_SECTION
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
//  year1=year1_prior(4);
//  age1=age1_prior(4);

// a<-0:19
//round(580662.6*exp(-.2)^a)

//      year1(1)= 580663; 
//      year1(2)= 475406; 
//      year1(3)= 389230; 
//      year1(4)= 318674; 
//      year1(5)= 260909; 
//      year1(6)= 213614; 
//      year1(7)= 174892; 
//      year1(8)= 143190; 
//      year1(9)= 117234; 
//      year1(10)= 95983; 
//      year1(11)= 78584; 
//      year1(12)= 64339; 
//      year1(13)= 52677; 
//      year1(14)= 43128; 
//      year1(15)= 35310; 
//      year1(16)= 28909; 
//      year1(17)= 23669; 
//      year1(18)= 19379; 
//      year1(19)= 15866; 
//      year1(20)= 12990; 

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


//cut out the observed data for the selected period  


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

PROCEDURE_SECTION
  ofv=0;
  age1_last_year_pdf=0;
  initialize_model();
//  cout<<"OK"<<endl;
  Population_Dynamics();
//  cout<<"OK"<<endl;

  ctotal_Likelihood();
  rvtotal_Likelihood();
  Cpage_Likelihood();
  RVpage_Likelihood();

  rvspring_total_Likelihood();
  RVspring_page_Likelihood();

  Age1_2007_pdf();

//  ofv+=ctotal_like*like_weight(1);
//  ofv+=rvtotal_like*like_weight(2);
//  ofv+=Cpage_like*like_weight(3);
//  ofv+=RVpage_like*like_weight(4);

  ofv+=age1_last_year_pdf;

  ofv+=ctotal_like;
  ofv+=rvtotal_like;
  ofv+=Cpage_like;
  ofv+=RVpage_like;

  ofv+=rvspring_total_like;
  ofv+=RVspring_page_like;

  ssb2003=ssb(2003);


  if(mceval_phase()) MCWrite();
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FUNCTION initialize_model  

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
//          RV_select(year) = Make_select_dhn(first_age,last_age,RV_Sfullold,RV_log_varLold,RV_log_varRold);
          select(year) = Make_select_logistic(first_age,last_age,C_Sfullold,C_log_varLold);
          }    

         for (year=1998;year<=last_year;year++)      
          {
          RVspring_select(year) = Make_select_logistic(first_age,last_age,RVspring_Sfull,RVspring_log_varL);
          }    




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION MCWrite
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

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FUNCTION Population_Dynamics

   for (year=first_year;year<=1987;year++)		// no catch
     {
    for (age=first_age;age<=last_age;age++)
       {
        F_full(year)=0;
        F(year,age)=0;
//   cout<<F(year,10)<<endl;
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


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//######################################################################################################################
//======================================================
FUNCTION dvariable LogNormal_LogLike(dvar_vector observed, dvar_vector predicted, double len,double sigma, double type)
// I changed observed to type dvar_vector.  This allowed me to reassign indicies.
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
//  sumlogr = sum(log(observed));
//  sumsqrlogratio = sum( square(  log(observed)-log(predicted) )  );  
//  neglogLike=sumsqrlogratio;
     if(type==2)
     {
     neglogLike = 1/(2*square(sigma))*sumsqrlogratio   + sumlogr   + len2*log(sigma)*sqrt(2*pi);  //specified sigma used
     }

    if(type==1 &last_phase())
    {
    neglogLike = 0.5*len2*log( (2*pi/len2)*sumsqrlogratio) + sumlogr + len2/2;  //est. sigma used
    }
  
  return(neglogLike);


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

//======================================================
FUNCTION dvariable Multinomial_LogLike(dvector observed, dvar_vector predicted, double n)
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
  
//   cout<<"term1= "<<term1<<endl;
//   cout<<"term2= "<<term2<<endl;
//   cout<<"term3= "<<term3<<endl;



//THIS ALLOWS YOU TO REM OUT VALUES IN THE DAT FILE BY MAKING THEM NEGATIVE:
 
  
//  neglogLike= -( gammln(n+1)-sum(gammln(observed)) + sum(elem_prod(observed,log(predicted+0.001)))   );     
 
  neglogLike= -(term1-term2+term3);

 // cout<<neglogLike<<endl;

  return(neglogLike);
 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FUNCTION ctotal_Likelihood
//this reassigns indices

  dvar_vector temp_comcaa(1,n_all_comcaa_years);
  dvar_vector obs1(1,n_all_comcaa_years);

  for(i=1;i<=n_all_comcaa_years;i++)
     {
     temp_comcaa(i)=pred_total_catch(all_comcaa_years(i));   
     obs1(i)=obs_total_catch(all_comcaa_years(i));   
     }
    ctotal_like=LogNormal_LogLike(obs1+0.1,temp_comcaa+0.1,n_all_comcaa_years,sigma1,2);
 //   ctotal_sigma=calc_sigma(obs1+0.01,temp_comcaa+0.1,n_all_comcaa_years);
 
//   ctotal_like = sum(pow(log(obs1+.1) - log(temp_comcaa+.1),2));


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION Cpage_Likelihood    

    for(year=1988;year<=last_year;year++)  
     {  
      for (age=first_age;age<=last_age;age++)
        {
        pred_Cpage(year,age)=Cage(year,age)/pred_total_catch(year);
        }
     }

    Cpage_like=-sum(elem_prod(obs_Cpage,log(pred_Cpage+.001)));
    Cpage_di=-2*elem_prod(obs_Cpage,log(pred_Cpage+.001));


//  double n;
//  dvector observed(first_age,last_age);
//  dvar_vector predicted(first_age,last_age);

//  Cpage_like=0;

//   for(year=first_year;year<=last_year;year++)  
//     {
//    for (age=first_age; age<=last_age;age++)    
//     {
//      observed(age)=all_obs_comcaa(year,age); 
//      predicted(age)=Cage(year,age)/pred_total_catch(year);
//     }

//      n=sum(observed);
//      Cpage_like+=Multinomial_LogLike(observed,predicted,n);
//     }

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FUNCTION rvtotal_Likelihood    
  dvar_vector temp_rvcaa(1,n_all_rvcaa_years);
  dvar_vector obs1(1,n_all_rvcaa_years);

  for(i=1;i<=n_all_rvcaa_years;i++)
     {
     temp_rvcaa(i)=pred_rvtotal_catch(all_rvcaa_years(i));   
     obs1(i)=obs_rvtotal_catch(all_rvcaa_years(i));   
     }
    rvtotal_like=LogNormal_LogLike(obs1+0.1,temp_rvcaa+0.1,n_all_rvcaa_years,sigma2,2);
//    rvtotal_sigma=calc_sigma(obs1+0.01,temp_rvcaa+0.1,n_all_rvcaa_years);

//     rvtotal_like = sum(pow(log(obs1+.1) - log(temp_rvcaa+.1),2));


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION RVpage_Likelihood    

     for(year=first_year;year<=last_year;year++)  
     {  
      for (age=first_age;age<=last_age;age++)
        {
        pred_RVpage(year,age)=RVage(year,age)/pred_rvtotal_catch(year);
        }
     }
       RVpage_like=-sum(elem_prod(obs_RVpage+.001,log(pred_RVpage+.001)));
       RVpage_di=-2*elem_prod(obs_RVpage+.001,log(pred_RVpage+.001));


//  double n;
//  dvector observed(first_age,last_age);
//  dvar_vector predicted(first_age,last_age);

//  RVpage_like=0;

//   for(year=first_year;year<=last_year;year++)  
//     {
//    for (age=first_age; age<=last_age;age++)    
//     {
//      observed(age)=all_obs_rvcaa(year,age); 
//      predicted(age)=RVage(year,age)/pred_rvtotal_catch(year);
//     }

//      n=sum(observed);
//      RVpage_like+=Multinomial_LogLike(observed,predicted,n);
//     }



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FUNCTION rvspring_total_Likelihood    
  dvar_vector temp_rvcaa(1,n_all_rvspring_years);
  dvar_vector obs1(1,n_all_rvspring_years);

  for(i=1;i<=n_all_rvspring_years;i++)
     {
     temp_rvcaa(i)=pred_rvspring_total_catch(all_rvspring_years(i));   
     obs1(i)=obs_rvspring_total_catch(all_rvspring_years(i));   
     }
//    rvspring_total_like=LogNormal_LogLike(obs1+0.1,temp_rvcaa+0.1,n_all_rvcaa_years,.2,2);
//    rvtotal_sigma=calc_sigma(obs1+0.01,temp_rvcaa+0.1,n_all_rvcaa_years);

     rvspring_total_like = sum(pow(log(obs1+.1) - log(temp_rvcaa+.01),2));

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION RVspring_page_Likelihood    

     for(year=1998;year<=last_year;year++)  
     {  
      for (age=first_age;age<=last_age;age++)
        {
        pred_RVspring_page(year,age)=RVspring_age(year,age)/pred_rvspring_total_catch(year);
        }
     }
       RVspring_page_like=-sum(elem_prod(obs_RVspring_page+.001,log(pred_RVspring_page+.001)));
       RVspring_page_di=-2*elem_prod(obs_RVspring_page+.001,log(pred_RVspring_page+.001));


//======================================================
FUNCTION dvariable Normal(dvariable observed,dvariable mean,double nsigma)

   dvariable normal;
   double pi=3.1415926535897;
   dvariable neglogpdf;
   //pdf=(1/(sqrt(2*pi)*nsigma)*exp(-(pow(observed-mean,2)/(2*pow(nsigma,2)))));
   //neglogpdf=-log(pdf);
    neglogpdf=log(sqrt(2*pi)*nsigma)+pow(observed-mean,2)/(2*pow(nsigma,2));
   return(neglogpdf);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FUNCTION Age1_2007_pdf
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
 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

REPORT_SECTION  
   REPORT(first_year);
   REPORT(last_year);
   REPORT(first_age);
   REPORT(last_age);
   REPORT(n_landings);
   REPORT(landing_year);
   REPORT(landings);
   REPORT(mat);
   REPORT(obs_comcaa);
   REPORT(obs_rvcaa);
   REPORT(obs_rvwaa);
   REPORT(obs_rvspring);
   REPORT(Nage);
   REPORT(M);
   REPORT(Cage);
   REPORT(RVage);
   REPORT(RVspring_age);
   REPORT(obs_rvspring_total_catch);
   REPORT(pred_rvspring_total_catch);
   REPORT(obs_Cpage);
   REPORT(pred_Cpage);
   REPORT(obs_RVpage);
   REPORT(pred_RVpage);
   REPORT(Cpage_di);
   REPORT(RVpage_di);
   REPORT(RVspring_page_di);
   REPORT(obs_RVspring_page);
   REPORT(pred_RVspring_page);
   REPORT(F);
   REPORT(F_full);
   REPORT(ssb);
   REPORT(C_Sfullold);
   REPORT(C_log_varLold);
   REPORT(C_log_varRold);
   REPORT(select);
   REPORT(RV_Sfullold);
   REPORT(RV_log_varLold);
   REPORT(RV_log_varRold);
   REPORT(RVspring_Sfull);
   REPORT(RVspring_log_varL);
   REPORT(RVspring_log_varR);
   REPORT(RV_log_qold);
   REPORT(RV_select);
   REPORT(obs_total_catch);
   REPORT(pred_total_catch);
   REPORT(year1);
   REPORT(RVwt_age);
   REPORT(mean_rvwaa);
   REPORT(Wtage);
   REPORT(ctotal_like);
   REPORT(Cpage_like);
   REPORT(rvtotal_like);
   REPORT(RVpage_like);
   REPORT(rvspring_total_like);
   REPORT(RVspring_page_like);
   REPORT(ctotal_sigma);
   REPORT(rvtotal_sigma);
   REPORT(ofv);


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


