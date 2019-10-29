/*******************************************************************
EPI 202 - Crude & Stratified Analyses for Person-time Data 
Created by: Sam Molsberry
Date: 9/1/2018
Dataset: Evans County example data

Reminders before beginning:
 1) If you would like to recreate the example using the Evans County data done in the documentation,
    you do not need to change any of this code aside from the file path where you have saved the csv fileT
 2) If you would like to change this code to do analyses on a different dataset, please adhere to the following:
	a) The outcome variables must be binary. To change this variable, please change the 
	   variable names and values indicating case and noncase as necessary 
	b) You must have a variable representing the amount of person-time that each subject contributes under
	   under follow-up. Please change the name of your person-time variable as necessary.
	c) The stratification variable(s) must be categorical. To change these variables, please change the 
	   variables names as necessary
	d) You do not need to change the names of any of the outputted objects, but if you do so  make sure to
	   also change any subsequent references to those objects
 3) Please see example documentation for a description of the variables in the example dataset
 4) Make sure not to uncomment any of the notes throughout (which appear in green text) or the code will not run
*******************************************************************/

/*Call in data - make sure to change file path to the location where 
your data is saved*/
proc import out= WORK.dat DATAFILE= "P:\Epi 202 2018\evans_example_dat.csv" 
            DBMS=csv REPLACE;
     GETNAMES=YES;
RUN;

/************************************/
/**Crude Analysis: Person-time data**/
/************************************/

data dat;
set dat;
if HTN^=.; *Change HTN to exposure variable and . to value indicating missing value as necessary;

/************************************************************************************************/
/*Next four lines create variables for each type of case and person-time (exposed or unexposed). 
To use a different data sets change 'HTN' to name of you exposure variable and 'CHD' to name of your 
outcome variable. Additionally change to values of 1 and 0 as necessary to be consistent with the 
values indicating, exposed v. unexposed and case v. noncase in your data. Do not change the value
of 1 that exposed_case and exposed_noncasee are set to.

Similarly, when creating the indicators for exposed person-time, change the name and value of your
exposure variable as well as the name of your person-time variable as necessary. In the example here,
the person-time variable is simply called 'person_time'.

After changing these variables and values as necessary, you do not need to make any further changes
to the crude analysis section of the code in order to complete the crude analyses*/
/************************************************************************************************/ 
if HTN=1 and CHD=1 then exposed_case=1;
if HTN=0 and CHD=1 then unexposed_case=1;

if HTN=1 then exposed_PT=person_time;
if HTN=0 then unexposed_PT=person_time;
run;

proc means data=dat sum maxdec=2;
ods output Summary=crude_open;
var exposed_case unexposed_case exposed_PT unexposed_PT;
run;

data crude_open;
set crude_open;
rename exposed_case_Sum=a
	unexposed_case_Sum=b
	exposed_PT_Sum=N1
	unexposed_PT_Sum=N0;
drop VNAME_exposed_case VNAME_unexposed_case VNAME_exposed_PT VNAME_unexposed_PT;
run; 

data crude_open;
set crude_open;
label a="a"
	b="b"
	N1="N1"
	N0="N0";
run;

data crude_open;
set crude_open;
/*Specify count related variables*/
M1=a+b;
T=N1+N0;

/*Calculate IRR & 95% CI*/
IRR=(a/N1)/(b/N0);
Variance_IRR_confint=(1/a)+(1/b);
IRR_confint_lower=round(exp(log(IRR)-1.96*sqrt(Variance_IRR_confint)),0.001);
IRR_confint_upper=round(exp(log(IRR)+1.96*sqrt(Variance_IRR_confint)),0.001);
IRR_confint=catx(',', IRR_confint_lower, IRR_confint_upper);

/*Calculate IRD & 95% CI*/
IRD=(a/N1)-(b/N0);
Variance_IRD_confint=(a/(N1**2))+(b/(N0**2));
IRD_confint_lower=round(IRD-1.96*sqrt(Variance_IRD_confint), 0.001);
IRD_confint_upper=round(IRD+1.96*sqrt(Variance_IRD_confint), 0.001);
IRD_confint=catx(',', IRD_confint_lower, IRD_confint_upper);

/*Test of Association*/
Expected_null=(N1*M1)/T;
Variance_testassoc=(M1*N1*N0)/(T**2);
Z_squared=((a-Expected_null)**2)/Variance_testassoc;
p_testassoc=1-probchi(Z_squared, 1);
run;

proc print data=crude_open;
var IRR IRR_confint IRD IRD_confint;
run;

proc print data=crude_open;
var Z_squared p_testassoc;
run;

/***************************************/
/*Stratified Analysis: Person-time data*/
/***************************************/

/************************************************************************************************/
/* In order conduct the stratified analyses, you must first run the initial data step of the crude
analysis section if you have not already done so in order to produce the dataset named 'dat'. 
Additionally, if you would like to adapt this code to be used on a different dataset, the only 
thing that must be changed now is to change the variable name 'SMK' to a list of your stratification 
factor(s) in the 'by' line of the proc sort and proc means steps below. Once you have done that, 
you do not need to change any of the following code in order to get the complete results of the 
stratified analyses.*/
/************************************************************************************************/

proc sort data=dat;
by SMK; *Change SMK to name(s) of your stratification factor(s);
run;

proc means data=dat sum maxdec=2;
ods output Summary=stratified_open;
var exposed_case unexposed_case exposed_PT unexposed_PT;
by SMK; *Change SMK to name(s) of your stratification factor(s);
run;

data stratified_open;
set stratified_open;
rename exposed_case_Sum=a_i
	unexposed_case_Sum=b_i
	exposed_PT_Sum=N1_i
	unexposed_PT_Sum=N0_i;
drop VNAME_exposed_case VNAME_unexposed_case VNAME_exposed_PT VNAME_unexposed_PT;
run; 

data stratified_open;
set stratified_open;
label a="a_i"
	b="b_i"
	N1="N1_i"
	N0="N0_i";
run; 

data stratified_open;
set stratified_open;
/*Specify count related variables*/
M1_i=a_i+b_i;
T_i=N1_i+N0_i;

/*Calculate IRR & 95% CI*/
IRR_i=(a_i/N1_i)/(b_i/N0_i);
Variance_IRR_confint_i=(1/a_i)+(1/b_i);
IRR_confint_lower_i=round(exp(log(IRR_i)-1.96*sqrt(Variance_IRR_confint_i)),0.001);
IRR_confint_upper_i=round(exp(log(IRR_i)+1.96*sqrt(Variance_IRR_confint_i)),0.001);
IRR_confint_i=catx(',', IRR_confint_lower_i, IRR_confint_upper_i);

/*Summary: Calculate MH IRR & 95% CI*/
MH_w_i=(b_i*N1_i)/T_i;
MH_numerator_term_i=(a_i*N0_i)/T_i;
MH_variance_numerator_i=(M1_i*N1_i*N0_i)/(T_i**2);
MH_variance_denom1_i=(a_i*N0_i)/T_i;
MH_variance_denom2_i=(b_i*N1_i)/T_i;

/*Calculate IRD & 95% CI*/
IRD_i=(a_i/N1_i)-(b_i/N0_i);
Variance_IRD_confint_i=(a_i/(N1_i**2))+(b_i/(N0_i**2));
IRD_confint_lower_i=round(IRD_i-1.96*sqrt(Variance_IRD_confint_i), 0.001);
IRD_confint_upper_i=round(IRD_i+1.96*sqrt(Variance_IRD_confint_i), 0.001);
IRD_confint_i=catx(',', IRD_confint_lower_i, IRD_confint_upper_i);

/*Summary: Calculate Inverse Variance IRD & 95% CI*/
IV_w_i=((N1_i**2)*(N0_i**2))/((a_i*(N0_i**2))+(b_i*(N1_i**2)));
IV_numerator_term_i=IV_w_i*IRD_i;
drop IRR_confint_lower_i IRR_confint_upper_i IRR_confint_lower_i IRR_confint_upper_i;

/*Test of No association*/
expected_value_i=(N1_i*M1_i)/T_i;
variance_term_i=(M1_i*N1_i*N0_i)/(T_i**2);
run;
/*Calculate MH IRR & 95% CI*/
proc sql;
title 'MH IRR';
select sum(MH_numerator_term_i)/sum(MH_w_i) as MH_IRR,
		exp(log(sum(MH_numerator_term_i)/sum(MH_w_i))-1.96*sqrt(sum(MH_variance_numerator_i)/(sum(MH_variance_denom1_i)*sum(MH_variance_denom2_i)))) as MH_IRR_lower,
		exp(log(sum(MH_numerator_term_i)/sum(MH_w_i))+1.96*sqrt(sum(MH_variance_numerator_i)/(sum(MH_variance_denom1_i)*sum(MH_variance_denom2_i)))) as MH_IRR_upper
from stratified_open;
quit;

/*Calculate IV IRD & 95% CI*/
proc sql;
title 'Inverse Variance IRD';
select sum(IV_numerator_term_i)/sum(IV_w_i) as IV_IRD,
	sum(IV_numerator_term_i)/sum(IV_w_i)-1.96*sqrt(1/sum(IV_w_i)) as IV_IRD_lower,
	sum(IV_numerator_term_i)/sum(IV_w_i)+1.96*sqrt(1/sum(IV_w_i)) as IV_IRD_upper
from stratified_open;
quit;
/*Test of Association*/
proc sql;
title 'Test of No Exposure-Disease Association';
select ((sum(a_i)-sum(expected_value_i))**2)/sum(variance_term_i) as Z_squared,
	1-probchi(((sum(a_i)-sum(expected_value_i))**2)/sum(variance_term_i),1) as p,
	1	as DF
from stratified_open;
quit;
/*Test of Heterogeneity: IRR*/
/*Store MH IRR for H statistic calculation*/
proc sql noprint;
select sum(MH_numerator_term_i)/sum(MH_w_i) into :MH_IRR 
from stratified_open;
quit;
proc sql;
title 'Tests of Homogeneity: IRR';
select sum(((log(IRR_i)-log(&MH_IRR))**2)/Variance_IRR_confint_i) as H_IRR,
	1-probchi(sum(((log(IRR_i)-log(&MH_IRR))**2)/Variance_IRR_confint_i), (count(*)-1)) as p,
	count(*)-1 as DF
from stratified_open;
quit;
/*Test of Heterogeneity: IRD*/
/*Store IV IRD for H statistic calculation*/
proc sql noprint;
select sum(IV_numerator_term_i)/sum(IV_w_i) into :IV_IRD 
from stratified_open;
quit;
proc sql;
title 'Tests of Homogeneity: IRD';
select sum(((IRD_i-&IV_IRD)**2)/Variance_IRD_confint_i) as H_IRD,
	1-probchi(sum(((IRD_i-&IV_IRD)**2)/Variance_IRD_confint_i), (count(*)-1)) as p,
	count(*)-1 as DF
from stratified_open;
quit;




