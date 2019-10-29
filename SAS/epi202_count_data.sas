/*******************************************************************
EPI 202 - Crude & Stratified Analyses for Count Data 
Created by: Sam Molsberry
Date: 9/1/2018
Dataset: Evans County example data

Reminders before beginning:
 1) If you would like to recreate the example using the Evans County data done in the documentation,
    you do not need to change any of this code aside from the file path where you have saved the csv fileT
 2) If you would like to change this code to do analyses on a different dataset, please adhere to the following:
	a) The exposure and outcome variables must be binary. To change these variables, please change the 
	   variables names and values indicating exposed/case and unexposed/noncase as necessary 
	b) The stratification variable(s) must be categorical. To change these variables, please change the 
	   variables names as necessary
	c) You do not need to change the names of any of the outputted objects, but if you do so  make sure to
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
/*Crude Analysis: Closed cohort data*/
/************************************/

data dat;
set dat;
if HTN^=.; *Change HTN to exposure variable and . to value indicating missing value as necessary;

/************************************************************************************************/
/*Next four lines create variables for each distinct combination of exposure and outcome. To use 
a different data sets change 'HTN' to name of you exposure variable and 'CHD' to name of your 
outcome variable. Additionally change to values of 1 and 0 as necessary to be consistent with the 
values indicating, exposed v. unexposed and case v. noncase in your data. Do not change the value
of 1 that exposed_case, exposed_noncase, unexposed_case, and unexposed_noncase are set to.

After changing these variables and values as necessary, you do not need to make any further changes
to the crude analysis section of the code in order to complete the crude analyses*/
/************************************************************************************************/
if HTN=1 and CHD=1 then exposed_case=1; 
if HTN=1 and CHD=0 then exposed_noncase=1;
if HTN=0 and CHD=1 then unexposed_case=1;
if HTN=0 and CHD=0 then unexposed_noncase=1;
run;


proc means data=dat sum maxdec=2;
ods output Summary=crude_closed;
var exposed_case exposed_noncase unexposed_case unexposed_noncase;
run;

data crude_closed;
set crude_closed;
rename exposed_case_Sum=a
	exposed_noncase_Sum=c
	unexposed_case_Sum=b
	unexposed_noncase_Sum=d;
drop VNAME_exposed_case VNAME_unexposed_case VNAME_exposed_noncase VNAME_unexposed_noncase;
run; 

data crude_closed;
set crude_closed;
label a="a"
	b="b"
	c="c"
	d="d";
run;

data crude_closed;
set crude_closed;
/*Specify count related variables*/
N1=a+c;
N0=b+d;
M1=a+b;
M0=c+d;
T=N1+N0;

/*Calculate CIR & 95% CI*/
CIR=(a/N1)/(b/N0);
Variance_CIR_confint=(c/(a*N1))+(d/(b*N0));
CIR_confint_lower=round(exp(log(CIR)-1.96*sqrt(Variance_CIR_confint)),0.001);
CIR_confint_upper=round(exp(log(CIR)+1.96*sqrt(Variance_CIR_confint)),0.001);
CIR_confint=catx(',', CIR_confint_lower, CIR_confint_upper);

/*Calculate CID & 95% CI*/
CID=(a/N1)-(b/N0);
Variance_CID_confint=((a*c)/(N1**3))+((b*d)/(N0**3));
CID_confint_lower=round(CID-1.96*sqrt(Variance_CID_confint), 0.001);
CID_confint_upper=round(CID+1.96*sqrt(Variance_CID_confint), 0.001);
CID_confint=catx(',', CID_confint_lower, CID_confint_upper);

/*Test of Association*/
Expected_null=(N1*M1)/T;
Variance_testassoc=(M1*M0*N1*N0)/(T**3);
Z_squared=((a-Expected_null)**2)/Variance_testassoc;
p_testassoc=1-probchi(Z_squared, 1);
run;

proc print data=crude_closed;
var CIR CIR_confint CID CID_confint;
run;

proc print data=crude_closed;
var Z_squared p_testassoc;
run;

/*****************************************/
/*Stratified Analysis: Closed cohort data*/
/*****************************************/

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
by SMK; *Change SMK to the list of your stratification factor(s);
run;

proc means data=dat sum maxdec=2;
ods output Summary=stratified_closed;
var exposed_case exposed_noncase unexposed_case unexposed_noncase;
by SMK; *Change SMK to the list of your stratification factor(s);
run;

data stratified_closed;
set stratified_closed;
rename exposed_case_Sum=a_i
	exposed_noncase_Sum=c_i
	unexposed_case_Sum=b_i
	unexposed_noncase_Sum=d_i;
drop VNAME_exposed_case VNAME_unexposed_case VNAME_exposed_noncase VNAME_unexposed_noncase;
run; 

data crude_open;
set crude_open;
label a_i="a_i"
	b_i="b_i"
	c_i="c_i"
	d_i="d_i";
run;

data stratified_closed;
set stratified_closed;
/*Specify count related variables*/
N1_i=a_i+c_i;
N0_i=b_i+d_i;
M1_i=a_i+b_i;
M0_i=c_i+d_i;
T_i=N1_i+N0_i;

/*Stratum-specific: Calculate CIR & 95% CI*/
CIR_i=(a_i/N1_i)/(b_i/N0_i);
Variance_CIR_confint_i=(c_i/(a_i*N1_i))+(d_i/(b_i*N0_i));
CIR_confint_lower_i=round(exp(log(CIR_i)-1.96*sqrt(Variance_CIR_confint_i)),0.001);
CIR_confint_upper_i=round(exp(log(CIR_i)+1.96*sqrt(Variance_CIR_confint_i)),0.001);
CIR_confint_i=catx(',', CIR_confint_lower_i, CIR_confint_upper_i);

/*Summary: Calculate MH CIR & 95% CI*/
MH_w_i=(b_i*N1_i)/T_i;
MH_numerator_term_i=(a_i*N0_i)/T_i;
MH_variance_numerator_i=((M1_i*N1_i*N0_i)-(a_i*b_i*T_i))/(T_i**2);
MH_variance_denom1_i=(a_i*N0_i)/T_i;
MH_variance_denom2_i=(b_i*N1_i)/T_i;

/*Stratum-specific: Calculate CID & 95% CI*/
CID_i=(a_i/N1_i)-(b_i/N0_i);
Variance_CID_confint_i=((a_i*c_i)/(N1_i**3))+((b_i*d_i)/(N0_i**3));
CID_confint_lower_i=round(CID_i-1.96*sqrt(Variance_CID_confint_i), 0.001);
CID_confint_upper_i=round(CID_i+1.96*sqrt(Variance_CID_confint_i), 0.001);
CID_confint_i=catx(',', CID_confint_lower_i, CID_confint_upper_i);

/*Summary: Calculate Inverse Variance CID & 95% CI*/
IV_w_i=((N1_i**3)*(N0_i**3))/(((N0_i**3)*a_i*c_i)+((N1_i**3)*b_i*d_i));
IV_numerator_term_i=IV_w_i*CID_i;
drop CIR_confint_lower_i CIR_confint_upper_i CID_confint_lower_i CID_confint_upper_i;

/*Test of No association*/
expected_value_i=(N1_i*M1_i)/T_i;
variance_term_i=(M1_i*M0_i*N1_i*N0_i)/(T_i**3);
run;
/*Calculate MH CIR & 95% CI*/
proc sql;
title 'MH CIR';
select sum(MH_numerator_term_i)/sum(MH_w_i) as MH_CIR,
		exp(log(sum(MH_numerator_term_i)/sum(MH_w_i))-1.96*sqrt(sum(MH_variance_numerator_i)/(sum(MH_variance_denom1_i)*sum(MH_variance_denom2_i)))) as MH_CIR_lower,
		exp(log(sum(MH_numerator_term_i)/sum(MH_w_i))+1.96*sqrt(sum(MH_variance_numerator_i)/(sum(MH_variance_denom1_i)*sum(MH_variance_denom2_i)))) as MH_CIR_upper
from stratified_closed;
quit;

/*Calculate IV CID & 95% CI*/
proc sql;
title 'Inverse Variance CID';
select sum(IV_numerator_term_i)/sum(IV_w_i) as IV_CID,
	sum(IV_numerator_term_i)/sum(IV_w_i)-1.96*sqrt(1/sum(IV_w_i)) as IV_CID_lower,
	sum(IV_numerator_term_i)/sum(IV_w_i)+1.96*sqrt(1/sum(IV_w_i)) as IV_CID_lower
from stratified_closed;
quit;
/*Test of Association*/
proc sql;
title 'Test of No Exposure-Disease Association';
select ((sum(a_i)-sum(expected_value_i))**2)/sum(variance_term_i) as Z_squared,
	1-probchi(((sum(a_i)-sum(expected_value_i))**2)/sum(variance_term_i),1) as p,
	1	as DF
from stratified_closed;
quit;
/*Test of Heterogeneity: CIR*/
/*Store MH CIR for H statistic calculation*/
proc sql noprint;
select sum(MH_numerator_term_i)/sum(MH_w_i) into :MH_CIR 
from stratified_closed;
quit;
proc sql;
title 'Tests of Homogeneity: CIR';
select sum(((log(CIR_i)-log(&MH_CIR))**2)/Variance_CIR_confint_i) as H_CIR,
	1-probchi(sum(((log(CIR_i)-log(&MH_CIR))**2)/Variance_CIR_confint_i), (count(*)-1)) as p,
	count(*)-1 as DF
from stratified_closed;
quit;
/*Test of Heterogeneity: CID*/
/*Store IV CID for H statistic calculation*/
proc sql noprint;
select sum(IV_numerator_term_i)/sum(IV_w_i) into :IV_CID 
from stratified_closed;
quit;
proc sql;
title 'Tests of Homogeneity: CID';
select sum(((CID_i-&IV_CID)**2)/Variance_CID_confint_i) as H_CID,
	1-probchi(sum(((CID_i-&IV_CID)**2)/Variance_CID_confint_i), (count(*)-1)) as p,
	count(*)-1 as DF
from stratified_closed;
quit;





