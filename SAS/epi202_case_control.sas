/*******************************************************************
EPI 202 - Crude & Stratified Analyses for Case-control Data 
Created by: Sam Molsberry
Date: 9/1/2018
Dataset: Evans County example data

Reminders before beginning:
 1) If you would like to recreate the example using the Evans County data done in the documentation,
    you do not need to change any of this code aside from the file path where you have saved the csv fileT
 2) If you would like to change this code to do analyses on a different dataset, please adhere to the following:
	a) The exposure and case-control status variables must be binary. To change these variables, please change the 
	   variables names and values indicating exposed/case and unexposed/control as necessary 
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

/*************************************/
/**Crude Analysis: Case-control data**/
/*************************************/

data dat;
set dat;
if caco^='NA';*Restricts to observations selected for case-control study;
if HTN^=.;*Restricts to observations with non-missing exposure data;

/************************************************************************************************/
/*Next four lines create variables for each distinct combination of exposure and case status. To use 
a different data sets change 'HTN' to name of you exposure variable and 'caco' to name of your 
outcome variable. Additionally change to values of 1/case and 0/control as necessary to be consistent with the 
values indicating, exposed v. unexposed and case v. noncase in your data. Do not change the value
of 1 that exposed_case, exposed_control, unexposed_case, and unexposed_control are set to.

After changing these variables and values as necessary, you do not need to make any further changes
to the crude analysis section of the code in order to complete the crude analyses*/
/************************************************************************************************/

if HTN=1 and caco='case' then exposed_case=1;
if HTN=0 and caco='case' then unexposed_case=1;
if HTN=1 and caco='control' then exposed_control=1;
if HTN=0 and caco='control' then unexposed_control=1;
run;

proc means data=dat sum maxdec=2;
ods output Summary=crude_cc;
var exposed_case unexposed_case exposed_control unexposed_control;
run;

data crude_cc;
set crude_cc;
rename exposed_case_Sum=a
	unexposed_case_Sum=b
	exposed_control_Sum=c
	unexposed_control_Sum=d;
drop VNAME_exposed_case VNAME_unexposed_case VNAME_exposed_control VNAME_unexposed_control;
run; 

data crude_cc;
set crude_cc;
label a="a"
	b="b"
	c="c"
	d="d";
run;

data crude_cc;
set crude_cc;
/*Specify count related variables*/
M1=a+b;
M0=c+d;
N1=a+c;
N0=b+d;
T=N1+N0;

/*Calculate OR & 95% CI*/
OR=(a*d)/(b*c);
Variance_OR_confint=(1/a)+(1/b)+(1/c)+(1/d);
OR_confint_lower=round(exp(log(OR)-1.96*sqrt(Variance_OR_confint)),0.001);
OR_confint_upper=round(exp(log(OR)+1.96*sqrt(Variance_OR_confint)),0.001);
OR_confint=catx(',', OR_confint_lower, OR_confint_upper);

/*Test of Association*/
Expected_null=(N1*M1)/T;
Variance_testassoc=(M1*M0*N1*N0)/((T**2)*(T-1));
Z_squared=((a-Expected_null)**2)/Variance_testassoc;
p_testassoc=1-probchi(Z_squared, 1);
run;

proc print data=crude_cc;
var OR OR_confint;
run;

proc print data=crude_cc;
var Z_squared p_testassoc;
run;

/****************************************/
/*Stratified Analysis: Case-control data*/
/****************************************/

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
ods output Summary=stratified_cc;
var exposed_case unexposed_case exposed_control unexposed_control;
by SMK; *Change SMK to the list of your stratification factor(s);
run;

data stratified_cc;
set stratified_cc;
rename exposed_case_Sum=a_i
	unexposed_case_Sum=b_i
	exposed_control_Sum=c_i
	unexposed_control_Sum=d_i;
drop VNAME_exposed_case VNAME_unexposed_case VNAME_exposed_control VNAME_unexposed_control;
run; 

data stratified_cc;
set stratified_cc;
label a_i="a_i"
	b_i="b_i"
	c_i="c_i"
	d_i="d_i";
run;

data stratified_cc;
set stratified_cc;
/*Specify count related variables*/
M1_i=a_i+b_i;
M0_i=c_i+d_i;
N1_i=a_i+c_i;
N0_i=b_i+d_i;
T_i=N1_i+N0_i;

/*Calculate OR & 95% CI*/
OR_i=(a_i*d_i)/(b_i*c_i);
Variance_OR_confint_i=(1/a_i)+(1/b_i)+(1/c_i)+(1/d_i);
OR_confint_lower_i=round(exp(log(OR_i)-1.96*sqrt(Variance_OR_confint_i)),0.001);
OR_confint_upper_i=round(exp(log(OR_i)+1.96*sqrt(Variance_OR_confint_i)),0.001);
OR_confint_i=catx(',', OR_confint_lower_i, OR_confint_upper_i);

/*Summary: Calculate MH IRR & 95% CI*/
MH_w_i=(b_i*c_i)/T_i;
MH_numerator_term_i=(a_i*d_i)/T_i;
/*Set up RGB Variance calculation*/
RGB_A_i=(a_i*d_i)/T_i;
RGB_B_i=(b_i*c_i)/T_i;
RGB_C_i=(a_i+d_i)/T_i;
RGB_D_i=(b_i+c_i)/T_i;

drop OR_confint_lower_i OR_confint_upper_i;

/*Test of No association*/
expected_value_i=(N1_i*M1_i)/T_i;
variance_term_i=(M1_i*M0_i*N1_i*N0_i)/((T_i**2)*(T_i-1));
run;
/*Calculate MH OR & 95% CI*/
/*Calculate RGB Variance*/
proc sql noprint;
select sum(RGB_A_i*RGB_C_i)/(2*sum(RGB_A_i)*sum(RGB_A_i))+(sum((RGB_A_i*RGB_D_i)+(RGB_B_i*RGB_C_i))/(2*(sum(RGB_A_i)*sum(RGB_B_i))))+(sum(RGB_B_i*RGB_D_i)/(2*(sum(RGB_B_i)**2))) into :RGB 
from stratified_cc;
quit;

proc sql;
title 'MH OR';
select sum(MH_numerator_term_i)/sum(MH_w_i) as MH_OR,
		exp(log(sum(MH_numerator_term_i)/sum(MH_w_i))-1.96*sqrt(&RGB)) as MH_OR_lower,
		exp(log(sum(MH_numerator_term_i)/sum(MH_w_i))+1.96*sqrt(&RGB)) as MH_OR_upper
from stratified_cc;
quit;

/*Test of Association*/
proc sql;
title 'Test of No Exposure-Disease Association';
select ((sum(a_i)-sum(expected_value_i))**2)/sum(variance_term_i) as Z_squared,
	1-probchi(((sum(a_i)-sum(expected_value_i))**2)/sum(variance_term_i),1) as p,
	1	as DF
from stratified_cc;
quit;
/*Test of Heterogeneity: IRR*/
/*Store MH IRR for H statistic calculation*/
proc sql noprint;
select sum(MH_numerator_term_i)/sum(MH_w_i) into :MH_OR 
from stratified_cc;
quit;
proc sql;
title 'Tests of Homogeneity: OR';
select sum(((log(OR_i)-log(&MH_OR))**2)/Variance_OR_confint_i) as H_OR,
	1-probchi(sum(((log(OR_i)-log(&MH_OR))**2)/Variance_OR_confint_i), (count(*)-1)) as p,
	count(*)-1 as DF
from stratified_cc;
quit;





