
/*Dupilcate Original dataset & time varaible*/
%macro TIMEAUCBOOT(DT=, Timepoint = , adj_var =  , 
					cls_var =  , Time2event = , 
					event = , Nboot = 10, method = km, seed = 123);
%let time2event_2 = &time2event._2;
Data orig;
set &DT;
&time2event_2 = &time2event;
run;
/*Run cox model on original sample*/
ods graphics on;
ods exclude all;
/*AUC(t) in original sample*/
proc phreg data=&DT plots=roc 
rocoptions(at= &Timepoint outroc = Res_original(rename = (_Sensitivity_ = Sens)) method = &method ) ;
class &cls_var/ param=reference;
model &time2event*&event(0)= &adj_var/ rl ;
run;
ods graphics off;
/*options nonotes nosource nosource2 errors=0;*/
/*ods exclude all;*/
proc sql;
select count(*) into: n_orig from Res_original; /*How many marker if use original data*/
/*Calculate original AUC*/
proc sort data=Res_original
(keep= _Specificity_ Sens);
by descending _Specificity_ sens;
run;
data area1;
set Res_original end=last;
x = 1 - _Specificity_;
xprev=lag(x);
yprev=lag(sens);
output;
if last then do;
xprev=x;
yprev=sens;
x=1;
sens=1;
output;
end;
run; 
data _null_;
retain area 0;
set area1(firstobs=2) end=last;
area=area+(sens+yprev)*(x-xprev)/2;
if last then call symput('ROC_orig',put(area, best20.));
run; 
/*Bootstrap original sample*/
proc surveyselect data=&DT NOPRINT seed=&seed outhits 
     out=Boot
     method=urs              /* resample with replacement */
     samprate=1              /* each bootstrap sample has N observations */
     reps=&Nboot; 	/* generate NumSamples bootstrap resamples */
run;
/*Cox model on Boot&n data*/
proc phreg data=Boot noprint;
by Replicate;
class &cls_var/ param=reference;
model &time2event*&event(0) = &adj_var/rl ridging=absolute;
/*Original data below*/
baseline out=orig_z xbeta=betaz covariates=orig timelist=&Timepoint;
run;
Data output;
N = .;
ROC_orig = .;
Roc_Boot = .;
Bias = .;
Model_obs = .;
run;
%Macro Boot_roc;
%do b = 1 %to &Nboot;
data orig_z1;
set orig_z(where =(Replicate = &b));
run;
proc sql;
select count(*) into: n from orig_z1;  /*First count how many betaz:n*/
select max(betaz) format = best20. into:m from orig_z1;
select betaz format = best20. into: c1 - : %sysfunc(compress(c&n)) from orig_z1
order by betaz;	/*Create n markers*/
quit;
/*Shat(t)*/
proc lifetest data=orig_z1 method=KM outsurv=shat noprint;
time &time2event_2*&event(0);
run;
proc sort data=shat;
by &time2event_2;
run;
data _null_;
set shat end=eof;
where &time2event_2 <= &Timepoint and not missing(survival);
if eof then call symput('shat',put(survival, best20.));
run;
/*The following Macro calculate SP and SN at each Markers for Boot'b'*/
proc datasets lib=work;
delete res: ;
run;
data res;
c = .;
boot_sens = 1;
boot_spec = 0;
run;
%macro Sen_Spe;
/*We calculate Sen and Spe at each M value*/
%let i = 1;
%do %until(%sysevalf(&&c&i. >= &m,boolean) = 1);
proc lifetest data=orig_z1 method=KM outsurv=shat_c noprint;
where betaz > &&c&i..; 
time &time2event_2*&event(0);
run;
proc sort data=shat_c;
by &time2event_2;
run;
data _null_;
set shat_c end=eof;
where &time2event_2 <= &Timepoint and not missing(survival);
if eof then call symput('shat_c',put(survival, best20.));
run;
proc sql;
select count(*) into: n_c from orig_z1 where betaz <= &&c&i.;
quit;
data Out;
c = &&c&i.;
boot_sens = (1 - &shat_c)*(1 - &n_c/&n)/(1 - &shat);
boot_spec = 1 - (&shat_c*(1 - &n_c/&n))/(&shat);
run;
proc append data=out base=res;
run;
%let i = %eval(&i+1);
%end;
/*%else %do;*/
data Out;
c = &&c&i.;
boot_sens = 0;
boot_spec = 1 ;
run;
proc append data=out base=res;
run;
/*%end;*/
/*%end;*/
%mend;
%Sen_Spe;
/*Step 5 Calculate AUC*/
proc sort data=res(rename=(Boot_sens = sens));
by descending boot_spec Sens;
run;
data area1;
set res end=last;
x = 1 - boot_spec;
xprev=lag(x);
yprev=lag(sens);
output;
if last then do;
xprev=x;
yprev=sens;
x=1;
sens=1;
 output;
end;
run; 
data _null_;
retain area 0;
set area1(firstobs=2) end=last;
area=area+(sens+yprev)*(x-xprev)/2;
if last then call symput('ROC_Boot',put(area, best20.));
run; 
/*Optimism with ROC Original*/
data out;
N = &b;
ROC_orig = &ROC_orig;
ROC_boot = &ROC_Boot;
Bias = ROC_boot - Roc_ORIG;
Model_obs = &n;
run;
proc append base=output data=out;
run;
%end;
%mend;
%Boot_roc;
proc sql;
select sum(Bias)/&Nboot format best20. into: mean from output;
quit;
/*proc datasets lib=work;*/
/*delete orig orig_z area1 boot shat shat_c;*/
/*run;*/
Data output;
set output;
if missing(N) then delete;
run;
Data conv;
set output;
if Model_obs < &n_orig then flag = 0; 
N1 = &n_orig;
if flag = 1 then output;
label N = "Replicate of Bootstrap" 
N1 = "Obs of Orginal Model" flag = "Model Convergence" 
Model_obs = "Obs of Bootstrap Model" Roc_Boot = "Bootstrap C-statistics";
drop ROC_orig Bias;
run;
proc sql;
select count(*) into: Num from conv;
quit;
Data out;
Method = "&method";
Time = &Timepoint;
AUC = &ROC_orig;
bias = &mean;
ADJ_AUC = AUC - abs(bias);
label AUC = "AUC original sample" bias = "Over-optimism" ADJ_AUC = "Bias-corrected AUC" Time = "Time point";
run;
ods exclude none;
ods graphics on;
options notes source source2 errors=20;
%if &num ^= 0 %then %do;
	proc print data=Conv label noobs;
	title "There are %sysfunc(compress(&num.)) of boostrap sample does not converge";
	var flag N Model_obs N1 Roc_Boot;
	run;
%end;
proc sgplot data=output noautolegend;
scatter x=n y=Bias;
series x=n y=bias;
refline 0;
title "Internal Estimate of Optimism";
run;
proc sgplot data=output ;
histogram Roc_Boot;
title "Internal Validation C-Statistics";
run;
proc print data=out label noobs;
title "Internal Validation";
run;
%mend;
