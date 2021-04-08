# TIMEAUCBOOT
Title: A SAS Macro to Perform Bootstrapping for Internal Validation of Time-Dependent Area Under the Curve

# Instruction for using TIMEAUCBOOT macro 
The %TIMEAUCBOOT macro parameters are shown below:
DATA=	(Required) SAS data set name.
TIMEPOINT=	(Required) Specific time point for AUC. 
EVENT= 	(Required) The outcome of interest. This has to be coded as 0 or 1, where 0 means censored and 1 means having the outcome of interest. 
TIME2EVENT=	(Required) Time to the event (or censoring).
ADJ_VAR=	All covariates adjusted for in the model. 
CLS_VAR=	(Required if categorical variable exist in ADJ_VAR) All categorical covariates from ADJ_VAR. The reference can be specified in this statement. 
NBOOT =	The number of bootstrap samples. The default is 10.  
SEED=	Random control seed. The default is 123. 
METHOD=KM	Methods of estimating time-dependent AUC. Currently only the Kaplan-Meier method is available. 

# Macro Code example
The following code used TRAC.csv file.
%TIMEAUCBOOT(
	DT		=	TRAC, 	
	TIMEPOINT 	= 	365, 
	ADJ_VAR 	=  	AGE_GRP REC_SF36_AVG REC_ISMALE, 
	CLS_VAR 	= 	AGE_GRP REC_ISMALE(REF = “1”), 
	TIME2EVENT 	= 	FOLLOW_TIME, 
	EVENT 	= 	EVENT, 
	NBOOT 	= 	100, 
	METHOD 	= 	KM,
	SEED 		= 	123
);

# Contact
Your comments and questions are valued and encouraged. Contact the author at:
Jialin Han 
Stanford University
Colin-Jialin.han@stanford.edu or Jialin.han@uconn.edu

