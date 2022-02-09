
/********************************************************************************************************************/
/******************SAS SYNTAX TO EXTRACT AND COMBINE NHIS 1997-2014 and Mortality 2015 DATA *************************/
/********************************************************************************************************************/

/*Location of individual data files*/
libname lib1997 "C:\...\NHIS Raw Data\1997";
libname lib1998 "C:\...\NHIS Raw Data\1998";
libname lib1999 "C:\...\NHIS Raw Data\1999";
libname lib2000 "C:\...\NHIS Raw Data\2000";
libname lib2001 "C:\...\NHIS Raw Data\2001";
libname lib2002 "C:\...\NHIS Raw Data\2002";
libname lib2003 "C:\...\NHIS Raw Data\2003";
libname lib2004 "C:\...\NHIS Raw Data\2004";
libname lib2005 "C:\...\NHIS Raw Data\2005";
libname lib2006 "C:\...\NHIS Raw Data\2006";
libname lib2007 "C:\...\NHIS Raw Data\2007";
libname lib2008 "C:\...\NHIS Raw Data\2008";
libname lib2009 "C:\...\NHIS Raw Data\2009";
libname lib2010 "C:\...\NHIS Raw Data\2010";
libname lib2011 "C:\...\NHIS Raw Data\2011";
libname lib2012 "C:\...\NHIS Raw Data\2012";
libname lib2013 "C:\...\NHIS Raw Data\2013";
libname lib2014 "C:\...\NHIS Raw Data\2014";
libname mort    "C:\...\Mortality Data\1997-2014";

/*Location of combined, cleaned data file */
libname combined "C:\..\NHIS clean\";


OPTIONS nofmterr;    /*needed since we're not loading the format for each NHIS data file*/
options nonotes;	 /*suppress notes in log - too many notes about missing formats*/





/*********************************************************************************************************************/
/****************************PART 1: EXTRACT ALL OF THE DATA FROM PUBLIC USE DATA FILES*******************************/
/********************************************************************************************************************/

/************ NHIS 1997 ***************/

/*Sample Adult File*/
data samadult1997; set lib1997.samadult;
	/* Create unique identifier*/
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX)); 
	
	/*select variables to keep*/
	keep PUBLICID SRVY_YR INTV_QRT HHX FMX PX WTFA_SA STRATUM PSU ALC1YR ALCLIFE ALC12MYR ALCAMT ALC5UPYR	
	SMKSTAT2 BMI HYPEV HYPDIFV DIBEV VIGFREQW VIGMIN MODFREQW MODMIN SAD NERVOUS RESTLESS HOPELESS 
	EFFORT WORTHLS;
run;

/*Person File*/
data personsx1997; set lib1997.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));	
	keep PUBLICID SEX AGE_P R_MARITL DOINGLW WHYNOWRK ORIGIN RACEREC MRACE_P EDUC RAT_CAT PHSTAT USBORN_P; 	
run;


		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis1997 as
				select L.*, R.*
				from samadult1997 as L
				LEFT JOIN personsx1997 as R
				on L.PUBLICID = R.PUBLICID;
			quit;



/************ NHIS 1998 ***************/

data samadult1998; set lib1998.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));					
	keep PUBLICID SRVY_YR INTV_QRT HHX FMX PX WTFA_SA STRATUM PSU ALC1YR ALCLIFE ALC12MYR ALCAMT ALC5UPYR	
		SMKSTAT2 BMI HYPEV HYPDIFV DIBEV VIGFREQW VIGMIN MODFREQW MODMIN 
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
	run;

data personsx1998; set lib1998.personsx;
	length PUBLICID $14;PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX)); 						
	keep PUBLICID SEX AGE_P R_MARITL DOINGLW WHYNOWRK HISPCODE EDUC MRACE_P RAT_CAT PHSTAT USBORN_P;	 
run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
			create table nhis1998 as
			select L.*, R.*
			from samadult1998 as L
			LEFT JOIN personsx1998 as R
			on L.PUBLICID = R.PUBLICID;
		quit;





/************ NHIS 1999 ***************/

data samadult1999; set lib1999.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));  				 
	keep PUBLICID SRVY_YR INTV_QRT HHX FMX PX WTFA_SA STRATUM PSU ALC1YR ALCLIFE ALC12MYR ALCAMT ALC5UPYR  	 	
		SMKSTAT2 BMI HYPEV HYPDIFV DIBEV VIGFREQW VIGMIN MODFREQW MODMIN 
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS MHDSADWY MHDINTWY MHDSAD2W; 
	run;
	

data personsx1999; set lib1999.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));  							
	keep PUBLICID SEX AGE_P R_MARITL DOINGLW WHYNOWRK HISPCODR MRACBR_P EDUC RAT_CAT PHSTAT USBORN_P;  
		run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis1999 as
				select L.*, R.*
				from samadult1999 as L
				LEFT JOIN personsx1999 as R
				on L.PUBLICID = R.PUBLICID;
			quit;



/************ NHIS 2000 ***************/

data samadult2000; set lib2000.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));  				 
	keep PUBLICID SRVY_YR INTV_QRT HHX FMX PX WTFA_SA STRATUM PSU ALC1YR ALCLIFE ALC12MYR ALCAMT ALC5UPYR  	 	
		SMKSTAT2 BMI HYPEV HYPDIFV DIBEV VIGFREQW VIGMIN MODFREQW MODMIN 
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
run;

data personsx2000; set lib2000.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));  							
	keep PUBLICID SEX AGE_P R_MARITL DOINGLW WHYNOWRK HISCOD_I MRACBP_I EDUC RAT_CAT PHSTAT USBRTH_P;  
run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2000 as
				select L.*, R.*
				from samadult2000 as L
				LEFT JOIN personsx2000 as R
				on L.PUBLICID = R.PUBLICID;
		quit;



/************ NHIS 2001 ***************/

data samadult2001; set lib2001.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));  				 
	keep PUBLICID SRVY_YR INTV_QRT HHX FMX PX WTFA_SA STRATUM PSU ALC1YR ALCLIFE ALC12MYR ALCAMT ALC5UPYR  	 	
		SMKSTAT2 BMI HYPEV HYPDIFV DIBEV VIGFREQW VIGMIN MODFREQW MODMIN 
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
	run;

data personsx2001; set lib2001.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));  							
	keep PUBLICID SEX AGE_P R_MARITL DOINGLW1 WHYNOWK1 HISCOD_I MRACBP_I EDUC RAT_CAT PHSTAT USBRTH_P;  
run;


		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2001 as
				select L.*, R.*
				from samadult2001 as L
				LEFT JOIN personsx2001 as R
				on L.PUBLICID = R.PUBLICID;
			quit;


/************ NHIS 2002 ***************/

data samadult2002; set lib2002.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));  				 
	keep PUBLICID SRVY_YR INTV_QRT HHX FMX PX WTFA_SA STRATUM PSU ALC1YR ALCLIFE ALC12MYR ALCAMT ALC5UPYR  	 	 
		SMKSTAT2 BMI HYPEV HYPDIFV DIBEV VIGFREQW VIGMIN MODFREQW MODMIN 
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
run;

data personsx2002; set lib2002.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));  							
	keep PUBLICID SEX AGE_P R_MARITL DOINGLW1 WHYNOWK1 HISCOD_I MRACBP_I EDUC RAT_CAT PHSTAT GEOBRTH;  
run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2002 as
				select L.*, R.*
				from samadult2002 as L
				LEFT JOIN personsx2002 as R
				on L.PUBLICID = R.PUBLICID;
			quit;


/************ NHIS 2003 ***************/

data samadult2003; set lib2003.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));  				 
	keep PUBLICID SRVY_YR INTV_QRT HHX FMX PX WTFA_SA STRATUM PSU ALC1YR ALCLIFE ALC12MYR ALCAMT ALC5UPYR   
		 SMKSTAT2 BMI HYPEV HYPDIFV DIBEV VIGFREQW VIGMIN MODFREQW MODMIN 
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
run;
		
data personsx2003; set lib2003.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||PX));  							
	keep PUBLICID SEX AGE_P R_MARITL DOINGLW1 WHYNOWK1 HISCODI2 MRACBPI2 EDUC RAT_CAT PHSTAT GEOBRTH;  
run;


		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2003 as
				select L.*, R.*
				from samadult2003 as L
				LEFT JOIN personsx2003 as R
				on L.PUBLICID = R.PUBLICID;
			quit;




/************ NHIS 2004 ***************/

/*Sample Adult File*/
data samadult2004; set lib2004.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));   
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR HHX FMX FPX WTFA_SA ALC1YR ALCLIFE ALC12MYR ALCAMT ALC5UPYR 
		SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
run;

/*Person File*/
data personsx2004; set lib2004.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));  						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI2 MRACBPI2 EDUC1 STRATUM PSU PHSTAT GEOBRTH;   /*in 2004 the PSU and STRATUM variables were only available on the person file*/ 
run;

/*Household File*/
data househld2004; set lib2004.househld;   /*in 2004 the Quater that the survey was completed in was only available on the household file*/ 
	keep HHX INTV_QRT;
run;

/*Family File*/
data familyxx2004; set lib2004.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT;
run;


		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2004_pre1 as
				select L.*, R.*
				from samadult2004 as L
				LEFT JOIN personsx2004 as R
				on L.PUBLICID = R.PUBLICID;
			quit;

		Proc sql;
				create table nhis2004_pre2 as
				select L.*, R.*
				from nhis2004_pre1 as L
				LEFT JOIN househld2004 as R
				on L.HHX = R.HHX;
			quit;

		Proc sql;
				create table nhis2004 as
				select L.*, R.*
				from nhis2004_pre2 as L
				LEFT JOIN familyxx2004 as R
				on L.FAM_ID = R.FAM_ID;
			quit;



/************ NHIS 2005 ***************/

data samadult2005; set lib2005.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));   
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR INTV_QRT HHX FMX FPX WTFA_SA STRATUM PSU ALC1YR ALCLIFE ALC12MYR ALCAMT 
		ALC5UPYR SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
run;

data personsx2005; set lib2005.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));  						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI2 MRACBPI2 EDUC1 PHSTAT GEOBRTH;  
run;

data familyxx2005; set lib2005.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT;
run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2005_pre1 as
				select L.*, R.*
				from samadult2005 as L
				LEFT JOIN personsx2005 as R
				on L.PUBLICID = R.PUBLICID;
			quit;

		Proc sql;
				create table nhis2005 as
				select L.*, R.*
				from nhis2005_pre1 as L
				LEFT JOIN familyxx2005 as R
				on L.FAM_ID = R.FAM_ID;
			quit;


/************ NHIS 2006 ***************/

data samadult2006; set lib2006.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX)); 				 
	length FAM_ID $14;     FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR INTV_QRT HHX FMX FPX WTFA_SA STRAT_P PSU_P ALC1YR ALCLIFE ALC12MYR ALCAMT 
		ALC5UPYR SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
run;

data personsx2006; set lib2006.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));  						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI3 MRACBPI2 EDUC1 PHSTAT GEOBRTH;  
run;

data familyxx2006; set lib2006.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT;
run;


		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2006_pre1 as
				select L.*, R.*
				from samadult2006 as L
				LEFT JOIN personsx2006 as R
				on L.PUBLICID = R.PUBLICID;
			quit;

		Proc sql;
				create table nhis2006 as
				select L.*, R.*
				from nhis2006_pre1 as L
				LEFT JOIN familyxx2006 as R
				on L.FAM_ID = R.FAM_ID;
			quit;


/************ NHIS 2007 ***************/

data samadult2007; set lib2007.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX)); 				 
	length FAM_ID $14;     FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR INTV_QRT HHX FMX FPX WTFA_SA STRAT_P PSU_P ALC1YR ALCLIFE ALC12MYR ALCAMT 
		ALC5UPYR SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN SAD 
		NERVOUS RESTLESS HOPELESS EFFORT WORTHLS DEPYR;
run;

data personsx2007; set lib2007.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));  						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI3 MRACBPI2 EDUC1 PHSTAT GEOBRTH;  
run;

data familyxx2007; set lib2007.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT3;
run;


		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2007_pre1 as
				select L.*, R.*
				from samadult2007 as L
				LEFT JOIN personsx2007 as R
				on L.PUBLICID = R.PUBLICID;
			quit;
			Proc sql;
				create table nhis2007 as
				select L.*, R.*
				from nhis2007_pre1 as L
				LEFT JOIN familyxx2007 as R
				on L.FAM_ID = R.FAM_ID;
			quit;


/************ NHIS 2008 ***************/

data samadult2008; set lib2008.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX)); 				 
	length FAM_ID $14;     FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR INTV_QRT HHX FMX FPX WTFA_SA STRAT_P PSU_P ALC1YR ALCLIFE ALC12MYR ALCAMT 
		ALC5UPYR SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS DEPRESS;
	run;

data personsx2008; set lib2008.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));  						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI3 MRACBPI2 EDUC1 PHSTAT GEOBRTH;  
run;

data familyxx2008; set lib2008.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT3;
run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2008_pre1 as
				select L.*, R.*
				from samadult2008 as L
				LEFT JOIN personsx2008 as R
				on L.PUBLICID = R.PUBLICID;
			quit;
		Proc sql;
				create table nhis2008 as
				select L.*, R.*
				from nhis2008_pre1 as L
				LEFT JOIN familyxx2008 as R
				on L.FAM_ID = R.FAM_ID;
			quit;


/************ NHIS 2009 ***************/

data samadult2009; set lib2009.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX)); 				 
	length FAM_ID $14;     FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR INTV_QRT HHX FMX FPX WTFA_SA STRAT_P PSU_P ALC1YR ALCLIFE ALC12MYR ALCAMT 
		ALC5UPYR SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
run;

data personsx2009; set lib2009.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));  						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI3 MRACBPI2 EDUC1 PHSTAT GEOBRTH;  
run;

data familyxx2009; set lib2009.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT3;
run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2009_pre1 as
				select L.*, R.*
				from samadult2009 as L
				LEFT JOIN personsx2009 as R
				on L.PUBLICID = R.PUBLICID;
			quit;
		Proc sql;
				create table nhis2009 as
				select L.*, R.*
				from nhis2009_pre1 as L
				LEFT JOIN familyxx2009 as R
				on L.FAM_ID = R.FAM_ID;
			quit;


/************ NHIS 2010 ***************/

data samadult2010; set lib2010.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX)); 				 
	length FAM_ID $14;     FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR INTV_QRT HHX FMX FPX WTFA_SA STRAT_P PSU_P ALC1YR ALCLIFE ALC12MYR ALCAMT 
		ALC5UPYR SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
run;

data personsx2010; set lib2010.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));  						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI3 MRACBPI2 EDUC1 PHSTAT GEOBRTH;  
run;

data familyxx2010; set lib2010.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT3;
run;


		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2010_pre1 as
				select L.*, R.*
				from samadult2010 as L
				LEFT JOIN personsx2010 as R
				on L.PUBLICID = R.PUBLICID;
			quit;
		Proc sql;
				create table nhis2010 as
				select L.*, R.*
				from nhis2010_pre1 as L
				LEFT JOIN familyxx2010 as R
				on L.FAM_ID = R.FAM_ID;
			quit;




/************ NHIS 2011 ***************/

data samadult2011; set lib2011.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX)); 				 
	length FAM_ID $14;     FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR INTV_QRT HHX FMX FPX WTFA_SA STRAT_P PSU_P ALC1YR ALCLIFE ALC12MYR ALCAMT 
		ALC5UPYR SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS;
run;

data personsx2011; set lib2011.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI3 MRACBPI2 EDUC1 PHSTAT GEOBRTH;  
run;

data familyxx2011; set lib2011.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT3;
run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2011_pre1 as
				select L.*, R.*
				from samadult2011 as L
				LEFT JOIN personsx2011 as R
				on L.PUBLICID = R.PUBLICID;
			quit ;

		Proc sql;
				create table nhis2011 as
				select L.*, R.*
				from nhis2011_pre1 as L
				LEFT JOIN familyxx2011 as R
				on L.FAM_ID = R.FAM_ID;
			quit;



/************ NHIS 2012 ***************/

data samadult2012; set lib2012.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX)); 				 
	length FAM_ID $14;     FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR INTV_QRT HHX FMX FPX WTFA_SA STRAT_P PSU_P ALC1YR ALCLIFE ALC12MYR ALCAMT 
		ALC5UPYR SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN
		SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS ADEPRSEV ADEPRSYR;
	run;

data personsx2012; set lib2012.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));  						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI3 MRACBPI2 EDUC1 PHSTAT GEOBRTH;  
run;

data familyxx2012; set lib2012.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT3;
run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2012_pre1 as
				select L.*, R.*
				from samadult2012 as L
				LEFT JOIN personsx2012 as R
				on L.PUBLICID = R.PUBLICID;
			quit;

		Proc sql;
				create table nhis2012 as
				select L.*, R.*
				from nhis2012_pre1 as L
				LEFT JOIN familyxx2012 as R
				on L.FAM_ID = R.FAM_ID;
			quit;




/************ NHIS 2013 ***************/

data samadult2013; set lib2013.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX)); 							 
	length FAM_ID $14;     FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR INTV_QRT HHX FMX FPX WTFA_SA STRAT_P PSU_P ALC1YR ALCLIFE ALC12MYR ALCAMT 
		ALC5UPYR SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN
		ASISAD ASINERV ASIRSTLS ASIHOPLS ASIEFFRT ASIWTHLS ; 
run;

data personsx2013; set lib2013.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX));  						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI3 MRACBPI2 EDUC1 PHSTAT GEOBRTH;		  
run;
	
data familyxx2013; set lib2013.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT3;
run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
			create table nhis2013_pre1 as
			select L.*, R.*
			from samadult2013 as L
			LEFT JOIN personsx2013 as R
			on L.PUBLICID = R.PUBLICID;
		quit;

		Proc sql;
				create table nhis2013 as
				select L.*, R.*
				from nhis2013_pre1 as L
				LEFT JOIN familyxx2013 as R
				on L.FAM_ID = R.FAM_ID;
			quit;



/************ NHIS 2014 ***************/

data samadult2014; set lib2014.samadult;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX)); 							 
	length FAM_ID $14;     FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep PUBLICID FAM_ID SRVY_YR INTV_QRT HHX FMX FPX WTFA_SA STRAT_P PSU_P ALC1YR ALCLIFE ALC12MYR ALCAMT 
		ALC5UPY1 SMKSTAT2 BMI HYPEV HYPDIFV DIBEV DOINGLWA WHYNOWKA VIGFREQW VIGMIN MODFREQW MODMIN
		ASISAD ASINERV ASIRSTLS ASIHOPLS ASIEFFRT ASIWTHLS ;
run;

data personsx2014; set lib2014.personsx;
	length PUBLICID $14; PUBLICID = trim(left(SRVY_YR||HHX||FMX||FPX)); 						 
	keep PUBLICID SEX AGE_P R_MARITL HISCODI3 MRACBPI2 EDUC1 PHSTAT GEOBRTH; 		 
run;

data familyxx2014; set lib2014.familyxx;  
	length FAM_ID $14; FAM_ID = trim(left(SRVY_YR||HHX||FMX));   
	keep FAM_ID RAT_CAT5;
run;

		/* Merge NHIS Data, keep those with data on the sample adult file */
		Proc sql;
				create table nhis2014_pre1 as
				select L.*, R.*
				from samadult2014 as L
				LEFT JOIN personsx2014 as R
				on L.PUBLICID = R.PUBLICID;
			quit;

		Proc sql;
				create table nhis2014 as
				select L.*, R.*
				from nhis2014_pre1 as L
				LEFT JOIN familyxx2014 as R
				on L.FAM_ID = R.FAM_ID;
			quit;


options notes;  /*Stop suppressing the notes in the log file */



/***************************************************/
/************Combine all NHIS data******************/
/***************************************************/
data nhis_all; 
	set nhis1997 nhis1998 nhis1999 nhis2000 nhis2001 nhis2002 nhis2003 nhis2004 nhis2005 
		nhis2006 nhis2007 nhis2008 nhis2009 nhis2010 nhis2011 nhis2012 nhis2013 nhis2014;
run;

proc contents data=nhis_all; run;


/*********************************************************************************************************************/
/*******************************************PART 2: RECODE NHIS DATA**************************************************/
/*********************************************************************************************************************/

data nhis_all; set nhis_all;
	
/*ALCOHOL USE - 1997-2014*/
	/*Data is consistent from 1997-2014, except that heavy drinking has a different variable name starting in 2014; make data consistent*/
	if SRVY_YR in (2014) then ALC5UPYR = ALC5UPY1;

	/*Remove 'refused', 'not ascertained' and 'don't know' */
	if ALC12MYR in (997,998,999) then ALC12MYR = .;  
	if ALCAMT in (97,98,99) then ALCAMT = .;
	if ALC5UPYR in (997,998,999) then ALC5UPYR =.;

	/* Input missing data - e.g. # days drank only asked to those who have had >12 drinks in any one year or their lifetime*/
	if ALC12MYR = . and (ALC1YR = 2 or ALCLIFE = 2) then ALC12MYR = 0;  
    /* input days drank in past year, for those who have had <12 drinks in any one year or their lifetime*/

	/* Create new variables - non-drinkers and heavy-drinkers*/
	none_drinker = 1;
		if ALC1YR = 1 or ALCLIFE = 1 then none_drinker = 0; 
		if ALC1YR in (7,8,9) and ALCLIFE in (7,8,9) then none_drinker = .;

	heavy_drinker = .;
		if ALC5UPYR >= 12 then heavy_drinker = 1;
		if .< ALC5UPYR <12 then heavy_drinker = 0;
		if none_drinker = 1 then heavy_drinker = 0;
		if ALCAMT >=5 then heavy_drinker = 1;

	drink_hist = .;
		if ALC12MYR = 0 and ALC1YR =2 then drink_hist = 0;	/* Never Drinkers*/
		if ALC12MYR = 0 and ALC1YR = 1 then drink_hist= 1;	/* Former Drinkers */ 
		if ALC12MYR >0 then drink_hist= 2;				/* Current Drinkers */ 
		if drink_hist = . and ALCAMT > 0 then drink_hist= 2;  /* Current Drinkers, though don't know how many drank days in past year*/ 
		if drink_hist = . and ALC12MYR=0 and ALCLIFE = 2 then drink_hist = 0; /* Never Drinkers*/
		if drink_hist = . and ALC12MYR=0 and ALCLIFE = 1 then drink_hist = 1; /* Former Drinkers */ 
			

	/******* Calculate daily grams of alcohol; assume 14 grams per drink*****/
	alc_daily_grams = ALC12MYR / 365 * ALCAMT * 14;   
			/* (days drank in past year) / 365 * (freq drinks per drinking day) * (14 grams per drink) */

			/* input alcohol_grams for those who didn’t drink in past year */				
			if ALC12MYR = 0 then alc_daily_grams = 0; 		 			
					
			/*input alcohol_grams if avg # drinks per drinking day is missing but # drinks a year is very small */
			if 0 < ALC12MYR <12 and ALCAMT= . and (ALC5UPYR=. or ALC5UPYR=0 or ALC5UPYR=1) then alc_daily_grams = 1;	
				
			
			/**** Add grams alcohol from heavy drinking days, among those whose avg drinks per day was <5 drinks; assume 5 drinks on heavy drinking days****/
				alc_daily_grams_heavy = ALC5UPYR / 365 * 5 * 14;  	
			  /* (heavy drinking days in past year) / 365 * (assume 5 drinks per day) * (14 grams per drink) */
					
				
			/* missing alc_grams_heavy set to 0 so that future calculations can still be calculated*/
			if alc_daily_grams_heavy = . then alc_daily_grams_heavy = 0;    
				
				/*if avg drinks per drinking data <5, then add heavy drinking grams to overall grams*/
			if .< ALCAMT <5 then alcohol_daily_grams = alc_daily_grams + alc_daily_grams_heavy;  
				
			/*if avg drinks per drinking data >=5, then do not add heavy drinking grams to overall grams*/	
			if alcohol_daily_grams =. then alcohol_daily_grams = alc_daily_grams;				

	/*****Create the Alcohol Use Categorical Variable*****/
		alcohol5 = .;
			/*Females*/
			if sex = 2 and alcohol_daily_grams = 0 then alcohol5 = 1; 		  /*abstinence*/
			if sex = 2 and 0 < alcohol_daily_grams <=20 then alcohol5 = 2;    /*low risk*/
			if sex = 2 and 20< alcohol_daily_grams <=40 then alcohol5 = 3; 	  /*medium risk*/
			if sex = 2 and 40< alcohol_daily_grams <=60 then alcohol5 = 4; 	  /*high risk*/
			if sex = 2 and alcohol_daily_grams >60 then alcohol5 = 5;   	  /*very high risk*/
			/*Males*/
			if sex = 1 and alcohol_daily_grams = 0 then alcohol5 = 1;  		  /*abstinence*/
			if sex = 1 and 0< alcohol_daily_grams <=40 then alcohol5 = 2;     /*low risk*/
			if sex = 1 and 40< alcohol_daily_grams <=60 then alcohol5 = 3;    /*medium risk*/
			if sex = 1 and 60< alcohol_daily_grams <=100 then alcohol5 = 4;   /*high risk*/
			if sex = 1 and alcohol_daily_grams >100 then alcohol5 = 5;   	  /*very high risk*/

	/*****Create the Alcohol Use Categorical Variable; separate abstainers into former and never drinkers*****/
		alcohol6 = .;
			if alcohol5 = 5 then alcohol6 = 6; /*very high risk*/
			if alcohol5 = 4 then alcohol6 = 5; /*high risk*/
			if alcohol5 = 3 then alcohol6 = 4; /*medium risk*/
			if alcohol5 = 2 then alcohol6 = 3; /*low risk*/
			if alcohol5 = 1 and drink_hist = 1 then alcohol6 = 2; /*Former drinker*/
			if alcohol5 = 1 and drink_hist = 0 then alcohol6 = 1; /*Never drinker*/
			


	/******* Calculate frequency of heavy drinking days *****/
	if ALC5UPYR = 0 then HED = 1;   		/* No HED (heavy episodic drinking)*/
		if ALC1YR = 2 then HED = 1;				/* no HED among abstainers */
		if ALCLIFE = 2 then HED = 1;			/* no HED among abstainers */
		if ALC12MYR = 0 then HED = 1;			/* no HED among abstainers */
	if 1 <= ALC5UPYR < 12 then HED = 2;   	/* HED less than once a month”, */
	if 12 <= ALC5UPYR < 52 then HED = 3;   	/* HED more than once a month but less than once a week */
	if ALC5UPYR >=52 then HED = 4;   		/* HED once a week or more */
	
		if HED = . and ALCAMT < 5 and ALC12MYR = 1 then HED = 1;
		if (HED = . or ALC5UPYR = 0) and ALCAMT>=5 and 1 <=ALC12MYR < 12 then HED = 2;   
		if (HED = . or ALC5UPYR = 0) and ALCAMT>=5 and 12 <=ALC12MYR < 52 then HED = 3;   
		if (HED = . or ALC5UPYR = 0) and ALCAMT>=5 and ALC12MYR >= 52 then HED = 4; 


/*AGE - 1997-2014*/
age = age_p;

/*SEX - 1997-2014*/
if sex = 2 then female = 1;
	if sex = 1 then female = 0;

/*MARITAL STATUS - 1997-2014 */
	if R_MARITL in (1,2,3,8) then married = 1;          /* Married/living together*/ 
		if R_MARITL in (4,5,6,7) then married = 0;    /* Not married nor living together*/ 
		if R_MARITL in (9) then married = .;

/* SMOKING - 1997-2014 */
	if SMKSTAT2 = 4 then smoking = 0;  	   /* never smoker */
		if SMKSTAT2 = 3 then smoking = 1;  /* former smoker */
		if SMKSTAT2 = 2 then smoking = 2;  /* Current some day smoker */
		if SMKSTAT2 = 1 then smoking = 3;  /* Current every day smoker */
		if SMKSTAT2 in (5,9) then smoking = .;  

/*HYPERTENSION - 1997-2014 */
	if HYPDIFV = 1 then hypertension = 1;				 /* Hypertension at least 2 times*/ 
		if HYPDIFV = 2 then hypertension = 0;			 /* No hypertension */ 
		if HYPEV = 2 then hypertension = 0;				 /* No hypertension */ 

/*DIABETES - 1997-2014*/
	if SRVY_YR in (2016, 2017, 2018) then DIBEV = DIBEV1; /* variable was renamed in 2016, but is similarly coded */
	if DIBEV = 1 then diabet = 2; 	 /* diabetes */ 
		if DIBEV = 3 then diabet = 1;  /* Borderline */ 
		if DIBEV = 2 then diabet = 0;  /* No diabetes */ 


/*EMPLOYMENT*/
	/* employment - 1997-2000*/ 
	if DOINGLW in (1,2) then employed = 1;         /* paid employment, student, or retired */
		if DOINGLW in (3,4) then employed = 0;   /* Not */
		if DOINGLW in (7,8,9) then employed = .;	
		if WHYNOWRK in(2,3) then employed = 1;	
		if WHYNOWRK in(1,4,5,6,7) then employed = 0; 

	/* employment - 2001-2003*/ 
	if DOINGLW1 in (1,2) then employed = 1;         	/* paid employment, student, or retired */
		if DOINGLW1 in (3,4,5) then employed = 0; 	/* Not */
		if DOINGLW1 in (7,8,9) then employed = .;	
		if WHYNOWK1 in(2,3,4,5) then employed = 1;	
		if WHYNOWK1 in(1,6,7,8,9,10) then employed = 0; 

	/* employment - 2004-2014*/ 
	if DOINGLWA in (1,2) then employed = 1;         	/* paid employment, student, or retired */
		if DOINGLWA in (3,4,5) then employed = 0; 	/* Not */
		if DOINGLWA in (7,8,9) then employed = .;	
		if WHYNOWKA in(2,3,4,5) then employed = 1;	
		if WHYNOWKA in(1,6,7,8,9,10) then employed = 0; 

/*BMI - 1997-2014*/
	if bmi= 99.99 then bmi = .;				/* remove 'unknown' category */
	if bmi < 18.5 then bmi_cat = 1;            	/*underweight*/
		if 18.5 =< bmi < 25 then bmi_cat = 2;  	/*healthy weight*/
		if 25 =< bmi < 30 then bmi_cat = 3;    	/*overweight*/
		if bmi >= 30 then bmi_cat = 4;   	   	/*obese*/
		if bmi = . then bmi_cat = .;
	
/* US Born */
	if USBORN_P = 1 then US_BORN = 1; 			/* US Born (50 states or DC; territories excluded) */
		if USBORN_P = 2 then US_BORN = 0;   	/* NOT US Born */
	if USBRTH_P = 1 then US_BORN = 1; 			/* US Born */
		if USBRTH_P = 2 then US_BORN = 0;   	/* NOT US Born */
	if GEOBRTH = 1 then US_BORN = 1; 			/* US Born */
		if GEOBRTH in (2,3) then US_BORN = 0;   /* NOT US Born */


/* ETHNICITY */
/* ethnicity - 2006-2014*/ 
	if HISCODI3 = 2 then ethnicity = 1; 			/*non-hispanic whites */
		if HISCODI3 = 3 then ethnicity = 2; 		/*non-hispanic blacks */
		if HISCODI3 = 1 then ethnicity = 3;  		/* hispanic*/ 
		if HISCODI3 in (4,5) then ethnicity = 4;  	/* non-hispanic other*/ 
		
	/* ethnicity - 2003-2005*/
	if HISCODI2 = 2 then ethnicity = 1; 			/*non-hispanic whites */
		if HISCODI2 = 3 then ethnicity = 2; 		/*non-hispanic blacks */
		if HISCODI2 = 1 then ethnicity = 3;    	    /* hispanic */ 
		if HISCODI2 = 4 then ethnicity = 4;  		/* non-hispanic other*/ 

	/* ethnicity - 2000-2002*/
	if HISCOD_I = 2 then ethnicity = 1; 			/*non-hispanic whites */
		if HISCOD_I = 3 then ethnicity = 2; 		/*non-hispanic blacks */
		if HISCOD_I = 1 then ethnicity = 3;    		/* hispanic*/ 
		if HISCOD_I = 4 then ethnicity = 4;    		/* non-hispanic other*/ 

	/* ethnicity - 1999*/
	if HISPCODR = 2 then ethnicity = 1; 			/*non-hispanic whites */
		if HISPCODR = 3 then ethnicity = 2; 		/*non-hispanic blacks */
		if HISPCODR = 1 then ethnicity = 3;   		/* hispanic*/ 
		if HISPCODR = 4 then ethnicity = 4;   		/* non-hispanic others*/ 

	/* ethnicity - 1998*/
	if HISPCODE = 2 then ethnicity = 1; 			/*non-hispanic whites */
		if HISPCODE = 3 then ethnicity = 2; 		/*non-hispanic blacks */
		if HISPCODE = 1 then ethnicity = 3;   		/* hispanic */ 
		if HISPCODE = 4 then ethnicity = 4;  		/* non-hispanic others*/ 

	/* ethnicity - 1997*/
	if ORIGIN = 2 and RACEREC = 1 then ethnicity = 1; 			/*non-hispanic whites */
		if ORIGIN = 2 and RACEREC = 2 then ethnicity = 2; 		/*non-hispanic blacks */
		if ORIGIN = 2 and RACEREC = 3 then ethnicity = 4; 		/*non-hispanic other */
		if ORIGIN = 1 then ethnicity = 3;   				    /* hispanic*/ 
	
	/* More detailed version of ethnicity*/
		/*Remove 'refused', 'not ascertained' and 'don't know' */
		if MRACE_P  in (97,98,99) then MRACE_P = .; 
		if MRACBR_P in (97,98,99) then MRACBR_P = .; 
		if MRACBP_I in (97,98,99) then MRACBP_I = .; 
		if MRACBPI2 in (97,98,99) then MRACBPI2 = .; 

	ethnicity_detail = ethnicity;
	if ethnicity_detail = 4 then ethnicity_detail=999; 
	if ethnicity_detail = 999 and (MRACE_P=3 or MRACBR_P=3 or MRACBP_I=3 or MRACBPI2=3) then ethnicity_detail = 4; /* American Indian/Alaska Native (AIAN)*/
	if ethnicity_detail = 999 and (MRACE_P in(6,7,12,15) or MRACBR_P in(6,7,12,15) or MRACBP_I in(6,7,12,15) or MRACBPI2 in(6,7,12,15)) then ethnicity_detail = 5; /* Asian Pacific Islander (API) */
	if ethnicity_detail = 999 and (MRACE_P in (16,17) or MRACBR_P in (16,17) or MRACBP_I in (16,17) or MRACBPI2 in (16,17)) then ethnicity_detail = 6; /* Other, including multiple*/
	if ethnicity_detail = 999 then ethnicity_detail = 6; /* Other */ 

/* ETHNICITY2 - differentiate US and not-US born Hispanics */
ethnicity2 = ethnicity; 
	if US_BORN = 1 and ethnicity = 3 then ethnicity2= 3;   /* US Borm Hispanic*/
	if US_BORN = 0 and ethnicity = 3 then ethnicity2= 5;   /* NOT US Borm Hispanic*/
	 

/* EDUCATION*/
	/* education -2004-2014 */
	if .< EDUC1 =<14 then edu = 1;   /* High school diploma or less */ 
		if EDUC1 in (15,16,17) then edu = 2;   /* Some college but no bachelor's degree */ 
		if EDUC1 in (18,19,20,21) then edu = 3;   /* bachelor's degree or more */ 
		if EDUC1 in (96,97,98,99) then edu = .; 

	/* education - 1997-2003 */
	if .< EDUC =<14 then edu = 1;   /* High school diploma or less */ 
		if EDUC in (15,16,17) then edu = 2;   /* Some college but no bachelor's degree */ 
		if EDUC in (18,19,20,21) then edu = 3;   /* bachelor's degree or more */ 
		if EDUC in (96,97,98,99) then edu = .; 

/* PHYSICAL ACTIVITY 1997-2014*/
		/*Remove 'refused', 'not ascertained' and 'don't know'; 'unable to' or 'never' coded as 0 */
		if VIGFREQW in (97, 98, 99) then VIGFREQW = .;  
			if VIGFREQW in (95, 96) then VIGFREQW = 0;
			if VIGMIN in (997,998,999) then VIGMIN = . ; 
			if VIGFREQW in (0, 95, 96) and VIGMIN=. then VIGMIN = 0;

		if MODFREQW in (97, 98, 99) then MODFREQW = .;  
			if MODFREQW in (95, 96) then MODFREQW = 0;  
			if MODMIN in (997,998,999) then MODMIN = . ; 
			if MODFREQW in (0, 95, 96) and MODMIN=. then MODMIN = 0;
	
		/* Calculate weekly minutes spent on physical activity*/
		min_wkly_vig_act = VIGFREQW * VIGMIN;
		min_wkly_mod_act = MODFREQW * MODMIN;

		min_wkly_mod_combined = (min_wkly_vig_act * 2) + min_wkly_mod_act; 

		/* Create categorical variable for physical activity*/
		if min_wkly_mod_combined = 0 then phy_act3 = 1;  				/*sedentary*/
			if 0 < min_wkly_mod_combined < 150  then phy_act3 = 2;  		/*Below recommendations*/
			if min_wkly_mod_combined >=150 then phy_act3 = 3;  			/*Meets recommendations*/
			if min_wkly_mod_combined = . then phy_act3 = .;
			if min_wkly_vig_act = . and min_wkly_mod_act >=150 then phy_act3 = 3;
			if min_wkly_vig_act >=75 and min_wkly_mod_act =. then phy_act3 = 3; 

/* Family income - ratio of family income to the poverty threshold*/
	/* income ratio -1997-2006 */
	if RAT_CAT in (1,2,3) then income = 1;   	  	/* poor: <100% of poverty threshold */
		if RAT_CAT in (4,5,6,7) then income = 2;  	/* near poor: 100-199% of poverty threshold */
		if RAT_CAT in (8,9,10,11) then income = 3;  /* middle income: 200-399% of poverty threshold */
		if RAT_CAT in (12,13,14) then income = 4;  	/* higher income: >=400% of poverty threshold */

	/* income ratio -2007-2013 */
	if RAT_CAT3 in (1,2,3,15) then income = 1;   	  	/* poor: <100% of poverty threshold */
		if RAT_CAT3 in (4,5,6,7,16) then income = 2;  	/* near poor: 100-199% of poverty threshold */
		if RAT_CAT3 in (8,9,10,11,17) then income = 3; 	/* middle income: 200-399% of poverty threshold */
		if RAT_CAT3 in (12,13,14,18) then income = 4;  	/* higher income: >=400% of poverty threshold */

	/* income ratio 2014 */
	if RAT_CAT5 in (1,2,3,15) then income = 1;   	  	/* poor: <100% of poverty threshold */
		if RAT_CAT5 in (4,5,6,7,16) then income = 2;  	/* near poor: 100-199% of poverty threshold */
		if RAT_CAT5 in (8,9,10,11,17) then income = 3; 	/* middle income: 200-399% of poverty threshold */
		if RAT_CAT5 in (12,13,14,18) then income = 4;  	/* higher income: >=400% of poverty threshold */

	if income = . then income = 5;						/* no income data*/


/*PSYCHOLOGICAL DISTRESS  - 1997-2014*/
	/*Data is consistent from 1997-2012, except that FEELINGS have a different variable name in 2013-2018; make data consistent*/
	if SRVY_YR in (2013, 2014) then SAD = ASISAD;
	if SRVY_YR in (2013, 2014) then EFFORT = ASIEFFRT;
	if SRVY_YR in (2013, 2014) then HOPELESS = ASIHOPLS;
	if SRVY_YR in (2013, 2014) then NERVOUS = ASINERV;
	if SRVY_YR in (2013, 2014) then RESTLESS = ASIRSTLS;
	if SRVY_YR in (2013, 2014) then WORTHLS = ASIWTHLS;
	

		/*Remove 'refused', 'not ascertained' and 'don't know' */
		If SAD in (7, 8, 9) then SAD=.;
		If EFFORT in (7, 8, 9) then EFFORT =.;
		If HOPELESS in (7, 8, 9) then HOPELESS =.;
		If NERVOUS in (7, 8, 9) then NERVOUS =.;
		If RESTLESS in (7, 8, 9) then RESTLESS =.;
		If WORTHLS in (7, 8, 9) then WORTHLS =.;
		
		/* Calculate psych distress score*/
		effort = 5 - effort;
		hopeless = 5 - hopeless;
		nervous = 5 - nervous;
		restless = 5 - restless;
		sad = 5 - sad;
		worthls = 5 - worthls;

		K6scale = hopeless + effort  + sad + worthls + nervous + restless;
			if hopeless = . or effort= . or sad= . or worthls = . or nervous= . or restless= . then K6scale = .;
		if K6scale <5 then K6scale3 = 1;   /*none to low psychological distress*/
		if 5=< K6scale <13 then K6scale3 = 2;  /*moderate psychological distress*/
		if K6scale >=13 then K6scale3 = 3;     /*severe psychological distress */
		if K6scale = . then K6scale3 =.; 


	/* health status, remove Refused, Not ascertained, Don't know*/
	if PHSTAT in (7,8,9) then PHSTAT = .;
	health_stat = PHSTAT; 

run;


			/* Check that the data were re-coded properly */
			/*alcohol use */ 
				proc freq data = nhis_all; tables alcohol5*ALC12MYR*ALCAMT*ALC1YR*ALCLIFE*ALC5UPYR*heavy_drinker /list missprint; run;
					proc freq data = nhis_all; by srvy_yr; tables alcohol5 /list missprint; run;
					/*export data to check*/
					data alc_use; set nhis_all; keep sex alcohol5 alcohol_daily_grams alc_daily_grams ALC12MYR ALCAMT ALC1YR ALCLIFE heavy_drinker ALC5UPYR alc_daily_grams_heavy; run;
					proc export data=alc_use dbms=xlsx replace outfile="C:\Users\klajd\Documents\CAMH\SIMAH\Public Use Data File\check_alcohol_use.xlsx"; run;
					proc freq data=nhis_all; tables drink_hist*ALC12MYR*ALC1YR*ALCLIFE / list missprint; run;
					proc freq data=nhis_all; tables alcohol5*alcohol6*drink_hist*ALC12MYR*ALCAMT*ALC1YR*ALCLIFE*ALC5UPYR / list missprint; run;
					proc freq data=nhis_all; tables HED*ALC12MYR*ALCAMT*ALC1YR*ALCLIFE*ALC5UPYR*alcohol5/list missprint; run;
			/* age */ proc means nmiss n mean std min max data=nhis_all; by srvy_yr; var age; run;
			/* sex */ proc freq data = nhis_all; by srvy_yr; tables female /list missprint; run;	
			/* married */ proc freq data = nhis_all; by srvy_yr; tables married /list missprint; run;	
						  proc freq data=nhis_all;  tables married*R_MARITL /list missprint; run;
			/* employed */ proc freq data = nhis_all; by srvy_yr; tables employed /list missprint; run;	
			/* smoking */ proc freq data = nhis_all; by srvy_yr; tables smoking /list missprint; run;				
						  proc freq data=nhis_all;  tables smoking*SMKSTAT2 /list missprint; run;
			/* hypertension */ proc freq data = nhis_all; by srvy_yr; tables hypertension /list missprint; run;	
				   		       proc freq data=nhis_all;  tables hypertension*HYPDIFV*HYPEV /list missprint; run;
			/* diabetes */ proc freq data = nhis_all; by srvy_yr; tables diabet /list missprint; run;	
						   proc freq data=nhis_all;  tables diabet*DIBEV /list missprint; run;
			/* BMI */ proc means nmiss n mean std min max data=nhis_all; by srvy_yr; var bmi; run;
			/* BMI_cat */ proc freq data = nhis_all; by srvy_yr; tables bmi_cat /list missprint; run;	
			/* ethnicity */ proc freq data = nhis_all; by srvy_yr; tables ethnicity /list missprint; run;	
			/* education */ proc freq data = nhis_all; by srvy_yr; tables edu /list missprint; run;	
						    proc freq data=nhis_all; tables edu*EDUC1*EDUC /list missprint; run;
			/* income*/ proc freq data=nhis_all; tables income; run;
			/* phy activity*/ 
			/* psych distress */ proc freq data=nhis_all; by srvy_yr; tables K6scale3; run;
									proc means nmiss n mean std min max data=nhis_all; by srvy_yr; var K6scale; run;


/***** Discard extra variables; specify variables to keep******/
data nhis_clean; set nhis_all;
	keep PUBLICID SRVY_YR INTV_QRT HHX FMX PX FPX WTFA_SA STRATUM PSU STRAT_P PSU_P PSTRAT PPSU
		 alcohol_daily_grams alcohol5 alcohol6 heavy_drinker drink_hist hed
		 female age ethnicity ethnicity_detail ethnicity2 married edu employed bmi bmi_cat smoking hypertension diabet phy_act3 income US_BORN health_stat
		 SAD NERVOUS RESTLESS HOPELESS EFFORT WORTHLS MHDSADWY MHDINTWY MHDSAD2W DEPYR DEPRESS
		 ADEPRSEV ADEPRSYR K6SCALE K6SCALE3;
run;
	
proc contents data=nhis_clean; run;




/*********************************************************************************************************************/
/*********************************PART 3: COMBINE & LABEL NHIS AND MORTALITY DATA*************************************/
/*********************************************************************************************************************/

proc sort data=nhis_clean; by PUBLICID; run;
proc sort data=mort.nhis_1997_to_2014_mort_2015; by PUBLICID; run;

Proc sql;
		create table nhis_mort_clean as
		select L.*, R.*
		from nhis_clean as L
		LEFT JOIN mort.nhis_1997_to_2014_mort_2015 as R
		on L.PUBLICID = R.PUBLICID;
	quit;

proc contents data=nhis_mort_clean; run;


/*Calculate follow-up time and sample weights*/
data nhis_mort_clean; set nhis_mort_clean;
	if MORTSTAT = 0 then DODYEAR = 2015;  	/*assign last possible follow-up to those alive */
	if MORTSTAT = 0 then DODQTR = 4;		/*assign last possible follow-up to those alive */
	yrs_followup = (DODYEAR - SRVY_YR) + ((DODQTR - INTV_QRT)*.25);
	
	/*for the public use file; some noise is introduced so follow-up time may be negative*/
	if (DODYEAR = SRVY_YR) and (DODQTR <= INTV_QRT) then yrs_followup = 0.1;  

	/* Sample weights and design variables*/ 
	new_weight = sa_wgt_new / 22;  /*divided by 22 because 22 survey years were pooled form NHIS*/
	new_psu = PSU;
	if new_psu = . then new_psu = psu_p;
	if new_psu = . then new_psu = PPSU ;
	new_stratum = STRATUM + 1000;
	if new_stratum = . then new_stratum = STRAT_P + 2000;
	if new_stratum = . then new_stratum = PSTRAT  + 3000;
	drop WTFA_SA sa_wgt_new PSU psu_p PPSU STRATUM STRAT_P PSTRAT;
run;  

/*Check data*/
proc print data=nhis_mort_clean (obs=20);
	var mortstat new_weight new_psu new_stratum;
run;



/**** Label the data ***/ 
/*First create a value library */
proc format library=combined.format;
		/*Values from NHIS data*/
		Value alc5_use 1 = "Abstinence"		2 = "Low risk"			3 = "Medium risk"
					  4 = "High risk"	      5 = "Very high risk";
		Value alc6_use 1 = "Never Drinker"		2 = "Former drinker" 	  3 = "Low risk"
					   4 = "Medium risk" 	5 = "High risk"	        6 = "Very high risk";	
		Value heavy_drink 0 = "Not occasional heavy drinker"
					1 = "Occasional heavy drinker";
		value hed 1 = "No HED"  2 = "HED < 1/month"  3="HED < 1 / week"   4 = "HED >= 1/week";
		value drink_hist 0="Never Drinker" 1="Former Drinker" 2="Current Drinker";
		value female 0 = "Male" 1="Female";
		value age 85="85+";
		value ethn 1 = "Non-Hispanic whites" 2="Non-Hispanic blacks" 3="Hispanic" 4="Other";
		value ethn_b 1 = "Non-Hispanic whites" 2="Non-Hispanic blacks" 3="Hispanic, US Born" 4="Other" 5="Hispanic, NOT US Born" ;
		value ethn_c 1 = "Non-Hispanic whites" 2="Non-Hispanic blacks" 3="Hispanic" 4="AI/AN" 5="API" 6="Other, including multiple";
		value married 1="Married or living together" 0 ="Not married/living togeter";
		value edu 1="Highschool or less" 2 = "Some college" 3= "Bachelor's or more";
		value employ 1 = "Paid employment, student or retired" 0 = "Not employed";
		value bmi 1 = "Underweight" 2 = "Healthy weight" 3="Overweight" 4="Obese";
		value smoke 0 = "Never smoker" 1 ="Former smoker" 2 ="Current some day smoker" 3="Current everyday smoker"; 
		value hyperten 1="Hypertension in past year" 0 = "No hypertension";
		value diabet 0 = "No" 1 ="Borderline" 2= "Yes";
		value phys_act 1 = "Sendentary" 2="Below recommendations" 3="Meets recommendations";
		value income 5 = "Missing" 1 = "Poor"  2="Near Poor"  3="Middle income"  4="Higher income";
		value health 1 = "Excellent" 2 = "Very good" 3 = "Good" 4="Fair" 5="Poor";
		value psydistress 1 = "None or low"  2="moderate"  3="severe";

		/*Values from Mortality data*/
		VALUE ELIGFMT  1 = "Eligible"		    
   					   2 = "Under age 18, not available for public release"		    
  					   3 = "Ineligible" ;
		VALUE MORTFMT  0 = "Assumed alive"	 1 = "Assumed deceased"	   . = "Ineligible or under age 18";
		VALUE FLAGFMT  0 = "No - Condition not listed as a multiple cause of death"
		    		   1 = "Yes - Condition listed as a multiple cause of death"  
		               . = "Assumed alive, under age 18,ineligible for mortality follow-up, or MCOD not available";
		VALUE QRTFMT   1 = "January-March"	 2 = "April-June"    3 = "July-September"	 4 = "October-December" 
		  	         . = "Ineligible, under age 18, or assumed alive";
		VALUE DODYFMT . = "Ineligible, under age 18, or assumed alive";
		VALUE $UCODFMT
			"001" = "Diseases of heart (I00-I09, I11, I13, I20-I51)"
			"002" = "Malignant neoplasms (C00-C97)"
			"003" = "Chronic lower respiratory diseases (J40-J47)"
			"004" = "Accidents (unintentional injuries) (V01-X59, Y85-Y86)"
			"005" = "Cerebrovascular diseases (I60-I69)"
			"006" = "Alzheimer's disease (G30)"
			"007" = "Diabetes mellitus (E10-E14)"
			"008" = "Influenza and pneumonia (J09-J18)"
			"009" = "Nephritis, nephrotic syndrome and nephrosis (N00-N07, N17-N19, N25-N27)"
			"010" = "All other causes (residual)" 
			"   " = "Ineligible, under age 18, assumed alive, or no cause of death data" ;
run;

/* Second, assigned the label libary to the data*/
data combined.nhis_mort_clean; set nhis_mort_clean;
OPTIONS FMTSEARCH=(combined.format);
Format 	alcohol5 alc5_use. 		alcohol6 alc6_use. 			heavy_drinker heavy_drink. 	
		hed hed.				drink_hist drink_hist.		female female.
		age age. 				ethnicity ethn.				ethnicity_detail ethn_c.
		ethnicity2 ethn_b. 		married married.			edu edu.					
		employed employ.		bmi_cat bmi.				smoking smoke.
		hypertension hyperten.	diabet diabet.				phy_act3 phys_act.
		income income.			health_stat health.			K6SCALE3 psydistress.
		ELIGSTAT ELIGFMT.			MORTSTAT MORTFMT.			UCOD_LEADING UCODFMT.
		DIABETES FLAGFMT.   		HYPERTEN FLAGFMT. 		DODQTR QRTFMT.           
		DODYEAR DODYFMT.;
run;

	
/*check*/
	proc freq data=combined.nhis_mort_clean; 
			tables alcohol5 alcohol6 heavy_drinker drink_hist HED female age ethnicity ethnicity_detail ethnicity2 married edu employed 
			        bmi_cat smoking hypertension diabet phy_act3 income health_stat us_born ELIGSTAT MORTSTAT K6SCALE3;
	run;
