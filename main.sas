/* Luka Šimek */
/* Projektni zadatak iz kolegija Računarska statistika */
/* januar 2025. */

ods trace off;
%include '/home/u64055593/sasuser.v94/projekt/fleishman.sas';
%include '/home/u64055593/sasuser.v94/projekt/odsoffon.sas';

%let seed=112025;
%let mu=0;
%let sigma=5;
%let nreps=500;

%let n_list = 10 15 20 50 100; %let n_len = 5;
%let gamma1_list = -2 0 2; %let gamma1_len = 3;
%let gamma2_list = 0 6 11; %let gamma2_len = 3;
/* gamma2=-3 je generalno nemoguce pa je izuzeto iz liste */
/* gamma2 je excess kurtosis */

data fl;
	gamma1=-999;
	gamma2=0;a=0;b=0;c=0;d=0;skewness=0;kurtosis=0;
	output;
run;

/* stvaramo tablicu gdje svakom paru (gamma1, gamma2) */
/* pridruzimo Fleishmanove koeficijente a, b, c, d */
%macro abcd;
%do j=1 %to &gamma1_len;
	%do k=1 %to &gamma2_len;
		%let gamma1 = %scan(&gamma1_list, &j, ' ');
		%let gamma2 = %scan(&gamma2_list, &k, ' ');
		
		data skewkurt;
			skewness=&gamma1;
			kurtosis=&gamma2+3;
			output;
		run;
		
		%odsoff;
		%fleishman;
		%odson;
		
		data fl_novo;
			gamma1=&gamma1;
			gamma2=&gamma2;
			merge fleishman;
			a=-c;
		run;
		
		proc append data=fl_novo base=fl; run;
	%end;
%end;

data fl;
	set fl;
	where gamma1 ne -999;
	keep gamma1 gamma2 a b c d;
run;
%mend;

%abcd;
/* proc print data=fl;run; */

data sim;
	gamma1=-999;
	gamma2=0;n=0;x=0;sample_id=0;i=0;
	output;
run;

/* generiranje podataka */
%macro sim;
%do i=1 %to &n_len;
	%do j=1 %to &gamma1_len;
		%do k=1 %to &gamma2_len;
			%let n = %scan(&n_list, &i, ' ');
			%let gamma1 = %scan(&gamma1_list, &j, ' ');
			%let gamma2 = %scan(&gamma2_list, &k, ' ');
			
			data fl_curr;
				set fl;
				where gamma1 eq &gamma1 and gamma2 eq &gamma2;
				keep a b c d;
			run;
			
			data sim_novo(keep=gamma1 gamma2 n sample_id i x);
				call streaminit(&seed+&n*&gamma1*&gamma2);
				merge fl_curr;
				gamma1=&gamma1;
				gamma2=&gamma2;
				n=&n;
				
				do sample_id=1 to &nreps;
					do i=1 to &n;
						x = rand('normal');
						x = a + b*x + c*x**2 + d*x**3;
						x = x*&sigma + &mu;
						output;
					end;
				end;
			run;
			
			proc append data=sim_novo base=sim; run;
			
		%end;
	%end;
%end;

data sim;
	set sim;
	where gamma1 ne -999
		and gamma1**2-2 le gamma2**2
		and 1.7022*gamma1**2-1.15 lt gamma2;
run;
%mend;

%sim;
/* proc print data=sim(obs=250);run; */

/* data sim_skraceni; */
/* 	set sim; */
/* 	where n=10 and gamma1=-2 and gamma2=6 and sample_id <= 5; */
/* run; */

/* racunanje pouzdanih intervala i pripadnosti */
%odsoff;
ods output basicintervals=pi;
proc univariate data=sim cibasic(alpha=0.1);
	var x;
	by n gamma1 gamma2 sample_id;
run;

ods output basicintervals=pi2;
proc univariate data=sim cibasic(alpha=0.05);
	var x;
	by n gamma1 gamma2 sample_id;
run;
%odson;

data pi(rename=(lowercl=l90 uppercl=u90));
	set pi;
	where parameter='Mean';
run;

data pi2(rename=(lowercl=l95 uppercl=u95));
	set pi2;
	where parameter='Mean';
run;

data pi;
	set pi;
	merge pi2;
	in90 = l90 lt &mu and &mu lt u90;
	in95 = l95 lt &mu and &mu lt u95;
	keep n gamma1 gamma2 sample_id l90 u90 l95 u95 in90 in95;
run;

/* proc print data=pi(obs=100);run; */

/* rezultati */
%odsoff;
ods output moments=rez;
proc univariate data=pi;
	var in90 in95;
	by n gamma1 gamma2;
run;
%odson;

data rez;
	set rez;
	where label1='Mean';
	if varname='in90' then nominalna=0.90; else nominalna=0.95;
	keep n gamma1 gamma2 nominalna nvalue1;
	rename nvalue1=empirijska;
run;

proc sort data=rez;
	by nominalna;
run;

ods listing gpath='/home/u64055593/sasuser.v94/projekt/';
options printerpath=png nodate papersize=('4.75in','10.75in');
ods _all_ close;
ods printer file="/home/u64055593/sasuser.v94/projekt/tablica90.png";
proc print data=rez(where=(nominalna=0.9));run;

ods printer file="/home/u64055593/sasuser.v94/projekt/tablica95.png";
proc print data=rez(where=(nominalna=0.95));run;
ods printer close;
ods listing;

data rez;
	set rez;
	gamma12 = catx('/', gamma1, gamma2);
	ngamma1 = catx('/', n, gamma1);
	ngamma2 = catx('/', n, gamma2);
run;

ods listing gpath='/home/u64055593/sasuser.v94/projekt/';
ods graphics / reset imagename="grafn90" imagefmt=png;
title 'ovisnost o n, 90%-tni interval';
proc sgplot data=rez(where=(nominalna=0.90));
series x=n y=empirijska / group=gamma12 lineattrs=(thickness=2 pattern=solid);
refline 0.90 / axis=y lineattrs=(thickness=1 color=black pattern=dash);
run; 

ods graphics / reset imagename="grafn95" imagefmt=png;
title 'ovisnost o n, 95%-tni interval';
proc sgplot data=rez(where=(nominalna=0.95));
series x=n y=empirijska / group=gamma12 lineattrs=(thickness=2 pattern=solid);
refline 0.95 / axis=y lineattrs=(thickness=1 color=black pattern=dash);
run;

ods graphics / reset imagename="graf190" imagefmt=png;
title 'ovisnost o gamma1, 90%-tni interval';
proc sgplot data=rez(where=(nominalna=0.90));
series x=gamma1 y=empirijska / group=ngamma2 lineattrs=(thickness=2 pattern=solid);
refline 0.90 / axis=y lineattrs=(thickness=1 color=black pattern=dash);
run;

ods graphics / reset imagename="graf195" imagefmt=png;
title 'ovisnost o gamma1, 95%-tni interval';
proc sgplot data=rez(where=(nominalna=0.95));
series x=gamma1 y=empirijska / group=ngamma2 lineattrs=(thickness=2 pattern=solid);
refline 0.95 / axis=y lineattrs=(thickness=1 color=black pattern=dash);
run;

ods graphics / reset imagename="graf290" imagefmt=png;
title 'ovisnost o gamma2, 90%-tni interval';
proc sgplot data=rez(where=(nominalna=0.90));
series x=gamma2 y=empirijska / group=ngamma1 lineattrs=(thickness=2 pattern=solid);
refline 0.90 / axis=y lineattrs=(thickness=1 color=black pattern=dash);
run;

ods graphics / reset imagename="graf295" imagefmt=png;
title 'ovisnost o gamma2, 95%-tni interval';
proc sgplot data=rez(where=(nominalna=0.95));
series x=gamma2 y=empirijska / group=ngamma1 lineattrs=(thickness=2 pattern=solid);
refline 0.95 / axis=y lineattrs=(thickness=1 color=black pattern=dash);
run;
