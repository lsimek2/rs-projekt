%macro ODSOFF;
	ods graphics off;
	ods exclude all;
	ods noresults;
%mend;

%macro ODSOn;
	ods graphics on;
	ods exclude none;
	ods results;
%mend;