options nonotes;
libname markers 'C:\Users\Brent.Hulke\Downloads\datafile';
ods select none;
%macro freq (range);
dm 'log;clear;';

proc freq data=markers.bialw nlevels;
tables &range;
ods output NLevels=levels;
run;

data levels (keep=TableVar);
set levels;
if NLevels = 3 then; else delete;
run;

%let i=1;

data _null_;
set levels nobs=n;
call symputx('nob', n);
stop;
run;

%do %until (&i>&nob);
	data _null_;
	set levels (obs=&i firstobs=&i);
	call symput('var', TableVar);
	stop;
	run;

	%put &var;

	proc freq data=markers.bialw ;
	tables &var / chisq(testp=(0.484375 0.03125 0.484375));
	output out=pre_chisq_test pchi;
	run;

	data pre_chisq_test;
	set pre_chisq_test;
	format Variant $16.;
	Variant="&var";
	run;

	proc append base=Varlist data=pre_chisq_test force;
	run;

	%let i=%eval(&i+1);
%end;
%mend;
%freq(_Variant_80001--_Variant_120000); /*_Variant_289678 */
%freq(_Variant_120001--_Variant_160000);
%freq(_Variant_160001--_Variant_200000);
%freq(_Variant_200001--_Variant_240000);
%freq(_Variant_240001--_Variant_280000);
%freq(_Variant_280001--_Variant_289678);

data Varlist_final (keep=P_PCHI Variant ID2);
set Varlist;
if P_PCHI < 0.10 then delete;
format ID2 BEST12.;
ID2=scan(Variant, 2, '_');
run;

proc sql;
create table work.bialf as (select * from work.Varlist_final as b
left join markers.bial as a  
on b.ID2 = a.ID)
;
run;

data bialf;
set bialf (drop=P_PCHI Variant ID2);
run;
quit;

