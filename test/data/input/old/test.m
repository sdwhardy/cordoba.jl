function mpc = case9
%CASE9    Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from Joe H. Chows book, p. 70.

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%%Buses type 3= swing, type2=PV refer to Matpower 6.0 Users Manual for description of fields
%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
1	3	0	0	0	0	1	1	0	220	1	1.1	0.9;
2	2	0	0	0	0	1	1	0	220	1	1.1	0.9;
3	2	0	0	0	0	1	1	0	220	1	1.1	0.9;
4	2	0	0	0	0	1	1	0	220	1	1.1	0.9;
];

%%Generators apf=0(wind generator),apf=1(onshore market ie load and genrator created) refer to Matpower 6.0 Users Manual for description of fields
%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
1	0		0	0	0	1	100	1	835	0	0	0	0	0	0	0	0	0	0	0	1;
2	0		0	0	0	1	100	1	835	0	0	0	0	0	0	0	0	0	0	0	1;
3	0		0	0	0	1	100	1	835	0	0	0	0	0	0	0	0	0	0	0	1;
4	0		0	0	0	1	100	1	835	0	0	0	0	0	0	0	0	0	0	0	0;
];

%%AC connections here (dummy only) refer to Matpower 6.0 Users Manual for description of fields
%% branch data
%	fbus	tbus	r	       x	    b	  rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1   2   0.0040   0.0400   0.00   0  0  0  0  0  0 -60  60;
];
%% Candidate AC connections here (dummy only) refer to Matpower 6.0 Users Manual for description of fields
%candidate branch data
%column_names%	f_bus	t_bus	br_r	br_x	br_b	rate_a	rate_b	rate_c	tap	shift	br_status	angmin	angmax	mva	construction_cost
mpc.ne_branch = [
1   2   0.0040   0.0400   0.00   0  0  0  0  0  0 -60  60	100.0	100000;
];

%%-----  OPF Data  -----%%
%% Cost of generators here - refer to Matpower 6.0 Users Manual for description of fields
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
2	 0	 0	 2	  10 0
2	 0	 0	 2	  10 0
2	 0	 0	 2	  10 0
2	 0	 0	 2	  10 0
];
%% All Candidate equipment has equivalent existing infrastructure field as well ie  mpc.branchdc_ne <->  mpc.branchdc similar to AC
%% Candidate DC buses here - refer to MATACDC 1.0 User’s Manual for description of fields
%% candidate dc bus data
%column_names%   busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin  Cdc
mpc.busdc_ne = [
	1              1       0       1       300         1.1     0.9     0;
	2              2       0       1       300         1.1     0.9     0;
	3              3       0       1       300         1.1     0.9     0;
	4              4       0       1       300         1.1     0.9     0;
];

%% Candidate Branches here - refer to MATACDC 1.0 Users Manual for description of fields
 %% candidate branches
 %column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB   rateC status cost
 mpc.branchdc_ne = [
 1	2	0.0	0.00	 0.00  1	 		1	 -1	 1.0	0;
 1	3	0.0	0.00	 0.00  1	 		1	 -1	 1.0	0;
 2	3	0.0	0.00	 0.00  1	 		1	 -1	 1.0	0;
 4	1	0.0	0.00	 0.00  1	 		1	 1	 1.0	0;
 4	2	0.0	0.00	 0.00  1	 		1	 1	 1.0	0;
 4	3	0.0	0.00	 0.00  1	 		1	 1	 1.0	0;
 ];

%% Candidate Converters here - refer to MATACDC 1.0 Users Manual for description of fields
%% All costs are lumped into candidate cable
%trans, filter, reactor, vmin vmax same as conv
%% candidate converters
%column_names%   busdc_i busac_i type_dc type_ac P_g   Q_g  islcc  Vtar    rtf   xtf  transformer tm   bf 	filter    rc     xc  reactor   basekVac Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  droop     Pdcset    Vdcset  dVdcset Pacmax Pacmin Qacmax Qacmin cost
mpc.convdc_ne = [
                1       1      2       3    4000    0   	0     1.0   0.001  0.1       0 			 1 	0.08 		0 		0.001   0.09 		0  				220    1.1     0.9    100000      1     	 0     0        0       0      0.0050    -52.7     1.0079     0  		4000  -4000   4000    -4000   0;
								2       2      3       2    4000    0   	0     1.0   0.001  0.1       0 			 1 	0.08 		0 		0.001   0.09 		0  				220    1.1     0.9    100000      1     	 0     0        0       0      0.0050    -52.7     1.0079     0  		4000  -4000   4000    -4000   0;
								3       3      3       2    4000    0   	0     1.0   0.001  0.1       0 			 1 	0.08 		0 		0.001   0.09 		0  				220    1.1     0.9    100000      1     	 0     0        0       0      0.0050    -52.7     1.0079     0  		4000  -4000   4000    -4000   0;
								4       4      3       2    4000    0   	0     1.0   0.001  0.1       0 			 1 	0.08 		0 		0.001   0.09 		0  				220    1.1     0.9    100000      1     	 0     0        0       0      0.0050    -52.7     1.0079     0  		4000  -4000   4000    -4000   0;
];
