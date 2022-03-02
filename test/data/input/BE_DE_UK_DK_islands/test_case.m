%TNEP optimizor input file

function mpc = owpp_tnep_map
mpc.version = '2';
mpc.baseMVA = 100;

%% bus data
	%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
	mpc.bus = [
1 3 0 0 0 0 1 1 0 220 1 1.1 0.9;
2 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
3 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
4 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
5 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
6 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
7 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
8 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
9 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
10 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
11 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
12 2 0 0 0 0 1 1 0 220 1 1.1 0.9;
];
%% generator data
	%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
	mpc.gen = [
1 0 0 0 0 1 100 1 220 0 0 0 0 0 0 0 0 0 0 0 1;
2 0 0 0 0 1 100 1 220 0 0 0 0 0 0 0 0 0 0 0 1;
3 0 0 0 0 1 100 1 220 0 0 0 0 0 0 0 0 0 0 0 1;
4 0 0 0 0 1 100 1 220 0 0 0 0 0 0 0 0 0 0 0 1;
5 0 0 0 0 1 100 1 220 0 0 0 0 0 0 0 0 0 0 0 1;
6 0 0 0 0 1 100 1 220 0 0 0 0 0 0 0 0 0 0 0 1;
7 0 0 0 0 1 100 1 220 0 0 0 0 0 0 0 0 0 0 0 1;
8 0 0 0 0 1 100 1 220 0 0 0 0 0 0 0 0 0 0 0 1;
9 0 0 0 0 1 100 1 20 0 0 0 0 0 0 0 0 0 0 0 1;
10 0 0 0 0 1 100 1 40 0 0 0 0 0 0 0 0 0 0 0 1;
11 0 0 0 0 1 100 1 40 0 0 0 0 0 0 0 0 0 0 0 1;
12 0 0 0 0 1 100 1 120 0 0 0 0 0 0 0 0 0 0 0 1;
];
%%Additional generator fields, type=0(wind generator), type=1(onshore market ie load and generator created)
	%column_names% 		type
	mpc.gen_type = [
1
1
1
1
1
1
1
1
0
0
0
0
];

	%% branch data
	%	fbus	tbus	r	       x	    b	  rateA	rateB	rateC	ratio	angle	status	angmin	angmax
	mpc.branch = [
	];

	%candidate branch data
	%column_names%	f_bus	t_bus	br_r	br_x	br_b	rate_a	rate_b	rate_c	tap	shift	br_status	angmin	angmax	mva	construction_cost
	mpc.ne_branch = [
1 2 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
1 3 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
1 4 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
1 5 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
1 6 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
1 7 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
1 8 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
1 9 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
1 10 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
1 11 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
1 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
2 3 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
2 4 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
2 5 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
2 6 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
2 7 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
2 8 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
2 9 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
2 10 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
2 11 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
2 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
3 4 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
3 5 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
3 6 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
3 7 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
3 8 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
3 9 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
3 10 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
3 11 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
3 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
4 5 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
4 6 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
4 7 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
4 8 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
4 9 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
4 10 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
4 11 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
4 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
5 6 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
5 7 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
5 8 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
5 9 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
5 10 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
5 11 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
5 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
6 7 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
6 8 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
6 9 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
6 10 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
6 11 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
6 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
7 8 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
7 9 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
7 10 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
7 11 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
7 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
8 9 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
8 10 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
8 11 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
8 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
9 10 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
9 11 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
9 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
10 11 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
10 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
11 12 0.004 0.04 0 1000 1000 1000 0 0 1 -60 60 100 10;
];
%%-----  OPF Data  -----%%
	%% generator cost data
	%	1	startup	shutdown	n	x1	y1	...	xn	yn
	%	2	startup	shutdown	n	c(n-1)	...	c0
	mpc.gencost = [
2 0 0 2 220 0
2 0 0 2 220 0
2 0 0 2 220 0
2 0 0 2 220 0
2 0 0 2 220 0
2 0 0 2 220 0
2 0 0 2 220 0
2 0 0 2 220 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
2 0 0 2 0 0
];
%% existing dc bus data
	%column_names%   busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin Cdc
	mpc.busdc = [
1 1 0 1 300 1.1 0.9 0;
2 2 0 1 300 1.1 0.9 0;
3 3 0 1 300 1.1 0.9 0;
4 4 0 1 300 1.1 0.9 0;
5 5 0 1 300 1.1 0.9 0;
6 6 0 1 300 1.1 0.9 0;
7 7 0 1 300 1.1 0.9 0;
8 8 0 1 300 1.1 0.9 0;
9 9 0 1 300 1.1 0.9 0;
10 10 0 1 300 1.1 0.9 0;
11 11 0 1 300 1.1 0.9 0;
12 12 0 1 300 1.1 0.9 0;
];
%% Candidate Branches here - refer to MATACDC 1.0 Users Manual for description of fields
	 %% candidate branches
	 %column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB   rateC status cost
	 mpc.branchdc_ne = [
1 2 0 0 0 10000 10000 10000 1 0;
1 3 0 0 0 10000 10000 10000 1 0;
1 4 0 0 0 10000 10000 10000 1 0;
1 5 0 0 0 10000 10000 10000 1 0;
1 6 0 0 0 10000 10000 10000 1 0;
1 7 0 0 0 10000 10000 10000 1 0;
1 8 0 0 0 10000 10000 10000 1 0;
1 9 0 0 0 10000 10000 10000 1 0;
1 10 0 0 0 10000 10000 10000 1 0;
1 11 0 0 0 10000 10000 10000 1 0;
1 12 0 0 0 10000 10000 10000 1 0;
2 3 0 0 0 10000 10000 10000 1 0;
2 4 0 0 0 10000 10000 10000 1 0;
2 5 0 0 0 10000 10000 10000 1 0;
2 6 0 0 0 10000 10000 10000 1 0;
2 7 0 0 0 10000 10000 10000 1 0;
2 8 0 0 0 10000 10000 10000 1 0;
2 9 0 0 0 10000 10000 10000 1 0;
2 10 0 0 0 10000 10000 10000 1 0;
2 11 0 0 0 10000 10000 10000 1 0;
2 12 0 0 0 10000 10000 10000 1 0;
3 4 0 0 0 10000 10000 10000 1 0;
3 5 0 0 0 10000 10000 10000 1 0;
3 6 0 0 0 10000 10000 10000 1 0;
3 7 0 0 0 10000 10000 10000 1 0;
3 8 0 0 0 10000 10000 10000 1 0;
3 9 0 0 0 10000 10000 10000 1 0;
3 10 0 0 0 10000 10000 10000 1 0;
3 11 0 0 0 10000 10000 10000 1 0;
3 12 0 0 0 10000 10000 10000 1 0;
4 5 0 0 0 10000 10000 10000 1 0;
4 6 0 0 0 10000 10000 10000 1 0;
4 7 0 0 0 10000 10000 10000 1 0;
4 8 0 0 0 10000 10000 10000 1 0;
4 9 0 0 0 10000 10000 10000 1 0;
4 10 0 0 0 10000 10000 10000 1 0;
4 11 0 0 0 10000 10000 10000 1 0;
4 12 0 0 0 10000 10000 10000 1 0;
5 6 0 0 0 10000 10000 10000 1 0;
5 7 0 0 0 10000 10000 10000 1 0;
5 8 0 0 0 10000 10000 10000 1 0;
5 9 0 0 0 10000 10000 10000 1 0;
5 10 0 0 0 10000 10000 10000 1 0;
5 11 0 0 0 10000 10000 10000 1 0;
5 12 0 0 0 10000 10000 10000 1 0;
6 7 0 0 0 10000 10000 10000 1 0;
6 8 0 0 0 10000 10000 10000 1 0;
6 9 0 0 0 10000 10000 10000 1 0;
6 10 0 0 0 10000 10000 10000 1 0;
6 11 0 0 0 10000 10000 10000 1 0;
6 12 0 0 0 10000 10000 10000 1 0;
7 8 0 0 0 10000 10000 10000 1 0;
7 9 0 0 0 10000 10000 10000 1 0;
7 10 0 0 0 10000 10000 10000 1 0;
7 11 0 0 0 10000 10000 10000 1 0;
7 12 0 0 0 10000 10000 10000 1 0;
8 9 0 0 0 10000 10000 10000 1 0;
8 10 0 0 0 10000 10000 10000 1 0;
8 11 0 0 0 10000 10000 10000 1 0;
8 12 0 0 0 10000 10000 10000 1 0;
9 10 0 0 0 10000 10000 10000 1 0;
9 11 0 0 0 10000 10000 10000 1 0;
9 12 0 0 0 10000 10000 10000 1 0;
10 11 0 0 0 10000 10000 10000 1 0;
10 12 0 0 0 10000 10000 10000 1 0;
11 12 0 0 0 10000 10000 10000 1 0;
];

	%% existing converters
	%column_names%   busdc_i busac_i type_dc type_ac P_g   Q_g  islcc  	Vtar 		rtf  xtf  transformer tm   bf  filter    rc      xc  reactor   basekVac Vmmax   Vmmin   Imax     status   LossA  LossB  LossCrec LossCinv   droop   Pdcset     Vdcset  dVdcset Pacmax Pacmin Qacmax   Qacmin cost
	mpc.convdc = [
1 1 2 3 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 19.25;
2 2 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 19.25;
3 3 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 19.25;
4 4 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 19.25;
5 5 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 19.25;
6 6 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 19.25;
7 7 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 19.25;
8 8 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 19.25;
9 9 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 51.5;
10 10 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 51.5;
11 11 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 51.5;
12 12 3 2 400000 0 0 1 0.001 0.1 0 1 0.08 0 0.001 0.09 0 220 1.1 0.9 100000 1 0 0 0 0 0.005 -52.7 1.0079 0 4000 -4000 4000 -4000 51.5;
];

	%trans, filter, reactor, vmin vmax same as conv
	%% candidate converters
	%column_names%   busdc_i busac_i type_dc type_ac P_g   Q_g  islcc  Vtar    rtf   xtf  transformer tm   bf 	filter    rc     xc  reactor   basekVac Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  droop     Pdcset    Vdcset  dVdcset Pacmax Pacmin Qacmax Qacmin cost
	mpc.convdc_ne = [
	                1000       3      1       1    400000    0   	0     1.0   0.001  0.1       0 			 1 	0.08 		0 		0.001   0.09 		0  				220    1.1     0.9    100000      1     	 0     0        0       0      0.0050    -52.7     1.0079     0  		1  -1   1    -1   1000000;
	];