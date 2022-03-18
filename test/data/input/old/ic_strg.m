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

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
1	3	0	0	0	0	1	1	0	220	1	1.1	0.9;
2	2	0	0	0	0	1	1	0	220	1	1.1	0.9;
3	2	0	0	0	0	1	1	0	220	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
1	0		0	0	0	1	100	1	835	0	0	0	0	0	0	0	0	0	0	0	0;
1	0		0	0	0	1	100	1	-835	-835	0	0	0	0	0	0	0	0	0	0	0;
2	0		0	0	0	1	100	1	835	0	0	0	0	0	0	0	0	0	0	0	0;
2	0		0	0	0	1	100	1	-835	-835	0	0	0	0	0	0	0	0	0	0	0;
3	0		0	0	0	1	100	1	100	100	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	       x	    b	  rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
];

%candidate branch data
%column_names%	f_bus	t_bus	br_r	br_x	br_b	rate_a	rate_b	rate_c	tap	shift	br_status	angmin	angmax	mva	construction_cost
mpc.ne_branch = [
3   1   0.0040   0.0400   0.00   500  500  500  0  0  1 -60  60	100.0	0;
3   2   0.0040   0.0400   0.00   500  500  500  0  0  1 -60  60	100.0	0;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
2	 0	 0	 2	  10 0
2	 0	 0	 2	  10 0
2	 0	 0	 2	  100 0
2	 0	 0	 2	  100 0
2	 0	 0	 2	  0 0
];

%% existing dc bus data
%column_names%   busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin Cdc
mpc.busdc = [
1              1       0       1       300         1.1     0.9     0;
2              2       0       1       300         1.1     0.9     0;
];

%% candidate dc bus data
%column_names%   busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin  Cdc
mpc.busdc_ne = [
	3              3       0       1       300         1.1     0.9     0;
];

%% existing dc branches
%column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB   rateC status
mpc.branchdc = [
1	2	0.0	0.00	 0.00  0	 0	 0	 0.0;
];

 %% candidate branches
 %column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB   rateC status cost
 mpc.branchdc_ne = [
 1	3	0.0	0.00	 0.00  1500	 		1500	 1500	 1.0	0;
 2	3	0.0	0.00	 0.00  1500	 		1500	 1500	 1.0	0;
 ];

 %% existing converters
 %column_names%   busdc_i busac_i type_dc type_ac P_g   Q_g  islcc  Vtar rtf    xtf  transformer tm   bf        filter    rc      xc    reactor   basekVac Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  droop       Pdcset    Vdcset  dVdcset Pacmax Pacmin Qacmax   Qacmin
 mpc.convdc = [
                     1       1       2       3    1000    0   	0   1.0  0.001  0.1       0 	 1   0.08 	  0      0.001   0.09      0  	      220    1.1     0.9    1000      1       0     0        0       0      0.0050    -52.7       1.0079     0     3500  -3500   3500    -3500;
                     2       2       3       2    1000    0   	0   1.0  0.001  0.1       0 	 1   0.08 	  0 	 0.001   0.09 	   0  	      220    1.1     0.9    1000      1       0     0        0       0      0.0050    -52.7       1.0079     0     3500  -3500   3500    -3500;
];

%trans, filter, reactor, vmin vmax same as conv
%% candidate converters
%column_names%   busdc_i busac_i type_dc type_ac P_g   Q_g  islcc  Vtar    rtf   xtf  transformer tm   bf 	filter    rc     xc  reactor   basekVac Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  droop     Pdcset    Vdcset  dVdcset Pacmax Pacmin Qacmax Qacmin cost
mpc.convdc_ne = [
                3       3      1       1    1000    0   	0     1.0   0.001  0.1       0 			 1 	0.08 		0 		0.001   0.09 		0  				220    1.1     0.9    1000      1     	 0     0        0       0      0.0050    -52.7     1.0079     0  		3500  -3500   3500    -3500   0;
];

% hours
mpc.time_elapsed = 1.0

%% storage data
%column_names%   storage_bus ps 	qs 	energy  energy_rating charge_rating  discharge_rating  charge_efficiency  discharge_efficiency  thermal_rating  qmin  	qmax  	r  		x  p_loss  	q_loss  status eq_cost inst_cost co2_cost 	max_energy_absorption 	stationary_energy_inflow 	stationary_energy_outflow 	self_discharge_rate	 	cost_abs 	cost_inj
mpc.ne_storage = [
											3	 		0.0	 0.0	 0.0	 			200.0	 			40.0	 						50.0	 					0.9	 									0.9	 								500	 			-100.0	140.0	   0.1	 0.0	 0.0	 		0.0	  	1 			35 		  	0 					0 					1000000									0 													0 													1e-4							0     0;
											3	 		0.0	 0.0	 0.0	 			100.0	 			20.0	 						25.0	 					0.9	 									0.9	 								250	 			-50.0	  70.0	   0.1	 0.0	 0.0	 		0.0	  	1 			17.5 			0 					0 					500000									0 													0 													1e-4							0     0;
											3	 		0.0	 0.0	 0.0	 			50.0	 			10.0	 						12.5	 					0.9	 									0.9	 								125	 			-25.0	  30.0	   0.1	 0.0	 0.0	 		0.0	  	1 			8.75 			0		 			  0 					250000									0 													0 													1e-4							0     0;
											3	 		0.0	 0.0	 0.0	 			25.0	 			5.0	 						  6.25	 					0.9	 									0.9	 								62.5	 		-12.5	  15.0	   0.1	 0.0	 0.0	 		0.0	  	1 			4.375 		0		 			  0 					125000									0 													0 													1e-4							0     0;
											3	 		0.0	 0.0	 0.0	 			12.5	 			2.5	 						  3.125	 					0.9	 									0.9	 								31.25	 		-6.25	  7.5 	   0.1	 0.0	 0.0	 		0.0	  	1 			2.1875 		0		 			  0 					62500						  			0 													0 													1e-4							0     0;
											3	 		0.0	 0.0	 0.0	 			6.25	 			1.25	 					  1.5625	 				0.9	 									0.9	 								15.625	 	3.125	  3.75 	   0.1	 0.0	 0.0	 		0.0	  	1 			1.09375 		0		 			  0 					31250							  		0 													0 													1e-4							0     0;
];
%3	 		0.0	 0.0	 0.0	 			800.0	 			160.0	 						200.0	 					0.9	 									0.9	 							 2000	 			-400.0	560.0	 	 0.1	 0.0	 0.0	 		0.0	  	1 			140 		  0 					0 					4000000									0 													0 													1e-4							0     0;
%3	 		0.0	 0.0	 0.0	 			400.0	 			80.0	 						100.0	 					0.9	 									0.9	 							 1000	 			-200.0	280.0	   0.1	 0.0	 0.0	 		0.0	  	1 			70 		  	0 					0 					2000000									0 													0 													1e-4							0     0;
