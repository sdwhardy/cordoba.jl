
function ppf_mainAC2mfile(r_df,ac_df)
	base_mva=100
	matfile = open("./test/binaries/test_case.m","w")#open the .mat file
	ppf_header(matfile,base_mva)#print top function data
	ppf_Buss(matfile,r_df)#prints the bus data
	ppf_Gens(matfile,r_df)#prints all generator (OWPP) data
	ppf_acbranches(matfile,ac_df)
	ppf_costs(matfile,r_df)
	ppf_BussDC(matfile,r_df)#prints the bus data
	ppf_dc_blanks(matfile)#prints the bus data
	#ppf_dcBuss(matfile,r_df)#prints the bus data
	#ppf_dcbranches(matfile,ac_df)
	close(matfile)#close the .mat file
end

function ppf_mainACDC2mfile(r_df,ac_df)
	base_mva=100
	matfile = open("./test/binaries/test_case.m","w")#open the .mat file
	ppf_header(matfile,base_mva)#print top function data
	ppf_Buss(matfile,r_df)#prints the bus data
	ppf_Gens(matfile,r_df)#prints all generator (OWPP) data
	ppf_acbranches(matfile,ac_df)
	ppf_costs(matfile,r_df)
	ppf_BussDC(matfile,r_df)#prints the bus data
	#ppf_dc_blanks(matfile)#prints the bus data
	#ppf_dcBuss(matfile,r_df)#prints the bus data
	#ppf_dcbranches(matfile,ac_df)
	close(matfile)#close the .mat file
end


function ppf_dc_blanks(mf)
	println(mf,"%% candidate dc bus data
	%column_names%   busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin  Cdc
	mpc.busdc_ne = [
		1000              3       0       1       300         1.1     0.9     0;
	];

	%% existing dc branches
	%column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB   rateC status
	mpc.branchdc = [
	1	2	0.0	0.00	 0.00  0	 0	 0	 0.0;
	];

	 %% candidate branches
	 %column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB   rateC status cost
	 mpc.branchdc_ne = [
	 1 2 0 0 0 0 0 0 0 0;
	 ];

	 %% existing converters
	 %column_names%   busdc_i busac_i type_dc type_ac P_g   Q_g  islcc  Vtar rtf    xtf  transformer tm   bf        filter    rc      xc    reactor   basekVac Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  droop       Pdcset    Vdcset  dVdcset Pacmax Pacmin Qacmax   Qacmin
	 mpc.convdc = [
	                     1       1       2       3    400000    0   	0   1.0  0.001  0.1       0 	 1   0.08 	  0      0.001   0.09      0  	      220    1.1     0.9    100000      1       0     0        0       0      0.0050    -52.7       1.0079     0     4000  -4000   4000    -4000;
	];

	%trans, filter, reactor, vmin vmax same as conv
	%% candidate converters
	%column_names%   busdc_i busac_i type_dc type_ac P_g   Q_g  islcc  Vtar    rtf   xtf  transformer tm   bf 	filter    rc     xc  reactor   basekVac Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  droop     Pdcset    Vdcset  dVdcset Pacmax Pacmin Qacmax Qacmin cost
	mpc.convdc_ne = [
	                1000       3      1       1    400000    0   	0     1.0   0.001  0.1       0 			 1 	0.08 		0 		0.001   0.09 		0  				220    1.1     0.9    100000      1     	 0     0        0       0      0.0050    -52.7     1.0079     0  		1  -1   1    -1   1000000;
	];")
end

function ppf_BussDC(mf,r_df)
	println(mf,"%% existing dc bus data
	%column_names%   busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin Cdc
	mpc.busdc = [")
	for r in r_df[:bus_dc]
		println(mf,r)
	end
	println(mf,"];")
end

function ppf_acbranches(mf,ac_df)
	println(mf,"
	%% branch data
	%	fbus	tbus	r	       x	    b	  rateA	rateB	rateC	ratio	angle	status	angmin	angmax
	mpc.branch = [
	];

	%candidate branch data
	%column_names%	f_bus	t_bus	br_r	br_x	br_b	rate_a	rate_b	rate_c	tap	shift	br_status	angmin	angmax	mva	construction_cost
	mpc.ne_branch = [")

	clms=filter!(x->x!="Row",names(ac_df))
	cables=[]
	for c in clms
		cables=vcat(cables,ac_df[Symbol(c)])
	end
	cables=filter!(x->!ismissing(x),cables)
	for r in cables
		println(mf,r)
	end
	println(mf,"];")
end

function ppf_Gens(mf,r_df)
	println(mf,"%% generator data
	%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
	mpc.gen = [")

	load_gen=filter!(x->!ismissing(x),vcat(r_df[:load],r_df[:gen]))
	load_gen_cost=filter!(x->!ismissing(x),vcat(r_df[:load_cost],r_df[:gen_cost]))
	for r in load_gen
		println(mf,r)
	end
	println(mf,"];")
end

function ppf_costs(mf,r_df)
	println(mf,"%%-----  OPF Data  -----%%
	%% generator cost data
	%	1	startup	shutdown	n	x1	y1	...	xn	yn
	%	2	startup	shutdown	n	c(n-1)	...	c0
	mpc.gencost = [")


	load_gen_cost=filter!(x->!ismissing(x),vcat(r_df[:load_cost],r_df[:gen_cost]))
	for r in load_gen_cost
		println(mf,r)
	end
	println(mf,"];")
end

function ppf_Buss(mf,r_df)
	println(mf,"%% bus data
	%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
	mpc.bus = [")
	for r in r_df[:bus]
		println(mf,r)
	end
	println(mf,"];")
end
function ppf_header(mf,base_mva)
	println(mf, "%TNEP optimizor input file")
	println(mf, "")
	println(mf, "function mpc = owpp_tnep_map")
	println(mf, "mpc.version = '2';")
	print(mf, "mpc.baseMVA = ")
	print(mf, base_mva)
	println(mf, ";")
	println(mf, "")
end
