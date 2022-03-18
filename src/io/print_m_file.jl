function topology_df(rt_ex, relax, _ac)
	nodes = DataFrame(XLSX.readtable(rt_ex*"input.xlsx", "node_generation")...)
	edges = DataFrame(XLSX.readtable(rt_ex*"input.xlsx", "connections_acdc")...)
	ac_cable_df = DataFrame(XLSX.readtable(rt_ex*"input.xlsx", "CABLES_AC_SET_UP")...)
	dc_cable_df = DataFrame(XLSX.readtable(rt_ex*"input.xlsx", "CABLES_DC_SET_UP")...)
	rem_df = DataFrame(XLSX.readtable(rt_ex*"input.xlsx", "REMAINDER")...)
	ppf_mainACDCStorage2mfile(rem_df,ac_cable_df,dc_cable_df,rt_ex, relax, _ac)
	file = rt_ex*"topology.m"
	data = PowerModels.parse_file(file)
	data,ics_ac=filter_AClines(data,edges,nodes)
	data,ics_dc=filter_DClines(data,edges,nodes)
	return data, ics_ac, ics_dc, nodes
end

function AC_cable_options(data,candidate_ics_ac,ics_ac,pu)
    z_base_ac=(data["bus"]["1"]["base_kv"])^2/pu
    data=additional_candidatesICS_AC(data,candidate_ics_ac,ics_ac)#adds additional candidates
    for (i,bac) in data["ne_branch"]
    data["ne_branch"][i]=candidateIC_cost_impedance_AC(bac,z_base_ac,pu);end
    data["ne_branch"]=unique_candidateIC_AC(data["ne_branch"])#keep only unique candidates
	temp_cables_ne=Dict{String,Any}()
	temp_cables=Dict{String,Any}()
	cable_pu_costs=Dict{String,Any}()
	cable_pu_r=Dict{String,Any}()
	cable_pu_x=Dict{String,Any}()
	for (i,acb) in enumerate(sort(OrderedCollections.OrderedDict(data["ne_branch"]), by=x->parse(Int64,x)))
		last(acb)["source_id"][2]=i
		push!(temp_cables_ne,string(i)=>last(acb))
		for (j,acb_con) in enumerate(sort(OrderedCollections.OrderedDict(data["branch"]), by=x->parse(Int64,x)))
			if (last(acb)["f_bus"]==last(acb_con)["f_bus"] && last(acb)["t_bus"]==last(acb_con)["t_bus"])
				last(acb_con)["source_id"][2]=j
				push!(temp_cables,string(j)=>last(acb_con))
				if (haskey(cable_pu_costs,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])))
					push!(cable_pu_costs[string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])],last(acb)["construction_cost"]/(last(acb)["rate_a"]))
					push!(cable_pu_r[string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])],last(acb)["br_r"])
					push!(cable_pu_x[string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])],last(acb)["br_x"])
				else
					push!(cable_pu_costs,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])=>[last(acb)["construction_cost"]/(last(acb)["rate_a"])])
					push!(cable_pu_r,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])=>[last(acb)["br_r"]])
					push!(cable_pu_x,string(last(acb_con)["f_bus"])*"_"*string(last(acb_con)["t_bus"])=>[last(acb)["br_x"]])
				end
				break
			end
		end
	end
	println(cable_pu_costs)
	data["ne_branch"]=temp_cables_ne
	temp_cables2=Dict{String,Any}()
	for (i,acb) in enumerate(sort(OrderedCollections.OrderedDict(temp_cables), by=x->parse(Int64,x)))
		trms=length(cable_pu_costs[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])
		trms_total=sum(b for b in cable_pu_costs[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])

		trms_R=length(cable_pu_r[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])
		trms_total_R=sum(b for b in cable_pu_r[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])

		trms_X=length(cable_pu_x[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])
		trms_total_X=sum(b for b in cable_pu_x[string(last(acb)["f_bus"])*"_"*string(last(acb)["t_bus"])])

		last(acb)["br_r"]=(trms_total_R/trms_R)
		last(acb)["br_x"]=(trms_total_X/trms_X)
		last(acb)["cost"]=(trms_total/trms)
		last(acb)["source_id"][2]=i
		push!(temp_cables2,string(i)=>last(acb))
	end
	data["branch"]=temp_cables2
    return data
end

function DC_cable_options(data,candidate_ics_dc,ics_dc,pu)
    z_base_dc=(data["busdc"]["1"]["basekVdc"])^2/pu
    data=additional_candidatesICS_DC(data,candidate_ics_dc,ics_dc)#adds additional candidates
    for (i,bdc) in data["branchdc_ne"]
    data["branchdc_ne"][i]=candidateIC_cost_impedance_DC(bdc,z_base_dc);end
    data["branchdc_ne"]=unique_candidateIC_DC(data["branchdc_ne"])#keep only unique candidates
	temp_cables_ne=Dict{String,Any}()
	temp_cables=Dict{String,Any}()
	cable_pu_costs=Dict{String,Any}()
	cable_pu_r=Dict{String,Any}()
	for (i,acb) in enumerate(sort(OrderedCollections.OrderedDict(data["branchdc_ne"]), by=x->parse(Int64,x)))
		last(acb)["source_id"][2]=i
		push!(temp_cables_ne,string(i)=>last(acb))
		for (j,acb_con) in enumerate(sort(OrderedCollections.OrderedDict(data["branchdc"]), by=x->parse(Int64,x)))
			if (last(acb)["fbusdc"]==last(acb_con)["fbusdc"] && last(acb)["tbusdc"]==last(acb_con)["tbusdc"])
				last(acb_con)["source_id"][2]=j
				push!(temp_cables,string(j)=>last(acb_con))
				break
			end
		end
		if (haskey(cable_pu_costs,string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])))
			push!(cable_pu_costs[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])],last(acb)["cost"]/(last(acb)["rateA"]/pu))
			push!(cable_pu_r[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])],last(acb)["r"])
		else
			push!(cable_pu_costs,string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])=>[last(acb)["cost"]/(last(acb)["rateA"]/pu)])
			push!(cable_pu_r,string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])=>[last(acb)["r"]])
		end
	end
	data["branchdc_ne"]=temp_cables_ne
	temp_cables2=Dict{String,Any}()
	for (i,acb) in enumerate(sort(OrderedCollections.OrderedDict(temp_cables), by=x->parse(Int64,x)))
		trms=length(cable_pu_costs[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		trms_total=sum(b for b in cable_pu_costs[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		#trms_max=minimum(b for b in cable_pu_costs[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		trms_R=length(cable_pu_r[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		trms_total_R=sum(r for r in cable_pu_r[string(last(acb)["fbusdc"])*"_"*string(last(acb)["tbusdc"])])
		last(acb)["cost"]=trms_total/trms
		#last(acb)["cost"]=trms_max
		last(acb)["r"]=trms_total_R/trms_R
		last(acb)["source_id"][2]=i
		push!(temp_cables2,string(i)=>last(acb))
	end
	data["branchdc"]=temp_cables2
    return data
end

function filter_DClines(data,edges,nodes)
    #size and length
    ics_dc=Tuple{Int64,Int64}[]
    dcc=filter(x->!ismissing(x),edges["DC_mva"])
    for (k, s) in enumerate(dcc)
        from_xy=utm_gps2xy((nodes["lat"][edges["DC_from"][k]],nodes["long"][edges["DC_from"][k]]))
        to_xy=utm_gps2xy((nodes["lat"][edges["DC_to"][k]],nodes["long"][edges["DC_to"][k]]))
        push!(ics_dc,(s,round(Int64,Geodesy.euclidean_distance(from_xy, to_xy, 31, true, Geodesy.wgs84)/1000*1.25)))
    end

    #filter dc connections
    dccbles2keep_ne=Dict[]
	dccbles2keep=Dict[]
    for (r,s) in enumerate(dcc)
        for (k,b) in data["branchdc_ne"]
            if (edges["DC_from"][r]==b["fbusdc"] && edges["DC_to"][r]==b["tbusdc"])
                push!(dccbles2keep_ne,deepcopy(b))
                break;
            end
        end
		for (k,b) in data["branchdc"]
            if (edges["DC_from"][r]==b["fbusdc"] && edges["DC_to"][r]==b["tbusdc"])
                push!(dccbles2keep,deepcopy(b))
                break;
            end
        end
    end
    data["branchdc_ne"]=Dict{String,Any}()
    for (k,c) in enumerate(dccbles2keep_ne)
        c["source_id"][2]=k
        push!(data["branchdc_ne"],string(k)=>c)
    end
	data["branchdc"]=Dict{String,Any}()
    for (k,c) in enumerate(dccbles2keep)
        c["source_id"][2]=k
        push!(data["branchdc"],string(k)=>c)
    end
    return data, ics_dc
end

function filter_AClines(data,edges,nodes)
    #size and length
    ics_ac=Tuple{Int64,Int64}[]
    acc=filter(x->!ismissing(x),edges["AC_mva"])
    for (k, s) in enumerate(acc)
        from_xy=utm_gps2xy((nodes["lat"][edges["AC_from"][k]],nodes["long"][edges["AC_from"][k]]))
        to_xy=utm_gps2xy((nodes["lat"][edges["AC_to"][k]],nodes["long"][edges["AC_to"][k]]))
        push!(ics_ac,(s,round(Int64,Geodesy.euclidean_distance(from_xy, to_xy, 31, true, Geodesy.wgs84)/1000*1.25)))
    end

    #filter ac connections
    accbles2keep_ne=Dict[];accbles2keep=Dict[]
    acc=filter!(x->!ismissing(x),edges["AC_mva"])
    for (r,s) in enumerate(acc)
        for (k,b) in data["ne_branch"]
            if (edges["AC_from"][r]==b["f_bus"] && edges["AC_to"][r]==b["t_bus"])
                push!(accbles2keep_ne,deepcopy(b))
                break;
            end
        end
    end
    data["ne_branch"]=Dict{String,Any}()
    for (k,c) in enumerate(accbles2keep_ne)
        c["source_id"][2]=k
        push!(data["ne_branch"],string(k)=>c)
    end

	for (r,s) in enumerate(acc)
        for (k,b) in data["branch"]
            if (edges["AC_from"][r]==b["f_bus"] && edges["AC_to"][r]==b["t_bus"])
                push!(accbles2keep,deepcopy(b))
                break;
            end
        end
    end
    data["branch"]=Dict{String,Any}()
    for (k,c) in enumerate(accbles2keep)
        c["source_id"][2]=k
        push!(data["branch"],string(k)=>c)
    end
    return data, ics_ac
end

#utm coordinates from gps
function utm_gps2xy(lla,north_south::Bool=true,zone_utm::Int64=31)
    utm_desired = Geodesy.UTMfromLLA(zone_utm, north_south, Geodesy.wgs84)#sets UTM zone
    utm = utm_desired(Geodesy.LLA(first(lla),last(lla)))#coverts to cartesian
    return utm
end
#=
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
=#
function ppf_mainACDCStorage2mfile(r_df,ac_df, dc_df, cdir,relax, _ac)
	base_mva=100
	matfile = open(cdir*"topology.m","w")#open the .mat file
	ppf_header(matfile,base_mva)#print top function data
	ppf_Buss(matfile,r_df)#prints the bus data
	ppf_Gens(matfile,r_df)#prints all generator (OWPP) data
	ppf_acbranches(matfile,ac_df,relax, _ac)
	ppf_costs(matfile,r_df)
	ppf_BussDC(matfile,r_df)#prints the bus data
	ppf_BussDC_ne_blank(matfile)#prints the bus data
	#ppf_dcbranches_blank(matfile)
	ppf_dcbranches(matfile,dc_df, relax)
	ppf_convs(matfile,r_df)
	ppf_dc_blank_conv_ne(matfile)
	ppf_storage(matfile,r_df)
	close(matfile)#close the .mat file
end

function ppf_storage(mf,r_df)
	println(mf, " % hours
	 mpc.time_elapsed = 1.0

	 %% storage data
	 % storage_bus   ps   qs energy energy_rating charge_rating discharge_rating charge_efficiency discharge_efficiency thermal_rating    qmin   qmax    r    x p_loss q_loss status
	 mpc.storage = [")
	 	for r in r_df[:storage]
 	   		println(mf,r)
    	end
	println(mf,"];")



	println(mf, "%% cost 36.6
	 %% storage additional data
	 %column_names% max_energy_absorption stationary_energy_inflow stationary_energy_outflow self_discharge_rate cost
	 mpc.storage_extra = [")
		 for r in r_df[:storage_extra]
			 println(mf,r)
		 end
	 println(mf,"];")

	 println(mf, "%% storage data
	 %column_names%   storage_bus ps 	qs 	energy  energy_rating charge_rating  discharge_rating  charge_efficiency  discharge_efficiency  thermal_rating  qmin  	qmax  	r  		x  p_loss  	q_loss  status eq_cost inst_cost co2_cost 	max_energy_absorption 	stationary_energy_inflow 	stationary_energy_outflow 	self_discharge_rate	 	cost_abs 	cost_inj on_off
	 mpc.ne_storage = [
	 1	 0.0 0.0 0.0 6400.0	1280.0 1600.0 0.9 0.9 16000	-3200.0	4480.0 0.1 0.0 0.0	0.0	 0 	100000 	100000	100000 	32000000 0 	0 1e-4 0 0 1;
	 											 											];")
end

function ppf_BussDC_ne_blank(mf)
	println(mf,"%% All Candidate equipment has equivalent existing infrastructure field as well ie  mpc.branchdc_ne <->  mpc.branchdc similar to AC
	%% Candidate DC buses here - refer to MATACDC 1.0 User’s Manual for description of fields
	%% candidate dc bus data
	%column_names%   busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin  Cdc
	mpc.busdc_ne = [
	5              4       0       1       300         1.1     0.9     0;
	];")
end

function ppf_convs(mf,r_df)
	println(mf,"
	%% existing converters
	%column_names%   busdc_i busac_i type_dc type_ac P_g   Q_g  islcc  	Vtar 		rtf  xtf  transformer tm   bf  filter    rc      xc  reactor   basekVac Vmmax   Vmmin   Imax     status   LossA  LossB  LossCrec LossCinv   droop   Pdcset     Vdcset  dVdcset Pacmax Pacmin Qacmax   Qacmin cost
	mpc.convdc = [")

	for r in r_df[:conv_dc]
		println(mf,r)
	end
	println(mf,"];")
end

function ppf_dcbranches_blank(mf)
	println(mf,"%% existing dc branches
	%column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB rateC cost status
	mpc.branchdc = [
	1	2	0.0	0.00	 0.00  0	 		0	 0 68.75	 0;
	];")
end

function ppf_dcbranches(mf,dc_df, relax)
	clms=filter!(x->x!="Row",names(dc_df))
	cables=[]
	for c in clms
		cables=vcat(cables,dc_df[Symbol(c)])
	end
	cables=filter!(x->!ismissing(x),cables)
	println(mf,"%% existing dc branches
	%column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB rateC cost status
	mpc.branchdc = [")
	for r in cables
		rv=split(r," ")
		filter!(x->x!="",rv)
		if (relax)
			println(mf,rv[1]*"	"*rv[2]*"	"*"0.0	0.00	 0.00  0	 		0	 0 0.125	 1;")
		else
			println(mf,rv[1]*"	"*rv[2]*"	"*"0.0	0.00	 0.00  0	 		0	 0 0.125	 0;")
		end
	end
	println(mf,"];")
	println(mf,"%% Candidate Branches here - refer to MATACDC 1.0 Users Manual for description of fields
	 %% candidate branches
	 %column_names%   fbusdc  tbusdc  r      l        c   rateA   rateB   rateC status cost
	 mpc.branchdc_ne = [")

	for r in cables
		println(mf,r)
	end
	println(mf,"];")
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
	];")
	ppf_dc_blank_conv_ne(mf)
end

function ppf_dc_blank_conv_ne(mf)
	println(mf,"
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

function ppf_acbranches(mf,ac_df, relax, _AC)
	clms=filter!(x->x!="Row",names(ac_df))
	cables=[]
	for c in clms
		cables=vcat(cables,ac_df[Symbol(c)])
	end
	cables=filter!(x->!ismissing(x),cables)
	println(mf,"
	%% branch data
	%	fbus	tbus	r	       x	    b	  rateA	rateB	rateC	ratio	angle	status	angmin	angmax
	mpc.branch = [")
	for r in cables
		rv=split(r," ")
		filter!(x->x!="",rv)
		if (relax)
			println(mf,rv[1]*"	"*rv[2]*"	"*"0.0040   0.0400   0.00   0  0  0  0  0  "*_AC*" -60  60;")
		else
			println(mf,rv[1]*"	"*rv[2]*"	"*"0.0040   0.0400   0.00   0  0  0  0  0  0 -60  60;")
		end
	end
	println(mf,"];")

	println(mf,"%%Additional generator fields, type=0(wind generator), type=1(onshore market ie load and generator created)
	%column_names% 		cost
	mpc.branch_cost = [")

	for r in cables
		println(mf,"1.14;")
	end
	println(mf,"];")

	println(mf,"
	%candidate branch data
	%column_names%	f_bus	t_bus	br_r	br_x	br_b	rate_a	rate_b	rate_c	tap	shift	br_status	angmin	angmax	mva	construction_cost
	mpc.ne_branch = [")
	#=for r in cables
		println(mf,r)
	end=#
	for r in cables
		rv=split(r," ")
		filter!(x->x!="",rv)
		println(rv)
		if (relax)
			println(mf,rv[1]*"	"*rv[2]*"	"*rv[3]*"	"*rv[4]*"	"*rv[5]*"	"*rv[6]*"	"*rv[7]*"	"*rv[8]*"	"*rv[9]*"	"*rv[10]*"	"*"0"*"	"*rv[12]*"	"*rv[13]*"	"*rv[14]*"	"*rv[15])
			#println(rv[1]*"	"*rv[2]*"	"*rv[3]*"	"*rv[4]*"	"*rv[5]*"	"*rv[6]*"	"*rv[7]*"	"*rv[8]*"	"*rv[9]*"	"*rv[10]*"	"*"0"*"	"*rv[12]*"	"*rv[13]*"	"*rv[14]*"	"*rv[15])
		else
			println(mf,rv[1]*"	"*rv[2]*"	"*rv[3]*"	"*rv[4]*"	"*rv[5]*"	"*rv[6]*"	"*rv[7]*"	"*rv[8]*"	"*rv[9]*"	"*rv[10]*"	"*_AC*"	"*rv[12]*"	"*rv[13]*"	"*rv[14]*"	"*rv[15])
			#println(rv[1]*"	"*rv[2]*"	"*rv[3]*"	"*rv[4]*"	"*rv[5]*"	"*rv[6]*"	"*rv[7]*"	"*rv[8]*"	"*rv[9]*"	"*rv[10]*"	"*"0"*"	"*rv[12]*"	"*rv[13]*"	"*rv[14]*"	"*rv[15])
		end
	end
	println(mf,"];")
end

function ppf_Gens(mf,r_df)
	println(mf,"%% generator data
	%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
	mpc.gen = [")

	load_gen=filter!(x->!ismissing(x),vcat(r_df[:load],r_df[:gen]))
	#load_gen_cost=filter!(x->!ismissing(x),vcat(r_df[:load_cost],r_df[:gen_cost]))
	for r in load_gen
		println(mf,r)
	end
	println(mf,"];")

	println(mf,"%%Additional generator fields, type=0(wind generator), type=1(onshore market ie load and generator created)
	%column_names% 		type invest
	mpc.gen_type = [")

	for r in r_df[:gen_type]
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
