
#load Time series data
function load_time_series(rt_ex, argz)
    #scenario_data=FileIO.load(rt_ex*"time_series_k"*string(argz["k"])*".jld2")
    scenario_data=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\time_series_k"*string(argz["k"])*".jld2")
    #keep only specified scenarios
    d_keys=keys(scenario_data);for k in d_keys;if !(issubset([string(k)],argz["scenario_names"]));delete!(scenario_data,k);else;y_keys=keys(scenario_data[k]);for y in y_keys;if !(issubset([string(y)],argz["scenario_years"]));delete!(scenario_data[k],y);end; end;end;end
    if (haskey(argz, "test") && argz["test"]==true)
        for k0 in d_keys; for k1 in keys(scenario_data[k0]); scenario_data[k0][k1]=scenario_data[k0][k1][1:2,:];end;end
    end
    ##################### Find minimum length scenario and Make all scenarios the same length
    ls=[];for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
    push!(ls,length(scenario_data[_sc][_yr].time_stamp))
    end;end;ls=minimum(ls)

    for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
    scenario_data[_yr][_sc]=scenario_data[_yr][_sc][1:ls,:]
    end;end
    return scenario_data, ls
end


############################################################################
#keeps k 24 clusters and formats RES time series 
#Test input:
#=
s = Dict(
"rt_ex"=>pwd()*"\\test\\data\\input\\UK_BE_DE_DK\\",#folder path
"scenario_data_file"=>"C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBEDEDK.jld2",
################# temperal parameters #################
"test"=>true,#if true smallest (2 hour) problem variation is built for testing
"scenario_planning_horizon"=>30,
"scenario_names"=>["NT","DE","GA"],#["NT","DE","GA"]
"k"=>6,#number of representative days modelled (24 hours per day)//#best for maintaining mean/max is k=6 2014, 2015
"res_years"=>["2014","2015"],#Options: ["2012","2013","2014","2015","2016"]//#best for maintaining mean/max is k=6 2014, 2015
"scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
"dr"=>0.04,#discount rate
"yearly_investment"=>1000000,
################ electrical parameters ################
"AC"=>"1",#0=false, 1=true
"owpp_mva"=>[2000,4000,4000],#mva of wf in MVA
"conv_lim_onshore"=>3000,#Max Converter size in MVA
"conv_lim_offshore"=>4000,#Max Converter size in MVA
"strg_lim_offshore"=>0.2,
"strg_lim_onshore"=>10,
"candidate_ics_ac"=>[1,4/5,3/5,2/5],#AC Candidate Cable sizes (fraction of full MVA)
"candidate_ics_dc"=>[1,4/5,3/5,2/5],#DC Candidate Cable sizes (fraction of full MVA)
################## optimization/solver setup options ###################
"output" => Dict("branch_flows" => false),
"eps"=>0.0001,#admm residual (100kW)
"beta"=>5.5,
"relax_problem" => false,
"conv_losses_mp" => true,
"process_data_internally" => false,
"corridor_limit" => true)=#

#load Time series data
#***#
function load_time_series_gentypes(s, scenario_data)
	#keeps data from specified scenarios only
    scenario_data["Generation"]["Scenarios"]=reduce_to_scenario_list(scenario_data["Generation"]["Scenarios"],s);
    scenario_data["Demand"]=reduce_to_scenario_list(scenario_data["Demand"],s);
    #Keep only specified markets
	countries=unique(vcat(s["map_gen_types"]["markets"][1],s["map_gen_types"]["markets"][2]));
    #scenario_data=reduce_to_market_list(scenario_data,countries);
    #keep only specified weather years
    scenario_data=reduce_to_weather_year_list(scenario_data,s);
    #keep only k specified days
    scenario_data["Generation"]["RES"], tss2keep = reduce_RES_to_k_days(scenario_data["Generation"]["RES"],s);
    scenario_data["Demand"]=reduce_DEMAND_to_k_days(scenario_data["Demand"],countries,tss2keep);
    #record number of hours
    country=first(keys(scenario_data["Generation"]["RES"]["Offshore Wind"]))
    year=first(keys(scenario_data["Generation"]["RES"]["Offshore Wind"][country]))
    s["hours_length"] = length(scenario_data["Generation"]["RES"]["Offshore Wind"][country][year].time_stamp)
	return scenario_data
end
#keys(scenario_data["Generation"]["RES"]["Solar PV"]["DE05"])
############## better copies
############################################################################################
#=keep only k specified RES days
#TEST INPUT 
_countries_in=["BE"]
_tss2keep_in=[Dates.DateTime("2020-01-01T00:00:00","yyyy-mm-ddTHH:MM:SS"),
 Dates.DateTime("2020-01-01T01:00:00","yyyy-mm-ddTHH:MM:SS")]
_d=Dict("2021"=>DataFrames.DataFrame(time_stamp=[Dates.DateTime("2021-01-01T00:00:00","yyyy-mm-ddTHH:MM:SS"),
Dates.DateTime("2021-01-01T01:00:00","yyyy-mm-ddTHH:MM:SS")],
                BE_MWh=[0.025,0.14855]))
_dict_in=Dict("Base"=>_d)
_dict_in["Base"]["2021"]
  
#TEST RUN
  _dict_out=reduce_DEMAND_to_k_days(_dict_in,_countries_in,_tss2keep_in)

#TEST OUTPUT
  print(_dict_out["Base"]["2021"])=#
#  Row │ time_stamp           BE_MWh  
#     │ DateTime             Float64 
#─────┼──────────────────────────────
#   1 │ 2020-01-01T00:00:00  0.025
#   2 │ 2020-01-01T01:00:00  0.14855
#***#
function reduce_DEMAND_to_k_days(_dict,_countries,_tss2keep);
    clmns2keep=[count for count in _countries];
    push!(clmns2keep,"time_stamp")
    for scenario in keys(_dict);
        #for yr in keys(_dict[scenario]);
            offset2020=Dates.Year(Dates.year(_dict[scenario][!,:time_stamp][1])-2020)
            _dict[scenario][!,:time_stamp]=_dict[scenario][!,:time_stamp].-offset2020;
            filter!(:time_stamp=>x->issubset([x],_tss2keep),_dict[scenario]);
            cs2delete=[name for name in names(_dict[scenario]) if !(issubset([name],clmns2keep))]
            DataFrames.select!(_dict[scenario],DataFrames.Not(cs2delete))
        #end;
    end
    return _dict
end
############################################################################################
#=keep only k specified RES days
#TEST INPUT 
_s_in=Dict("k"=>365,"test"=>true)
_d=Dict("2021"=>DataFrames.DataFrame(time_stamp=[Dates.DateTime("2021-01-01T00:00:00","yyyy-mm-ddTHH:MM:SS"),
Dates.DateTime("2021-01-01T01:00:00","yyyy-mm-ddTHH:MM:SS")],
                BE_MWh=[0.025,0.14855]))
_dict_in=Dict(
  "Offshore Wind" => Dict("BE"=>deepcopy(_d)),
  "Solar PV"      => Dict("BE"=>deepcopy(_d)),
  "Onshore Wind"  => Dict("BE"=>deepcopy(_d)))

  
#TEST RUN
  _dict_out,ts_out=reduce_to_k_days(_dict_in,_s_in)

#TEST OUTPUT
  print(_dict_out["Solar PV"]["BE"]["2021"])=#
#  Row │ time_stamp           BE_MWh  
#     │ DateTime             Float64 
#─────┼──────────────────────────────
#   1 │ 2020-01-01T00:00:00  0.025
#   2 │ 2020-01-01T01:00:00  0.14855
#_dict=scenario_data["Generation"]["RES"]
#_dict["Offshore Wind"]["FR01"]["2014"]
#country="DE05"
#year="2014"
#res_type="Onshore Wind"
#***#
function reduce_RES_to_k_days(_dict,_s)
    #keep only k specified days, shift year to 2020 (base year for all simmulations)
    ks=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\yearly_cluster_4EU.jld2")
    tss2keep=[]
    #RES sources
    for country in keys(_dict["Offshore Wind"])
        for year in unique(vcat(collect(keys(_dict["Offshore Wind"][country])),collect(keys(_dict["Onshore Wind"][country])),collect(keys(_dict["Solar PV"][country]))))
            if haskey(_dict["Offshore Wind"][country],year)
            ts2keep=_s["k"]<365 ? ks[year][_s["k"]] : _dict["Offshore Wind"][country][year].time_stamp;
            elseif haskey(_dict["Onshore Wind"][country],year)
                ts2keep=_s["k"]<365 ? ks[year][_s["k"]] : _dict["Onshore Wind"][country][year].time_stamp;
            else haskey(_dict["Solar PV"][country],year)
                ts2keep=_s["k"]<365 ? ks[year][_s["k"]] : _dict["Solar PV"][country][year].time_stamp;end
            if (haskey(_s, "test") && _s["test"]==true)
                ts2keep=ts2keep[1:2]
            end
            for res_type in keys(_dict)
                if haskey(_dict[res_type][country],year)
                filter!(:time_stamp=>x->issubset([x],ts2keep),_dict[res_type][country][year])
                #place common timestamp in 2020
			    offset2020=Dates.Year(2020-parse(Int64,year))
                _dict[res_type][country][year][!,:time_stamp]=_dict[res_type][country][year][!,:time_stamp].+offset2020
			    tss2keep=vcat(tss2keep,ts2keep.+offset2020);end
            end
        end
    end    
    return _dict,unique!(tss2keep)
end
############################################################################################
#keep only specified weather years
#TEST INPUT
#=_s_in=Dict("res_years"=>["2020"])
_dict_in=Dict("Generation" => Dict(
    "RES" => Dict(
        "Offshore Wind" => Dict("BE" => Dict("2020"=>"keep","2021"=>"don't keep"),"UK" => Dict("2020"=>"keep","2021"=>"don't keep")),
        "Solar PV" => Dict("BE" => Dict("2020"=>"keep","2021"=>"don't keep"),"UK" => Dict("2020"=>"keep","2021"=>"don't keep")),
        "Onshore Wind" => Dict("BE" => Dict("2020"=>"keep","2021"=>"don't keep"),"UK" => Dict("2020"=>"keep","2021"=>"don't keep")))))
#run 
_dict_out=reduce_to_weather_year_list(_dict_in,_s_in)

#TEST OUTPUT
print(_dict_out)
Dict("Generation" => Dict(
    "RES" => Dict(
        "Offshore Wind" => Dict("BE" => Dict("2020" => "keep"),"UK" => Dict("2020" => "keep")),
        "Solar PV" => Dict("BE" => Dict("2020" => "keep"),"UK" => Dict("2020" => "keep")),
        "Onshore Wind" => Dict("BE" => Dict("2020" => "keep"),"UK" => Dict("2020" => "keep")))))=#
#***#
function reduce_to_weather_year_list(_dict,_s)
	#keep only specified weather years
        for country in unique(vcat(collect(keys(_dict["Generation"]["RES"]["Offshore Wind"])), collect(keys(_dict["Generation"]["RES"]["Onshore Wind"])), collect(keys(_dict["Generation"]["RES"]["Solar PV"]))))
            years=[];
            if (haskey(_dict["Generation"]["RES"]["Offshore Wind"], country)); years=vcat(years, collect(keys(_dict["Generation"]["RES"]["Offshore Wind"][country])));end
            if (haskey(_dict["Generation"]["RES"]["Onshore Wind"], country)); years=vcat(years, collect(keys(_dict["Generation"]["RES"]["Onshore Wind"][country])));end
            if (haskey(_dict["Generation"]["RES"]["Solar PV"], country)); years=vcat(years, collect(keys(_dict["Generation"]["RES"]["Solar PV"][country])));end
            for year in unique(years)
                if !(issubset([string(year)],_s["res_years"]));
                    if (haskey(_dict["Generation"]["RES"]["Offshore Wind"], country) && haskey(_dict["Generation"]["RES"]["Offshore Wind"][country], year))
                    delete!(_dict["Generation"]["RES"]["Offshore Wind"][country],year);end
                    if (haskey(_dict["Generation"]["RES"]["Onshore Wind"], country) && haskey(_dict["Generation"]["RES"]["Onshore Wind"][country], year))
                    delete!(_dict["Generation"]["RES"]["Onshore Wind"][country],year);end
                    if (haskey(_dict["Generation"]["RES"]["Solar PV"], country) && haskey(_dict["Generation"]["RES"]["Solar PV"][country], year))
                    delete!(_dict["Generation"]["RES"]["Solar PV"][country],year);end
                end;
            end;
        end
        return _dict
    end
############################################################################################
#Keep only specified markets
#TEST INPUT
#=_countries_in=["UK","BE"]
_dict_in=Dict("Generation"=>Dict(
    "RES"=>Dict(
        "Offshore Wind"=>Dict("BE"=>"keep","UK"=>"keep","DE"=>"don't keep"),
        "Onshore Wind"=>Dict("BE"=>"keep","UK"=>"keep","DE"=>"don't keep"),
        "Solar PV"=>Dict("BE"=>"keep","UK"=>"keep","DE"=>"don't keep")),
    "Scenarios"=>Dict(
        "Base"=>Dict("2020"=>Dict("BE"=>"keep","UK"=>"keep","DE"=>"don't keep")),
        "NT"=>Dict("2020"=>Dict("BE"=>"keep","UK"=>"keep","DE"=>"don't keep")))))

#run         
 _dict_out=reduce_to_market_list(_dict_in,_countries_in)

#TEST OUTPUT
print(_dict_out)   
Dict{String,Dict{String,Dict{String,V} where V}}
("Generation" => Dict(
    "Scenarios" => Dict(
        "Base" => Dict("2020" => Dict("BE" => "keep","UK" => "keep")),
        "NT" => Dict("2020" => Dict("BE" => "keep","UK" => "keep"))),
    "RES" => Dict(
        "Offshore Wind" => Dict("BE" => "keep","UK" => "keep"),
        "Solar PV" => Dict("BE" => "keep","UK" => "keep"),
        "Onshore Wind" => Dict("BE" => "keep","UK" => "keep")))) 
=#

function reduce_to_market_list(_dict,countries)
    #Keep only specified markets
    d_keys=keys(_dict["Generation"]["Scenarios"]);
    for k in d_keys;
        c_keys=keys(_dict["Generation"]["Scenarios"][k]);
        for c in c_keys;
            if !(issubset([string(c)],countries));
                delete!(_dict["Generation"]["Scenarios"][k],c);
            end;
        end;
    end
    for key in keys(_dict["Generation"]["RES"]["Offshore Wind"])
        if !(issubset([string(key)],countries));
            delete!(_dict["Generation"]["RES"]["Offshore Wind"],key);
            delete!(_dict["Generation"]["RES"]["Onshore Wind"],key);
            delete!(_dict["Generation"]["RES"]["Solar PV"],key);
        end;
    end
    return _dict
end
############################################################################################
#keeps data from specified scenarios only
#s_in=Dict("scenario_names"=>["NT","DE"])
#_dict_in=Dict("Base"=>"in output","NT"=>"in output","DE"=>"in output","GA"=>"NOT in output")
#_dict_out=reduce_to_scenario_list(_dict_in,s_in)   
#Dict{String,String} with 3 entries:
#  "Base" => "in output"
#  "NT"   => "in output"
#  "DE"   => "in output" 
#***#
function reduce_to_scenario_list(_dict,s)
    namestokeep=vcat(s["scenario_names"],"Base")
    d_keys=keys(_dict);
    for k in d_keys;
        if !(issubset([string(k)],namestokeep));
            delete!(_dict,k);
        end;
    end
    return _dict
end
############################################################################
########################## Duplicate copy pre-changes end
############################################################################


#multi period problem setup
function multi_period_setup(ls,scenario_data,data, markets, infinite_grid, argz, s)
    #################### Multi-period input parameters #######################
    all_scenario_data,data,scenario, dim = multi_period_stoch_year_setup(ls,argz["scenario_years"],argz["scenario_names"],scenario_data,data);
    scenario["planning_horizon"] = argz["scenario_planning_horizon"]; # in years, to scale generation cost
    extradata,data =create_profile_sets_mesh(dim, data, all_scenario_data, markets, infinite_grid, [data["baseMVA"] for wf in argz["owpp_mva"]])
    extradata = create_profile_sets_rest(dim, extradata, data)
    #########################################################################
    #################### Scale cost data
    scale_cost_data_hourly!(extradata, scenario)
    extradata=costs_datas_wREZ(extradata, s, argz)

    ######################################################
    xtradata = Dict{String,Any}()
    xtradata["dim"] = Dict{String,Any}()
    xtradata["dim"] = dim
    ######################################################
    # Create data dictionary where time series data is included at the right place
    mn_data = _PMACDC.multinetwork_data(data, xtradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))
    return mn_data, extradata
end

#Organizes nw numbers per scenario-year
function multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data)
    scenario = Dict{String, Any}("hours" => ls,"years" => length(scenario_years), "sc_names" => Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    all_scenario_data=DataFrames.DataFrame()
    #set problem dimension
    dim = scenario["hours"] * length(scenario_years) * length(scenario_names)
    for _nm in scenario_names; push!(scenario["sc_names"], _nm=> Dict{String, Any}());for _yr in scenario_years; push!(scenario["sc_names"][_nm], _yr=> []);end;end

    for (s,(_sc, data_by_sc)) in enumerate(scenario_data);
        data["scenario"][string(s)] = Dict()
        data["scenario_prob"][string(s)] = 1/(length(scenario_names))
        for (t,(_yr, data_by_yr)) in enumerate(data_by_sc);
            all_scenario_data=vcat(all_scenario_data,scenario_data[_sc][_yr])
            start_idx=(s-1)*scenario["hours"]*length(scenario_years)
            start_idx=start_idx+(t-1)*scenario["hours"]
            for h in 1 : scenario["hours"]
                network = start_idx + h
                h2=h+(t-1)*scenario["hours"]
                data["scenario"][string(s)]["$h2"] = network
                push!(scenario["sc_names"][_sc][_yr],network)
            end
        end;
    end
    return all_scenario_data,data,scenario, dim
end
#[println(string(i)*": "*string(data["gen"][string(i)]["gen_bus"])) for i=1:length(data["gen"])]
#println(keys(extradata["gen"]["225"]["gen_bus"]))
#println(keys(s["xd"]["gen"]["101"]))
#println(keys(all_gens["onshore"]["BE"]["Solar PV"]))
#multi period problem setup
#data["convdc"]
#extradata["gen"]["1"]
#***#
function multi_period_setup_wgen_type(scenario_data,data, all_gens, s)
    #################### Multi-period input parameters #######################
    data, s = multi_period_stoch_year_setup_wgen_type(s,data);
    extradata, data = create_profile_sets_mesh_wgen_type(data, all_gens, scenario_data, s)
    extradata = create_profile_sets_rest_wgen_type(extradata, data, s)
    #########################################################################
    #################### Scale cost data
    scale_cost_data_hourly!(extradata, s["scenario"])
    s["xd"]=costs_datas_wREZ(extradata, s)

    ######################################################
    xtradata = Dict{String,Any}()
    xtradata["dim"] = Dict{String,Any}()
    xtradata["dim"] = s["dim"]
    ######################################################
    # Create data dictionary where time series data is included at the right place
    mn_data = _PMACDC.multinetwork_data(data, xtradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))
    for (g,gen) in mn_data["nw"]
		gen["gen"]=Dict(gen["gen"]);end
    return mn_data, s
end

#Organizes nw numbers per scenario-year
function multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data)
    scenario = Dict{String, Any}("hours" => ls,"years" => length(scenario_years), "sc_names" => Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    all_scenario_data=DataFrames.DataFrame()
    #set problem dimension
    dim = scenario["hours"] * length(scenario_years) * length(scenario_names)
    for _nm in scenario_names; push!(scenario["sc_names"], _nm=> Dict{String, Any}());for _yr in scenario_years; push!(scenario["sc_names"][_nm], _yr=> []);end;end

    for (s,(_sc, data_by_sc)) in enumerate(scenario_data);
        data["scenario"][string(s)] = Dict()
        data["scenario_prob"][string(s)] = 1/(length(scenario_names))
        for (t,(_yr, data_by_yr)) in enumerate(data_by_sc);
            all_scenario_data=vcat(all_scenario_data,scenario_data[_sc][_yr])
            start_idx=(s-1)*scenario["hours"]*length(scenario_years)
            start_idx=start_idx+(t-1)*scenario["hours"]
            for h in 1 : scenario["hours"]
                network = start_idx + h
                h2=h+(t-1)*scenario["hours"]
                data["scenario"][string(s)]["$h2"] = network
                push!(scenario["sc_names"][_sc][_yr],network)
            end
        end;
    end
    return all_scenario_data,data,scenario, dim
end

#Organizes nw numbers per scenario-year
#***#
function multi_period_stoch_year_setup_wgen_type(s,data)
	scenario_names=unique([sn[1:2] for sn in s["scenario_names"]])
    scenario = Dict{String, Any}("hours" => s["hours_length"],"years" => length(s["scenario_years"]), "sc_names" => Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    #set problem dimension

    
    for _nm in scenario_names; for _by in s["res_years"]; push!(scenario["sc_names"], string(_nm)*string(_by)=> Dict{String, Any}());for _yr in s["scenario_years"]; push!(scenario["sc_names"][string(_nm)*string(_by)], _yr=> []);end;end;end
    dim = scenario["hours"] * length(s["scenario_years"]) * length(scenario["sc_names"])
    for (s,(k_sc,_sc)) in enumerate(scenario["sc_names"]);
		data["scenario"][string(s)] = Dict()
        data["scenario_prob"][string(s)] = 1/length(scenario["sc_names"])
        for (t,(k_yr,_yr)) in enumerate(sort(OrderedCollections.OrderedDict(_sc), by=x->parse(Int64,x)));
            start_idx=(s-1)*scenario["hours"]*length(_sc)
            start_idx=start_idx+(t-1)*scenario["hours"]
            for h in 1 : scenario["hours"]
                network = start_idx + h
                h2=h+(t-1)*scenario["hours"]
                data["scenario"][string(s)]["$h2"] = network
                push!(scenario["sc_names"][string(k_sc)][k_yr],network)
            end
        end;
    end
    s["scenario"]=scenario
    s["scenario"]["planning_horizon"] = s["scenario_planning_horizon"]; # in years, to scale generation cost
    s["dim"]=dim
    return data,s
end

#=function multi_period_stoch_year_setup_wgen_type(ls,res_years,scenario_years,scenario_names,data)
	scenario_names=unique([sn[1:2] for sn in scenario_names])
    scenario = Dict{String, Any}("hours" => ls,"years" => length(scenario_years), "sc_names" => Dict{String, Any}())
    data["scenario"] = Dict{String, Any}()
    data["scenario_prob"] = Dict{String, Any}()
    #set problem dimension
    dim = scenario["hours"] * length(scenario_years) * length(res_years) * length(scenario_names)
    for _nm in scenario_names; for _by in res_years; push!(scenario["sc_names"], string(_nm)*string(_by)=> Dict{String, Any}());for _yr in scenario_years; push!(scenario["sc_names"][string(_nm)*string(_by)], _yr=> []);end;end;end

    for (s,(k_sc,_sc)) in enumerate(scenario["sc_names"]);
		data["scenario"][string(s)] = Dict()
        data["scenario_prob"][string(s)] = 1/length(scenario["sc_names"])
        for (t,(k_yr,_yr)) in enumerate(sort(OrderedCollections.OrderedDict(_sc), by=x->parse(Int64,x)));
            start_idx=(s-1)*scenario["hours"]*length(_sc)
            start_idx=start_idx+(t-1)*scenario["hours"]
            for h in 1 : scenario["hours"]
                network = start_idx + h
                h2=h+(t-1)*scenario["hours"]
                data["scenario"][string(s)]["$h2"] = network
                push!(scenario["sc_names"][string(k_sc)][k_yr],network)
            end
        end;
    end
    return data,scenario, dim
end=#

#loads generator cost and profile time series In multi-period simulation
#=function create_profile_sets_mesh(number_of_hours, data_orig, zs_data, zs, inf_grid, owpp_mva)
    pu=data_orig["baseMVA"]
    e2me=1000000/pu#into ME/PU
    extradata = Dict{String,Any}()
    data=Dict{String,Any}();data["gen"]=Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
    data["gen"]=sort!(OrderedCollections.OrderedDict(data_orig["gen"]), by=x->parse(Int64,x))
    for (g, gen) in data["gen"]
        extradata["gen"][g] = Dict{String,Any}()
        extradata["gen"][g]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][g]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end
    for d in 1:number_of_hours
        #Day ahead BE
        #source generator
        for (g, gen) in data["gen"]
            if (gen["type"]>0)#market generator onshore
                extradata["gen"][g]["pmax"][1, d] = inf_grid/pu
                extradata["gen"][g]["pmin"][1, d] = 0
                extradata["gen"][g]["cost"][d] = [(zs_data[!,"EUR_da"*zs[1][gen["gen_bus"]]][d])/e2me,0]
            else#wind gen
                extradata["gen"][g]["pmax"][1, d] = (zs_data[!,"Wnd_MWh"*zs[2][gen["gen_bus"]-length(zs[1])]][d])*owpp_mva[gen["gen_bus"]-length(zs[1])]/pu

                #extradata["gen"][g]["pmax"][1, d]
                #zs_data[!,"Wnd_MWh"*zs[2][gen["gen_bus"]-length(zs[1])]][d]
                #owpp_mva[gen["gen_bus"]-length(zs[1])]/pu

                extradata["gen"][g]["pmin"][1, d] = 0
                extradata["gen"][g]["cost"][d] = [0,0]
            end
        end
    end
    #add loads
    loads=Dict{String,Any}()
    num_of_gens=length(data["gen"])
    for (g, gen) in sort!(OrderedCollections.OrderedDict(extradata["gen"]), by=x->parse(Int64,x))
        if (data["gen"][g]["type"]>0)#market generator onshore
            load=deepcopy(data["gen"][g])
            load["index"]=num_of_gens+1
            load["source_id"][2]=num_of_gens+1
            load["pmin"]=deepcopy(load["pmax"])*-1
            load["pmax"]=0
            push!(loads,string(num_of_gens+1)=>deepcopy(load))
            num_of_gens=num_of_gens+1
        end
    end
    for (l, load) in loads
        extradata["gen"][l] = Dict{String,Any}()
        extradata["gen"][l]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][l]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][l]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end

    for d in 1:number_of_hours
        #Day ahead BE
        #source generator
        for (l, load) in loads
            if (load["type"]>0)#market generator onshore
                extradata["gen"][l]["pmax"][1, d] = 0
                extradata["gen"][l]["pmin"][1, d] = (inf_grid/pu)*-1
                extradata["gen"][l]["cost"][d] = [(zs_data[!,"EUR_da"*zs[1][load["gen_bus"]]][d])/e2me,0]
                push!(data_orig["gen"],l=>load)
            else#wind gen
            end
        end
    end

    #set ["type"]
    for (g, gen) in data_orig["gen"]
        gen["type"]=0
    end
    return extradata,data_orig
end=#
#extradata["gen"]["1"]["cost"]
#***#
function create_profile_sets_mesh_wgen_type(data_orig, all_gens, scenario_data, s)
	genz=[];wfz=[]
    pu=data_orig["baseMVA"]
    e2me=1000000/pu#into ME/PU
    extradata = Dict{String,Any}()
    data=Dict{String,Any}();data["gen"]=Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = s["dim"]
    extradata["gen"] = Dict{String,Any}()
	#demand_curve = Dict{String,Any}()

	for (c,country) in all_gens["onshore"];
		for (fuel,type) in country;
			for (num,g) in type;
				push!(data["gen"],num=>g)
				;end;end;end

	for (c,country) in all_gens["offshore"];
		for (num,g) in country;
		push!(data["gen"],num=>g);end;end

    data["gen"]=sort!(OrderedCollections.OrderedDict(data["gen"]), by=x->parse(Int64,x))
	for (g, gen) in data["gen"]
        extradata["gen"][string(g)] = Dict{String,Any}()
        extradata["gen"][string(g)]["pmax"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["gen"][string(g)]["pmin"] = Array{Float64,2}(undef, 1, s["dim"])
		extradata["gen"][string(g)]["cost"] = [Vector{Float64}() for i=1:s["dim"]]
		if (gen["type"]==0)
		extradata["gen"][string(g)]["invest"] = Array{Float64,2}(undef, 1, s["dim"]);
        extradata["gen"][string(g)]["wf_pmax"] = Array{Float64,2}(undef, 1, s["dim"]);end
    end
    

    for (k_sc,sc) in s["scenario"]["sc_names"]
		for (k_yr,yr) in sc
			k_yr_sd=k_yr=="2020" ? "2025" : k_yr;
			k_sc_sd=k_yr=="2020" ? "NT" : k_sc[1:2]
			for (h,d) in enumerate(yr)
		        #Day ahead BE
		        #onshore generators
				for (xy, country) in all_gens["onshore"]
		        	for (fuel, type) in country
						for (g, gen) in type
		            	#market generator onshore
						S_row=filter(:Generation_Type=>x->x==fuel, scenario_data["Generation"]["Scenarios"][string(k_sc_sd)*string(k_yr_sd)][xy])[!,:Capacity]
						if isempty(S_row)
							extradata["gen"][string(g)]["pmax"][1, d] = 0
							extradata["gen"][string(g)]["cost"][d] = [0,0]
							extradata["gen"][string(g)]["pmin"][1, d] = 0
		#					extradata["gen"][string(g)]["gen_status"][1, d] = 0
						else
							if issubset([fuel],keys(scenario_data["Generation"]["RES"]))
								CF=scenario_data["Generation"]["RES"][fuel][xy][k_sc[3:6]][!,Symbol(xy)][h]
								extradata["gen"][string(g)]["pmax"][1, d] = CF*S_row[1]/pu
							else
								extradata["gen"][string(g)]["pmax"][1, d] = S_row[1]/pu
							end
							extradata["gen"][string(g)]["pmin"][1, d] = 0
			                extradata["gen"][string(g)]["cost"][d] = [scenario_data["Generation"]["costs"][fuel]/e2me,0]
							#if !(haskey(demand_curve,xy));push!(demand_curve,xy=>DataFrames.DataFrame(:generation=>[],:fuel_cost=>[]));end
							#push!(demand_curve[xy],[S_row[1]/pu,scenario_data["Generation"]["costs"][fuel]/e2me])
							push!(genz,(parse(Int64,g),S_row[1]/pu))
						end
                    end;
                end;
            end
            
				#Wind power Plants
				for (j,(xy, country)) in enumerate(all_gens["offshore"])
		        	for (i,(g, gen)) in enumerate(country)
					#wind gen
						CF=scenario_data["Generation"]["RES"]["Offshore Wind"][xy][k_sc[3:6]][!,Symbol(xy)][h]
		                extradata["gen"][string(g)]["pmax"][1, d] = CF
						extradata["gen"][string(g)]["pmin"][1, d] = 0
		                extradata["gen"][string(g)]["cost"][d] = [0,0]
						extradata["gen"][string(g)]["invest"][1, d] = gen["invest"]
                        #extradata["gen"][string(g)]["wf_pmax"][1, d] = s["owpp_mva"][j]/pu
                        extradata["gen"][string(g)]["wf_pmax"][1, d] = gen["pmax"]*pu
						#push!(wfz,(parse(Int64,g),s["owpp_mva"][j]/pu))
                        push!(wfz,(parse(Int64,g),gen["pmax"]*pu))
				    end;
                end
			end
        end
    end
	#set ["type"]
    for (g, gen) in data["gen"]
		if (isapprox(sum(extradata["gen"][g]["pmax"]),0,atol=10e-3) && isapprox(sum(extradata["gen"][g]["pmin"]),0,atol=10e-3))
			gen["gen_status"]=0;
		else
		end
    end
	#=xy="DE"
	if !(haskey(demand_curve,xy));push!(demand_curve,xy=>DataFrames.DataFrame(:generation=>[],:fuel_cost=>[]));end
	push!(demand_curve[xy],[1,3])=#
    #adjust demand curve
	#=demand_curve_reduced=Dict()
	for (cuntree, df) in demand_curve
		#cuntree_total=sum(df[!,:generation])
		unique_costs=unique(df[!,:fuel_cost])
		top=maximum(unique_costs)
		bottom=minimum(unique_costs)
		step=(top-bottom)/9
		if !(haskey(demand_curve_reduced,cuntree));push!(demand_curve_reduced,cuntree=>DataFrames.DataFrame(:generation=>[],:fuel_cost=>[]));end
		#for fuel_cost in unique_costs
		for i = 1:1:10
			#capacity=sum(filter(:fuel_cost=>x->x==fuel_cost,df)[!,:generation])
			#push!(demand_curve_reduced[cuntree],[capacity/cuntree_total,fuel_cost])
            #println(string(cuntree)*"|"*string(capacity/cuntree_total)*"|"*string(fuel_cost))
            #push!(demand_curve_reduced[cuntree],[1/length(unique_costs),fuel_cost])
            push!(demand_curve_reduced[cuntree],[1/10,top-(i-1)*step])
			#println(string(cuntree)*"|"*string(1/10)*"|"*string(top-(i-1)*step))

		end
	end=#
	#add loads
	#all_gens["onshore"]["UK"]["Gas"]["8"]=g

	num_of_gens=length(data["gen"])
	if !(haskey(s["map_gen_types"],"loads"));push!(s["map_gen_types"],"loads"=>Dict());end
	loads=Dict{String,Any}()
	for _row in eachrow(filter(:type=>x->x==1,s["nodes"]))
		if !(haskey(loads,_row[:country]));push!(loads,_row[:country]=>Dict());end
		load=deepcopy(first(first(all_gens["onshore"][_row[:country]])[2])[2])
		#for row in eachrow(df)
			#reset load parameters
	    	load["type"]=1
			load["gen_status"]=1
	        load["index"]=num_of_gens+1
	        load["source_id"][2]=num_of_gens+1
	        load["pmin"]=-1
            load["pmax"]=-1
	        #load["pmax"]=0
			#load["cost"]=[deepcopy(row[:fuel_cost]),0]
            load["cost"]=[0.015,0]
	        push!(loads[_row[:country]],string(num_of_gens+1)=>deepcopy(load))
	        num_of_gens=num_of_gens+1

			#map load to country
			if !(haskey(s["map_gen_types"]["loads"],_row[:country]));push!(s["map_gen_types"]["loads"],_row[:country]=>[]);end
			push!(s["map_gen_types"]["loads"][_row[:country]],num_of_gens)
	end;#end
    
	for (cuntree, dic) in loads
    	for (l, load) in dic
	        extradata["gen"][string(l)] = Dict{String,Any}()
	        extradata["gen"][string(l)]["pmax"] = Array{Float64,2}(undef, 1, s["dim"])
	        extradata["gen"][string(l)]["pmin"] = Array{Float64,2}(undef, 1, s["dim"])
	        extradata["gen"][string(l)]["cost"] = [Vector{Float64}() for i=1:s["dim"]]
    end;end

	#onshore loads
    #=k_sc="GA2014";s["scenario"]["sc_names"][k_sc]
    k_yr="2040";yr=sc[k_yr]
    h=1;d=yr[1]
    cuntree="CH00";dic=loads[cuntree]
    l="728";load=dic[l]=#
	for (k_sc,sc) in s["scenario"]["sc_names"]
		for (k_yr,yr) in sc
			k_yr_sd=k_yr=="2020" ? "2025" : k_yr;
			k_sc_sd=k_yr=="2020" ? "NT" : k_sc[1:2]
			for (h,d) in enumerate(yr)
		        #loads
				for (cuntree, dic) in loads
			        for (l, load) in sort!(OrderedCollections.OrderedDict(dic), by=x->parse(Int64,x))
                        #NOTE: Left here ts is off - 2014 time stamp vs 2020 BLNK time stamp would work but...
                        #Just this final load to be done, however, the rest still needs to be confirmed! 
                        if (cuntree=="BLNK")
						    ts=scenario_data["Generation"]["RES"]["Onshore Wind"][first(keys(scenario_data["Generation"]["RES"]["Onshore Wind"]))][k_sc[3:6]][!,:time_stamp][h]
                        else
                            ts=scenario_data["Generation"]["RES"]["Onshore Wind"][cuntree][k_sc[3:6]][!,:time_stamp][h]
                        end
                        #println(k_sc_sd*k_yr_sd)
						S_row=filter(:time_stamp=>x->x==ts,scenario_data["Demand"][k_sc_sd*k_yr_sd])[!,Symbol(cuntree)]
		                extradata["gen"][string(l)]["pmax"][1, d] = load["pmax"]*S_row[1]/pu
		                extradata["gen"][string(l)]["pmin"][1, d] = load["pmin"]*S_row[1]/pu
                        #println(string(extradata["gen"][string(l)]["pmax"][1, d])*" "*string(extradata["gen"][string(l)]["pmin"][1, d]))
		                extradata["gen"][string(l)]["cost"][d] = load["cost"]
		                push!(data["gen"],string(l)=>load)
						push!(genz,(parse(Int64,l),S_row[1]/pu))
			        end
				end
    end;end;end
	unique!(x->first(x),genz)
	unique!(x->first(x),wfz)
	push!(s,"genz"=>genz)
    push!(s,"wfz"=>wfz)
	data_orig["gen"]=data["gen"]
    return extradata, data_orig
end

#=function create_profile_sets_mesh_wgen_type(number_of_hours, scenario, data_orig, all_gens, scenario_data, markets, s, argz, map_gen_types)
	genz=[];wfz=[]

    pu=data_orig["baseMVA"]
    e2me=1000000/pu#into ME/PU
    extradata = Dict{String,Any}()
    data=Dict{String,Any}();data["gen"]=Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
	demand_curve = Dict{String,Any}()

	for (c,country) in all_gens["onshore"];
		for (fuel,type) in country;
			for (num,g) in type;
				push!(data["gen"],num=>g)
				;end;end;end

	for (c,country) in all_gens["offshore"];
		for (num,g) in country;
		push!(data["gen"],num=>g);end;end

    data["gen"]=sort!(OrderedCollections.OrderedDict(data["gen"]), by=x->parse(Int64,x))
	for (g, gen) in data["gen"]
        extradata["gen"][string(g)] = Dict{String,Any}()
        extradata["gen"][string(g)]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][string(g)]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
		extradata["gen"][string(g)]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
		if (gen["type"]==0)
		extradata["gen"][string(g)]["invest"] = Array{Float64,2}(undef, 1, number_of_hours);
        extradata["gen"][string(g)]["wf_pmax"] = Array{Float64,2}(undef, 1, number_of_hours);end
    end

	for (k_sc,sc) in scenario["sc_names"]
		for (k_yr,yr) in sc
			#k_yr_sd=k_yr=="2020" ? "2025" : k_yr;
			k_sc_sd=k_yr=="2020" ? "Base" : k_sc[1:2]
			for (h,d) in enumerate(yr)
		        #Day ahead BE
		        #onshore generators
				for (xy, country) in all_gens["onshore"]
		        	for (fuel, type) in country
						for (g, gen) in type
		            	#market generator onshore
						S_row=filter(:Generation_Type=>x->x==fuel, scenario_data["Generation"]["Scenarios"][k_sc_sd][k_yr][xy])[!,:Capacity]
						if isempty(S_row)
							extradata["gen"][string(g)]["pmax"][1, d] = 0
							extradata["gen"][string(g)]["cost"][d] = [0,0]
							extradata["gen"][string(g)]["pmin"][1, d] = 0
		#					extradata["gen"][string(g)]["gen_status"][1, d] = 0
						else
							if issubset([fuel],keys(scenario_data["Generation"]["RES"]))
								CF=scenario_data["Generation"]["RES"][fuel][xy][k_sc[3:6]][!,Symbol(xy*"_MWh")][h]
								extradata["gen"][string(g)]["pmax"][1, d] = CF*S_row[1]/pu
							else
								extradata["gen"][string(g)]["pmax"][1, d] = S_row[1]/pu
							end
							extradata["gen"][string(g)]["pmin"][1, d] = 0
			                extradata["gen"][string(g)]["cost"][d] = [scenario_data["Generation"]["costs"][fuel]/e2me,0]
							if !(haskey(demand_curve,xy));push!(demand_curve,xy=>DataFrames.DataFrame(:generation=>[],:fuel_cost=>[]));end
							push!(demand_curve[xy],[S_row[1]/pu,scenario_data["Generation"]["costs"][fuel]/e2me])
							push!(genz,(parse(Int64,g),S_row[1]/pu))
						end
				end;end;end
				#Wind power Plants
				for (j,(xy, country)) in enumerate(all_gens["offshore"])
		        	for (i,(g, gen)) in enumerate(country)
					#wind gen
						CF=scenario_data["Generation"]["RES"]["Offshore Wind"][xy][k_sc[3:6]][!,Symbol(xy*"_MWh")][h]
		                extradata["gen"][string(g)]["pmax"][1, d] = CF
						extradata["gen"][string(g)]["pmin"][1, d] = 0
		                extradata["gen"][string(g)]["cost"][d] = [0,0]
						extradata["gen"][string(g)]["invest"][1, d] = gen["invest"]
                        extradata["gen"][string(g)]["wf_pmax"][1, d] = argz["owpp_mva"][j]/pu
						push!(wfz,(parse(Int64,g),argz["owpp_mva"][j]/pu))
				end;end
			end
        end
    end
	#set ["type"]
    for (g, gen) in data["gen"]
		if (isapprox(sum(extradata["gen"][g]["pmax"]),0,atol=10e-3) && isapprox(sum(extradata["gen"][g]["pmin"]),0,atol=10e-3))
			gen["gen_status"]=0;
		else
		end
    end
	#=xy="DE"
	if !(haskey(demand_curve,xy));push!(demand_curve,xy=>DataFrames.DataFrame(:generation=>[],:fuel_cost=>[]));end
	push!(demand_curve[xy],[1,3])=#
    #adjust demand curve
	demand_curve_reduced=Dict()
	for (cuntree, df) in demand_curve
		#cuntree_total=sum(df[!,:generation])
		unique_costs=unique(df[!,:fuel_cost])
		top=maximum(unique_costs)
		bottom=minimum(unique_costs)
		step=(top-bottom)/9
		if !(haskey(demand_curve_reduced,cuntree));push!(demand_curve_reduced,cuntree=>DataFrames.DataFrame(:generation=>[],:fuel_cost=>[]));end
		#for fuel_cost in unique_costs
		for i = 1:1:10
			#capacity=sum(filter(:fuel_cost=>x->x==fuel_cost,df)[!,:generation])
			#push!(demand_curve_reduced[cuntree],[capacity/cuntree_total,fuel_cost])
            #println(string(cuntree)*"|"*string(capacity/cuntree_total)*"|"*string(fuel_cost))
            #push!(demand_curve_reduced[cuntree],[1/length(unique_costs),fuel_cost])
            push!(demand_curve_reduced[cuntree],[1/10,top-(i-1)*step])
			#println(string(cuntree)*"|"*string(1/10)*"|"*string(top-(i-1)*step))

		end
	end

	#add loads
	#all_gens["onshore"]["UK"]["Gas"]["8"]=g
	num_of_gens=length(data["gen"])
	if !(haskey(map_gen_types,"loads"));push!(map_gen_types,"loads"=>Dict());end
	loads=Dict{String,Any}()
	for (cuntree,df) in demand_curve_reduced
		if !(haskey(loads,cuntree));push!(loads,cuntree=>Dict());end
		load=deepcopy(first(first(all_gens["onshore"][cuntree])[2])[2])
		for row in eachrow(df)
			#reset load parameters
	    	load["type"]=1
			load["gen_status"]=1
	        load["index"]=num_of_gens+1
	        load["source_id"][2]=num_of_gens+1
	        load["pmin"]=deepcopy(row[:generation])*-1
	        load["pmax"]=deepcopy(row[:generation])*-1
			load["cost"]=[deepcopy(row[:fuel_cost]),0]
	        push!(loads[cuntree],string(num_of_gens+1)=>deepcopy(load))
	        num_of_gens=num_of_gens+1

			#map load to country
			if !(haskey(map_gen_types["loads"],cuntree));push!(map_gen_types["loads"],cuntree=>[]);end
			push!(map_gen_types["loads"][cuntree],num_of_gens)
	end;end
	for (cuntree, dic) in loads
    	for (l, load) in dic
	        extradata["gen"][string(l)] = Dict{String,Any}()
	        extradata["gen"][string(l)]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
	        extradata["gen"][string(l)]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
	        extradata["gen"][string(l)]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end;end

	#onshore loads
	for (k_sc,sc) in scenario["sc_names"]
		for (k_yr,yr) in sc
			#k_yr_sd=k_yr=="2020" ? "2025" : k_yr;
			k_sc_sd=k_yr=="2020" ? "Base" : k_sc[1:2]
			for (h,d) in enumerate(yr)
		        #loads
				for (cuntree, dic) in loads
			        for (l, load) in sort!(OrderedCollections.OrderedDict(dic), by=x->parse(Int64,x))
						ts=scenario_data["Generation"]["RES"]["Onshore Wind"][cuntree][k_sc[3:6]][!,:time_stamp][h]
						S_row=filter(:time_stamp=>x->x==ts,scenario_data["Demand"][k_sc_sd][k_yr])[!,Symbol(cuntree*"_MWh")]
		                extradata["gen"][string(l)]["pmax"][1, d] = 0
		                extradata["gen"][string(l)]["pmin"][1, d] = load["pmin"]*S_row[1]/pu
		                extradata["gen"][string(l)]["cost"][d] = load["cost"]
		                push!(data["gen"],string(l)=>load)
						push!(genz,(parse(Int64,l),S_row[1]/pu))
			        end
				end
    end;end;end
	unique!(x->first(x),genz)
	unique!(x->first(x),wfz)
	push!(s,"genz"=>genz)
    push!(s,"wfz"=>wfz)
	data_orig["gen"]=data["gen"]
    return extradata, data_orig, map_gen_types
end=#
#=
function create_profile_sets_mesh_wgen_type(number_of_hours, scenario, data_orig, all_gens, scenario_data, markets, s, argz, map_gen_types)
	genz=[];wfz=[]

    pu=data_orig["baseMVA"]
    e2me=1000000/pu#into ME/PU
    extradata = Dict{String,Any}()
    data=Dict{String,Any}();data["gen"]=Dict{String,Any}()
    extradata["dim"] = Dict{String,Any}()
    extradata["dim"] = number_of_hours
    extradata["gen"] = Dict{String,Any}()
	demand_curve = Dict{String,Any}()

	for (c,country) in all_gens["onshore"];
		for (fuel,type) in country;
			for (num,g) in type;
				push!(data["gen"],num=>g)
				;end;end;end

	for (c,country) in all_gens["offshore"];
		for (num,g) in country;
		push!(data["gen"],num=>g);end;end

    data["gen"]=sort!(OrderedCollections.OrderedDict(data["gen"]), by=x->parse(Int64,x))
	for (g, gen) in data["gen"]
        extradata["gen"][string(g)] = Dict{String,Any}()
        extradata["gen"][string(g)]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][string(g)]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
		extradata["gen"][string(g)]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
		if (gen["type"]==0)
		extradata["gen"][string(g)]["invest"] = Array{Float64,2}(undef, 1, number_of_hours);end
    end

	for (k_sc,sc) in scenario["sc_names"]
		for (k_yr,yr) in sc
			k_yr_sd=k_yr=="2020" ? "2025" : k_yr;
			k_sc_sd=k_yr=="2020" ? "NT" : k_sc[1:2]
			for (h,d) in enumerate(yr)
		        #Day ahead BE
		        #onshore generators
				for (xy, country) in all_gens["onshore"]
		        	for (fuel, type) in country
						for (g, gen) in type
		            	#market generator onshore
						S_row=filter(:Generation_Type=>x->x==fuel, scenario_data["Generation"]["Scenarios"][k_sc_sd*k_yr_sd][xy])[:Capacity]
						if isempty(S_row)
							extradata["gen"][string(g)]["pmax"][1, d] = 0
							extradata["gen"][string(g)]["cost"][d] = [0,0]
							extradata["gen"][string(g)]["pmin"][1, d] = 0
		#					extradata["gen"][string(g)]["gen_status"][1, d] = 0
						else
							if issubset([fuel],keys(scenario_data["Generation"]["RES"]))
								CF=scenario_data["Generation"]["RES"][fuel][xy][k_sc[3:6]][!,Symbol(xy*"_MWh")][h]
								extradata["gen"][string(g)]["pmax"][1, d] = CF*S_row[1]/pu
							else
								extradata["gen"][string(g)]["pmax"][1, d] = S_row[1]/pu
							end
							extradata["gen"][string(g)]["pmin"][1, d] = 0
			                extradata["gen"][string(g)]["cost"][d] = [scenario_data["Generation"]["costs"][fuel]/e2me,0]
							push!(demand_curve,string(g)=>(xy,scenario_data["Generation"]["costs"][fuel]/e2me))
							push!(genz,(parse(Int64,g),S_row[1]/pu))
						end
				end;end;end
				#Wind power Plants
				for (j,(xy, country)) in enumerate(all_gens["offshore"])
		        	for (i,(g, gen)) in enumerate(country)
					#wind gen
						CF=scenario_data["Generation"]["RES"]["Offshore Wind"][xy][k_sc[3:6]][!,Symbol(xy*"_MWh")][h]
		                extradata["gen"][string(g)]["pmax"][1, d] = CF
						extradata["gen"][string(g)]["pmin"][1, d] = 0
		                extradata["gen"][string(g)]["cost"][d] = [0,0]
						extradata["gen"][string(g)]["invest"][1, d] = gen["invest"]
						push!(wfz,(parse(Int64,g),argz["owpp_mva"][j]/pu))
				end;end
			end
        end
    end
	#set ["type"]
    for (g, gen) in data["gen"]
		if (isapprox(sum(extradata["gen"][g]["pmax"]),0,atol=10e-3) && isapprox(sum(extradata["gen"][g]["pmin"]),0,atol=10e-3))
			delete!(demand_curve,g)
			gen["gen_status"]=0;
		else
		end
    end
    #add loads
    loads=Dict{String,Any}()
    num_of_gens=length(data["gen"])
	if !(haskey(map_gen_types,"loads"));push!(map_gen_types,"loads"=>Dict());end
    for (g, gen) in sort!(OrderedCollections.OrderedDict(data["gen"]), by=x->parse(Int64,x))
	#for (g, gen) in sort!(OrderedCollections.OrderedDict(data_orig["gen"]))
        if (gen["type"]>0 && gen["gen_status"]==1)#market generator onshore
			#adjust demand curve
			push!(demand_curve,string(num_of_gens+1)=>deepcopy(demand_curve[g]))
			delete!(demand_curve,g)
			country=first(demand_curve[string(num_of_gens+1)])
			if !(haskey(demand_curve,country));push!(demand_curve,country=>0);end
			demand_curve[country]=demand_curve[country]+1

			#set load
            load=deepcopy(gen)
            load["index"]=num_of_gens+1
            load["source_id"][2]=num_of_gens+1
            load["pmin"]=deepcopy(gen["pmax"])*-1
            load["pmax"]=0
            push!(loads,string(num_of_gens+1)=>deepcopy(load))
            num_of_gens=num_of_gens+1

			#map load to country
			if !(haskey(map_gen_types["loads"],country));push!(map_gen_types["loads"],country=>[]);end
			push!(map_gen_types["loads"][country],num_of_gens)
        end
    end
    for (l, load) in loads
        extradata["gen"][string(l)] = Dict{String,Any}()
        extradata["gen"][string(l)]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][string(l)]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["gen"][string(l)]["cost"] = [Vector{Float64}() for i=1:number_of_hours]
    end
	println("number of loads: "*string(length(loads)))
	#onshore loads
	for (k_sc,sc) in scenario["sc_names"]
		for (k_yr,yr) in sc
			k_yr_sd=k_yr=="2020" ? "2025" : k_yr;
			k_sc_sd=k_yr=="2020" ? "NT" : k_sc[1:2]
			for (h,d) in enumerate(yr)
		        #loads
		        for (load_num,(l, load)) in enumerate(sort!(OrderedCollections.OrderedDict(loads), by=x->parse(Int64,x)))
		            if (load["type"]>0)#market generator onshore
						country=first(demand_curve[l])
						cst=last(demand_curve[l])
						country_gens=demand_curve[country]
						ts=scenario_data["Generation"]["RES"]["Onshore Wind"][country][k_sc[3:6]][!,:time_stamp][h]
						S_row=filter(:time_stamp=>x->x==ts,scenario_data["Demand"][k_sc_sd*k_yr_sd])[!,Symbol(country*"_MWh")]
		                extradata["gen"][string(l)]["pmax"][1, d] = 0
		                extradata["gen"][string(l)]["pmin"][1, d] = ((S_row[1]/pu)/country_gens)*-1
		                extradata["gen"][string(l)]["cost"][d] = [cst,0]
		                push!(data["gen"],string(l)=>load)
						push!(genz,(parse(Int64,l),S_row[1]/pu))
		            end
		        end
    end;end;end
	unique!(x->first(x),genz)
	unique!(x->first(x),wfz)
	push!(s,"genz"=>genz)
    push!(s,"wfz"=>wfz)
	data_orig["gen"]=data["gen"]
    return extradata, data_orig, map_gen_types
end
=#

######################## Scaling cost data ###########################
#scale investment to hourly cost spread over the year
#***#
function scale_cost_data_hourly!(data, scenario)
    #rescale_hourly = x -> (8760*scenario["planning_horizon"] / (scenario["hours"])) * x # scale hourly costs to the planning horizon
    rescale_hourly = x -> (8760*scenario["planning_horizon"] / (scenario["hours"]*scenario["years"])) * x # scale hourly costs to the planning horizon
    rescale_total  = x -> (                                1 / (scenario["hours"])) * x # scale total costs to the planning horizon

    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "cost", rescale_hourly)
    end
    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "invest", rescale_total)
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        _PM._apply_func!(branch, "construction_cost", rescale_total)
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "eq_cost", rescale_total)
        _PM._apply_func!(strg, "inst_cost", rescale_total)
        _PM._apply_func!(strg, "cost_abs", rescale_hourly)
        _PM._apply_func!(strg, "cost_inj", rescale_hourly)
    end
    for (b, branch) in get(data, "branchdc", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (b, branch) in get(data, "branch", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
    end
    for (s, strg) in get(data, "storage", Dict{String,Any}())
        _PM._apply_func!(strg, "cost", rescale_total)
    end
end

function scale_cost_data_2hourly!(data, scenario)
    rescale_hourly = x -> (8760*scenario["planning_horizon"] / (scenario["hours"])) * x # scale hourly costs to the planning horizon
    rescale_total  = x -> (                                1 / (scenario["hours"])) * x # scale total costs to the planning horizon

    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "cost", rescale_hourly)
    end
    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "invest", rescale_total)
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        _PM._apply_func!(branch, "construction_cost", rescale_total)
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "eq_cost", rescale_total)
        _PM._apply_func!(strg, "inst_cost", rescale_total)
        _PM._apply_func!(strg, "cost_abs", rescale_hourly)
        _PM._apply_func!(strg, "cost_inj", rescale_hourly)
    end
    for (b, branch) in get(data, "branchdc", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (b, branch) in get(data, "branch", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
    end
    for (s, strg) in get(data, "storage", Dict{String,Any}())
        _PM._apply_func!(strg, "cost", rescale_total)
    end
end


#scales data to hourly cost spread over horizon years
function scale_cost_data_2yearly!(data, scenario)
    rescale_hourly = x -> (8760*scenario["planning_horizon"] / (scenario["hours"]*scenario["years"])) * x # scale hourly costs to the planning horizon
    rescale_total  = x -> (                                1 / (scenario["hours"]*scenario["years"])) * x # scale total costs to the planning horizon

    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "cost", rescale_hourly)
    end
    for (g, gen) in data["gen"]
        _PM._apply_func!(gen, "invest", rescale_total)
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        _PM._apply_func!(branch, "construction_cost", rescale_total)
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "eq_cost", rescale_total)
        _PM._apply_func!(strg, "inst_cost", rescale_total)
        _PM._apply_func!(strg, "cost_abs", rescale_hourly)
        _PM._apply_func!(strg, "cost_inj", rescale_hourly)
    end
    for (b, branch) in get(data, "branchdc", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (b, branch) in get(data, "branch", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", rescale_total)
    end
    for (c, conv) in get(data, "convdc", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", rescale_total)
    end
    for (s, strg) in get(data, "storage", Dict{String,Any}())
        _PM._apply_func!(strg, "cost", rescale_total)
    end
end


##################### NPV calculations ###############################
#translates costs (array) to yearly NPV value
#***#
function npvs_costs_datas_wREZ(data, scenario, _yrs, _dr)
    _scs=data["scenario"]
    base_yr=parse(Int64,_yrs[1])
    _hrs=deepcopy(scenario["hours"])
    #_yrs=[k for k in keys(scenario["sc_names"][_scs[1]])]
    for (_sci,_sc) in _scs
        _sc_temp=sort!(OrderedCollections.OrderedDict(deepcopy(_sc)), by=x->parse(Int64,x))
        _sc_first=first(_sc_temp)[2]
        _yr=1
        for (_str,_num) in _sc_temp
            if (_num<=_sc_first+_yr*_hrs-1)
                #data["nw"][string(_num)]=deepcopy(npv_cost_data_wREZ(deepcopy(data["nw"][string(_num)]),base_yr,parse(Int64,_yrs[_yr]),scenario["planning_horizon"], _dr))
                data["nw"][string(_num)]=deepcopy(npv_cost_data_wREZ(data["nw"][string(_num)],base_yr,parse(Int64,_yrs[_yr]),scenario["planning_horizon"], _dr))
                if (_num==_sc_first+_yr*_hrs-1)
                    _yr=_yr+1
                end
            end
        end
    end
    return data
end

#=function calc_single_convdc_cost_npv(i, b_cost, nw)
    cost = 0.0
    cost += b_cost * _PM.var(pm,nw,:p_pacmax,i)
    return cost
end=#
#translates cost (single) to yearly NPV value
function npv_cost_data_wREZ(data,base_yr,current_yr,_ph,_dr::Float64=0.04)
    #println(((_ph-(current_yr-base_yr))/_ph))
    function npv_yearly(x)
        cost = (1 / (1+_dr)^(current_yr-base_yr)) * x *((_ph-(current_yr-base_yr))/_ph)# npv
        return deepcopy(cost)
    end
    function npv_hourly(x)
        cost = (1 / (1+_dr)^(current_yr-base_yr)) * x# npv
        return deepcopy(cost)
    end

    for (g, gen) in get(data, "gen", Dict{String,Any}())
        gen["cost"]=npv_hourly(gen["cost"])
    end
    for (g, gen) in get(data, "gen", Dict{String,Any}())
        gen["invest"]=npv_yearly(gen["invest"])
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        branch["construction_cost"]=npv_yearly(branch["construction_cost"])
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        branch["cost"]=npv_yearly(branch["cost"])
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        conv["cost"]=npv_yearly(conv["cost"])
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        strg["eq_cost"]=npv_yearly(strg["eq_cost"])
        strg["inst_cost"]=npv_yearly(strg["inst_cost"])
        strg["cost_abs"]=npv_hourly(strg["cost_abs"])
        strg["cost_inj"]=npv_hourly(strg["cost_inj"])
    end
    for (b, branch) in get(data, "branchdc", Dict{String,Any}())
        branch["cost"]=npv_yearly(branch["cost"])
    end
    for (b, branch) in get(data, "branch", Dict{String,Any}())
        branch["cost"]=npv_yearly(branch["cost"])
    end
    for (c, conv) in get(data, "convdc", Dict{String,Any}())
        data["convdc"][c]["cost"]=deepcopy(npv_yearly(conv["cost"]))
    end
    for (s, strg) in get(data, "storage", Dict{String,Any}())
        strg["cost"]=npv_yearly(strg["cost"])
    end
    return data
end
#=
function npv_cost_data_wREZ_wdeepcopy(data,base_yr,current_yr,_ph,_dr::Float64=0.04)
    #println(((_ph-(current_yr-base_yr))/_ph))
    npv_yearly = x -> (1 / (1+_dr)^(current_yr-base_yr)) * x *((_ph-(current_yr-base_yr))/_ph)# npv
    npv_hourly = x -> (1 / (1+_dr)^(current_yr-base_yr)) * x# npv
    for (g, gen) in get(data, "gen", Dict{String,Any}())
        _PM._apply_func!(gen, "cost", npv_hourly)
    end
    for (g, gen) in get(data, "gen", Dict{String,Any}())
        _PM._apply_func!(gen, "invest", npv_yearly)
    end
    for (b, branch) in get(data, "ne_branch", Dict{String,Any}())
        _PM._apply_func!(branch, "construction_cost", npv_yearly)
    end
    for (b, branch) in get(data, "branchdc_ne", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", npv_yearly)
    end
    for (c, conv) in get(data, "convdc_ne", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", npv_yearly)
    end
    for (s, strg) in get(data, "ne_storage", Dict{String,Any}())
        _PM._apply_func!(strg, "eq_cost", npv_yearly)
        _PM._apply_func!(strg, "inst_cost", npv_yearly)
        _PM._apply_func!(strg, "cost_abs", npv_hourly)
        _PM._apply_func!(strg, "cost_inj", npv_hourly)
    end
    for (b, branch) in get(data, "branchdc", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", npv_yearly)
    end
    for (b, branch) in get(data, "branch", Dict{String,Any}())
        _PM._apply_func!(branch, "cost", npv_yearly)
    end
    for (c, conv) in get(data, "convdc", Dict{String,Any}())
        _PM._apply_func!(conv, "cost", npv_yearly)
    end
    for (s, strg) in get(data, "storage", Dict{String,Any}())
        _PM._apply_func!(strg, "cost", npv_yearly)
    end
    return data
end=#

function create_profile_sets_rest(number_of_hours, extradata, data_orig)
    pu=data_orig["baseMVA"]
    e2me=1000000/pu#into ME/PU
    data=Dict{String,Any}()#;data["convdc"]=Dict{String,Any}()

    extradata["convdc"] = Dict{String,Any}()
    data["convdc"]=sort!(OrderedCollections.OrderedDict(data_orig["convdc"]), by=x->parse(Int64,x))
    for (c, cnv) in data["convdc"]
        extradata["convdc"][c] = Dict{String,Any}()
        extradata["convdc"][c]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (c, cnv) in data["convdc"]
                extradata["convdc"][c]["cost"][1, d] = cnv["cost"]
        end
    end

    data["gen"]=sort!(OrderedCollections.OrderedDict(data_orig["gen"]), by=x->parse(Int64,x))
    for (g, gen) in data["gen"]
        if !(haskey(extradata["gen"], g));extradata["gen"][g] = Dict{String,Any}();end
        extradata["gen"][g]["invest"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (g, gen) in data["gen"]
                extradata["gen"][g]["invest"][1, d] = gen["invest"]
        end
    end

    extradata["ne_branch"] = Dict{String,Any}()
    data["ne_branch"]=sort!(OrderedCollections.OrderedDict(data_orig["ne_branch"]), by=x->parse(Int64,x))
    for (b, br) in data["ne_branch"]
        extradata["ne_branch"][b] = Dict{String,Any}()
        extradata["ne_branch"][b]["construction_cost"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (b, br) in data["ne_branch"]
                extradata["ne_branch"][b]["construction_cost"][1, d] = br["construction_cost"]
        end
    end

    extradata["branchdc_ne"] = Dict{String,Any}()
    data["branchdc_ne"]=sort!(OrderedCollections.OrderedDict(data_orig["branchdc_ne"]), by=x->parse(Int64,x))
    for (b, br) in data["branchdc_ne"]
        extradata["branchdc_ne"][b] = Dict{String,Any}()
        extradata["branchdc_ne"][b]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (b, br) in data["branchdc_ne"]
                extradata["branchdc_ne"][b]["cost"][1, d] = br["cost"]
        end
    end

    extradata["convdc_ne"] = Dict{String,Any}()
    data["convdc_ne"]=sort!(OrderedCollections.OrderedDict(data_orig["convdc_ne"]), by=x->parse(Int64,x))
    for (c, cnv) in data["convdc_ne"]
        extradata["convdc_ne"][c] = Dict{String,Any}()
        extradata["convdc_ne"][c]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (c, cnv) in data["convdc_ne"]
                extradata["convdc_ne"][c]["cost"][1, d] = cnv["cost"]
        end
    end

    extradata["ne_storage"] = Dict{String,Any}()
    data["ne_storage"]=sort!(OrderedCollections.OrderedDict(data_orig["ne_storage"]), by=x->parse(Int64,x))
    for (s, stg) in data["ne_storage"]
        extradata["ne_storage"][s] = Dict{String,Any}()
        extradata["ne_storage"][s]["eq_cost"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][s]["inst_cost"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][s]["cost_abs"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][s]["cost_inj"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (s, stg) in data["ne_storage"]
                extradata["ne_storage"][s]["eq_cost"][1, d] = stg["eq_cost"]
                extradata["ne_storage"][s]["inst_cost"][1, d] = stg["inst_cost"]
                extradata["ne_storage"][s]["cost_abs"][1, d] = stg["cost_abs"]
                extradata["ne_storage"][s]["cost_inj"][1, d] = stg["cost_inj"]
        end
    end

    extradata["storage"] = Dict{String,Any}()
    data["storage"]=sort!(OrderedCollections.OrderedDict(data_orig["storage"]), by=x->parse(Int64,x))
    for (s, stg) in data["storage"]
        extradata["storage"][s] = Dict{String,Any}()
        extradata["storage"][s]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (s, stg) in data["storage"]
                extradata["storage"][s]["cost"][1, d] = stg["cost"]
        end
    end

    extradata["branchdc"] = Dict{String,Any}()
    data["branchdc"]=sort!(OrderedCollections.OrderedDict(data_orig["branchdc"]), by=x->parse(Int64,x))
    for (b, br) in data["branchdc"]
        extradata["branchdc"][b] = Dict{String,Any}()
        extradata["branchdc"][b]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (b, br) in data["branchdc"]
                extradata["branchdc"][b]["cost"][1, d] = br["cost"]
        end
    end

    extradata["branch"] = Dict{String,Any}()
    data["branch"]=sort!(OrderedCollections.OrderedDict(data_orig["branch"]), by=x->parse(Int64,x))
    for (b, br) in data["branch"]
        extradata["branch"][b] = Dict{String,Any}()
        extradata["branch"][b]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (b, br) in data["branch"]
                extradata["branch"][b]["cost"][1, d] = br["cost"]
        end
    end
    return extradata
end

#***#
function create_profile_sets_rest_wgen_type(extradata, data_orig, s)
    pu=data_orig["baseMVA"]
    e2me=1000000/pu#into ME/PU
    data=Dict{String,Any}()#;data["convdc"]=Dict{String,Any}()

    extradata["convdc"] = Dict{String,Any}()
    data["convdc"]=sort!(OrderedCollections.OrderedDict(data_orig["convdc"]), by=x->parse(Int64,x))
    for (c, cnv) in data["convdc"]
        extradata["convdc"][c] = Dict{String,Any}()
        extradata["convdc"][c]["cost"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["convdc"][c]["Pacmax"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["convdc"][c]["Pacmin"] = Array{Float64,2}(undef, 1, s["dim"])
    end
    for d in 1:s["dim"]
        for (c, cnv) in data["convdc"]
                extradata["convdc"][c]["cost"][1, d] = cnv["cost"]
                extradata["convdc"][c]["Pacmin"][1, d] = 0.0
                if issubset([parse(Int64,c)],s["onshore_nodes"])
                    #println(c*" onshore node")
                    extradata["convdc"][c]["Pacmax"][1, d] = s["conv_lim_onshore"]/pu
                elseif issubset([parse(Int64,c)],s["offshore_nodes"])
                    #println(c*" offshore node")
                    extradata["convdc"][c]["Pacmax"][1, d] = s["conv_lim_offshore"]/pu
                else
                    extradata["convdc"][c]["Pacmax"][1, d] = data_orig["convdc"][c]["Pacmax"]
                end
        end
    end

    extradata["ne_branch"] = Dict{String,Any}()
    data["ne_branch"]=sort!(OrderedCollections.OrderedDict(data_orig["ne_branch"]), by=x->parse(Int64,x))
    for (b, br) in data["ne_branch"]
        extradata["ne_branch"][b] = Dict{String,Any}()
        extradata["ne_branch"][b]["construction_cost"] = Array{Float64,2}(undef, 1, s["dim"])
    end
    for d in 1:s["dim"]
        for (b, br) in data["ne_branch"]
                extradata["ne_branch"][b]["construction_cost"][1, d] = br["construction_cost"]
        end
    end

    extradata["branchdc_ne"] = Dict{String,Any}()
    data["branchdc_ne"]=sort!(OrderedCollections.OrderedDict(data_orig["branchdc_ne"]), by=x->parse(Int64,x))
    for (b, br) in data["branchdc_ne"]
        extradata["branchdc_ne"][b] = Dict{String,Any}()
        extradata["branchdc_ne"][b]["cost"] = Array{Float64,2}(undef, 1, s["dim"])
    end
    for d in 1:s["dim"]
        for (b, br) in data["branchdc_ne"]
                extradata["branchdc_ne"][b]["cost"][1, d] = br["cost"]
        end
    end

    extradata["convdc_ne"] = Dict{String,Any}()
    data["convdc_ne"]=sort!(OrderedCollections.OrderedDict(data_orig["convdc_ne"]), by=x->parse(Int64,x))
    for (c, cnv) in data["convdc_ne"]
        extradata["convdc_ne"][c] = Dict{String,Any}()
        extradata["convdc_ne"][c]["cost"] = Array{Float64,2}(undef, 1, s["dim"])
    end
    for d in 1:s["dim"]
        for (c, cnv) in data["convdc_ne"]
                extradata["convdc_ne"][c]["cost"][1, d] = cnv["cost"]
        end
    end

    extradata["ne_storage"] = Dict{String,Any}()
    data["ne_storage"]=sort!(OrderedCollections.OrderedDict(data_orig["ne_storage"]), by=x->parse(Int64,x))
    for (b, stg) in data["ne_storage"]
        extradata["ne_storage"][b] = Dict{String,Any}()
        extradata["ne_storage"][b]["eq_cost"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["ne_storage"][b]["inst_cost"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["ne_storage"][b]["cost_abs"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["ne_storage"][b]["cost_inj"] = Array{Float64,2}(undef, 1, s["dim"])
    end
    for d in 1:s["dim"]
        for (b, stg) in data["ne_storage"]
                extradata["ne_storage"][b]["eq_cost"][1, d] = stg["eq_cost"]
                extradata["ne_storage"][b]["inst_cost"][1, d] = stg["inst_cost"]
                extradata["ne_storage"][b]["cost_abs"][1, d] = stg["cost_abs"]
                extradata["ne_storage"][b]["cost_inj"][1, d] = stg["cost_inj"]
        end
    end

    extradata["storage"] = Dict{String,Any}()
    data["storage"]=sort!(OrderedCollections.OrderedDict(data_orig["storage"]), by=x->parse(Int64,x))
    for (b, stg) in data["storage"]
        extradata["storage"][b] = Dict{String,Any}()
        extradata["storage"][b]["cost"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["storage"][b]["pmax"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["storage"][b]["pmin"] = Array{Float64,2}(undef, 1, s["dim"])
    end
    for d in 1:s["dim"]
        for (b, stg) in data["storage"]
                extradata["storage"][b]["cost"][1, d] = stg["cost"]
                if (issubset([stg["storage_bus"]],s["offshore_nodes"]))
                    extradata["storage"][b]["pmax"][1, d] = s["strg_lim_offshore"]
                else
                    extradata["storage"][b]["pmax"][1, d] = s["strg_lim_onshore"];end
                extradata["storage"][b]["pmin"][1, d] = 0.0
        end
    end
    

    extradata["branchdc"] = Dict{String,Any}()
    data["branchdc"]=sort!(OrderedCollections.OrderedDict(data_orig["branchdc"]), by=x->parse(Int64,x))
    for (b, br) in data["branchdc"]
        extradata["branchdc"][b] = Dict{String,Any}()
        extradata["branchdc"][b]["cost"] = Array{Float64,2}(undef, 1, s["dim"])
        ##################
        extradata["branchdc"][b]["rateA"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["branchdc"][b]["r"] = Array{Float64,2}(undef, 1, s["dim"]) 
        ##################
    end
    for d in 1:s["dim"]
        for (b, br) in data["branchdc"]
                extradata["branchdc"][b]["cost"][1, d] = br["cost"]
                #######################################################
                extradata["branchdc"][b]["rateA"][1, d] = br["rateA"]
                extradata["branchdc"][b]["r"][1, d] = br["r"]
                #######################################################
        end
    end

    extradata["branch"] = Dict{String,Any}()
    data["branch"]=sort!(OrderedCollections.OrderedDict(data_orig["branch"]), by=x->parse(Int64,x))
    for (b, br) in data["branch"]
        extradata["branch"][b] = Dict{String,Any}()
        extradata["branch"][b]["cost"] = Array{Float64,2}(undef, 1, s["dim"])
        ##################
        extradata["branch"][b]["rateA"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["branch"][b]["br_r"] = Array{Float64,2}(undef, 1, s["dim"])
        extradata["branch"][b]["br_x"] = Array{Float64,2}(undef, 1, s["dim"]) 
        ##################
    end
    for d in 1:s["dim"]
        for (b, br) in data["branch"]
                extradata["branch"][b]["cost"][1, d] = br["cost"]
                extradata["branch"][b]["rateA"][1, d] = br["rateA"]
                extradata["branch"][b]["br_r"][1, d] = br["br_r"]
                extradata["branch"][b]["br_x"][1, d] = br["br_x"]
        end
    end
    return extradata
end

#=
function create_profile_sets_rest_wgen_type(number_of_hours, extradata, data_orig, argz)
    pu=data_orig["baseMVA"]
    e2me=1000000/pu#into ME/PU
    data=Dict{String,Any}()#;data["convdc"]=Dict{String,Any}()

    extradata["convdc"] = Dict{String,Any}()
    data["convdc"]=sort!(OrderedCollections.OrderedDict(data_orig["convdc"]), by=x->parse(Int64,x))
    for (c, cnv) in data["convdc"]
        extradata["convdc"][c] = Dict{String,Any}()
        extradata["convdc"][c]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["convdc"][c]["Pacmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["convdc"][c]["Pacmin"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (c, cnv) in data["convdc"]
                extradata["convdc"][c]["cost"][1, d] = cnv["cost"]
                extradata["convdc"][c]["Pacmax"][1, d] = argz["conv_lim"]/pu
                extradata["convdc"][c]["Pacmin"][1, d] = 0.0
        end
    end

    extradata["ne_branch"] = Dict{String,Any}()
    data["ne_branch"]=sort!(OrderedCollections.OrderedDict(data_orig["ne_branch"]), by=x->parse(Int64,x))
    for (b, br) in data["ne_branch"]
        extradata["ne_branch"][b] = Dict{String,Any}()
        extradata["ne_branch"][b]["construction_cost"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (b, br) in data["ne_branch"]
                extradata["ne_branch"][b]["construction_cost"][1, d] = br["construction_cost"]
        end
    end

    extradata["branchdc_ne"] = Dict{String,Any}()
    data["branchdc_ne"]=sort!(OrderedCollections.OrderedDict(data_orig["branchdc_ne"]), by=x->parse(Int64,x))
    for (b, br) in data["branchdc_ne"]
        extradata["branchdc_ne"][b] = Dict{String,Any}()
        extradata["branchdc_ne"][b]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (b, br) in data["branchdc_ne"]
                extradata["branchdc_ne"][b]["cost"][1, d] = br["cost"]
        end
    end

    extradata["convdc_ne"] = Dict{String,Any}()
    data["convdc_ne"]=sort!(OrderedCollections.OrderedDict(data_orig["convdc_ne"]), by=x->parse(Int64,x))
    for (c, cnv) in data["convdc_ne"]
        extradata["convdc_ne"][c] = Dict{String,Any}()
        extradata["convdc_ne"][c]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (c, cnv) in data["convdc_ne"]
                extradata["convdc_ne"][c]["cost"][1, d] = cnv["cost"]
        end
    end

    extradata["ne_storage"] = Dict{String,Any}()
    data["ne_storage"]=sort!(OrderedCollections.OrderedDict(data_orig["ne_storage"]), by=x->parse(Int64,x))
    for (s, stg) in data["ne_storage"]
        extradata["ne_storage"][s] = Dict{String,Any}()
        extradata["ne_storage"][s]["eq_cost"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][s]["inst_cost"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][s]["cost_abs"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["ne_storage"][s]["cost_inj"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (s, stg) in data["ne_storage"]
                extradata["ne_storage"][s]["eq_cost"][1, d] = stg["eq_cost"]
                extradata["ne_storage"][s]["inst_cost"][1, d] = stg["inst_cost"]
                extradata["ne_storage"][s]["cost_abs"][1, d] = stg["cost_abs"]
                extradata["ne_storage"][s]["cost_inj"][1, d] = stg["cost_inj"]
        end
    end

    extradata["storage"] = Dict{String,Any}()
    data["storage"]=sort!(OrderedCollections.OrderedDict(data_orig["storage"]), by=x->parse(Int64,x))
    for (s, stg) in data["storage"]
        extradata["storage"][s] = Dict{String,Any}()
        extradata["storage"][s]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["storage"][s]["pmax"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["storage"][s]["pmin"] = Array{Float64,2}(undef, 1, number_of_hours)
    end
    for d in 1:number_of_hours
        for (s, stg) in data["storage"]
                extradata["storage"][s]["cost"][1, d] = stg["cost"]
                extradata["storage"][s]["pmax"][1, d] = argz["strg_lim_onshore"]#argz["strg_lim_offshore"]
                extradata["storage"][s]["pmin"][1, d] = 0.0
        end
    end
    

    extradata["branchdc"] = Dict{String,Any}()
    data["branchdc"]=sort!(OrderedCollections.OrderedDict(data_orig["branchdc"]), by=x->parse(Int64,x))
    for (b, br) in data["branchdc"]
        extradata["branchdc"][b] = Dict{String,Any}()
        extradata["branchdc"][b]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
        ##################
        extradata["branchdc"][b]["rateA"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["branchdc"][b]["r"] = Array{Float64,2}(undef, 1, number_of_hours) 
        ##################
    end
    for d in 1:number_of_hours
        for (b, br) in data["branchdc"]
                extradata["branchdc"][b]["cost"][1, d] = br["cost"]
                #######################################################
                extradata["branchdc"][b]["rateA"][1, d] = br["rateA"]
                extradata["branchdc"][b]["r"][1, d] = br["r"]
                #######################################################
        end
    end

    extradata["branch"] = Dict{String,Any}()
    data["branch"]=sort!(OrderedCollections.OrderedDict(data_orig["branch"]), by=x->parse(Int64,x))
    for (b, br) in data["branch"]
        extradata["branch"][b] = Dict{String,Any}()
        extradata["branch"][b]["cost"] = Array{Float64,2}(undef, 1, number_of_hours)
        ##################
        extradata["branch"][b]["rateA"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["branch"][b]["br_r"] = Array{Float64,2}(undef, 1, number_of_hours)
        extradata["branch"][b]["br_x"] = Array{Float64,2}(undef, 1, number_of_hours) 
        ##################
    end
    for d in 1:number_of_hours
        for (b, br) in data["branch"]
                extradata["branch"][b]["cost"][1, d] = br["cost"]
                extradata["branch"][b]["rateA"][1, d] = br["rateA"]
                extradata["branch"][b]["br_r"][1, d] = br["br_r"]
                extradata["branch"][b]["br_x"][1, d] = br["br_x"]
        end
    end
    return extradata
end=#

function costs_datas_wREZ(extradata, s)
    function npv_yearly(x,current_yr)
		#println(string(x)*" "*string(current_yr))
        cost = (1 / (1+s["dr"])^(current_yr-base_year)) * x *((s["scenario_planning_horizon"]-(current_yr-base_year))/s["scenario_planning_horizon"])# npv
        return deepcopy(cost)
    end
    function npv_hourly(x,current_yr)
        cost = (1 / (1+s["dr"])^(current_yr-base_year)) * x# npv
        return deepcopy(cost)
    end
    function mip(cost0, cost1)
        cost=cost0-cost1
        return deepcopy(cost)
    end
    base_year=parse(Int64,s["scenario_years"][1])
    sl=s["scenarios_length"]
    yl=s["years_length"]
    hl=s["hours_length"]
    if (haskey(extradata,"convdc"))
        for (c,cnv) in extradata["convdc"]
            for (n,cst) in enumerate(cnv["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["convdc"][c]["cost"][n]=npv_yearly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
        end
    end
    if (haskey(extradata,"gen"))
        for (g,gen) in extradata["gen"]
            if (haskey(gen,"cost"))
                for (n,cst) in enumerate(gen["cost"])
                    _sc=floor(Int64,(n-1)/(yl*hl))
                    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                    extradata["gen"][g]["cost"][n]=npv_hourly(cst,parse(Int64,s["scenario_years"][_yr]))
                end
            end
        end
        for (g,gen) in extradata["gen"]
            if (haskey(gen,"invest"))
                for (n,cst) in enumerate(gen["invest"])
                    _sc=floor(Int64,(n-1)/(yl*hl))
                    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
					#println(typeof(cst))
                    extradata["gen"][g]["invest"][n]=npv_yearly(cst,parse(Int64,s["scenario_years"][_yr]))
                end
            end
        end
    end
    if (haskey(extradata,"ne_branch"))
        for (b,br) in extradata["ne_branch"]
            for (n,cst) in enumerate(br["construction_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["ne_branch"][b]["construction_cost"][n]=npv_yearly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
            for (n,cst) in enumerate(br["construction_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                if (n<=(_sc)*yl*hl+(yl-1)*hl)
                    extradata["ne_branch"][b]["construction_cost"][n]=mip(cst,extradata["ne_branch"][b]["construction_cost"][n+hl])
                end
            end
        end
    end
    if (haskey(extradata,"branchdc_ne"))
        for (b,br) in extradata["branchdc_ne"]
            for (n,cst) in enumerate(br["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["branchdc_ne"][b]["cost"][n]=npv_yearly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
            for (n,cst) in enumerate(br["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                if (n<=(_sc)*yl*hl+(yl-1)*hl)
                    extradata["branchdc_ne"][b]["cost"][n]=mip(cst,extradata["branchdc_ne"][b]["cost"][n+hl])
                end
            end
        end
    end
    if (haskey(extradata,"convdc_ne"))
        for (c,cnv) in extradata["convdc_ne"]
            for (n,cst) in enumerate(cnv["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["convdc_ne"][c]["cost"][n]=npv_yearly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
            for (n,cst) in enumerate(cnv["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                if (n<=(_sc)*yl*hl+(yl-1)*hl)
                    extradata["convdc_ne"][c]["cost"][n]=mip(cst,extradata["convdc_ne"][c]["cost"][n+hl])
                end
            end
        end
    end
    if (haskey(extradata,"ne_storage"))
        for (b,stg) in extradata["ne_storage"]
            for (n,cst) in enumerate(stg["eq_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["ne_storage"][b]["eq_cost"][n]=npv_yearly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
            for (n,cst) in enumerate(stg["eq_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                if (n<=(_sc)*yl*hl+(yl-1)*hl)
                    extradata["ne_storage"][b]["eq_cost"][n]=mip(cst,extradata["ne_storage"][b]["eq_cost"][n+hl])
                end
            end
        end
        for (b,stg) in extradata["ne_storage"]
            for (n,cst) in enumerate(stg["inst_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["ne_storage"][b]["inst_cost"][n]=npv_yearly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
            for (n,cst) in enumerate(stg["inst_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                if (n<=(_sc)*yl*hl+(yl-1)*hl)
                    extradata["ne_storage"][b]["inst_cost"][n]=mip(cst,extradata["ne_storage"][b]["inst_cost"][n+hl])
                end
            end
        end
        for (b,stg) in extradata["ne_storage"]
            for (n,cst) in enumerate(stg["cost_abs"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["ne_storage"][b]["cost_abs"][n]=npv_hourly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
        end
        for (b,stg) in extradata["ne_storage"]
            for (n,cst) in enumerate(stg["cost_inj"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["ne_storage"][b]["cost_inj"][n]=npv_hourly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
        end
    end
    if (haskey(extradata,"branchdc"))
        for (b,br) in extradata["branchdc"]
            for (n,cst) in enumerate(br["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["branchdc"][b]["cost"][n]=npv_yearly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
        end
    end
    if (haskey(extradata,"branch"))
        for (b,br) in extradata["branch"]
            for (n,cst) in enumerate(br["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["branch"][b]["cost"][n]=npv_yearly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
        end
    end
    if (haskey(extradata,"storage"))
        for (b,stg) in extradata["storage"]
            for (n,cst) in enumerate(stg["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["storage"][b]["cost"][n]=npv_yearly(cst,parse(Int64,s["scenario_years"][_yr]))
            end
        end
    end
    return extradata
end

function costs_datas_wREZ(extradata, s, argz)
    function npv_yearly(x,current_yr)
		#println(string(x)*" "*string(current_yr))
        cost = (1 / (1+argz["dr"])^(current_yr-base_year)) * x *((argz["scenario_planning_horizon"]-(current_yr-base_year))/argz["scenario_planning_horizon"])# npv
        return deepcopy(cost)
    end
    function npv_hourly(x,current_yr)
        cost = (1 / (1+argz["dr"])^(current_yr-base_year)) * x# npv
        return deepcopy(cost)
    end
    function mip(cost0, cost1)
        cost=cost0-cost1
        return deepcopy(cost)
    end
    base_year=parse(Int64,argz["scenario_years"][1])
    sl=s["scenarios_length"]
    yl=s["years_length"]
    hl=s["hours_length"]
    if (haskey(extradata,"convdc"))
        for (c,cnv) in extradata["convdc"]
            for (n,cst) in enumerate(cnv["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["convdc"][c]["cost"][n]=npv_yearly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
        end
    end
    if (haskey(extradata,"gen"))
        for (g,gen) in extradata["gen"]
            if (haskey(gen,"cost"))
                for (n,cst) in enumerate(gen["cost"])
                    _sc=floor(Int64,(n-1)/(yl*hl))
                    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                    extradata["gen"][g]["cost"][n]=npv_hourly(cst,parse(Int64,argz["scenario_years"][_yr]))
                end
            end
        end
        for (g,gen) in extradata["gen"]
            if (haskey(gen,"invest"))
                for (n,cst) in enumerate(gen["invest"])
                    _sc=floor(Int64,(n-1)/(yl*hl))
                    _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
					#println(typeof(cst))
                    extradata["gen"][g]["invest"][n]=npv_yearly(cst,parse(Int64,argz["scenario_years"][_yr]))
                end
            end
        end
    end
    if (haskey(extradata,"ne_branch"))
        for (b,br) in extradata["ne_branch"]
            for (n,cst) in enumerate(br["construction_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["ne_branch"][b]["construction_cost"][n]=npv_yearly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
            for (n,cst) in enumerate(br["construction_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                if (n<=(_sc)*yl*hl+(yl-1)*hl)
                    extradata["ne_branch"][b]["construction_cost"][n]=mip(cst,extradata["ne_branch"][b]["construction_cost"][n+hl])
                end
            end
        end
    end
    if (haskey(extradata,"branchdc_ne"))
        for (b,br) in extradata["branchdc_ne"]
            for (n,cst) in enumerate(br["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["branchdc_ne"][b]["cost"][n]=npv_yearly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
            for (n,cst) in enumerate(br["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                if (n<=(_sc)*yl*hl+(yl-1)*hl)
                    extradata["branchdc_ne"][b]["cost"][n]=mip(cst,extradata["branchdc_ne"][b]["cost"][n+hl])
                end
            end
        end
    end
    if (haskey(extradata,"convdc_ne"))
        for (c,cnv) in extradata["convdc_ne"]
            for (n,cst) in enumerate(cnv["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["convdc_ne"][c]["cost"][n]=npv_yearly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
            for (n,cst) in enumerate(cnv["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                if (n<=(_sc)*yl*hl+(yl-1)*hl)
                    extradata["convdc_ne"][c]["cost"][n]=mip(cst,extradata["convdc_ne"][c]["cost"][n+hl])
                end
            end
        end
    end
    if (haskey(extradata,"ne_storage"))
        for (s,stg) in extradata["ne_storage"]
            for (n,cst) in enumerate(stg["eq_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["ne_storage"][s]["eq_cost"][n]=npv_yearly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
            for (n,cst) in enumerate(stg["eq_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                if (n<=(_sc)*yl*hl+(yl-1)*hl)
                    extradata["ne_storage"][s]["eq_cost"][n]=mip(cst,extradata["ne_storage"][s]["eq_cost"][n+hl])
                end
            end
        end
        for (s,stg) in extradata["ne_storage"]
            for (n,cst) in enumerate(stg["inst_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["ne_storage"][s]["inst_cost"][n]=npv_yearly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
            for (n,cst) in enumerate(stg["inst_cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                if (n<=(_sc)*yl*hl+(yl-1)*hl)
                    extradata["ne_storage"][s]["inst_cost"][n]=mip(cst,extradata["ne_storage"][s]["inst_cost"][n+hl])
                end
            end
        end
        for (s,stg) in extradata["ne_storage"]
            for (n,cst) in enumerate(stg["cost_abs"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["ne_storage"][s]["cost_abs"][n]=npv_hourly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
        end
        for (s,stg) in extradata["ne_storage"]
            for (n,cst) in enumerate(stg["cost_inj"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["ne_storage"][s]["cost_inj"][n]=npv_hourly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
        end
    end
    if (haskey(extradata,"branchdc"))
        for (b,br) in extradata["branchdc"]
            for (n,cst) in enumerate(br["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["branchdc"][b]["cost"][n]=npv_yearly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
        end
    end
    if (haskey(extradata,"branch"))
        for (b,br) in extradata["branch"]
            for (n,cst) in enumerate(br["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["branch"][b]["cost"][n]=npv_yearly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
        end
    end
    if (haskey(extradata,"storage"))
        for (s,stg) in extradata["storage"]
            for (n,cst) in enumerate(stg["cost"])
                _sc=floor(Int64,(n-1)/(yl*hl))
                _yr=ceil(Int64,(n-_sc*(yl*hl))/(hl))
                extradata["storage"][s]["cost"][n]=npv_yearly(cst,parse(Int64,argz["scenario_years"][_yr]))
            end
        end
    end
    return extradata
end

#ensures binary candidates (array) costs sum to proper NPV value over the number of years
function npvs_costs_datas_4mip(data, scenario, _yrs, _dr)
    _scs=data["scenario"]
    _hrs=deepcopy(scenario["hours"])
    for (_sci,_sc) in _scs
        _sc=sort!(OrderedCollections.OrderedDict(_sc), by=x->parse(Int64,x))
        _sc_first=first(_sc)[2]
        for (_str,_num) in _sc
            if (_num<=_sc_first+(length(_yrs)-1)*_hrs-1)
                data["nw"][string(_num)]=deepcopy(npv_cost_data_4mip(deepcopy(data["nw"][string(_num)]),data["nw"][string(_num+_hrs)]))
            end
        end
    end
    return data
end

#ensures binary candidate (single) cost sum to proper NPV value over the number of years
function npv_cost_data_4mip(data0,data1)
    function mip(cost0, cost1)
        return cost0-cost1
    end

    for (b, branch) in get(data0, "ne_branch", Dict{String,Any}())
        data0["ne_branch"][b]["construction_cost"] = mip(data0["ne_branch"][b]["construction_cost"],data1["ne_branch"][b]["construction_cost"])
    end
    for (b, branch) in get(data0, "branchdc_ne", Dict{String,Any}())
        data0["branchdc_ne"][b]["cost"] = mip(data0["branchdc_ne"][b]["cost"],data1["branchdc_ne"][b]["cost"])
    end
    for (c, conv) in get(data0, "convdc_ne", Dict{String,Any}())
        data0["convdc_ne"][c]["cost"] = mip(data0["convdc_ne"][c]["cost"],data1["convdc_ne"][c]["cost"])
    end
    return data0
end

#returns vector of NPV investment limits per year
#***#
function max_invest_per_year(argz)
    max_invest=Float64[]
    for _yr in argz["scenario_years"]
        push!(max_invest,argz["yearly_investment"]/((1+argz["dr"])^(parse(Int64,_yr)-parse(Int64,argz["scenario_years"][1]))));end
    return max_invest
end
######################################################################
