################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, OrderedCollections
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels

##################### Input parameters #################################
rt_ex=pwd()*"\\test\\data\\input\\UK_DE_DK\\"#folder path
argz = Dict(
"test"=>false,#if true smallest (2 hour) problem variation is built for testing
"scenario_planning_horizon"=>30,
"scenario_names"=>["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"],#Options: ["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
"k"=>2,#number of representative days modelled (24 hours per day)//Must add clustered time series for each k
"scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
"owpp_mva"=>[4000],#mva of wf in MVA
"conv_lim"=>4000,#Max Converter size in MVA
"candidate_ics_ac"=>[1/5],#AC Candidate Cable sizes (fraction of full MVA)
"candidate_ics_dc"=>[1,4/5,3/5,1/2],#DC Candidate Cable sizes (fraction of full MVA)
"dr"=>0.04,#discount rate
"yearly_investment"=>1000000)

################## optimization/solver setup options ###################
s = Dict("output" => Dict("branch_flows" => false),
"home_market"=>[],#nodes within Zonal market
"balancing_reserve"=>0.3,#zonal market must be defined to have any effect
"AC"=>"0",#0=false, 1=true
"eps"=>0.0001,#admm residual (100kW)
"relax_problem" => false,
"conv_losses_mp" => false,
"process_data_internally" => false,
"corridor_limit" => true,
"strg_lim_offshore"=>0.2,
"strg_lim_onshore"=>10,
"max_invest_per_year"=>_CBD.max_invest_per_year(argz))
########################################################################
scenario_data, ls = _CBD.load_time_series(rt_ex,argz)
#file = rt_ex*"topology.m"
#data = PowerModels.parse_file(file)
################## Run MIP Formulation ###################
mn_data, data, argz, s = main_ACDC_chandra(rt_ex,argz, s);#Build data structure for given options
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
#with OBZ
s["home_market"]=[]#nodes within Zonal market
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print solution


mn_data["nw"]["1"]
year1_keys=[];year2_keys=[];year3_keys=[]
hours_per_timeline=length(mn_data["scenario"]["1"])
2*hours_per_timeline/3
for _s in sort(OrderedCollections.OrderedDict(mn_data["scenario"]), by=x->parse(Int64,x))
    length(first(_s))
    for _k in sort(OrderedCollections.OrderedDict(last(_s)), by=x->parse(Int64,x))
        if (parse(Int64,first(_k))<=hours_per_timeline/3)
            push!(year1_keys, last(_k))
        elseif (parse(Int64,first(_k))<=2*hours_per_timeline/3)
            push!(year2_keys, last(_k))
        else
            push!(year3_keys, last(_k))
        end
    end
end
println(year3_keys)

benefits_y1=0
benefits_y2=0
benefits_y3=0
wbenefits_y1=0
wbenefits_y2=0
wbenefits_y3=0
w_y1=0
w_y2=0
w_y3=0

result_mip["solution"]["nw"]["1"]["gen"]
for (nw_i,nw) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]), by=x->parse(Int64,x))
   
    if (issubset([nw_i],string.(year1_keys))) 
        for (g_i,g) in nw["gen"]
            benefits_y1=benefits_y1+g["pg"]*mn_data["nw"][nw_i]["gen"][g_i]["cost"][1]
        end
        wbenefits_y1=wbenefits_y1+nw["gen"]["4"]["pg"]*minimum([mn_data["nw"][nw_i]["gen"]["1"]["cost"][1],mn_data["nw"][nw_i]["gen"]["2"]["cost"][1],mn_data["nw"][nw_i]["gen"]["3"]["cost"][1]])
        w_y1=w_y1+nw["gen"]["4"]["pg"]
    elseif (issubset([nw_i],string.(year2_keys)))
        for (g_i,g) in nw["gen"]
            benefits_y2=benefits_y2+g["pg"]*mn_data["nw"][nw_i]["gen"][g_i]["cost"][1]
        end
        wbenefits_y2=wbenefits_y2+nw["gen"]["4"]["pg"]*minimum([mn_data["nw"][nw_i]["gen"]["1"]["cost"][1],mn_data["nw"][nw_i]["gen"]["2"]["cost"][1],mn_data["nw"][nw_i]["gen"]["3"]["cost"][1]])
        w_y2=w_y2+nw["gen"]["4"]["pg"]
    else
        for (g_i,g) in nw["gen"]
            benefits_y3=benefits_y3+g["pg"]*mn_data["nw"][nw_i]["gen"][g_i]["cost"][1]
        end
        wbenefits_y3=wbenefits_y3+nw["gen"]["4"]["pg"]*minimum([mn_data["nw"][nw_i]["gen"]["1"]["cost"][1],mn_data["nw"][nw_i]["gen"]["2"]["cost"][1],mn_data["nw"][nw_i]["gen"]["3"]["cost"][1]])
        w_y3=w_y3+nw["gen"]["4"]["pg"]
    end 
end

benefits_y1/24
benefits_y2/12
benefits_y3/12
benefits_y1/24+benefits_y2/12+benefits_y3/12
w_y1/12*100*8760/48*5
w_y2/12*100*8760/48*10
w_y3/12*100*8760/48*10 

wbenefits_y1/24
wbenefits_y2/12
wbenefits_y3/12
wbenefits_y1/24+wbenefits_y2/12+wbenefits_y3/12
 
    
function main_ACDC_chandra(rt_ex,argz, s)
    ################# Load topology files ###################################
    _CBD.topology_df(rt_ex, s["relax_problem"], s["AC"])#creates .m file
    data, ics_ac, ics_dc, nodes = _CBD.filter_mfile_cables(rt_ex)#loads resulting topology and filters for candidate cables
    for (_is, _s) in data["storage"]; _s["cost"]=_s["cost"]*100;end
    
    #ics_dc=[(4000,550),(4000,760),(4000,250),(4000,443),(4000,212),(4000,311)]
    ############### defines size and market of genz and wfs ###################
    infinite_grid, genz, wfz, markets = _CBD.genz_n_wfs(argz["owpp_mva"],nodes,data["baseMVA"])
    push!(argz,"genz"=>genz)
    push!(argz,"wfz"=>wfz)
    #################### Calculates cable options for AC lines
    data=_CBD.AC_cable_options(data,argz["candidate_ics_ac"],ics_ac,data["baseMVA"])
    _CBD.print_topology_data_AC(data,markets)#print to verify
    #################### Calculates cable options for DC lines
    data=_CBD.DC_cable_options(data,argz["candidate_ics_dc"],ics_dc,data["baseMVA"])
    _CBD.additional_params_PMACDC(data)
    _CBD.print_topology_data_DC(data,markets)#print to verify
    ##################### load time series data ##############################
    scenario_data, ls = _CBD.load_time_series(rt_ex,argz)
    push!(argz,"ls"=>ls)
    ##################### multi period setup #################################
    mn_data = multi_period_setup(ls, scenario_data, data, markets, infinite_grid, argz)
    s=_CBD.update_settings(s, argz, data)
#   if (haskey(s,"home_market") && length(s["home_market"])>0)
#      mn_data=zonal_adjust(mn_data, s);end
    return  mn_data, data, argz, s
end

function multi_period_setup(ls,scenario_data,data, markets, infinite_grid, argz)
    #################### Multi-period input parameters #######################
    all_scenario_data,data,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,argz["scenario_years"],argz["scenario_names"],scenario_data,data)
    scenario["planning_horizon"] = argz["scenario_planning_horizon"] # in years, to scale generation cost
    extradata,data =_CBD.create_profile_sets_mesh(dim, data, all_scenario_data, markets, infinite_grid, [data["baseMVA"] for wf in argz["owpp_mva"]])
    #########################################################################
    #################### Scale cost data
    #[println(b*" "*string(br["construction_cost"])) for (b,br) in data["ne_branch"]];println()
    _CBD.scale_cost_data_2hourly!(data, scenario)#infrastructure investments
    #[println(b*" "*string(br["construction_cost"])) for (b,br) in data["ne_branch"]];println()
    _CBD.scale_cost_data_2yearly!(extradata, scenario)#energy cost benefits
    #[println(b*" "*string(br["construction_cost"])) for (b,br) in data["ne_branch"]];println()

    # Create data dictionary where time series data is included at the right place
    mn_data = _PMACDC.multinetwork_data(data, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))
    #[println(k*" "*b*" "*string(br["construction_cost"])) for (k,nw) in mn_data["nw"] for (b,br) in nw["ne_branch"]];println()
    # scale all to NPV
    #mn_data_mip= _CBD.npvs_costs_datas(mn_data_mip, scenario, scenario_years)#sum of years must equal total
   #mn_data = npvs_costs_datas_wREZ(mn_data, scenario, argz["scenario_years"], argz["dr"])#sum of years must equal total
    #[println(k*" "*b*" "*string(br["construction_cost"])) for (k,nw) in mn_data["nw"] for (b,br) in nw["ne_branch"]];println()
    mn_data = _CBD.npvs_costs_datas_4mip(mn_data, scenario, argz["scenario_years"], argz["dr"])#future investment at y scaled to year y=0
    return mn_data
end




########################################
################## Run Convex Formulation ################
s["relax_problem"]=true
mn_data, data, argz, s = _CBD.main_ACDC_wstrg(rt_ex,argz, s)#Build data structure for given options
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print convex solution
mn_data, data, s = _CBD.convex2mip(result_mip, data, mn_data, s)#Convert convex to candicate cables
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap" => 1e-4)#In large problems a larger MIPgap (1-5%?) may be desirable
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
_CBD.print_solution_data(result_mip, data, argz)#print MIP solution

################## Run ADMM Formulation ################
s["relax_problem"]=true
mn_data, data, argz, s = _CBD.main_ACDC_wstrg(rt_ex,argz, s)#Build data structure for given options
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 0)
result_mip = _CBD.admm_4_AjAwAgAuAo_main(mn_data, gurobi, s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print Convex (ADMM) solution
mn_data, data, s = _CBD.convex2mip(result_mip, data, mn_data, s)#Convert convex to candicate cables
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap" => 1e-4)#In large problems a larger MIPgap (1-5%?) may be desirable
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print MIP solution
