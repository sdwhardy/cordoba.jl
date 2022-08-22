################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, OrderedCollections, FileIO
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels

##################### Input parameters #################################
rt_ex=pwd()*"\\test\\data\\input\\UK_BE_DK4000\\"#folder path
argz = Dict(
"test"=>false,#if true smallest (2 hour) problem variation is built for testing
"scenario_planning_horizon"=>30,
"scenario_names"=>["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"],#Options: ["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
"k"=>5,#number of representative days modelled (24 hours per day)//Must add clustered time series for each k
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
#scenario_data, ls = _CBD.load_time_series(rt_ex,argz)


#file = rt_ex*"topology.m"
#data = PowerModels.parse_file(file)
################## Run MIP Formulation ###################
mn_data, data, argz, s = _CBD.main_ACDC_chandra(rt_ex,argz, s);#Build data structure for given options

result_mip=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\chandra\\result_ukdedk_k54000.jld2")
_CBD.print_solution_data(result_mip, data, argz)#print solution


year1_keys=[];year2_keys=[];year3_keys=[]
hours_per_timeline=length(mn_data["scenario"]["1"])

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

benefits_y1=0
benefits_y2=0
benefits_y3=0
wbenefits_y1=0
wbenefits_y2=0
wbenefits_y3=0
w_y1=0
w_y2=0
w_y3=0

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
w_y1/12*100*8760/(24*argz["k"])*5
w_y2/12*100*8760/(24*argz["k"])*10
w_y3/12*100*8760/(24*argz["k"])*10 
wbenefits_y1/24
wbenefits_y2/12
wbenefits_y3/12
wbenefits_y1/24+wbenefits_y2/12+wbenefits_y3/12
println(string(-1*benefits_y1/24)*" "*string(wbenefits_y1/24)*" "*string(w_y1/12*100*8760/(24*argz["k"])*5)) 
println(string(-1*benefits_y2/12)*" "*string(wbenefits_y2/12)*" "*string(w_y2/12*100*8760/(24*argz["k"])*5)) 
println(string(-1*benefits_y3/12)*" "*string(wbenefits_y3/12)*" "*string(w_y3/12*100*8760/(24*argz["k"])*5)) 
