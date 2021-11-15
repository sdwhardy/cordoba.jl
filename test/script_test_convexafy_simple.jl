################## loads external packages ##############################
using Ipopt, Juniper, JuMP, Cbc, Gurobi, Cbc, XLSX, DataFrames, Dates, CSV, FileIO
import cordoba; const _CBD = cordoba#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
import FlexPlan; const _FP = FlexPlan
import InfrastructureModels; const _IM = InfrastructureModels
include("../test/data/conv_spec.jl")
rez=[]
for k=10:10:250
    ################### ENTSO-E scenario ####################################
    scenario_names=["EU20"]
    scenario_years=["2020"]
    ################### Input data for script_test.jl ########################
    owpp_mva=4000#mva of wf (in de)
    #interconnectors format: (mva,(km,type)) type: -1:both ends onshore, 0:both ends offshore, 1:onshore/offshore
    #must be in same order as .m file ie uk(1)-de(2), uk(1)-dk(3), de(2)-dk(3), wf(4)-uk(1), wf(4)-de(2), wf(4)-dk(3)
    ics=[(2000,(550,-1))];
    #location of nodes
    markets=["UK","DE"]#must be in same order as .m file gens
    #candidate_ics=[1,4/5,3/5,1/2]#Candidate Cables
    infinite_grid=sum(first.(ics))+owpp_mva#ensures enough generation and consumption in all markets

    ################## reads .m input file name ######################
    casename = "test_convexafy_simple"
    file = "./test/data/input/$casename.m"
    data = _PM.parse_file(file)#load data in PM format
    _PMACDC.process_additional_data!(data)#add extra DC model data

    ##################### load time series data ##############################
    #scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series
    scenario_data = _CBD.get_scenario_year_tss_convexafy(scenario_names,scenario_years)#Retrieve the scenario time series

    daily_ts=daily_tss(DateTime.(scenario_data["EU20"]["2020"]["time_stamp"]))
    daily_wnd=daily_tss(Float64.(scenario_data["EU20"]["2020"]["Wnd_MWhDE"]))
    daily_de=daily_tss(Float64.(scenario_data["EU20"]["2020"]["EUR_daDE"]))
    daily_uk=daily_tss(Float64.(scenario_data["EU20"]["2020"]["EUR_daUK"]))


    #=kmedoids_clusters_deuk=kmedoids(pairwise(Euclidean(), sqrt.((daily_uk.^2).+(daily_de.^2)), sqrt.((daily_uk.^2).+(daily_de.^2))), k)
    #kmedoids_clusters_wnd=kmedoids(pairwise(Euclidean(), daily_wnd, daily_wnd), k)
    wnd=Vector{Float64}(); for medoid in kmedoids_clusters_deuk.medoids; wnd=vcat(wnd,daily_wnd[:,medoid]); end
    #kmedoids_clusters_de=kmedoids(pairwise(Euclidean(), daily_de, daily_de), k)
    de=Vector{Float64}(); for medoid in kmedoids_clusters_deuk.medoids; de=vcat(de,daily_de[:,medoid]); end
    #kmedoids_clusters_uk=kmedoids(pairwise(Euclidean(), daily_uk, daily_uk), k)
    uk=Vector{Float64}(); for medoid in kmedoids_clusters_deuk.medoids; uk=vcat(uk,daily_uk[:,medoid]); end=#



    kshape_clusters_deuk=ks.kshape(ks.zscore(sqrt.((daily_uk.^2).+(daily_de.^2))',axis=1), k)
    ts=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); ts=vcat(ts,daily_ts[:,rand(clusters)+1]); end;end
    #=wnd=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); wnd=vcat(wnd,daily_wnd[:,rand(clusters)+1]); end;end
    de=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); de=vcat(de,daily_de[:,rand(clusters)+1]); end;end
    uk=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); uk=vcat(uk,daily_uk[:,rand(clusters)+1]); end;end=#

    #=kshape_clusters_indexes=ks.kshape(ks.zscore(sqrt.((daily_uk.^2).+(daily_de.^2))',axis=1), k)
    kmedoids_clusters=[]
    for k_shape_cluster_indexes in  last.(kshape_clusters_indexes)
        kshape_clusters=Vector{Float64}();
        for d in k_shape_cluster_indexes;
            if (length(kshape_clusters)==0)
                kshape_clusters=daily_wnd[:,d+1]
            else
            kshape_clusters=hcat(kshape_clusters,daily_wnd[:,d+1]); end
        end
        try
            D=pairwise(Euclidean(), kshape_clusters, kshape_clusters)
            push!(kmedoids_clusters,kmedoids(D, 1).medoids[1])
        catch
            println("error caught")
            println(length(k_shape_cluster_indexes))
        end
    end
    wnd=Vector{Float64}(); for medoid in kmedoids_clusters; wnd=vcat(wnd,daily_wnd[:,medoid]); end
    de=Vector{Float64}(); for medoid in kmedoids_clusters; de=vcat(de,daily_de[:,medoid]); end
    uk=Vector{Float64}(); for medoid in kmedoids_clusters; uk=vcat(uk,daily_uk[:,medoid]); end=#

    #=kmeans_clusters_wnd=kmeans(daily_wnd, k)
    wnd=Vector{Float64}(); for mean in kmeans_clusters_wnd.centers; wnd=vcat(wnd,mean); end
    kmeans_clusters_de=kmeans(daily_de, k)
    de=Vector{Float64}(); for mean in kmeans_clusters_de.centers; de=vcat(de,mean); end
    kmeans_clusters_uk=kmeans(daily_uk, k)
    uk=Vector{Float64}(); for mean in kmeans_clusters_uk.centers; uk=vcat(uk,mean); end=#

    #filter!(row -> row.time_stamp in [hour for day in days for hour in day], scenario_data["EU20"]["2020"])

    filter!(row -> row.time_stamp in ts, scenario_data["EU20"]["2020"])

    #scenario_data_old=deepcopy(scenario_data)
    ##################### Sample scenarios
    ls=[];for (_sc, data_by_sc) in scenario_data; for (_yr, data_by_yr) in data_by_sc;
    #scenario_data[_sc][_yr]=_CBD.set_of_hours(data_by_yr,[i for i=1:1:24],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31],[1,2,3,4,5,6,7,8,9,10,11,12])
    #scenario_data[_sc][_yr]=_CBD.set_of_hours(data_by_yr,[i for i=1:12:24],[2],[1])
    push!(ls,length(scenario_data[_sc][_yr].time_stamp))
    end;end

    ##################### Make all scenarios the same length
    for (_sc, data_by_sc) in scenario_data; for (_yr, data_by_yr) in data_by_sc;
    scenario_data[_sc][_yr]=scenario_data[_sc][_yr][1:minimum(ls),:]
    end;end

    #################### Multi-period input parameters #######################
    #n=number_of_hours = 8760 # Number of time points in DE
    #all_scenario_data,data,scenario, dim = _CBD.multi_period_setup(minimum(ls),scenario_years,scenario_names,scenario_data,data)
    all_scenario_data,data,scenario, dim = _CBD.multi_period_stoch_year_setup(minimum(ls),scenario_years,scenario_names,scenario_data,data)

    scenario["planning_horizon"] = 30/length(scenario_years) # in years, to scale generation cost

    extradata,data = _CBD.create_profile_sets_mesh(dim, data, all_scenario_data, markets, infinite_grid, owpp_mva)

    # Scale cost data
    _CBD.scale_cost_data_cordoba_convexafy!(data, scenario)
    _CBD.scale_cost_data_cordoba!(extradata, scenario)

    # Create data dictionary where time series data is included at the right place
    mn_data = _PMACDC.multinetwork_data(data, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "per_unit", "source_version"]))

    #select solver
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)

    #optimation settings
    scenes_years=[];for (i,(k,_scene)) in enumerate(scenario["sc_names"]);push!(scenes_years,[]);for (j, _yr) in _scene; push!(scenes_years[i],_yr);end;end
    s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => minimum(ls))

    #run optimization
    #result = _CBD.cordoba_mp_acdctnepopf_stoch(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
    result = _CBD.cordoba_mp_acdctnepopf_convexafy(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
    push!(rez,(k,result["objective"]))
    println("k: "*string(k))
end
result["solution"]["nw"]["1"]
result["solution"]["nw"]["2"]

for (k,n) in result["solution"]["nw"]
    println(k*": pg: "*string(n["gen"]["4"]["pg"]))
end


result["solution"]["nw"]["1"]
result["solution"]["nw"]["50160"]

#save("../../package_results/cordoba/result_ukdedk_stochFULL.jld2","result",result)
#save("../../package_results/cordoba/mn_data_ukdedk_stochFULL.jld2","mn_data",mn_data)

#resultACDC = _CBD.cordoba_mp_acdctnepopf(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)

profit_1to2=0
profit_1to3=0
profit_2to3=0
profit_1to4=0
profit_2to4=0
profit_3to4=0
wnd_energy=0
wnd_profit=0
capex=0
mn_data["nw"]["1"]["gen"]
for (i,n) in result["solution"]["nw"]
    #if (parse(Int,i)>5*8360 && parse(Int,i)<6*8360)
        wnd_energy=wnd_energy+n["gen"]["4"]["pg"]
        wnd_profit=wnd_profit+n["gen"]["4"]["pg"]*max(mn_data["nw"][i]["gen"]["1"]["cost"][1],mn_data["nw"][i]["gen"]["2"]["cost"][1],mn_data["nw"][i]["gen"]["3"]["cost"][1])
        #profit_1to2=profit_1to2+n["branchdc_ne"]["3"]["pt"]*(mn_data["nw"][i]["gen"]["2"]["cost"][1]-mn_data["nw"][i]["gen"]["1"]["cost"][1])
        profit_1to3=profit_1to3+n["branchdc_ne"]["5"]["pt"]*(mn_data["nw"][i]["gen"]["3"]["cost"][1]-mn_data["nw"][i]["gen"]["1"]["cost"][1])
        profit_2to3=profit_2to3+n["branchdc_ne"]["9"]["pt"]*(mn_data["nw"][i]["gen"]["3"]["cost"][1]-mn_data["nw"][i]["gen"]["2"]["cost"][1])
        profit_1to4=profit_1to4+n["branchdc_ne"]["13"]["pt"]*(mn_data["nw"][i]["gen"]["1"]["cost"][1]-mn_data["nw"][i]["gen"]["4"]["cost"][1])
        profit_2to4=profit_2to4+n["branchdc_ne"]["17"]["pt"]*(mn_data["nw"][i]["gen"]["2"]["cost"][1]-mn_data["nw"][i]["gen"]["4"]["cost"][1])
        #profit_3to4=profit_3to4+n["branchdc_ne"]["23"]["pt"]*(mn_data["nw"][i]["gen"]["3"]["cost"][1]-mn_data["nw"][i]["gen"]["4"]["cost"][1])
        capex=capex+(mn_data["nw"][i]["branchdc_ne"]["5"]["cost"]+mn_data["nw"][i]["branchdc_ne"]["9"]["cost"]+mn_data["nw"][i]["branchdc_ne"]["13"]["cost"]+mn_data["nw"][i]["branchdc_ne"]["17"]["cost"])#+mn_data["nw"][i]["branchdc_ne"]["23"]["cost"])
    #end
end
println(string(profit_1to2)*" "*
string(profit_1to3)*" "*
string(profit_2to3)*" "*
string(profit_1to4)*" "*
string(profit_2to4)*" "*
string(profit_3to4)*" "*
string(capex)*" "*
string(result["objective"])*" "*
string(wnd_energy)*" "*
string(wnd_profit)*" "*
string(8360)*" "*
string(scenario["planning_horizon"]))




uk2dkwind=deepcopy(resultACDC)
uk2dkwind_dta=deepcopy(mn_data)


for (key,nw) in resultACDC["solution"]["nw"]
    for (key2,b) in nw["ne_storage"]
    #println(key*" -  "*key2*" dr:  "*string(b["sd_ne"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["discharge_rating"]))
    #println(key*" -  "*key2*" cr: "*string(b["sc_ne"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["charge_rating"]))
    #println(key*" -  "*key2*" cr: "*string(b["sc_ne"]*mn_data["nw"][key]["ne_storage"][key2]["cost_abs"])*" "*" dr:  "*string(b["sd_ne"]*mn_data["nw"][key]["ne_storage"][key2]["cost_inj"]))
    println(key*" -  "*key2*" cst: "*string(mn_data["nw"][key]["ne_storage"][key2]["inst_cost"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["eq_cost"]))
    #println(key*" -  "*key2*" cst: "*string(mn_data["nw"][key]["ne_storage"][key2]["cost_abs"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["cost_inj"]))
end;end



plotly()
width=750
height=500
p=plot(size = (width, height),xaxis = ("Clusters", font(20, "Courier")),yaxis = ("Normalized Objective", font(20, "Courier")))
plot!(p,first.(rez),last.(rez)./(-14517.108066419452),seriestype=:scatter,color=:red, label="Kmedoids",size = (width, height))
plot!(p,first.(rez),[1 for i in last.(rez)],color=:black,size = (width, height))
gui()

res_full_med=deepcopy(result["objective"])
rez_kmedoid=deepcopy(rez)

res_full_med_wnd=deepcopy(result["objective"])
rez_kmedoid_wnd=deepcopy(rez)

res_full_means_wnd=deepcopy(result["objective"])
rez_kmeans_wnd=deepcopy(rez)

res_full_means_de=deepcopy(result["objective"])
rez_kmeans_de=deepcopy(rez)

res_full_means_uk=deepcopy(result["objective"])
rez_kmeans_uk=deepcopy(rez)

res_full_means_beuk=deepcopy(result["objective"])
rez_kmeans_beuk=deepcopy(rez)
