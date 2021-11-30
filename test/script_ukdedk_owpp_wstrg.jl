using Ipopt, Juniper, JuMP, Cbc, Gurobi, Cbc, XLSX, DataFrames, Dates, CSV, Distances
using MultivariateStats, Clustering
import cordoba; const _CBD = cordoba
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
import FlexPlan; const _FP = FlexPlan
import InfrastructureModels; const _IM = InfrastructureModels
include("../test/data/conv_spec.jl")
using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
#Conda.pip("install","dtw-python")
#dtw = pyimport_conda("dtw", "*")

#################### Multi-period INPUT PARAMETERS #######################
scenario_names=["STRG20"]
#scenario_names=["DG20"]

scenario_years=["2020"]
#scenario_years=["2020","2030"]
#scenario_years=["2020"]

scenario_planning_horizon=30

################### Input data for script_test.jl ########################
owpp_mva=[4000]#mva of wf (in de)
#owpp_mva=[0]#mva of wf (in de)

#interconnectors format: (mva,km)
#ics=[(4000,550),(4000,760),(4000,250),(4000,470),(4000,145),(4000,246)];
ics=[(4000,145)];
conv_lim=4000

#location of nodes
#markets_wfs=[["UK","DE","DK"],["DE"]]#must be in same order as .m file gens
markets_wfs=[["DE"],["DE"]]
infinite_grid=sum(first.(ics))+sum(owpp_mva)#ensures enough generation and consumption in all markets

#800MWh with ic_mva=500;owpp_mva=1000#mva,ic_length=600;#km,owpp_km=300#km,candidates=[[0.5,0.4],[0.5,0.4]],x=3 in storage_costs
##############################
rez=[]
for k=5:1:300
    casename = "ic_ukdedk_owpp_wstrg"
    file = "./test/data/input/$casename.m"
    data = _PM.parse_file(file)#load data in PM format
    candidate_ics=[1,4/5,3/5,1/2]#Candidate Cable sizes
    z_base_dc=(data["busdc"]["1"]["basekVdc"])^2/data["baseMVA"]
    data=_CBD.additional_candidatesICS(data,candidate_ics,ics)
    for (i,bdc) in data["branchdc_ne"]
    data["branchdc_ne"][i]=_CBD.candidateIC_cost_impedance(bdc,z_base_dc);end
    data["branchdc_ne"]=_CBD.unique_candidateIC(data["branchdc_ne"])#keep only unique candidates

    #delete!(data["branchdc_ne"],"4");delete!(data["branchdc_ne"],"3");delete!(data["branchdc_ne"],"2")
    #delete!(data["branchdc_ne"],"5");delete!(data["branchdc_ne"],"6")
    _FP.add_storage_data!(data) # Add addtional storage data model
    data=_CBD.storage_costs(data,2030)#adjust storage parameters - default: 2021 also 2030 included

    #for (i,bdc) in data["branchdc_ne"]#onshore to onshore connections
    #data["branchdc_ne"][i]=_CBD.candidateIC_cost(bdc);end

    #data["branchdc_ne"]=_CBD.unique_candidateIC(data["branchdc_ne"])#keep only unique candidates


    _PMACDC.process_additional_data!(data)#add extra DC model data
    ##################### Network cost/elec PARAMETERS ######################

    _CBD.converter_parameters_rxb(data)#adjust converter parameters
    ##########################################################################
    ##################### load time series data ##############################
    #k=5
    scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series
    #k=365
    ##################### Cluster time series data ###########################
    for (sc,yrs_ts) in scenario_data
        for (yr,ts) in yrs_ts
            daily_ts=_CBD.half_daily_tss(DateTime.(ts[!,"time_stamp"]),6)
            #daily_dk=_CBD.daily_tss(Float64.(ts[!,"Wnd_MWhDK"]))
            daily_de=_CBD.half_daily_tss(Float64.(ts[!,"EUR_daDE"]),6)
            #daily_uk=_CBD.daily_tss(Float64.(ts[!,"EUR_daUK"]))

            #=kshape_clusters_deuk=ks.kshape(ks.zscore(sqrt.((daily_de.^2))',axis=1), k)
            ts=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); ts=vcat(ts,daily_ts[:,rand(clusters)+1]); end;end
=#
            kmedoids_clusters=Clustering.kmedoids(Distances.pairwise(Distances.Euclidean(), daily_de, daily_de), k)
            ts=Vector{Float64}(); for cluster in kmedoids_clusters.medoids; ts=vcat(ts,daily_ts[:,cluster]); end

            filter!(row -> row.time_stamp in ts, scenario_data[sc][yr])
        end
    end


    #scenario_data=load("./test/data/input/EUSTDG_TS_k"*string(k)*".jld2")

    ##################### Find minimum length scenario
    #scenario_data=load("./test/data/input/EUSTDG_TS_k5.jld2")
    ls=[];for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
    #scenario_data[_sc][_yr]=scenario_data[_sc][_yr][1:2,:]
    push!(ls,length(scenario_data[_sc][_yr].time_stamp))
    end;end;ls=minimum(ls)

    ##################### Make all scenarios the same length
    for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
    scenario_data[_yr][_sc]=scenario_data[_yr][_sc][1:ls,:]
    end;end
    #set problem dimension
    #dim = scenario["hours"]

     # create a dictionary to pass time series data to data dictionary
     all_scenario_data,data,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data)
     scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost

    extradata,data = _CBD.create_profile_sets_mesh(dim, data, all_scenario_data, markets_wfs, infinite_grid, owpp_mva)
    extradata,data = _CBD.add_storage_profile(dim, data, extradata, all_scenario_data, markets_wfs, ls)

    #################### Scale cost data

    # Scale cost data
    _CBD.scale_cost_data_per_scenario!(data, scenario)
    _CBD.scale_cost_data_per_scenario!(extradata, scenario)


    # Create data dictionary where time series data is included at the right place
    mn_data = _PMACDC.multinetwork_data(data, extradata, Set{String}(["source_type","scenario","name", "source_version", "per_unit"]))
    [println(string(i)*", er: "*string(s["energy_rating"])*", mea: "*string(s["max_energy_absorption"])*", dr: "*string(s["discharge_rating"])*", cr: "*string(s["charge_rating"])*", ca: "*string(s["cost_abs"])*", ci: "*string(s["cost_inj"])) for (i,s) in mn_data["nw"]["1"]["ne_storage"]]

    #select solver
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)

    #settings
    s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false)

    #run optimization
    resultACDC =_CBD.cordoba_strg_tnep(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
    [println(i) for (i,res) in resultACDC["solution"]["nw"]["1"]["ne_storage"] if (res["isbuilt"]==1)]
    println();println()
    [println(i) for (i,res) in resultACDC["solution"]["nw"]["1"]["ne_storage"] if (res["isbuilt"]==1)]
    push!(rez,(k,resultACDC["objective"]))
    println("k: "*string(k))
end
-13360.189595985426
    #resultACDC = _PMACDC.run_tnepopf(data, _PM.DCPPowerModel, gurobi, setting = s)

    return resultACDC, mn_data
end
display(mn_data["gen"])

pu=100
e2me=1000000/pu
pwr=(8760*30*pu)/number_of_hours
profit_122=0
profit_522=0
profit_223=0
profit_321=0
wnd_energy=0
cbl_profit=0
bat_profit=0
for (i,n) in resultACDC["solution"]["nw"]
    wnd_energy=wnd_energy+n["gen"]["2"]["pg"]*(pwr)
    cbl_profit=cbl_profit+n["branchdc_ne"]["1"]["pt"]*mn_data["nw"][i]["gen"]["1"]["cost"][1]
    #profit_122=profit_125+abs(n["branchdc_ne"]["2"]["pt"])*(pwr)*abs(uk2dkwind_dta["nw"][i]["gen"]["1"]["cost"][1]-uk2dkwind_dta["nw"][i]["gen"]["2"]["cost"][1])*e2me
    #profit_522=profit_522+abs(n["branchdc_ne"]["14"]["pt"])*abs(uk2dkFull_dta["nw"][i]["gen"]["5"]["cost"][1]-uk2dkFull_dta["nw"][i]["gen"]["2"]["cost"][1])
    #profit_223=profit_223+abs(n["branchdc_ne"]["8"]["pt"])*abs(uk2dkFull_dta["nw"][i]["gen"]["2"]["cost"][1]-uk2dkFull_dta["nw"][i]["gen"]["3"]["cost"][1])
    #profit_321=profit_321+abs(n["branchdc_ne"]["4"]["pt"])*abs(uk2dkFull_dta["nw"][i]["gen"]["3"]["cost"][1]-uk2dkFull_dta["nw"][i]["gen"]["1"]["cost"][1])

    #profit_uk2dk=profit_uk2dk+abs(n["branchdc_ne"]["6"]["pt"])*abs(mn_data["nw"][i]["gen"]["1"]["cost"][1]-mn_data["nw"][i]["gen"]["3"]["cost"][1])
    #profit_de2dk=profit_de2dk+abs(n["branchdc_ne"]["10"]["pt"])*abs(mn_data["nw"][i]["gen"]["2"]["cost"][1]-mn_data["nw"][i]["gen"]["3"]["cost"][1])
    end
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



r=deepcopy(resultACDC)
m=deepcopy(mn_data)
d=deepcopy(data)
z=deepcopy(zs_data)
#Row  │ time_stamp          │ Wnd_MWhDE   │ EUR_daDE │ EUR_idDE │ MWh_upDE │ EUR_upDE │ MWh_dwnDE │ EUR_dwnDE

bats=Dict{String,Any}()
for (i,c) in r["solution"]["nw"]["1"]["ne_storage"]
    push!(bats,i=>DataFrame("sc_ne"=>[],"sd_ne"=>[],"e_abs_ne"=>[],"sc_price"=>[],"sd_price"=>[]))
end
for i= 1:length(r["solution"]["nw"])
    for (j,b) in r["solution"]["nw"][string(i)]["ne_storage"]
        push!(bats[j].sc_ne,b["sc_ne"])
        push!(bats[j].sd_ne,b["sd_ne"])
        push!(bats[j].e_abs_ne,b["se_ne"])
        push!(bats[j].sc_price,m["nw"][string(i)]["ne_storage"][j]["cost_abs"])
        push!(bats[j].sd_price,m["nw"][string(i)]["ne_storage"][j]["cost_inj"])
    end
end

mx_up=max.(bats["1"].sd_price,bats["2"].sd_price,bats["3"].sd_price,bats["4"].sd_price,bats["5"].sd_price,bats["6"].sd_price,bats["7"].sd_price,bats["8"].sd_price,bats["9"].sd_price)
mx_dwn=max.(bats["1"].sc_price,bats["2"].sc_price,bats["3"].sc_price,bats["4"].sc_price,bats["5"].sc_price,bats["6"].sc_price,bats["7"].sc_price,bats["8"].sc_price,bats["9"].sc_price)
mx_up_mwh=max.(bats["1"].sd_price,bats["2"].sd_price,bats["3"].sd_price,bats["4"].sd_price,bats["5"].sd_price,bats["6"].sd_price,bats["7"].sd_price,bats["8"].sd_price,bats["9"].sd_price)
mx_dwn_mwh=max.(bats["1"].sc_price,bats["2"].sc_price,bats["3"].sc_price,bats["4"].sc_price,bats["5"].sc_price,bats["6"].sc_price,bats["7"].sc_price,bats["8"].sc_price,bats["9"].sc_price)

up=sum.(bats["1"].sd_price,bats["2"].sd_price,bats["3"].sd_price,bats["4"].sd_price,bats["5"].sd_price,bats["6"].sd_price,bats["7"].sd_price,bats["8"].sd_price,bats["9"].sd_price)
dwn=sum([bats["1"].sc_ne,bats["2"].sc_ne,bats["3"].sc_ne,bats["4"].sc_ne,bats["5"].sc_ne,bats["6"].sc_ne,bats["7"].sc_ne,bats["8"].sc_ne,bats["9"].sc_ne])
up=sum([bats["1"].sd_ne,bats["2"].sd_ne,bats["3"].sd_ne,bats["4"].sd_ne,bats["5"].sd_ne,bats["6"].sd_ne,bats["7"].sd_ne,bats["8"].sd_ne,bats["9"].sd_ne])
sumry=DataFrame("ts"=>[i for i in 1:length(up)],"up"=>up,"dwn"=>dwn,"mx_up"=>mx_up,"mx_dwn"=>mx_dwn,"mwh_up"=>z["MWh_upDE"],"mwh_dwn"=>z["MWh_dwnDE"])
filter(sumry

dwns=filter(row -> row.dwn > 0, sumry)
bz=DataFrame("ts"=>[i for i in 1:length(z["MWh_upDE"])],"ffr"=>z["EUR_upDE"],"dc"=>bats["9"].sd_ne,"nt"=>bats["9"].e_abs_ne)
s=3110;e=3150
bz=bz[s:e,:]
ups=filter(row -> row.dc != 0, bz)

#cz=DataFrame("ts"=>[i for i in 1:length(z["MWh_dwnDE"])],"ffr"=>min.(z["EUR_dwnDE"],z["EUR_idDE"]),"dc"=>bats["9"].sc_ne,"nt"=>bats["9"].e_abs_ne)
cz=DataFrame("ts"=>[i for i in 1:length(z["MWh_dwnDE"])],"ffr"=>z["EUR_dwnDE"],"dc"=>bats["9"].sc_ne,"nt"=>bats["9"].e_abs_ne)

cz=cz[s:e,:]
dwns=filter(row -> row.dc != 0, cz)
using Plots
plotly()
width=750
height=500
#p=plot(size = (width, height),xaxis = ("Nodes", font(40, "Courier")),yaxis = ("Depth [m]", font(40, "Courier")))
p=plot(size = (width, height),markersize=2,seriestype=:scatter,xaxis = ("T.S.", font(20, "Courier")),yaxis = ("Norm Cost of Energy", font(20, "Courier")))
#p=plot(size = (width, height),xlims=(590,910),ylims=(-4,3),xaxis = ("T.S", font(20, "Courier")),yaxis = ("", font(20, "Courier")))
#plot!(p,1:1:length(mx_up),mx_up,color = :red,markersize=3, markershape = :circle, label="",size = (width, height))
plot!(p,first.(cs[1]),last.(cs[1]),color = :red,seriestype=:scatter,markersize=3, markershape = :circle, label="",size = (width, height))
plot!(p,first.(cs[2]),last.(cs[2]),color = :blue,seriestype=:scatter,markersize=3, markershape = :circle, label="",size = (width, height))
plot!(p,first.(cs[3]),last.(cs[3]),color = :black,seriestype=:scatter,markersize=3, markershape = :circle, label="",size = (width, height))
plot!(p,first.(cs[4]),last.(cs[4]),color = :green,seriestype=:scatter,markersize=3, markershape = :circle, label="",size = (width, height))
plot!(p,first.(cs[5]),last.(cs[5]),color = :yellow,seriestype=:scatter,markersize=3, markershape = :circle, label="",size = (width, height))

gui()
#plot!(p,dwns.ts,dwns.mwh_dwn,mx_dwn,color = :blue,markersize=3, markershape = :uptriangle, label="",size = (width, height))
#gui()

#plotly()
#width=750
#height=500
#p=plot(size = (width, height),xaxis = ("Nodes", font(40, "Courier")),yaxis = ("Depth [m]", font(40, "Courier")))
#p=plot(size = (width, height),markersize=2,seriestype=:scatter,xaxis = ("T.S.", font(20, "Courier")),yaxis = ("Norm Charge", font(20, "Courier")))
#p=plot(size = (width, height),xlims=(590,910),ylims=(-4,3),xaxis = ("T.S", font(20, "Courier")),yaxis = ("", font(20, "Courier")))
#plot!(p,1:1:length(mx_up),mx_up,color = :red,markersize=3, markershape = :circle, label="",size = (width, height))
plot!(p,(bz.ts),(((bz.nt)./maximum(bz.nt)).*2).-1,markersize=8,color = :black, label="",size = (width, height))


#plot!(p,dwns.ts,dwns.mwh_dwn,mx_dwn,color = :blue,markersize=3, markershape = :uptriangle, label="",size = (width, height))
gui()


data["scenario"] = Dict{String, Any}()
data["scenario_prob"] = Dict{String, Any}()
windgenprofile_beuk=[];gencostid_beuk=[];gencost_beuk=[];ic_nflow_losses=[];
zs_data=[]
for z in zs
    df=CSV.read("./test/data/input/"*z*"data.csv", DataFrames.DataFrame)
    colnames = ["time_stamp","Wnd_MWh"*z,"EUR_da"*z,"EUR_id"*z,"MWh_up"*z,"EUR_up"*z,"MWh_dwn"*z,"EUR_dwn"*z]
    names!(df, Symbol.(colnames))
    push!(zs_data,df)
end
zsd=zs_data[1];for z in zs_data[2:end];
zsd=innerjoin(zsd,z, makeunique=true,on=:time_stamp);end
Ytr=_CBD.PCA_cluster(zsd[2:7])
Ytr = convert(Matrix, zsd[1:1:end,2:end])
cluster=_CBD.kmeans_cluster(Ytr)
cs=[[],[],[],[],[],[],[]]

for (i,ass) in enumerate(cluster.assignments); push!(cs[ass],(eachrow(zsd)[i].Wnd_MWhDE,eachrow(zsd)[i].EUR_daDE));end
#cluster=kmedoid_cluster(Ytr)
n_samples=_CBD.n_samps(cluster,n)
zsd=zsd[n_samples,:];


using Plots
plotly()
width=750
height=500
p=plot(size = (width, height),yticks = 0.6:0.3:1.5,xlims=(0,300),ylims=(0.4,1.7),xaxis = ("Clusters", font(20, "Courier")),yaxis = ("Normalized Objective", font(20, "Courier")))
plot!(p,first.(rez),last.(rez)./(-13360.189595985426),seriestype=:scatter,color=:red, label="Kmedoids",size = (width, height))
plot!(p,first.(rez),[1 for i in last.(rez)],color=:black,size = (width, height))
gui()


rez300=deepcopy(rez)
