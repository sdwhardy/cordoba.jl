# Import packages and create short names
import DataFrames; const _DF = DataFrames
using FileIO, JLD2
import CSV, Dates
import ExcelFiles; const _EF = ExcelFiles
import JuMP
import Gurobi
import Feather
import PowerModels; const _PM = PowerModels
import JSON
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0)

# Add auxiliary functions to construct input and scenario data dictionary
include("get_value.jl")
include("data.jl")
include("load_data.jl")
include("get_grid_data.jl")


Generation_Demand=Dict{}("Generation"=>Dict(),"Demand"=>Dict(),"ntcs"=>Dict(),"nodes"=>Dict())
push!(Generation_Demand["Generation"],"Scenarios"=>Dict())
gen_types=[]
gen_costs=Dict()
#Demand
for scenario in ["DE2030","DE2040","GA2030","GA2040","NT2025","NT2030","NT2040"]
    push!(Generation_Demand["Demand"],scenario=>_DF.DataFrame())
    ntcs, nodes, arcs, capacity, demand, gen_types, gen_costs, node_positions = get_grid_data(scenario)
    demand=demand[1:8760,:]
    if issubset([scenario],["NT2025","NT2030"])
    time_stamp=[Dates.DateTime(parse(Int64,scenario[3:end]),parse(Int64,t[1][2:end]),parse(Int64,t[2][2:end]),parse(Int64,t[3][2:end])) for t in split.(demand[!,:PATTERN],",")]
    else
    time_stamp=Dates.DateTime.(filter(x->!ismissing(x),demand[!,:YEAR]),filter(x->!ismissing(x),demand[!,:MONTH]),filter(x->!ismissing(x),demand[!,:DAY]),filter(x->!ismissing(x),demand[!,:Period]))
    end
    Generation_Demand["Demand"][scenario][!,:time_stamp]=time_stamp
    for zone in names(demand)[5:end]
        c=Float64.(filter(x->!ismissing(x),demand[!,Symbol(zone)]))
        if !(isempty(c))
        Generation_Demand["Demand"][scenario][!,Symbol(string(zone))]=c;end;end
        Generation_Demand["Demand"][scenario] =Generation_Demand["Demand"][scenario][_DF.completecases(Generation_Demand["Demand"][scenario]), :]
        _DF.disallowmissing!(Generation_Demand["Demand"][scenario])
end

for scenario in [("NT","National Trends"),("DE","Distributed Energy"),("GA","Global Ambition")]
    for year in [2025,2030,2040]
        if (year==2025 && !(first(scenario)=="NT"))
        else
            push!(Generation_Demand["Generation"]["Scenarios"],first(scenario)*string(year)=>Dict())
            ntcs, nodes, arcs, capacity, demand, gen_types, gen_costs, node_positions = get_grid_data(first(scenario)*string(year))
            push!(Generation_Demand["Generation"]["Scenarios"],first(scenario)*string(year)=>Dict())
            _C=filter!(Symbol("Scenario")=>x->x==last(scenario),capacity)
            _C=filter!(Symbol("Year")=>x->x==year,_C)
            _C=filter!(Symbol("Climate Year")=>x->x==2007,_C)
            zones=unique(_C[!,Symbol("Node/Line")])
            for zone in zones
                C=filter(Symbol("Node/Line")=>x->x==zone,_C)
                if !(isempty(C))
                    push!(Generation_Demand["Generation"]["Scenarios"][first(scenario)*string(year)],zone=>_DF.DataFrame("Generation_Type"=>C[!,Symbol("Generator_ID")],"Capacity"=>C[!,Symbol("Value")]))
end;end;end;end;end

#Find set of keys with full data
demand_keys=collect(names(Generation_Demand["Demand"]["NT2025"]))[2:end]
gen_keys=collect(keys(Generation_Demand["Generation"]["Scenarios"]["NT2025"]))
demand_scenarios=collect(keys(Generation_Demand["Demand"]))
for scenario in demand_scenarios;intersect!(demand_keys,names(Generation_Demand["Demand"][scenario]));end
gen_scenarios=collect(keys(Generation_Demand["Generation"]["Scenarios"]))
for scenario in gen_scenarios;intersect!(gen_keys,keys(Generation_Demand["Generation"]["Scenarios"][scenario]));end
intersect!(gen_keys,demand_keys)

#reduce data set
for scenario in gen_scenarios;
    for country in collect(keys(Generation_Demand["Generation"]["Scenarios"][scenario]));
    if issubset([country],gen_keys)
    else
        delete!(Generation_Demand["Generation"]["Scenarios"][scenario],country)
    end;end;end

for scenario in demand_scenarios;
    Generation_Demand["Demand"][scenario]=Generation_Demand["Demand"][scenario][!,Symbol.(vcat("time_stamp",gen_keys))]
end    


push!(Generation_Demand["Generation"],"RES"=>Dict())
push!(Generation_Demand["Generation"]["RES"],"Solar PV"=>Dict())
push!(Generation_Demand["Generation"]["RES"],"Offshore Wind"=>Dict())
push!(Generation_Demand["Generation"]["RES"],"Onshore Wind"=>Dict())
pv_all, wind_onshore_all, wind_offshore_all = load_res_data()

for country in unique(vcat(pv_all[!,:area], wind_onshore_all[!,:area], wind_offshore_all[!,:area]))
    pv=filter(Symbol("area")=>x->x==country,pv_all)
    wind_onshore=filter(Symbol("area")=>x->x==country,wind_onshore_all)
    wind_offshore=filter(Symbol("area")=>x->x==country,wind_offshore_all)
    push!(Generation_Demand["Generation"]["RES"]["Solar PV"],country=>Dict())
    push!(Generation_Demand["Generation"]["RES"]["Offshore Wind"],country=>Dict())
    push!(Generation_Demand["Generation"]["RES"]["Onshore Wind"],country=>Dict())
    for year in ["2012","2013","2014","2015","2016"]
        time_stamp_pv=Dates.DateTime.(parse(Float64,year),filter(x->!ismissing(x),pv[!,:month]),filter(x->!ismissing(x),pv[!,:day]),filter(x->!ismissing(x),pv[!,:hour]))
        time_stamp_on=Dates.DateTime.(parse(Float64,year),filter(x->!ismissing(x),wind_onshore[!,:month]),filter(x->!ismissing(x),wind_onshore[!,:day]),filter(x->!ismissing(x),wind_onshore[!,:hour]))
        time_stamp_off=Dates.DateTime.(parse(Float64,year),filter(x->!ismissing(x),wind_offshore[!,:month]),filter(x->!ismissing(x),wind_offshore[!,:day]),filter(x->!ismissing(x),wind_offshore[!,:hour]))
        _pv=Float64.(filter(x->!ismissing(x),pv[!,Symbol(year)]))
        _offwnd=Float64.(filter(x->!ismissing(x),wind_offshore[!,Symbol(year)]))
        _onwnd=Float64.(filter(x->!ismissing(x),wind_onshore[!,Symbol(year)]))
        if (length(time_stamp_pv)==length(_pv) && length(_pv)==8760);push!(Generation_Demand["Generation"]["RES"]["Solar PV"][country],year=>_DF.DataFrame("time_stamp"=>time_stamp_pv,country=>_pv));end
        if (length(time_stamp_off)==length(_offwnd)&&length(_offwnd)==8760);push!(Generation_Demand["Generation"]["RES"]["Offshore Wind"][country],year=>_DF.DataFrame("time_stamp"=>time_stamp_off,country=>_offwnd));end
        if (length(time_stamp_on)==length(_onwnd)&&length(_onwnd)==8760);push!(Generation_Demand["Generation"]["RES"]["Onshore Wind"][country],year=>_DF.DataFrame("time_stamp"=>time_stamp_on,country=>_onwnd));end
    end
end


ns=deepcopy(nodes[!,:node_id])
for _id in ns
    if (issubset([_id],gen_keys))
    else
        println(_id)
        nodes=nodes[(nodes.node_id .!= _id), :]
    end
end
delete!(Generation_Demand,"nodes")
push!(Generation_Demand["Generation"],"keys"=>gen_types)
push!(Generation_Demand["Generation"],"costs"=>gen_costs)
push!(Generation_Demand["Generation"],"ntcs"=>ntcs)
push!(Generation_Demand["Generation"],"nodes"=>nodes)

file = "C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_4EU.jld2"
FileIO.save(file, Generation_Demand)
datas=FileIO.load(file)

blank=deepcopy(datas["Generation"]["Scenarios"]["NT2025"]["BE00"][1:1,[:Generation_Type,:Capacity]])
push!(blank,["VOLL",0.0])
blank=blank[2:2,[:Generation_Type,:Capacity]]
for k in keys(datas["Generation"]["Scenarios"])
    push!(datas["Generation"]["Scenarios"][k],"BLNK"=>blank)
end

col_names=["time_stamp","BLNK"]
blank=deepcopy(datas["Generation"]["RES"]["Solar PV"]["BE00"])
for k in keys(blank)
    _DF.rename!(blank[k], Symbol.(col_names))
    blank[k][!,Symbol("BLNK")]=blank[k][!,Symbol("BLNK")].*0.0
end
k="Offshore Wind"
for k in keys(datas["Generation"]["RES"])
    push!(datas["Generation"]["RES"][k],"BLNK"=>blank)
end

zerodemand=datas["Demand"]["NT2025"][!,:BE00].*0
for k in keys(datas["Demand"])
    datas["Demand"][k][!,:BLNK]=zerodemand
end

save(file, datas)
########################################
##################### Cluster time series data ###########################
using FileIO, PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
scenario_data=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_4EU.jld2")

#divides time series into 24 hour groups
function daily_tss(ts)
    ts_daily=ts[1:24]
    for i=25:24:length(ts)
        ts_daily=hcat(ts_daily,ts[i:i+23])
    end
    return ts_daily
end
#count="BE";country=scenario_data["Generation"]["RES"]["Offshore Wind"][count]
#year="2014";ts=country[year]
daily_tss_data=Dict()
for (gen,res) in scenario_data["Generation"]["RES"]
    for (count,country) in res
            for (year,ts) in country
                if !(haskey(daily_tss_data,year));push!(daily_tss_data,year=>Dict());end
                if !(haskey(daily_tss_data[year],country));push!(daily_tss_data[year],"time_stamp"=>daily_tss(Dates.DateTime.(ts[!,"time_stamp"])));end
                if !(haskey(daily_tss_data[year],country));push!(daily_tss_data[year],"ts"=>daily_tss(Float64.(ts[!,count])).^2);
                else; daily_tss_data[year]["ts"]=daily_tss_data[year]["ts"].+daily_tss(Float64.(ts[!,count])).^2;end
            end
        end
    end
k=2;year="2014";tss=daily_tss_data[year]
yearly_cluster=Dict()
for k in 2:10
    for (year,tss) in daily_tss_data
        if !(haskey(yearly_cluster,year));push!(yearly_cluster,year=>Dict());end
        l2_norm=sqrt.(tss["ts"])
        kshape_clusters=ks.kshape(ks.zscore(sqrt.(l2_norm)',axis=1), k)
        ts=Vector{Float64}(); for clusters in last.(kshape_clusters); if (length(clusters)>0); ts=vcat(ts,tss["time_stamp"][:,rand(clusters)+1]); end;end
        if !(haskey(yearly_cluster[year],k));push!(yearly_cluster[year],k=>ts);end
    end
end

file = "C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\yearly_cluster_4EU.jld2"
save(file, yearly_cluster)

scenario_data["Generation"]["RES"]["GA2040"]["BLNK"]
#=yc=FileIO.load(file)
file = "C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\yearly_cluster_4UKFRBENLDEDKNO.jld2"
save(file, yc)=#

##################################
#=market=Dict()
market["base_years"]=Dict()
market["scenarios"]=Dict()
for year in ["2017","2018","2019","2020"]
dfs=[]
for country in ["BE","DK","UK","DE"]
market_data=CSV.read("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\time_series_data_w_BE\\Based_on_"*year*"_prices\\"*country*"_0.csv")
unique!(market_data, Symbol("MTU (CET)"))
_DF.rename!(market_data, Symbol("MTU (CET)")=>"time_stamp", Symbol("Day-ahead Price [EUR/MWh]")=>"EUR_da"*country)
market_data=market_data[:,1:2]
market_data[!,:time_stamp]=[Dates.DateTime(parse(Int64,ts[26:29]),parse(Int64,ts[23:24]),parse(Int64,ts[20:21]),parse(Int64,ts[31:32]),parse(Int64,ts[34:35])) for ts in market_data[!,:time_stamp]]

_DF.filter!(Symbol("EUR_da"*country) => x -> !(ismissing(x) || isnothing(x)), market_data)
_DF.filter!(Symbol("EUR_da"*country) => x -> x!="N/A", market_data)
_DF.filter!(Symbol("EUR_da"*country) => x -> x!="#VALUE!", market_data)
#if (country=="UK" && year =="2019")
if (typeof(market_data[!,Symbol("EUR_da"*country)][1])==typeof("String"))
    println(year)
    market_data[!,Symbol("EUR_da"*country)]=parse.(Float64,market_data[!,Symbol("EUR_da"*country)])
else
    println(country)
    println(year)
    market_data[!,Symbol("EUR_da"*country)]=Float64.(market_data[!,Symbol("EUR_da"*country)])
end


push!(dfs,market_data)
end
market_data=dfs[1];for df in dfs[2:end]; market_data=_DF.innerjoin(market_data,df,on=:time_stamp) ;end
push!(market["base_years"],year=>market_data);end

#ST->NT, DG->DE, GCA->GA
#BE, DE, DK, UK
#from Scenario Building 2018 Outputs.xlsx
costs=[[85.3,83.6,83.8,82.6],[46.0,45.2,46.0,42.1],[69.8,66.9,67.4,68.3],[65.7,65.3,66.8,64.9],[69.1,68.4,78.5,68.6],[50.9,50.6,50.5,53.1]]
market["scenarios"]=Dict()
for (i,scenario) in enumerate(["NT2030","NT2040","DE2030","DE2040","GA2030","GA2040"])
    push!(market["scenarios"],scenario=>Dict())
    for (j,country) in enumerate(["BE","DE","DK","UK"])
        push!(market["scenarios"][scenario],country=>costs[i][j])
    end
end




Generation_Demand["Market"]=market

wind=Dict()
for year in ["2017","2018","2019","2020"]
dfs=[]
for country in ["BE","DK","UK","DE"]
wind_data=CSV.read("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\time_series_data_w_BE\\Based_on_"*year*"_prices\\"*country*"_wnd.csv")
unique!(wind_data, Symbol("MTU"))
_DF.rename!(wind_data, Symbol("MTU")=>"time_stamp", Symbol("Wind Offshore  - Actual Aggregated [MW]")=>"Wnd_WWh"*country)
_DF.select!(wind_data, _DF.Not(Symbol("Area")))
wind_data[!,:time_stamp]=[Dates.DateTime(parse(Int64,ts[26:29]),parse(Int64,ts[23:24]),parse(Int64,ts[20:21]),parse(Int64,ts[31:32]),parse(Int64,ts[34:35])) for ts in wind_data[!,:time_stamp]]

_DF.filter!(Symbol("Wnd_WWh"*country) => x -> !(ismissing(x) || isnothing(x)), wind_data)
_DF.filter!(Symbol("Wnd_WWh"*country) => x -> x!="N/A", wind_data)
if (country=="UK" && year !="2020")
    println(year)
    wind_data[!,Symbol("Wnd_WWh"*country)]=parse.(Float64,wind_data[!,Symbol("Wnd_WWh"*country)])
else
    wind_data[!,Symbol("Wnd_WWh"*country)]=Float64.(wind_data[!,Symbol("Wnd_WWh"*country)])
end

push!(dfs,wind_data)
end
wind_data=dfs[1];for df in dfs[2:end]; wind_data=_DF.innerjoin(wind_data,df,on=:time_stamp) ;end
push!(wind,year=>wind_data);end
Generation_Demand["Wind"]=wind=#
