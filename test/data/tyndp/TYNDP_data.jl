# Import packages and create short names
include("C:\\Users\\shardy\\Documents\\julia\\packages\\AC-DC-CBA\\su-s_data.jl")


Generation_Demand=Dict{}("Generation"=>Dict(),"Demand"=>Dict())
push!(Generation_Demand["Generation"],"Scenarios"=>Dict())
gen_types=[]
gen_costs=Dict()
ntcs=_DF.DataFrame() 
nodes=_DF.DataFrame()
#Demand
for scenario in ["DE2030","DE2040","GA2030","GA2040","NT2025","NT2030","NT2040"]
    push!(Generation_Demand["Demand"],scenario=>_DF.DataFrame())
    push!(Generation_Demand["Generation"]["Scenarios"],scenario=>Dict())
    ntcs, nodes, arcs, capacity, demand, gen_types, gen_costs, node_positions = get_grid_data(scenario)
    
    #demand setup
    demand=demand[1:8760,:]
    if issubset([scenario],["NT2025","NT2030"])
    time_stamp=[Dates.DateTime(parse(Int64,scenario[3:end]),parse(Int64,t[1][2:end]),parse(Int64,t[2][2:end]),parse(Int64,t[3][2:end])) for t in split.(demand[!,:PATTERN],",")]
    demand=demand[:,2:end]
    else
    time_stamp=Dates.DateTime.(filter(x->!ismissing(x),demand[!,:YEAR]),filter(x->!ismissing(x),demand[!,:MONTH]),filter(x->!ismissing(x),demand[!,:DAY]),filter(x->!ismissing(x),demand[!,:Period]))
    demand=demand[:,5:end]
    end
    Generation_Demand["Demand"][scenario][!,:time_stamp]=time_stamp
    for zone in names(demand)
        c=Float64.(filter(x->!ismissing(x),demand[!,Symbol(zone)]))
        if !(isempty(c))
        Generation_Demand["Demand"][scenario][!,Symbol(string(zone))]=c;
        end;
    end
    Generation_Demand["Demand"][scenario] =Generation_Demand["Demand"][scenario][_DF.completecases(Generation_Demand["Demand"][scenario]), :]
    _DF.disallowmissing!(Generation_Demand["Demand"][scenario])

    #generation setup 
    scenario_full = "National Trends"
    scenario_full = scenario[1:2]=="GA" ? "Global Ambition" : scenario_full
    scenario_full = scenario[1:2]=="DE" ? "Distributed Energy" : scenario_full
    _C=filter!(Symbol("Scenario")=>x->x==scenario_full,capacity)
    _C=filter!(Symbol("Year")=>x->x==parse(Float64,scenario[3:end]),_C)
    _C=filter!(Symbol("Climate Year")=>x->x==2007,_C)
    zones=unique(_C[!,Symbol("Node/Line")])
    for zone in zones
        C=filter(Symbol("Node/Line")=>x->x==zone,_C)
        if !(isempty(C))
            push!(Generation_Demand["Generation"]["Scenarios"][scenario],zone=>_DF.DataFrame("Generation_Type"=>C[!,Symbol("Generator_ID")],"Capacity"=>C[!,Symbol("Value")]))
        end;
    end;  
end

#Find set of keys with full data
demand_keys=collect(names(Generation_Demand["Demand"]["NT2025"]))[2:end]
gen_keys=collect(keys(Generation_Demand["Generation"]["Scenarios"]["NT2025"]))

demand_scenarios=collect(keys(Generation_Demand["Demand"]))
for scenario in demand_scenarios;intersect!(demand_keys,names(Generation_Demand["Demand"][scenario]));end
gen_scenarios=collect(keys(Generation_Demand["Generation"]["Scenarios"]))
for scenario in gen_scenarios;intersect!(gen_keys,keys(Generation_Demand["Generation"]["Scenarios"][scenario]));end
intersect!(gen_keys,demand_keys)
intersect!(demand_keys,gen_keys)

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

push!(Generation_Demand["Generation"],"keys"=>gen_types)
push!(Generation_Demand["Generation"],"costs"=>gen_costs)
push!(Generation_Demand["Generation"],"ntcs"=>ntcs)
push!(Generation_Demand["Generation"],"nodes"=>nodes)

file = "C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_4EU.jld2"
FileIO.save(file, Generation_Demand)

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
for k in 1:1
    for (year,tss) in daily_tss_data
        if !(haskey(yearly_cluster,year));push!(yearly_cluster,year=>Dict());end
        l2_norm=sqrt.(tss["ts"])
        kshape_clusters=ks.kshape(ks.zscore(sqrt.(l2_norm)',axis=1), k)
        ts=Vector{Float64}(); for clusters in last.(kshape_clusters); if (length(clusters)>0); ts=vcat(ts,tss["time_stamp"][:,rand(clusters)+1]); end;end
        if !(haskey(yearly_cluster[year],k));push!(yearly_cluster[year],k=>ts);end
    end
end

file = "C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\yearly_cluster_4EU.jld2"

save(file, scenario_data)


scenario_data=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\yearly_cluster_4EU.jld2")
for (k_sc,_sc) in scenario_data
    scenario_data[k_sc][1]=yearly_cluster[k_sc][1]
end

scenario_data["2012"][1]