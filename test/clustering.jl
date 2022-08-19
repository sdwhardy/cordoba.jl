using Gurobi, JuMP, DataFrames, FileIO, JLD2, Dates, OrderedCollections, CSV
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels


scenario_data=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBEDEDK.jld2")
ts=scenario_data["Generation"]["RES"]["Offshore Wind"]["DE"]["2012"]
################### ENTSO-E scenario description ####################################
#scenario_names=["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
scenario_names=["EU17","EU18","EU19","EU20"]
#scenario_names=["EU19"]

#scenario_years=["2020","2030","2040"]
#scenario_years=["2020","2030"]
scenario_years=["2020"]

scenario_planning_horizon=30

################### Input data for script_test.jl ########################
owpp_mva=[4000]#mva of wf (in de)
#owpp_mva=[0]#mva of wf (in de)

#interconnectors format: (mva,km)
ics=[(4000,550),(4000,760),(4000,250),(4000,470),(4000,145),(4000,246)];
#ics=[(2000,550),(2000,760),(2000,250),(0,470),(0,145),(0,246)];
#ics=[(4000,145)];
conv_lim=4000

#location of nodes
markets_wfs=[["UK","DE","DK"],["DE"]]#must be in same order as .m file gens
#markets_wfs=[["DE"],["DE"]]
infinite_grid=sum(first.(ics))+sum(owpp_mva)#ensures enough generation and consumption in all markets

##################### load time series data ##############################
k=2
scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series

##################### Cluster time series data ###########################

for (sc,yrs_ts) in scenario_data
    for (yr,ts) in yrs_ts
        println(sc*" "*yr)

        filter!(row -> ismissing(row.time_stamp)==false, scenario_data[sc][yr])
        filter!(row -> ismissing(row.Wnd_MWhDE)==false, scenario_data[sc][yr])
        filter!(row -> ismissing(row.EUR_daDE)==false, scenario_data[sc][yr])
        filter!(row -> ismissing(row.EUR_daUK)==false, scenario_data[sc][yr])
        filter!(row -> ismissing(row.EUR_daDK)==false, scenario_data[sc][yr])
        daily_ts=_CBD.daily_tss(DateTime.(ts[!,"time_stamp"]))
        daily_dew=_CBD.daily_tss(Float64.(ts[!,"Wnd_MWhDE"]))
        daily_de=_CBD.daily_tss(Float64.(ts[!,"EUR_daDE"]))
        daily_uk=_CBD.daily_tss(Float64.(ts[!,"EUR_daUK"]))
        daily_dk=_CBD.daily_tss(Float64.(ts[!,"EUR_daDK"]))

        kshape_clusters_deuk=ks.kshape(ks.zscore(sqrt.((daily_uk.^2).+(daily_de.^2).+(daily_dk.^2).+(daily_dew.^2))',axis=1), k)
        ts=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); ts=vcat(ts,daily_ts[:,rand(clusters)+1]); end;end
        filter!(row -> row.time_stamp in ts, scenario_data[sc][yr])
    end
end
save("./test/data/input/EUSTDG17t020_TS_k"*string(k)*".jld2",scenario_data)=#
scenario_data=load("./test/data/input/EUSTDG17t020_TS_k"*string(k)*".jld2")
d_keys=keys(scenario_data);for k in d_keys;if !(issubset([string(k)],scenario_names));delete!(scenario_data,k);else;y_keys=keys(scenario_data[k]);for y in y_keys;if !(issubset([string(y)],scenario_years));delete!(scenario_data[k],y);end; end;end;end



#divides time series into 24 hour groups
function daily_tss(ts)
    ts_daily=ts[1:24]
    for i=25:24:length(ts)-24
        ts_daily=hcat(ts_daily,ts[i:i+23])
    end
    return ts_daily
end

Ytr = convert(Matrix, zsd[1:1:end,2:end])'
    cluster=kmeans_cluster(Ytr)