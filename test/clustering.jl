#=using Gurobi, JuMP, DataFrames, FileIO, JLD2, Dates, OrderedCollections, CSV, Clustering, Distances
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels


scenario_data=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBEDEDK.jld2")

k_exmaple=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\yearly_cluster_4UKBEDEDK.jld2")
    
bins=bin_khours(scenario_data)
##################### Cluster time series data ###########################

for res_yr in ["2012","2013","2014","2015","2016"]
    #for cunt in ["UK","BE","DE","DK"]
    tsde=scenario_data["Generation"]["RES"]["Offshore Wind"]["DE"][res_yr]
    tsdk=scenario_data["Generation"]["RES"]["Offshore Wind"]["DK"][res_yr]
    tsbe=scenario_data["Generation"]["RES"]["Offshore Wind"]["BE"][res_yr]
    daily_dew=daily_tss(Float64.(tsde[!,"DE_MWh"]))
    daily_dkw=daily_tss(Float64.(tsdk[!,"DK_MWh"]))
    daily_bew=daily_tss(Float64.(tsbe[!,"BE_MWh"]))
    ave_de=sum(daily_dew)/length(daily_dew)
    max_de=maximum(daily_dew)
    ave_dk=sum(daily_dkw)/length(daily_dkw)
    max_dk=maximum(daily_dkw)
    ave_be=sum(daily_bew)/length(daily_bew)
    max_be=maximum(daily_bew)
    for k=2:1:10
        for (k_bs, bs) bins[res_yr][string(k)]
            for (k_b,b) in bs
function bin_khours(scenario_data) 
    bins=Dict()
    for k=2:1:10
        for res_yr in ["2012","2013","2014","2015","2016"]
            #for cunt in ["UK","BE","DE","DK"]
            tsde=scenario_data["Generation"]["RES"]["Offshore Wind"]["DE"][res_yr]
            tsdk=scenario_data["Generation"]["RES"]["Offshore Wind"]["DK"][res_yr]
            tsbe=scenario_data["Generation"]["RES"]["Offshore Wind"]["BE"][res_yr]
            if !(haskey(bins,res_yr));push!(bins,res_yr=>Dict());end
            if !(haskey(bins[res_yr],string(k)));push!(bins[res_yr],string(k)=>Dict());end
            #if !(haskey(bins[res_yr],cunt));push!(bins[res_yr],cunt=>Dict());end
            
            push!(bins,"ts"=>daily_tss(DateTime.(tsde[!,"time_stamp"])))

            daily_dew=daily_tss(Float64.(tsde[!,"DE_MWh"]))
            daily_dkw=daily_tss(Float64.(tsdk[!,"DK_MWh"]))
            daily_bew=daily_tss(Float64.(tsbe[!,"BE_MWh"]))
            daily_w=sqrt.(daily_dew.^2+daily_dkw.^2+daily_bew.^2)
            D = pairwise(Euclidean(), daily_w, daily_w)
            rez=kmedoids(D, k)
            
            for (a,ass) in enumerate(rez.assignments)
            for (b,in_b) in enumerate(rez.counts)
                if (ass==b)
                    if !(haskey(bins[res_yr][string(k)],string(b)))
                        push!(bins[res_yr][string(k)],string(b)=>[a])
                    else
                        push!(bins[res_yr][string(k)][string(b)],a)
                    end
                end
            end;end
        end
    end
    return bins
end
for res_yr in ["2012","2013","2014","2015","2016"]
    for cunt1 in ["UK","BE","DE","DK"]
        for cunt2 in ["UK","BE","DE","DK"]
            if (cunt1!=cunt2)
                for (k1, s1) in bins[res_yr][cunt1]
                    for (k2, s2) in bins[res_yr][cunt1]
                        _g1=intersect(s1,s2)
_g1=intersect(bins["2014"]["BE"]["1"],bins["2014"]["DE"]["1"],bins["2014"]["DK"]["1"],bins["2014"]["UK"]["1"])
_g1=intersect(bins["2014"]["BE"]["2"],bins["2014"]["DE"]["2"],bins["2014"]["DK"]["2"],bins["2014"]["UK"]["2"])

intersect(bede,bins["2014"]["UK"]["1"])            
        
        Ytr = convert(Matrix, daily_dew)'
        
        kshape_clusters_deuk=ks.kshape(ks.zscore(sqrt.((daily_uk.^2).+(daily_de.^2).+(daily_dk.^2).+(daily_dew.^2))',axis=1), k)
        ts=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); ts=vcat(ts,daily_ts[:,rand(clusters)+1]); end;end
        filter!(row -> row.time_stamp in ts, scenario_data[sc][yr])
    end
end
save("./test/data/input/EUSTDG17t020_TS_k"*string(k)*".jld2",scenario_data)
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
    =#