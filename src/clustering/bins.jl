####################################### Clustering with bins ###################
function sample_cluster(tss_bins,tss,n)
    sample_tss=[]
    tinf=length(tss[1].time_stamp)
    loops=0
    while (length(sample_tss)<n && loops<10000)
        sample_i=rand(1:tinf)
        sample=tss[1][sample_i,:].time_stamp
        ok,tss_bins=check_sample(sample,tss,tss_bins)#working on this function
        if (ok)
            println(length(sample_tss))
            push!(sample_tss,sample)#store sample set
            loops=0
        end
        loops=loops+1
    end
    return sample_tss
end

function check_sample(sample,tss,tss_bins)
    ok=true;locations=[]
    for (key0,ts_bins) in tss_bins
        for (key1,ts_bin) in ts_bins
            if (length(filter(row -> row.time_stamp == sample, ts_bin["ts"]).time_stamp)>0)
                if (ts_bin["full"]==false)
                    push!(locations,(key0,key1))
                else
                    ok=false
                    @goto bad_sample
                end
            end
        end
    end
    if (length(locations)==length(tss_bins))
        for location in locations
            if (haskey(tss_bins[location[1]][location[2]], "chosen")==false);
                push!(tss_bins[location[1]][location[2]],"chosen"=>1);
            else
                tss_bins[location[1]][location[2]]["chosen"]=tss_bins[location[1]][location[2]]["chosen"]+1
            end
            if (tss_bins[location[1]][location[2]]["elements"]-tss_bins[location[1]][location[2]]["chosen"]==0)
                tss_bins[location[1]][location[2]]["full"]=true
            end
        end
    else
        ok=false
    end
    @label bad_sample
    return ok,tss_bins
end

function cluster_probs(tss_bins,tss,n)
    tinf=length(tss[1].time_stamp)
    for (k0,ts_bins) in tss_bins
        for (k1,ts_bin) in ts_bins
            elements=length(ts_bin["ts"].time_stamp)==0 ? 0 : ceil(Int64,n*length(ts_bin["ts"].time_stamp)/tinf)
            push!(ts_bin,"elements"=>elements)
            if (elements>0); push!(ts_bin,"full"=>false);else push!(ts_bin,"full"=>true);end
        end
    end
    return tss_bins
end
#tss is an array of time series arays, n is the desired sample size
function cluster_ts(tss,n)
    number_tss=length(tss)
    tss_bins=Dict{String,Any}()
    for (key,ts) in enumerate(tss)
        push!(tss_bins,string(key)=>ts_binify(ts))#sort into bins
    end
    tss_bins=cluster_probs(tss_bins,tss,n)
    return tss_bins
end

function ts_binify(ts)
    ts_bins=Dict{String,Any}()
    nrm_ts=normalize_ts(deepcopy(ts))
    nm_bins=10
    for b=1/nm_bins:1/nm_bins:1
        for p in eachrow(nrm_ts)
            if (p.price<=b && p.price>(b-1/nm_bins))
                if (haskey(ts_bins, string(b))==false);push!(ts_bins,string(b)=>Dict{String,Any}());push!(ts_bins[string(b)],"ts"=>DataFrame("time_stamp"=>[]));end
                push!(ts_bins[string(b)]["ts"],p[1:1])
            end
        end
    end
    return ts_bins
end
