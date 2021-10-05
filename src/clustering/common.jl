
function PCA_cluster(zsd)
    # split half to training set
    Xtr = convert(Matrix, zsd[1:1:end,2:end])'

    # train a PCA model, allowing up to n dimensions
    M = MultivariateStats.fit(MultivariateStats.PCA, Xtr; maxoutdim=length(Xtr[:,1]))

    # apply PCA model to testing set
    Ytr = MultivariateStats.transform(M, Xtr)

    return Ytr
end

function kmeans_cluster(Ytr)
    D = Distances.pairwise(Distances.Euclidean(), Ytr, Ytr)
    strt=3;nd=10

    clusters=[Clustering.kmeans(Ytr, k) for k = strt:1:nd];cluster_quality=[]
    for (i,cluster) in enumerate(clusters)
        push!(cluster_quality,(i,quality_of_cluster(cluster,D)))
    end
    return clusters[findmax(last.(last.(cluster_quality)))[2]]
end

function n_samps(cluster,n)
    bins=[]; ns=[]; n_samples=[]
    for count in cluster.counts; push!(bins,Vector{Int64}());end
    for (i,ass) in enumerate(cluster.assignments); push!(bins[ass],i);end
    for (i,bin) in enumerate(bins); for j=1:1:ceil(Int64,length(bin)/length(cluster.assignments)*n); push!(ns,i); end;end
    for n in ns
        push!(n_samples,rand(bins[n]))
    end
    return sort!(n_samples)
end

function kmedoid_cluster(Ytr)
    D = Distances.pairwise(Distances.Euclidean(), Ytr, Ytr)
    strt=3;nd=10

    clusters=[Clustering.kmedoids(D, k) for k = strt:1:nd];cluster_quality_m=[]
    for (i,cluster) in enumerate(clusters)
        push!(cluster_quality_m,(i,quality_of_cluster(cluster,D)))
    end
    return clusters[findmax(last.(last.(cluster_quality_m)))[2]]
end

function quality_of_cluster(cluster,D)
    sil = Clustering.silhouettes(cluster, D)
    return (MultivariateStats.mean(sil),MultivariateStats.var(sil),count(i->(i<0), sil),(length(sil)-count(i->(i<0), sil))/length(sil)+MultivariateStats.mean(sil)+(1-MultivariateStats.var(sil))+1e-3*length(cluster.counts))
end

function normalize_ts(nrm_ts)
    mn=minimum(nrm_ts.price)
    nrm_ts.price=nrm_ts.price.-mn#set minimum to zero
    mx=maximum(nrm_ts.price)
    nrm_ts.price=nrm_ts.price./mx#normalize
    return nrm_ts
end


function rand_set_days(ukbe_data,dn,wn,mn)
    sampled_set=DataFrame()
    selected_months=[]; for j=1:mn; months=[i for i=1:12 if !issubset(i,selected_months)]; push!(selected_months,rand(months));end;sort!(selected_months)
    for mnth in selected_months;#println("month: "*string(mnth))
        hours_in_month=filter(e->month(e.time_stamp)==mnth,ukbe_data)
        days_in_month=round(Int64,length(hours_in_month[!,:time_stamp])/24)
        weeks_in_month=ceil(Int64,days_in_month/7);selected_weeks=[];
        #for w=1:wn; weeks=[i for i=1:weeks_in_month if !issubset(i,selected_weeks)]; push!(selected_weeks,rand(weeks));end;sort!(selected_weeks);selected_days=[];
        selected_weeks=wn;selected_days=[];
        for wek in selected_weeks;#println(wek)
            for d=1:dn; days=[i for i in (wek-1)*7+1:min((wek)*7,days_in_month) if !issubset(i,selected_days)]; if (length(days)>0); push!(selected_days,rand(days));end;end;sort!(selected_days);end
        for dy in selected_days;
            hours_in_day=filter(e->day(e.time_stamp)==dy,hours_in_month)
            sampled_set=vcat(sampled_set,hours_in_day)
        end
    end
    sort!(sampled_set)
    return sampled_set
end
