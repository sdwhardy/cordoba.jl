#"connection,cost,time,ic_mva,ic_km,owpp_mva,owpp_km"
function store_data(wcfile,resultACDC, ic_mva, owpp_mva, ic_length, owpp_km)
    if (resultACDC["solution"]["nw"]["1"]["ne_branch"]["1"]["built"]>0)
        connection=-1
    elseif (resultACDC["solution"]["nw"]["1"]["ne_branch"]["2"]["built"]>0)
        connection=0
    else
        connection=1
    end
    println(wcfile,string(connection)*", "*string(resultACDC["objective"])*", "*string(resultACDC["solve_time"])*", "*string(ic_mva)*", "*string(ic_length)*", "*string(owpp_mva)*", "*string(owpp_km))
end


# Example code to print list of built HVDC branches and converters
function display_results(result)
    built_cv = []
    built_br = []
    built_ACbr = []
    for (c, conv) in result["solution"]["convdc_ne"]
        if isapprox(conv["isbuilt"] , 1; atol = 0.01)
            print("Conv: $c \n")
            push!(built_cv,c)
        end
    end
    for (b, branch) in result["solution"]["branchdc_ne"]
        if isapprox(branch["isbuilt"] , 1; atol = 0.01)
            print("DCBranch: $b \n")
            push!(built_br,b)
        end
    end
    for (b, branch) in result["solution"]["ne_branch"]
        if isapprox(branch["built"] , 1; atol = 0.01)
            print("ACBranch: $b \n")
            push!(built_ACbr,b)
        end
    end
    return built_cv, built_br, built_ACbr
end


#=
plotly()
width=750
height=500
#p=plot(size = (width, height),xaxis = ("Nodes", font(40, "Courier")),yaxis = ("Depth [m]", font(40, "Courier")))
p=plot(size = (width, height),xticks = 0:200:1000,xlims=(0,1000),markersize=2,seriestype=:scatter,xaxis = ("Clusters", font(20, "Courier")),yaxis = ("Silhouette Mean", font(20, "Courier")))
#p=plot(size = (width, height),xlims=(590,910),ylims=(-4,3),xaxis = ("T.S", font(20, "Courier")),yaxis = ("", font(20, "Courier")))
plot!(p,2:1:1000,first.(last.(cluster_quality)),color = :red,seriestype=:scatter,markersize=3, markershape = :circle, label="Kmeans",size = (width, height))
plot!(p,2:1:1000,first.(last.(cluster_quality_m)),color = :blue,seriestype=:scatter,markersize=3, markershape = :circle, label="Kmedoids",size = (width, height))
gui()

#BE-UK average energy price (NEMO commissioned start of 2019)
uk_be=[(2021, 9.87),(2020, 9.87),
(2019, 9.76),
(2019, 11.05),
(2018, 11.05),
(2017, 11.09)]

using Plots
plotly()
width=750
height=500
#p=plot(size = (width, height),xaxis = ("Nodes", font(40, "Courier")),yaxis = ("Depth [m]", font(40, "Courier")))
p=plot(size = (width, height),xlims=(2017,2021),xaxis = ("Year", font(16, "Courier")),yaxis = ("Mean Diff in Euro/MWh (UK-BE)", font(16, "Courier")))
plot!(p,first.(uk_be),last.(uk_be),color = :red, label="",size = (width, height))
gui()=#

#=
# train a PCA model, allowing up to n dimensions
for k in 1:10; push!(_pca,fit(PCA, Xtr; maxoutdim=k)); end

plotly()
width=750
height=500
#p=plot(size = (width, height),xaxis = ("Nodes", font(40, "Courier")),yaxis = ("Depth [m]", font(40, "Courier")))
p=plot(size = (width, height),xticks = 0:200:1000,xlims=(0,1000),markersize=2,seriestype=:scatter,xaxis = ("Clusters", font(20, "Courier")),yaxis = ("Silhouette Mean", font(20, "Courier")))
#p=plot(size = (width, height),xlims=(590,910),ylims=(-4,3),xaxis = ("T.S", font(20, "Courier")),yaxis = ("", font(20, "Courier")))
plot!(p,2:1:1000,first.(last.(cluster_quality)),color = :red,seriestype=:scatter,markersize=3, markershape = :circle, label="Kmeans",size = (width, height))
plot!(p,2:1:1000,first.(last.(cluster_quality_m)),color = :blue,seriestype=:scatter,markersize=3, markershape = :circle, label="Kmedoids",size = (width, height))
gui()

using PlotlyJS, CSV, DataFrames
D = pairwise(Euclidean(), Ytr, Ytr)
strt=2;nd=8

clusters=[Clustering.kmeans(Ytr, k) for k = strt:1:nd];cluster_quality=DataFrame
for (i,cluster) in enumerate(clusters)
    push!(cluster_quality,(i+1,Clustering.silhouettes(cluster, D)))
end
j=7
cluster_qualit=DataFrame(cluster=[first(cluster_quality[j]) for i=1:1:length(last(cluster_quality[j]))],mean=last(cluster_quality[j]))


color_vec = fill("lightslategray", 5)
color_vec[2] = "crimson"
sort!(cluster_qualit.mean)
for rw in cluster_qualit.mean
    println(rw)
end



###########################################
kkm = [(k,mean(Clustering.silhouettes(Clustering.kmedoids(D, k), D))) for k = 2:1:20]
kk = [(k,mean(Clustering.silhouettes(, D))) for k = 2:1:20]
clusters=10:10:1000;rez=[]
for cluster in clusters
    clustered_data = kmeans(Ytr, cluster; maxiter=500)#kmeans cluster the data
    push!(rez,(cluster,clustered_data.totalcost))
    println(cluster)
end
#=silhouettes(clustered_data.assignments,clustered_data.counts,clustered_data.counts)
bins=[Vector{}() for i=1:clusters]
[push!(bins[j],Xtr_labels[i]) for (i,j) in enumerate(clustered_data.assignments)]=#
plotly()
width=750
height=500
#p=plot(size = (width, height),xaxis = ("Nodes", font(40, "Courier")),yaxis = ("Depth [m]", font(40, "Courier")))
p=plot(size = (width, height),xticks = 0:10:100,xlims=(0,100),markersize=2,seriestype=:scatter,xaxis = ("Clusters", font(20, "Courier")),yaxis = ("Variance Silhouette", font(20, "Courier")))
#p=plot(size = (width, height),xlims=(590,910),ylims=(-4,3),xaxis = ("T.S", font(20, "Courier")),yaxis = ("", font(20, "Courier")))
plot!(p,first.(kk),last.(kk),color = :red,seriestype=:scatter,markersize=3, markershape = :square, label="KMEANS",size = (width, height))
plot!(p,first.(kkm),last.(kkm),color = :blue,seriestype=:scatter,markersize=3, markershape = :circle, label="KMEDOIDS",size = (width, height))
gui()

# visualize first 3 principal components in 3D interacive plot
p = scatter(setosa[1,:],setosa[2,:],setosa[3,:],marker=:circle,linewidth=0)
scatter!(versicolor[1,:],versicolor[2,:],versicolor[3,:],marker=:circle,linewidth=0)
scatter!(virginica[1,:],virginica[2,:],virginica[3,:],marker=:circle,linewidth=0)
plot!(p,xlabel="PC1",ylabel="PC2",zlabel="PC3")
gui()
=#
