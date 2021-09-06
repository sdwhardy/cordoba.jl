location="Belgium"
df = DataFrame(XLSX.readtable("/Users/shardy/Documents/GitHub/hv_offshore_topology_optimization/v1.0/layout/input_data/"*location*"/data.xlsx", "wind_data")...)
wnd="Norther"
cbl_data=get_220kV_cables()
cable_tables=Array{Array{cable,1},1}()
km=72.5
#ps=[50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0,550.0,600.0,650.0,700.0,750.0,800.0,850.0,900.0,950.0,1000.0,1050.0,1100.0,1150.0,1200.0,1250.0,1300.0,1350.0,1400.0,1450.0,1500.0]
#ps=[50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0,550.0,600.0,650.0,700.0]
#
ps=[1200.0]

osss=Vector{}()
cbls=Vector{}()
osss_costs=Vector{}()
cbls_costs=Vector{}()
ks=get_Cost_Data()
xfo_data=get_Xfo_Data()
xfo_on=transformer()
for p in ps
    oss=bus()
    xfo_off=transformer()
    xfo_off.relia=get_offshore_xfo_failure_data(xfo_off)
    xfo_off.wnd=wnd
    xfo_off.mva=p
    xfo_on=deepcopy(xfo_off)
    xfo_off=xfo_oss(xfo_off,ks,xfo_data)
    push!(oss.xfmrs,xfo_off)
    xfo_pcc(xfo_on,ks,xfo_data)

    plat=platform()
    plat.mva=p
    plat.wnd=wnd
    plat=cost_ac_platform(plat,ks)
    plat=adjust_base_ac_platform(plat,ks)
    push!(oss.plat,plat)
    push!(osss,oss)

cable()

    push!(cable_tables,solve_cable_interval((1,km,p),cbl_data,wnd,ks))
    cable_dictionary=make_dictionaries_cables([cable_tables])
    sg_qc=cable_dictionary["220.0"][string(p)][length(cable_dictionary["220.0"][string(p)])].costs.sg+cable_dictionary["220.0"][string(p)][length(cable_dictionary["220.0"][string(p)])].costs.qc
    cbl_cost=cable_dictionary["220.0"][string(p)][length(cable_dictionary["220.0"][string(p)])].costs.ttl-sg_qc-cable_dictionary["220.0"][string(p)][length(cable_dictionary["220.0"][string(p)])].costs.rlc-cable_dictionary["220.0"][string(p)][length(cable_dictionary["220.0"][string(p)])].costs.eens-cable_dictionary["220.0"][string(p)][length(cable_dictionary["220.0"][string(p)])].costs.cm
    push!(cbls,cable_dictionary["220.0"][string(p)][length(cable_dictionary["220.0"][string(p)])])
    println(sg_qc)
    println(sg_qc)
    println(sg_qc)
    push!(cbls_costs,cbl_cost/km)
    push!(osss_costs,(plat.costs.ttl-plat.costs.cm+xfo_off.costs.ttl-xfo_off.costs.cm-xfo_off.costs.eens-xfo_off.costs.tlc+sg_qc/2))
    push!(osss_costs,(plat.costs.ttl-plat.costs.cm+sg_qc/2))
    #push!(osss_costs,xfo_on.costs.ttl+sg_qc/2)
end

cbls_costs[1]*km

println(osss_costs)

plotly()
p=plot()
plot!(p,[ps],[cbls_costs],xticks = 0:200:800,xlims=(0,800),ylims=(0,3),color = :black,label="",xaxis = ("MW", font(20, "Courier")),yaxis = ("ME", font(20, "Courier")))
gui()=#

plotly()
p=plot()
plot!(p,[ps],[osss_costs],xticks = 0:100:800,xlims=(0,700),ylims=(0,250),color = :black,label="",xaxis = ("MW", font(20, "Courier")),yaxis = ("ME", font(20, "Courier")))
gui()

osss_costs_orig=deepcopy(osss_costs)



GMIN=minimum([wf.mva for wf in owpps])
A=range(1,stop=30)
base_b=zeros(Int8,1,length(owpps))
B=Array{Tuple{Array{Int8,2},Int64},1}()
combos=collect(Iterators.flatten([combinations(A, k) for k = 1:9]))
