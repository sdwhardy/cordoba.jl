    ################## loads external packages ##############################
    using Ipopt, Gurobi, JuMP, FileIO, JLD2, Dates, OrderedCollections, CSV, DataFrames, XLSX
    using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
    import cordoba; const _CBD = cordoba#Cordoba package backend - under development
    import PowerModelsACDC; const _PMACDC = PowerModelsACDC
    import PowerModels; const _PM = PowerModels
    rt=pwd()
    include(rt*"\\test\\binaries\\print_m_file.jl")
    include(rt*"\\aux\\post_process\\functions.jl")

    rt_ex=rt*"\\test\\data\\input\\BE_DE_UK_DK_islands\\"
    nodes = DataFrame(XLSX.readtable(rt_ex*"north_sea_intercons_nodes.xlsx", "node_generation")...)
    edges = DataFrame(XLSX.readtable(rt_ex*"north_sea_intercons_nodes.xlsx", "connections_acdc")...)
    ac_cable_df = DataFrame(XLSX.readtable(rt_ex*"north_sea_intercons_nodes.xlsx", "CABLES_AC_SET_UP")...)
    dc_cable_df = DataFrame(XLSX.readtable(rt_ex*"north_sea_intercons_nodes.xlsx", "CABLES_DC_SET_UP")...)
    rem_df = DataFrame(XLSX.readtable(rt_ex*"north_sea_intercons_nodes.xlsx", "REMAINDER")...)
    ppf_mainACDC2mfile(rem_df,ac_cable_df,dc_cable_df,rt_ex)

    casename = "test_case"
    file = rt_ex*"$casename.m"
    data = PowerModels.parse_file(file)
    data,ics_ac=_CBD.filter_AClines(data,edges,nodes)
    data,ics_dc=_CBD.filter_DClines(data,edges,nodes)


    pu=data["baseMVA"]
    owpp_mva=[2000,4000,4000,12000]./pu#mva of wf (in de)
    infinite_grid=sum(owpp_mva)
    markets_wfs=[String[],String[]]#UK,DE,DK must be in same order as .m file gens
    for (k,cunt) in enumerate(nodes["country"])
        if (nodes["type"][k]>0)
        push!(markets_wfs[1],cunt);else
        push!(markets_wfs[2],cunt);end
    end
    ######################## seperate gnz and wfs ###################
    genz=[];wfz=[]
    for i=1:1:length(markets_wfs[1]); push!(genz,(i,infinite_grid));end
    for i=1:1:length(markets_wfs[1]); push!(genz,(i+length(markets_wfs[1])+length(markets_wfs[2]),infinite_grid));end
    for i=1:1:length(markets_wfs[2]); push!(wfz,(i+length(markets_wfs[1]),owpp_mva[i]));end



    #################### Calculates cable options for AC lines
    candidate_ics_ac=[1]
    data=_CBD.additional_candidatesICS_AC(data,candidate_ics_ac,ics_ac)#adds additional candidates
    z_base_ac=(data["bus"]["1"]["base_kv"])^2/pu
    for (i,bac) in data["ne_branch"]
    data["ne_branch"][i]=_CBD.candidateIC_cost_impedance_AC(bac,z_base_ac,pu);end
    data["ne_branch"]=_CBD.unique_candidateIC_AC(data["ne_branch"])#keep only unique candidates
    print_topology_data_AC(data,markets_wfs)#print to verify

    #################### Calculates cable options for DC lines
    candidate_ics_dc=[1]#Candidate Cable sizes
    data=_CBD.additional_candidatesICS_DC(data,candidate_ics_dc,ics_dc)#adds additional candidates
    z_base_dc=(data["busdc"]["1"]["basekVdc"])^2/pu
    for (i,bdc) in data["branchdc_ne"]
    data["branchdc_ne"][i]=_CBD.candidateIC_cost_impedance_DC(bdc,z_base_dc);end
    data["branchdc_ne"]=_CBD.unique_candidateIC_DC(data["branchdc_ne"])#keep only unique candidates
    print_topology_data_DC(data,markets_wfs)#print to verify


    #######################
    PowerModelsACDC.process_additional_data!(data)
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer)
    s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false)
    _CBD.converter_parameters_rxb(data)

    print_topology_data_DC(data,[[],[]])
