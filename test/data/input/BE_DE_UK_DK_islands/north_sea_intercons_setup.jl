    ################## loads external packages ##############################
    using Ipopt, Gurobi, JuMP, FileIO, JLD2, Dates, OrderedCollections, CSV, DataFrames, XLSX
    using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
    import cordoba; const _CBD = cordoba#Cordoba package backend - under development
    import PowerModelsACDC; const _PMACDC = PowerModelsACDC
    import PowerModels; const _PM = PowerModels
    include("../../../binaries/print_m_file.jl")
    include("../../../../aux/post_process/functions.jl")

    rem_df = DataFrame(XLSX.readtable("./test/data/input/BE_DE_UK_DK_islands/north_sea_intercons_nodes.xlsx", "REMAINDER")...)
    ac_cable_df = DataFrame(XLSX.readtable("./test/data/input/BE_DE_UK_DK_islands/north_sea_intercons_nodes.xlsx", "CABLES_AC_SET_UP")...)
    dc_cable_df = DataFrame(XLSX.readtable("./test/data/input/BE_DE_UK_DK_islands/north_sea_intercons_nodes.xlsx", "CABLES_DC_SET_UP")...)
    ppf_main2mfile(rem_df,ac_cable_df)
    ics_ac=[(500,10),(475,15),(450,20),(425,25),(400,30),(375,35),(350,40),(325,45),(300,50),(275,55),(250,60),(550,65),(600,70),(650,75),(700,80),(750,65),(800,70),(900,75),(1000,80)];#UK,DE,DK

    casename = "test_case"
    file = "./test/data/input/BE_DE_UK_DK_islands/$casename.m"
    data = PowerModels.parse_file(file)
    candidate_ics_ac=[1,10/20,4/20]
    data=_CBD.additional_candidatesICS_AC(data,candidate_ics_ac,ics_ac)#adds additional candidates
    z_base_ac=(data["bus"]["1"]["base_kv"])^2/data["baseMVA"]
    for (i,bac) in data["ne_branch"]
    data["ne_branch"][i]=_CBD.candidateIC_cost_impedance_AC(bac,z_base_ac);end


    PowerModelsACDC.process_additional_data!(data)
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer)
    s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false)
    _CBD.converter_parameters_rxb(data)
    data["ne_branch"]=_CBD.unique_candidateIC_AC(data["ne_branch"])#keep only unique candidates
    print_topology_data_AC(data,[[],[]])


    ##

    resultACDC = _PMACDC.run_acdctnepopf(data, _PM.DCPPowerModel, gurobi, setting = s)

    print_solution_data_AC(resultACDC, data)
