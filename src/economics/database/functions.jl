#this file constructs the equipment lookup database for system powers, lengths
#that is stored in the ocean layout
################################################################################
################################# Cables #######################################
################################################################################
#=kbles=database["cables"]
km_mva=km_mva_set[16]
kbl66=kbles["66.0"][string(km_mva[2])][1]=#
#**
function set_cables_per_km_cost(kbles,km_mva_set)
    eps=10^-2
    ks=get_Cost_Data()
    for km_mva in km_mva_set
        for kbl66 in kbles["66.0"][string(km_mva[2])]
            c1=kbl66.costs.ttl
            c0=deepcopy(mvac_cable(kbl66.mva,kbl66.length-1,kbl66.wnd,kbles["66.0"][string(km_mva[2])],ks)).costs.ttl
            kbl66.costs.perkm_ttl=deepcopy(c1-c0)
            kbl66.mx_rng=kbl66.length+eps
        end
        for kbl220 in kbles["220.0"][string(km_mva[2])]
            c1=kbl220.costs.ttl
            c0=deepcopy(hvac_cable(kbl220.mva,kbl220.length-1,kbl220.wnd,kbles["220.0"][string(km_mva[2])],ks)).costs.ttl
            kbl220.costs.perkm_ttl=deepcopy(c1-c0)
            kbl220.mx_rng=kbl220.length+eps
        end
        for kbl400 in kbles["400.0"][string(km_mva[2])]
            c1=kbl400.costs.ttl
            c0=deepcopy(hvac_cable(kbl400.mva,kbl400.length-1,kbl400.wnd,kbles["400.0"][string(km_mva[2])],ks)).costs.ttl
            kbl400.costs.perkm_ttl=deepcopy(c1-c0)
            kbl400.mx_rng=kbl400.length+eps
        end
        for kbl300 in kbles["300.0"][string(km_mva[2])]
            c1=kbl300.costs.ttl
            c0=deepcopy(hvdc_cable(kbl300.mva,kbl300.length-1,kbl300.wnd,kbles["300.0"][string(km_mva[2])],ks)).costs.ttl
            kbl300.costs.perkm_ttl=deepcopy(c1-c0)
            kbl300.mx_rng=kbl300.length+eps
        end
    end
    return kbles
end

#**
function get_cable_table(km_mva_set,wnd)
    ks=get_Cost_Data()
    cable_tables=Array{Array{Array{cable,1},1},1}()
    #HVDC
    cbl_datas=[get_300kV_cables()]
    for cbl_data in cbl_datas
        cable_table=Array{Array{cable,1},1}()
        for km_mva in km_mva_set
            cbl=cable()
            cbl.mva=km_mva[2]
            cbl.wnd=wnd
            push!(cable_table,[optimal_hvdc_cable(cbl,cbl_data,km_mva[1],ks)])
        end
        push!(cable_tables,deepcopy(cable_table))
    end
    #HVAC
    cbl_datas=[get_220kV_cables(),get_400kV_cables()]
    for cbl_data in cbl_datas
        cable_table=Array{Array{cable,1},1}()
        for km_mva in km_mva_set
            cable_interval=Tuple
            cable_interval=(1,km_mva[1],km_mva[2])
            push!(cable_table,solve_cable_interval(cable_interval,cbl_data,wnd,ks))
        end
        push!(cable_tables,deepcopy(cable_table))
    end
    #MVAC
    cbl_data=get_66kV_cables()
    cable_table=Array{Array{cable,1},1}()
    for km_mva in km_mva_set
        cable_interval=Tuple
        cable_interval=(1,km_mva[1],km_mva[2])
        push!(cable_table,solve_cable_interval(cable_interval,cbl_data,wnd,ks))
    end
    push!(cable_tables,deepcopy(cable_table))
    cable_dictionary=make_dictionaries_cables(cable_tables)
    return cable_dictionary
end

#finds transistions in cable sizes per up to km at mva for kv
#**
function solve_cable_interval(cable_interval,cbl_data,wnd,ks)
    openQ=Array{Array{cable,1},1}()
    closedQ=Array{Array{cable,1},1}()
    cable_table=Array{cable,1}()

    #calculate cable at min of interval distance
    cbl=cable()
    km_min=cable_interval[1]
    cbl.mva=cable_interval[3]
    cbl.wnd=wnd
    cbl.length=km_min
    cbl_min=optimal_ac_cable(cbl,cbl_data,km_min,ks)

    #calculate cable at max of interval distance
    cbl=cable()
    km_max=cable_interval[2]
    cbl.mva=cable_interval[3]
    cbl.wnd=wnd
    cbl.length=km_max
    cbl_max=optimal_ac_cable(cbl,cbl_data,km_max,ks)

    if ((cbl_min.size==cbl_max.size) && (cbl_min.num==cbl_max.num))
        push!(cable_table,cbl_max)
    else
        push!(openQ,[cbl_min,cbl_max])
        while (length(openQ)>0)
            openQ=halfway_cable(openQ,cbl_data,ks)
            closedQ,openQ=update_Qs(closedQ,openQ)
        end
        cable_table_oversized=Array{cable,1}()
        for q in closedQ
            push!(cable_table_oversized,deepcopy(q[1]))
            push!(cable_table_oversized,deepcopy(q[2]))
        end
        push!(cable_table_oversized,cbl_max)
        for (index,cbl) in enumerate(cable_table_oversized[1:length(cable_table_oversized)-1])
            if ((cbl.size==cable_table_oversized[index+1].size) && (cbl.num==cable_table_oversized[index+1].num))
            else
                push!(cable_table,deepcopy(cbl))
            end
        end
        push!(cable_table,cbl_max)
    end
    return cable_table
end

#moves finished intervals form the openQ to the closedQ
#**
function update_Qs(closedQ,openQ)
    openQ_new=Array{Array{cable,1},1}()
    for q in openQ
        if ((q[1].size==q[2].size) && (q[1].num==q[2].num))
            push!(closedQ,deepcopy(q))
        elseif ((q[1].length+1>=q[2].length))
            push!(closedQ,deepcopy(q))
        else
            push!(openQ_new,deepcopy(q))
        end
    end
    return closedQ,openQ_new
end

#devides open q at mid point
#**
function halfway_cable(openQ,cbl_data,ks)
    openQ_new=Array{Array{cable,1},1}()
    for q in openQ
        cbl=cable()
        km_halfway=(q[1].length+q[2].length)/2
        cbl.mva=q[1].mva
        cbl.wnd=q[1].wnd
        cbl_halfway=optimal_ac_cable(cbl,cbl_data,km_halfway,ks)
        if ((cbl_halfway.size==q[2].size) && (cbl_halfway.num==q[2].num))
            push!(openQ_new,[q[1],cbl_halfway])
        elseif ((cbl_halfway.size==q[1].size) && (cbl_halfway.num==q[1].num))
            push!(openQ_new,[cbl_halfway,q[2]])
        else
            push!(openQ_new,[q[1],cbl_halfway])
            push!(openQ_new,[cbl_halfway,q[2]])
        end
    end
    return openQ_new
end

#finds the hvdc cable for the look up table
#**

function optimal_hvdc_cable(cbl0,cbl_data,km,ks)
    cbl=cable()
    cbl.costs.ttl=Inf
    cbl.elec.volt=cbl_data[1][1]
    cbl.length=km
    cbl.mva=cbl0.mva
    cbl0.costs.ttl=Inf

    for _cd in cbl_data
        num=2
        capacity=(_cd[1]*_cd[5])*10^-3
        while (cbl0.mva>num*capacity)
            num=num+2
        end
        fillOut_cable_struct_dc(cbl0,_cd,km,num)
        cbl0=cost_hvdc_cable(cbl0,ks)
        if (cbl0.costs.ttl<cbl.costs.ttl)
            cbl=deepcopy(cbl0)
        end
    end
    return cbl
end

############################## Change made here ########################### Temporary kana???
#returns the optimal cable for a hvac for the look up table
#**
function optimal_ac_cable(cbl0,cbl_data,km,ks)
    cbl=cable()
    cbl.costs.ttl=Inf
    cbl.elec.volt=cbl_data[1][1]
    cbl.length=km
    cbl.mva=cbl0.mva
    ##########################
    #=if (cbl0.elec.volt==66.0)
        cbl.mva=cbl0.mva*1.333
    else
        cbl.mva=cbl0.mva
    end=#
    ##########################
    cbl0.costs.ttl=Inf
    for cd in cbl_data
        #num=1
        ######################
        if (cbl0.elec.volt==66.0)
            num=6
        else
            num=1
        end
        ######################
        exception_66kv=false
        cap_at_km=get_newQ_Capacity(cbl0.elec.freq,km,cd[1],cd[4]*10^-9,cd[5])
        #max number of cables in parallel = 12
        max_in_parallel=12#if changing check economics main function mvac_cable(mva,km,wnd,cable_array,ks)
        GMAX=2000
        if (cap_at_km/cbl0.mva>1/max_in_parallel || (cbl.elec.volt==66.0 && cbl0.mva<=3000))
            while ((cbl0.mva>num*cap_at_km) && num<max_in_parallel)
                num=deepcopy(num+1)
            end
            if ((cbl.elec.volt==66.0 && cbl0.mva>num*cap_at_km))
                num=1
                exception_66kv=true
                cbl0.mva=cbl0.mva/2
                while ((cbl0.mva>num*cap_at_km) && num<max_in_parallel)
                    num=deepcopy(num+1)
                end
            end
            fillOut_cable_struct_ac(cbl0,cd,km,num)
            #cbl0=cost_hvac_cable_o2o(cbl0,ks)
            if (cbl0.elec.volt==66.0)
                cbl0=cost_mvac_cable(cbl0,ks)
                if (exception_66kv==true)
                    cbl0.num=cbl0.num*2
                    cbl0.mva=cbl0.mva*2
                    cbl0.costs.perkm_cpx=cbl0.costs.perkm_cpx*2
                    cbl0.costs.sg=cbl0.costs.sg*2
                    cbl0.costs.cpx_p=cbl0.costs.cpx_p*2
                    cbl0.costs.cpx_i=cbl0.costs.cpx_i*2
                    cbl0.costs.rlc=cbl0.costs.rlc*2
                    cbl0.costs.cm=cbl0.costs.cm*2
                    cbl0.costs.eens=cbl0.costs.eens*2
                    cbl0.costs.ttl=cbl0.costs.ttl*2
                    cbl0.costs.grand_ttl=cbl0.costs.ttl
                end
            elseif (cbl0.elec.volt==220.0 || cbl0.elec.volt==400.0)
                cbl0=cost_hvac_cable(cbl0,ks)
            else
            end
        end
        if (cbl0.costs.ttl<cbl.costs.ttl)
            cbl=deepcopy(cbl0)
        end
    end
    return cbl
end

#puts cables into dictionaries
#**
function make_dictionaries_cables(cable_table)
    all_cable_dict=Dict{String,Dict{String,Vector{cable}}}()
    for cable_sets in cable_table
        kv_cable_dict=Dict{String,Vector{cable}}()
        for cable_set in cable_sets
             push!(kv_cable_dict,(string(cable_set[1].mva)=>cable_set))
        end
        push!(all_cable_dict,(string(cable_sets[1][1].elec.volt)=>kv_cable_dict))
    end
    return all_cable_dict
end
################################################################################
######################## Transformers and converters ###########################
################################################################################
#**
function get_xfoConv_table(km_mva_set,wnd)
    ks=get_Cost_Data()
    #HVDC
    conv_table=Array{Array{converter,1},1}()
    for km_mva in km_mva_set
        conv_off=converter()
        conv_off.relia=get_offshore_conv_failure_data(conv_off)
        conv_off.wnd=wnd
        conv_off.mva=km_mva[2]
        conv_off.elec.mva=km_mva[2]
        conv_off.num=1
        conv_off=cost_hvdc_oss(conv_off,ks)
        conv_off=adjust_base_hvdc_offshore_converter(conv_off,ks)

        conv_on=converter()
        conv_on.relia=get_onshore_conv_failure_data(conv_on)
        conv_on.wnd=wnd
        conv_on.mva=km_mva[2]
        conv_on.elec.mva=km_mva[2]
        conv_on.num=1
        conv_on=cost_hvdc_pcc(conv_on,ks)
        conv_on=adjust_base_hvdc_onshore_converter(conv_on,ks)
        push!(conv_table,[conv_off,conv_on])
    end

    #HVAC
    xfo_table=Array{Array{transformer,1},1}()
    xfo_data=get_Xfo_Data()
    for km_mva in km_mva_set
        xfo_off=transformer()
        xfo_off.relia=get_offshore_xfo_failure_data(xfo_off)
        xfo_off.wnd=wnd
        xfo_off.mva=km_mva[2]
        xfo_off=xfo_oss(xfo_off,ks,xfo_data)

        xfo_on=transformer()
        xfo_on.relia=get_onshore_xfo_failure_data(xfo_on)
        xfo_on.wnd=wnd
        xfo_on.mva=km_mva[2]
        xfo_on=xfo_pcc(xfo_on,ks,xfo_data)
        push!(xfo_table,[xfo_off,xfo_on])
    end

    xfo_dictionary,conv_dictionary=make_dictionaries_xfo_conv(xfo_table,conv_table)
    return xfo_dictionary,conv_dictionary
end


#puts transformers and converters into dictionaries
#**
function make_dictionaries_xfo_conv(xfo_table,conv_table)
    all_xfo_dict=Dict{String,Dict{String,transformer}}()
    for xfo_sets in xfo_table
        xfo_dict=Dict{String,transformer}()
        push!(xfo_dict,("offshore"=>xfo_sets[1]))
        push!(xfo_dict,("onshore"=>xfo_sets[2]))
        push!(all_xfo_dict,(string(xfo_sets[1].mva)=>xfo_dict))
    end
    all_conv_dict=Dict{String,Dict{String,converter}}()
    for conv_sets in conv_table
        conv_dict=Dict{String,converter}()
        push!(conv_dict,("offshore"=>conv_sets[1]))
        push!(conv_dict,("onshore"=>conv_sets[2]))
        push!(all_conv_dict,(string(conv_sets[1].mva)=>conv_dict))
    end
    return all_xfo_dict,all_conv_dict
end
################################################################################
############################# Database Bits ####################################
################################################################################
#checks if HVDC is cheapest option for longest ranges at each power for a PCC connection
#**
function check4hvdc(equipment_database,km_mva_set,wnd,kv)
    hvdc_bit=false
    for km_mva in km_mva_set
        eqps=optimal_mog2pcc_wMOG(km_mva[1],km_mva[2],kv,wnd,equipment_database)
        cheapest=argmin([eqps["300.0"]["cost"],eqps["220.0"]["cost"],eqps["400.0"]["cost"]])
        cheapest_connection_2pcc=[eqps["300.0"],eqps["220.0"],eqps["400.0"]][cheapest]
        if (haskey(cheapest_connection_2pcc,"cable"))
            if (cheapest_connection_2pcc["cable"].elec.volt==300.0)
                hvdc_bit=true
            end
        end
    end
    return hvdc_bit
end

#checks if mid point compensation is cheapest option for longest ranges at each power for a PCC connection
#**
function check4mpc_ac(equipment_database,km_mva_set,wnd,kv)
    mpc_bit=false
    for km_mva in km_mva_set
        eqps=optimal_mog2pcc_wMOG(km_mva[1],km_mva[2],kv,wnd,equipment_database)
        cheapest=argmin([eqps["300.0"]["cost"],eqps["220.0"]["cost"],eqps["400.0"]["cost"]])
        cheapest_connection_2pcc=[eqps["300.0"],eqps["220.0"],eqps["400.0"]][cheapest]
        if (haskey(cheapest_connection_2pcc,"cable"))
            if (cheapest_connection_2pcc["cable"].mpc_ac==true)
                mpc_bit=true
            end
        end
    end
    return mpc_bit
end

#checks if hvac is needed if HVDC is possible
#**
function check4hvac(equipment_database,km_mva_set,wnd,kv)
    hvac_bit=false
    for km_mva in km_mva_set
        eqps=optimal_mog2pcc_wMOG(1,km_mva[2],kv,wnd,equipment_database)
        cheapest=argmin([eqps["300.0"]["cost"],eqps["220.0"]["cost"],eqps["400.0"]["cost"]])
        cheapest_connection_2pcc=[eqps["300.0"],eqps["220.0"],eqps["400.0"]][cheapest]
        if (haskey(cheapest_connection_2pcc,"plat_aft_ac"))
            if (cheapest_connection_2pcc["plat_aft_ac"].acdc=="ac")
                hvac_bit=true
            end
        end
    end
    return hvac_bit
end
