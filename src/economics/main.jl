using XLSX, DataFrames
include("functions.jl")
#This file contains the exported functions
################################################################################
##########################  Default wind profile ###############################
df = DataFrame(XLSX.readtable("../cordoba/src/economics/input_data/data.xlsx", "wind_data")...)
wnd=wndF_wndPrf([getproperty(df,Symbol("Norther"))])
wind_module.save_wind4_module(wnd,"Norther")
################################################################################
########################### Base equipment functions  ##########################
################################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#NOTE EENS is set to zero for flexlan !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Cost of HVAC cable of capacity mva and length km
function AC_cbl(mva,km)
    cbl=cable();cbl.mva=mva;cbl.wnd="Norther"
    ac_cbl=optimal_ac_cable(cbl,get_220kV_cables(),km,get_Cost_Data())
    if (ac_cbl.num<1); ac_cbl.costs.cpx_p=1e12; end
    return ac_cbl
end
#Cost of HVDC cable of capacity mva and length km

function DC_cbl(mva,km)
    cbl=cable();cbl.mva=mva;cbl.wnd="Norther"
    dc_cbl=optimal_hvdc_cable(cbl,get_300kV_cables(),km,get_Cost_Data())
    return dc_cbl
end
#Cost of AC platform of capacity mva and depth m
function AC_plat(mva, m::Float64=17.1)
    plat=platform();
    plat.mva=mva
    plat.wnd="Norther"
    plat=cost_ac_platform(plat,get_Cost_Data())
    plat=adjust_base_ac_platform(plat,get_Cost_Data())
    #adjusts costs based on depth default and neutral depth is 17.1m
    multiplier=0.0136*abs(m) + 0.7676
    plat.costs.cm=plat.costs.cm*multiplier
    plat.costs.cpx=plat.costs.cpx*multiplier
    plat.costs.ttl=plat.costs.cpx+plat.costs.cm
    return plat
end
#Cost of DC platform of capacity mva and depth m

function DC_plat(mva, m::Float64=17.1)
    plat=platform();
    plat.mva=mva
    plat.wnd="Norther"
    plat=cost_dc_platform(plat,get_Cost_Data())
    plat=adjust_base_dc_platform(plat,get_Cost_Data())
    #adjusts costs based on depth default and neutral depth is 17.1m
    multiplier=0.0136*abs(m) + 0.7676
    plat.costs.cm=plat.costs.cm*multiplier
    plat.costs.cpx=plat.costs.cpx*multiplier
    plat.costs.ttl=plat.costs.cpx+plat.costs.cm
    return plat
end
#offshore transformer of capacity mva
function off_xfm(mva)
    xfo_off=transformer()
    xfo_off.relia=get_offshore_xfo_failure_data(xfo_off)
    xfo_off.wnd="Norther"
    xfo_off.mva=mva
    xfo_off=xfo_oss(xfo_off,get_Cost_Data(),get_Xfo_Data())
    return xfo_off
end
#onshore transformer of capacity mva
function on_xfm(mva)
    xfo_on=transformer()
    xfo_on.relia=get_offshore_xfo_failure_data(xfo_on)
    xfo_on.wnd="Norther"
    xfo_on.mva=mva
    xfo_pcc(xfo_on,get_Cost_Data(),get_Xfo_Data())
    return xfo_on
end
#cost of offshore converter of size mva
function off_conv(mva)
    conv_off=converter()
    conv_off.relia=get_offshore_conv_failure_data(conv_off)
    conv_off.wnd="Norther"
    conv_off.mva=mva
    conv_off.elec.mva=mva
    conv_off.num=1
    conv_off=cost_hvdc_oss(conv_off,get_Cost_Data())
    conv_off=adjust_base_hvdc_offshore_converter(conv_off,get_Cost_Data())
    return conv_off
end
#cost of onshore converter of size mva
function on_conv(mva)
    conv_on=converter()
    conv_on.relia=get_onshore_conv_failure_data(conv_on)
    conv_on.wnd="Norther"
    conv_on.mva=mva
    conv_on.elec.mva=mva
    conv_on.num=1
    conv_on=cost_hvdc_pcc(conv_on,get_Cost_Data())
    conv_on=adjust_base_hvdc_onshore_converter(conv_on,get_Cost_Data())
    return conv_on
end
################################################################################
########################### Finding MV max distance ############################
################################################################################
#**
function find_max_mv_transmission(owpps,database)
    ks=get_Cost_Data()
    for owpp in owpps
        pnt0=0.001
        pnt1=deepcopy(database["cables"]["66.0"][string(owpp.mva)][length(database["cables"]["66.0"][string(owpp.mva)])].length)
        mv_cable_pnt0=deepcopy(mvac_cable(owpp.mva,pnt0,owpp.wnd,database["cables"]["66.0"][string(owpp.mva)],ks))
        mv_cable_pnt1=deepcopy(mvac_cable(owpp.mva,pnt1,owpp.wnd,database["cables"]["66.0"][string(owpp.mva)],ks))
        plat=platform();plat.costs.ttl=0.0
        plat.mva=owpp.mva
        plat.wnd=owpp.wnd
        plat=cost_ac_platform(plat,ks)
        plat=adjust_base_ac_platform(plat,ks)
        ####################################### Addition
        xfo_off=transformer()
        xfo_off.relia=get_offshore_xfo_failure_data(xfo_off)
        xfo_off.wnd=owpp.wnd
        xfo_off.mva=owpp.mva
        xfo_off=xfo_oss(xfo_off,ks,get_Xfo_Data())
        ########################################
        plat.costs.ttl=plat.costs.ttl+xfo_off.costs.ttl
        np=0
        if (database["bits"]["hvac"]==true || database["bits"]["mpc_ac"]==true)
            hv_cable_pnt0_220=deepcopy(hvac_cable(owpp.mva,pnt0,owpp.wnd,database["cables"]["220.0"][string(owpp.mva)],ks))
            hv_cable_pnt0_400=deepcopy(hvac_cable(owpp.mva,pnt0,owpp.wnd,database["cables"]["400.0"][string(owpp.mva)],ks))
            if (hv_cable_pnt0_220.costs.ttl<hv_cable_pnt0_400.costs.ttl)
                hv_cable_pnt0=hv_cable_pnt0_220
                hv_c="220.0"
            else
                hv_cable_pnt0=hv_cable_pnt0_400
                hv_c="400.0"
            end
        elseif (database["bits"]["hvdc"]==true)
            hv_c="300.0"
        end

        if (mv_cable_pnt0.costs.ttl<hv_cable_pnt0.costs.ttl+plat.costs.ttl)
            hv_cable_pnt1=deepcopy(hvac_cable(owpp.mva,pnt1,owpp.wnd,database["cables"][hv_c][string(owpp.mva)],ks))
            #if AC cable range is maxed out this backs it up to a feasible distance
            while (hv_cable_pnt1.costs.ttl==Inf)
                pnt1=pnt1-25
                hv_cable_pnt1=deepcopy(hvac_cable(owpp.mva,pnt1,owpp.wnd,database["cables"][hv_c][string(owpp.mva)],ks))
                mv_cable_pnt1=deepcopy(mvac_cable(owpp.mva,pnt1,owpp.wnd,database["cables"]["66.0"][string(owpp.mva)],ks))
            end
            np=pnt1
            if (mv_cable_pnt1.costs.ttl>hv_cable_pnt1.costs.ttl+plat.costs.ttl)
                while (pnt1-pnt0>1)
                    np=(pnt0+pnt1)/2
                    mv_cable_np=deepcopy(mvac_cable(owpp.mva,np,owpp.wnd,database["cables"]["66.0"][string(owpp.mva)],ks))
                    hv_cable_np=deepcopy(hvac_cable(owpp.mva,np,owpp.wnd,database["cables"][hv_c][string(owpp.mva)],ks))
                    if (mv_cable_np.costs.ttl<hv_cable_np.costs.ttl+plat.costs.ttl)
                        pnt0=deepcopy(np)
                    elseif (mv_cable_np.costs.ttl>hv_cable_np.costs.ttl+plat.costs.ttl)
                        pnt1=deepcopy(np)
                    else
                        break
                    end
                end
            end
        end
        ############### addition #########
        # Account for collector circuit###
        ###### 6MW/km square assumed######
        #### Diameter of area travelled###
        collector=2*sqrt((owpp.mva/6)/pi)
        np_minus_collector=np-collector
        if (np_minus_collector<0)
            #np_minus_collector=0
            np_minus_collector=0.1
        end
        owpp.mv_zone=deepcopy(np_minus_collector)
        ##################################
        #owpp.mv_zone=deepcopy(np)
        owpp.kV=66.0
    end
    return owpps
end

################################################################################
############################## Finding lookup table ############################
################################################################################
#**
function get_equipment_tables(km_mva_set,wnd,kv)
    cable_database=get_cable_table(km_mva_set,wnd)
    xfo_database,conv_database=get_xfoConv_table(km_mva_set,wnd)
    database=Dict{String, Dict{String,Any}}()
    push!(database,"cables"=>cable_database)
    push!(database,"transformers"=>xfo_database)
    push!(database,"converters"=>conv_database)

    #vheck if hvdc is to be included
    bits=Dict([("hvdc", true), ("mpc_ac", false), ("hvac", true)])
    push!(database,"bits"=>bits)
    hvdc_bit=check4hvdc(database,km_mva_set,wnd,kv)

    #check hvac connections
    hvac_bit=true
    if (hvdc_bit==true)
        hvac_bit=check4hvac(database,km_mva_set,wnd,kv)
    end
    mpc_bit=false

    #final bits sets
    bits=Dict([("hvdc", hvdc_bit), ("mpc_ac", mpc_bit), ("hvac", hvac_bit)])
    push!(database,"bits"=>bits)
    database["cables"]=set_cables_per_km_cost(database["cables"],km_mva_set)
    return database
end

################################################################################
######################## Optimal equipment using lookup table ##################
################################################################################
#finds the offshore transformer using the look up table
#when increasing meshed transformer xfo_oss_meshed is used
#**
function xfo_oss(xfo0,ks,xfo_data)
    xfo=transformer()
    xfo.costs.ttl=Inf
    xfo.mva=xfo0.mva
    for xd in xfo_data
        if ((xd+10<=xfo0.mva) && (6*xd>xfo0.mva))
            xfo0.num=1
            xfo0.elec.mva=xd
            while ((xfo0.num*xfo0.elec.mva)<xfo0.mva)
                xfo0.num=xfo0.num+1
            end
            xfo0=cost_xfo_oss(xfo0,ks)
            if (xfo0.costs.ttl<xfo.costs.ttl)
                xfo=deepcopy(xfo0)
            end
        end
    end
    return xfo
end

#finds the onshore transformer using the look up table
#**
function xfo_pcc(xfo0,ks,xfo_data)
    xfo=transformer()
    xfo.costs.ttl=Inf
    xfo.mva=xfo0.mva
    for xd in xfo_data
        if ((xd+10<=xfo0.mva) && ((3*xd>xfo0.mva) || ((xfo0.mva>2250) && (6*xd>xfo0.mva))))
            xfo0.num=1
            xfo0.elec.mva=xd
            while ((xfo0.num*xfo0.elec.mva)<xfo0.mva)
                xfo0.num=xfo0.num+1
            end
            xfo0=cost_xfo_pcc(xfo0,ks)
            if (xfo0.costs.ttl<xfo.costs.ttl)
                xfo=deepcopy(xfo0)
            end
        end
    end
    return xfo
end

#finds the hvdc cable using the look up table
#**
function hvdc_cable(mva,km,wnd,cable_array,ks)
    cbl=deepcopy(cable_array[1])
    cbl.wnd=wnd
    cbl.length=km
    cbl=cost_hvdc_cable(cbl,ks)
    return cbl
end

#finds the hvac cable using the look up table
#**
function hvac_cable(mva,km,wnd,cable_array,ks)
    cbl=cable()
    cbl.costs.ttl=Inf
    cbl.costs.grand_ttl=Inf
    for cbl0 in cable_array
        if (km<=cbl0.length+0.01)
            #+0.01 is just to be certain floating point error does not occur
            cbl=deepcopy(cbl0)
            break
        end
    end
    if (cbl.size!=0.0)
        cbl.wnd=wnd
        cbl.length=km
        cbl.elec.mva=get_newQ_Capacity(cbl.elec.freq,km,cbl.elec.volt,cbl.elec.farrad,cbl.elec.amp)
        cbl=get_cbl_failure_data(cbl)
        cbl=cost_hvac_cable(cbl,ks)
    end
    return cbl
end

#finds the mvac cable using the look up table
#**
function mvac_cable(mva,km,wnd,cable_array,ks)
    cbl=cable()
    cbl.costs.ttl=Inf
    cbl.costs.grand_ttl=Inf
    for cbl0 in cable_array
        if (km<=cbl0.length+0.001)#+0.001 is just to be certain floating point error does not occur
            cbl=deepcopy(cbl0)
            break
        end
    end
    max_in_parallel=12#if changing check economics/database/functions function optimal_ac_cable(cbl0,cbl_data,km,ks)
    if (cbl.num>max_in_parallel && cbl.size!=0.0)#66kv up to 2000MW eception
        cbl.mva=cbl.mva/2
        cbl.num=cbl.num/2
        cbl.costs.perkm_cpx=cbl.costs.perkm_cpx/2
        cbl.elec.mva=get_newQ_Capacity(cbl.elec.freq,km,cbl.elec.volt,cbl.elec.farrad,cbl.elec.amp)
        cbl.wnd=wnd
        cbl.length=km
        cbl=get_cbl_failure_data(cbl)
        cbl=cost_mvac_cable(cbl,ks)
        cbl.num=cbl.num*2
        cbl.mva=cbl.mva*2
        cbl.costs.perkm_cpx=cbl.costs.perkm_cpx*2
        cbl.costs.sg=cbl.costs.sg*2
        cbl.costs.cpx_p=cbl.costs.cpx_p*2
        cbl.costs.cpx_i=cbl.costs.cpx_i*2
        cbl.costs.rlc=cbl.costs.rlc*2
        cbl.costs.cm=cbl.costs.cm*2
        cbl.costs.eens=cbl.costs.eens*2
        cbl.costs.ttl=cbl.costs.ttl*2
        cbl.costs.grand_ttl=cbl.costs.ttl
    elseif (cbl.size!=0.0)
        cbl.elec.mva=get_newQ_Capacity(cbl.elec.freq,km,cbl.elec.volt,cbl.elec.farrad,cbl.elec.amp)
        cbl.wnd=wnd
        cbl.length=km
        cbl=get_cbl_failure_data(cbl)
        cbl=cost_mvac_cable(cbl,ks)
        if (cbl.costs.ttl==Inf)
            cbl.costs.perkm_ttl=10^9
            cbl.costs.grand_ttl=0.0
        else
            cbl.costs.grand_ttl=cbl.costs.ttl
        end
    end
    return cbl
end

#Converter
#calculates the cost of an offshore HVDC converter
#each line needing a converter shares the total before the fixed cost is added at the end
#**
function cost_hvdc_oss(conv,ks)
    #capex converter
    conv.costs.cpx=capex_hvdc(conv,ks)
    #cost of losses
    conv.costs.tlc=cost_tlc(conv,ks)
    #corrective maintenance
    conv.costs.cm=cost_cm(conv.costs.cpx,ks.opx_co)
    #eens calculation
    conv.costs.eens=cost_eens(conv,ks)
    #totals the xfo cost
    conv.costs.ttl=cost_conv_sum(conv)
    return conv
end

#once all incoming lines with converters are combined the fixed cost is added
#**
function adjust_base_hvdc_offshore_converter(conv,ks)
    #adjust to add base cost capex converter
    conv.costs.cpx=conv.costs.cpx+ks.conv_d
    #corrective maintenance
    conv.costs.cm=conv.costs.cm+ks.conv_d*ks.opx_co*npv_years()
    #totals the converter cost
    conv.costs.ttl=cost_conv_sum(conv)
    return conv
end

#calculates the cost of an onshore HVDC converter
#each line needing a converter shares the total before the fixed cost is added at the end
#**
function cost_hvdc_pcc(conv,ks)
    #capex converter
    conv.costs.cpx=capex_hvdc(conv,ks)
    #cost of losses
    conv.costs.tlc=cost_tlc(conv,ks)
    #corrective maintenance
    conv.costs.cm=cost_cm(conv.costs.cpx,ks.opx_cp)
    #eens calculation
    conv.costs.eens=cost_eens(conv,ks)
    #totals the xfo cost
    conv.costs.ttl=cost_conv_sum(conv)
    return conv
end

#once all incoming lines with converters are combined the fixed cost is added
#**
function adjust_base_hvdc_onshore_converter(conv,ks)
    #adjust to add base cost capex converter
    conv.costs.cpx=conv.costs.cpx+ks.conv_d
    #corrective maintenance
    conv.costs.cm=conv.costs.cm+ks.conv_d*ks.opx_cp*npv_years()
    #totals the converter cost
    conv.costs.ttl=cost_conv_sum(conv)
    return conv
end

################################### platform ###################################
#calculates incremental cost of an AC platform NOT including fixed cost
#**
function cost_ac_platform(plat,ks)
    #capex oss
    plat.costs.cpx=capex_plat_ac(plat,ks)
    #corrective maintenance
    plat.costs.cm=cost_cm(plat.costs.cpx,ks.opx_pl)
    #totals the plat cost
    plat.costs.ttl=plat.costs.cpx+plat.costs.cm
    return plat
end

#once all incoming lines established the fixed cost is added to the of incremental sum - AC
#**
function adjust_base_ac_platform(plat,ks)
    #adjust to add base cost capex converter
    plat.costs.cpx=plat.costs.cpx+ks.pac_f
    #corrective maintenance
    plat.costs.cm=plat.costs.cm+ks.pac_f*ks.opx_pl*npv_years()
    #totals the converter cost
    plat.costs.ttl=plat.costs.cpx+plat.costs.cm
    return plat
end

#calculates incremental cost of an DC platform NOT including fixed cost
#**
function cost_dc_platform(plat,ks)
    #capex oss
    plat.costs.cpx=capex_plat_dc(plat,ks)
    #corrective maintenance
    plat.costs.cm=cost_cm(plat.costs.cpx,ks.opx_pl)
    #totals the plat cost
    plat.costs.ttl=plat.costs.cpx+plat.costs.cm
    return plat
end

#once all incoming lines established the fixed cost is added to the of incremental sum - DC
#**
function adjust_base_dc_platform(plat,ks)
    #adjust to add base cost capex converter
    plat.costs.cpx=plat.costs.cpx+ks.pdc_h
    #corrective maintenance
    plat.costs.cm=plat.costs.cm+ks.pdc_h*ks.opx_pl*npv_years()
    #totals the converter cost
    plat.costs.ttl=plat.costs.cpx+plat.costs.cm
    return plat
end
