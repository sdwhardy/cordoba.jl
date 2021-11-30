include("struct.jl")#
include("input_data/functions.jl")#
include("wind/functions.jl")#
include("wind/wind_module.jl")#
include("eens/functions.jl")#
include("database/functions.jl")
#include("input_data/test_cases.jl")#

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#NOTE EENS is set to zero for flexlan !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

################################################################################
############################ Support Functions #################################
################################################################################
# Equipment Cost Calculations are from Techno-Economic Analysis of HVAC, HVDC and OFAC Offshore Wind Power Connections. S. Hardy, K. Van Brusselen, S. Hendrix
# and suppliment data is from North sea grid Annex

################################### Cables #####################################

#calculates the cost of a given hvac cost
#**
function cost_hvac_cable(cbl,ks)
    #cost of losses in the cable
    cbl.costs.rlc=cost_rlc(cbl,ks)
    #cost and size of cable compensation placed on OSS - divide by 2 for each OSS or PCC
    cbl.costs.qc,cbl.reactors=cost_qc_hvac(cbl,ks)
    #cost of switchgear placed on OSS
    cbl.costs.sg=cost_sg_hvac(cbl)
    #capex of cable
    cbl.costs.cpx_p,cbl.costs.cpx_i=capex_cable(cbl,ks.lac)
    #cost of corrective maintenance
    cbl.costs.cm=cost_cm(cbl.costs.cpx_p,ks.opx_c)
    #cost of expected energy not served
    cbl.costs.eens=cost_eens(cbl,ks)
    #totals the cable cost
    cbl.costs.ttl=cost_cbl_sum(cbl)
    cbl.costs.grand_ttl=cbl.costs.ttl
    return cbl
end

##MVAC
#calculates the cost of a given mvac cable
#**
function cost_mvac_cable(cbl,ks)
    #cost of losses in the cable
    cbl.costs.rlc=cost_rlc(cbl,ks)
    #cost of switchgear placed on OSS and OWPP
    cbl.costs.sg=cost_sg_hvac(cbl)
    #capex of cable
    cbl.costs.cpx_p,cbl.costs.cpx_i=capex_cable(cbl,ks.lmac)
    #cost of corrective maintenance
    cbl.costs.cm=cost_cm(cbl.costs.cpx_p,ks.opx_c)
    #cost of expected energy not served
    cbl.costs.eens=cost_eens(cbl,ks)
    #totals the cable cost
    cbl.costs.ttl=cost_cbl_sum(cbl)
    cbl.costs.grand_ttl=cbl.costs.ttl
    return cbl
end
#calculates the cost of a given hvdc cable
#**
function cost_hvdc_cable(cbl,ks)
    #cost of losses in the hvdc cable
    cbl.costs.rlc=cost_rlc_hvdc(cbl,ks)
    #capex of cable
    cbl.costs.cpx_p,cbl.costs.cpx_i=capex_cable(cbl,ks.ldc)
    #cost of corrective maintenance
    cbl.costs.cm=cost_cm(cbl.costs.cpx_p,ks.opx_c)
    #cost of expected energy not served
    cbl.costs.eens=cost_eens(cbl,ks)
    #totals the cable cost
    cbl.costs.ttl=cost_cbl_sum(cbl)
    cbl.costs.grand_ttl=cbl.costs.ttl
    return cbl
end

#CAPEX of cables
#**
function capex_cable(cbl,lay)
    cpx_p=cbl.length*cbl.num*cbl.costs.perkm_cpx#1.775 is installation cost per core per km ref: North Sea Grid Final Report Annex
    k=0.985^(cbl.num)
    cpx_i=cbl.length*cbl.num*lay*k#1.775=DC, 1.905=HVAC, 1.69=MVAC
    return cpx_p,cpx_i
end

#AC RLC Calculation
#**
function cost_rlc(cbl,ks)
    #AC cable loss cost **
    eta=0.994#assumed efficiency of transformers
    A=cbl.mva*eta#eta is the efficiency of a feeder transformer, s power to transmit
    B=cbl.elec.volt*sqrt(3)
    I=A/B#cable current
    R=(cbl.length*cbl.elec.ohm)/(cbl.num)#cable resistance
    #delta is related to wind profile, T_op lifetime hours and E_op cost of energy
    rlc=I^2*R*npv_hours()*ks.E_op*wind_module.wind_profs[cbl.wnd].delta#I^2R losses times cost factors
    return rlc
end

#AC RLC Calculation
#**
function cost_rlc_hvdc(cbl,ks)
    #AC cable loss cost **
    eta=0.994#assumed efficiency of transformers
    A=cbl.mva*eta#eta is the efficiency of a feeder transformer, s power to transmit
    B=cbl.elec.volt*2
    I=A/B#cable current
    R=(cbl.length*cbl.elec.ohm)/(cbl.num)#cable resistance
    #delta is related to wind profile, T_op lifetime hours and E_op cost of energy
    rlc=I^2*R*npv_hours()*ks.E_op*wind_module.wind_profs[cbl.wnd].delta#I^2R losses times cost factors
    return rlc
end

#HVAC compensation cost oss/PCC
#**
function cost_qc_hvac(cbl,ks)
    f=cbl.elec.freq
    A=cbl.elec.farrad*cbl.length*cbl.num
    Q=2*pi*f*cbl.elec.volt^2*A
    qc,reactors=get_reactors(Q/2)
    opx=(qc*0.0015)*2*npv_years()
    cbl.costs.cm=cbl.costs.cm+opx
    return 2*qc,reactors
end

#corrective maintenance of all equipment
#North Sea Grid OPEX values used
#**
function cost_cm(cpx,k)
    cm=cpx*k*npv_years()
    return cm
end

#Cost of switch gear is added to the cables directly
#**
function cost_sg_hvac(cbl)
    sg=0
    if (cbl.elec.volt==400.0)
        sg=2*4.545
    elseif (cbl.elec.volt==220.0)
        sg=2*3.09
    elseif (cbl.elec.volt==132.0)
        sg=2*2.4
    elseif (cbl.elec.volt==66.0)
        sg=2*1.8
    else
        println(cbl.size)
        println(cbl.num)
        println(cbl.mva)
        println(cbl.costs.ttl)
        println(cbl.elec.mva)
        println("Problem Placing Switchgear")
    end
    #should add OPEX
    opx=(sg*0.007)*npv_years()
    cbl.costs.cm=cbl.costs.cm+opx
    return sg
end

#sums all cable costs and returns the total**
#**
function cost_cbl_sum(cbl)
    ttl=cbl.costs.rlc+cbl.costs.cpx_i+cbl.costs.cpx_p+cbl.costs.cm+cbl.costs.qc+cbl.costs.eens+cbl.costs.sg
    if (ttl==NaN)
        ttl=Inf
    end
    cbl.costs.grand_ttl=ttl
    return ttl
end

################################# Converters ###################################
#totals the cost of HVDC converter
#**
function cost_conv_sum(conv)
    ttl=conv.costs.tlc+conv.costs.cpx+conv.costs.cm+conv.costs.eens
    return ttl
end
#capex of a converter
#**
function capex_hvdc(conv,ks)
    conv.costs.cpx = ks.conv_c*conv.mva
    return conv.costs.cpx
end

################################ Transformers ##################################
##Transformers
#calculates the cost of a given transformer on an OSS
#**
function cost_xfo_oss(xfo,ks)
    #capex oss
    xfo.costs.cpx_p,xfo.costs.cpx_i=capex_oss(xfo,ks)
    #cost of losses
    xfo.costs.tlc=cost_tlc(xfo,ks)
    #corrective maintenance
    xfo.costs.cm=cost_cm(xfo.costs.cpx_p,ks.opx_x)
    #eens calculation
    xfo.costs.eens=cost_eens(xfo,ks)
    #totals the xfo cost
    xfo.costs.ttl=cost_xfo_sum(xfo)
    return xfo
end

#calculates the cost of a given transformer at the PCC
#**
function cost_xfo_pcc(xfo,ks)
    #capex pcc
    xfo.costs.cpx_p,xfo.costs.cpx_i=capex_pcc(xfo)
    #cost of losses
    xfo.costs.tlc=cost_tlc(xfo,ks)
    #corrective maintenance
    xfo.costs.cm=cost_cm(xfo.costs.cpx_p,ks.opx_x)
    #eens calculation
    xfo.costs.eens=cost_eens(xfo,ks)
    #totals the xfo cost
    xfo.costs.ttl=cost_xfo_sum(xfo)
    return xfo
end

#Offshore transformer CAPEX
#**
function capex_oss(xfo,ks)
    A=(1+ks.dc*(xfo.num-2))
    B=(ks.f_ct+ks.p_ct)
    cpx=A*B*xfo.num*xfo.elec.mva
    return cpx/2,cpx/2#cost is split 50-50 over procurement and installation
end

#PCC transformer CAPEX
#**
function capex_pcc(xfo)
    cpx=xfo.num*0.0071*(xfo.elec.mva)^1.05
    return cpx*0.8,cpx*0.2#cost is split 80-20 over procurement and installation
end

#xfo Losses Calculation
#power factor of 1 assumed
#**
function cost_tlc(xfo,ks)
    pf=1
    tlc=xfo.mva*pf*(1-xfo.eta)*npv_hours()*ks.E_op*wind_module.wind_profs[xfo.wnd].delta
    return tlc
end

#sums all xfo costs and returns the total**
#**
function cost_xfo_sum(xfo)
    ttl=xfo.costs.tlc+xfo.costs.cpx_p+xfo.costs.cpx_i+xfo.costs.cm+xfo.costs.eens
    return ttl
end

################################### platform ###################################

#depth ref A review of foundations of offshore wind energy convertors: Current status and future perspectives
#AC platform incremental capex
#**
#depth scale factors from Europe's onshore and offshore wind energy potential An assessment of environmental and economic constraints
#EEA Technical report No 6/2009, ISSN 1725-2237
#10–20 20–30 30–40 40–50
#1.000 1.067 1.237 1.396
#multiple is approximated by 0.0136*depth_in_meters + 0.7676
function capex_plat_ac(plat,ks)
    plat.costs.cpx=ks.pac_e*plat.mva#*(plat.depth/15)#15m set as base depth - original with ks.pac_e=16
    #plat.costs.cpx=ks.pac_e*plat.mva^1.3#new to reflect better the base cost
    return plat.costs.cpx
end

#depth ref A review of foundations of offshore wind energy convertors: Current status and future perspectives
#DC platform incremental capex
#**
function capex_plat_dc(plat,ks)
    plat_base=platform()
    plat_base.mva=plat.mva
    plat_base.wnd=plat.wnd
    plat.costs.cpx=(ks.pdc_g*plat.mva)-capex_plat_ac(plat,ks)#*(plat.depth/15)#15m set as base depth
    return plat.costs.cpx
end

############################### Net Present value ##############################
#Changes to a net present value
#**
function npv_years()
    #1+1/(1.04)^1+1/(1.04)^2+...+1/(1.04)^25
    return 16.62
end

#Changes to a net present value
#**
function npv_hours()
    #16.62*365.25*24=145690.91
    return 145691
end
