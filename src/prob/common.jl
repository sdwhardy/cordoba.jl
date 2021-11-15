
#=function gen(data,mva)
    data["pmax"]=mva
    data["pmin"]=mva
    return data
end

function converterdc_ne(conv,nde,mbase)
    #NOTE change
    conv["cost"]=nde["xfm_off"].costs.cpx_p+nde["xfm_off"].costs.cpx_i+nde["conv_off"].costs.cpx+nde["plat"].costs.cpx+nde["cableWF"].costs.cpx_p+nde["cableWF"].costs.cpx_i;
    conv["Pacmax"]=nde["conv_off"].mva/mbase
    conv["Pacmin"]=-1*nde["conv_off"].mva/mbase
    conv["Qacmax"]=nde["conv_off"].mva/mbase
    conv["Qacmin"]=-1*nde["conv_off"].mva/mbase
    return conv
end

function convdc_ne(conv,nde,mbase)
    #NOTE change
    conv["cost"]=nde["dc"]["xfm"].costs.cpx_p+nde["dc"]["xfm"].costs.cpx_i+nde["dc"]["conv"].costs.cpx
    #conv["cost"]=0
    if (haskey(nde["dc"],"plat")); conv["cost"]=conv["cost"]+nde["dc"]["plat"].costs.cpx;end
    #if (haskey(nde["dc"],"plat")); conv["cost"]=0;end
    #########
    conv["Pacmax"]=nde["mva"]/mbase
    conv["Pacmin"]=-1*nde["mva"]/mbase
    conv["Qacmax"]=nde["mva"]/mbase
    conv["Qacmin"]=-1*nde["mva"]/mbase
    return conv
end

function branchdc_on(brc,dcss,mbase)
    brc["cost"]=brc["cost"]+dcss["xfm_on"].costs.cpx_p+dcss["xfm_on"].costs.cpx_i+dcss["conv_on"].costs.cpx
    return brc
end

function converterdc(conv,nde,mbase)
    conv["Pacmax"]=nde.mva/mbase
    conv["Pacmin"]=-1*nde.mva/mbase
    conv["Qacmax"]=nde.mva/mbase
    conv["Qacmin"]=-1*nde.mva/mbase
    return conv
end=#

#=function branchdc_ne(b_ne,dc_can,mbase)
    b_ne["rateA"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["rateB"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["rateC"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["r"]=dc_can.elec.ohm

    #NOTE change
    b_ne["cost"]=dc_can.costs.cpx_p+dc_can.costs.cpx_i
    #b_ne["cost"]=0
    #######
    return b_ne
end=#

#=function ne_branch(ne_b,ac_can,mbase)
    ne_b["mva"]=ac_can.num*ac_can.elec.mva
    ne_b["rate_a"]=ac_can.num*ac_can.elec.mva/mbase
    ne_b["rate_b"]=ac_can.num*ac_can.elec.mva/mbase
    ne_b["rate_c"]=ac_can.num*ac_can.elec.mva/mbase
    ne_b["br_r"]=ac_can.elec.ohm
    ne_b["br_x"]=ac_can.elec.xl
    #data["ne_branch"][datai]["construction_cost"]=hoa["candidate_lines"]["ac"][hoai].costs.cpx_p+hoa["candidate_lines"]["ac"][hoai].costs.cpx_i+hoa["candidate_lines"]["ac"][hoai].costs.qc+hoa["candidate_lines"]["ac"][hoai].costs.sg
    #ne_b["construction_cost"]=(ac_can.costs.cpx_p+ac_can.costs.cpx_i)/mbase

    #NOTE change
    ne_b["construction_cost"]=ac_can.costs.cpx_p+ac_can.costs.cpx_i
    #ne_b["construction_cost"]=0
    ##########
    return ne_b
end=#

#=function branch_xfo(ne_b,ne_xfo,mbase)
    #ne_b["construction_cost"]=ne_b["construction_cost"]+(ne_xfo.costs.cpx_p+ne_xfo.costs.cpx_i)/mbase

    #NOTE change
    ne_b["construction_cost"]=ne_b["construction_cost"]+ne_xfo.costs.cpx_p+ne_xfo.costs.cpx_i
    #ne_b["construction_cost"]=ne_b["construction_cost"]+0
    ##########
    return ne_b
end=#

#=function branch_plat(ne_b,ne_plat,mbase)
    #ne_b["construction_cost"]=ne_b["construction_cost"]+(ne_plat.costs.cpx)/mbase

    #NOTE change
    ne_b["construction_cost"]=ne_b["construction_cost"]+ne_plat.costs.cpx
    #ne_b["construction_cost"]=ne_b["construction_cost"]+0
    ##########
    return ne_b
end=#

#=function pm_dict_candidateICs(data,hoa,candidates)
    #WF AC candidate lines
    for (i,k) in enumerate(vcat(candidates[2],candidates[2]))
        if (i<=length(candidates[2]))
            data["ne_branch"][string(i)]=ne_branch(data["ne_branch"][string(i)], hoa["candidates"]["wf"]["ac"][string(k)]["cableBE"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_xfo(data["ne_branch"][string(i)],hoa["candidates"]["wf"]["ac"][string(k)]["xfm_on"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_xfo(data["ne_branch"][string(i)],hoa["candidates"]["wf"]["ac"][string(k)]["xfm_off"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_plat(data["ne_branch"][string(i)], hoa["candidates"]["wf"]["ac"][string(k)]["plat"], data["gen"]["1"]["mbase"])
        else
            data["ne_branch"][string(i)]=ne_branch(data["ne_branch"][string(i)], hoa["candidates"]["wf"]["ac"][string(k)]["cableUK"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_xfo(data["ne_branch"][string(i)],hoa["candidates"]["wf"]["ac"][string(k)]["xfm_on"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_xfo(data["ne_branch"][string(i)],hoa["candidates"]["wf"]["ac"][string(k)]["xfm_off"], data["gen"]["1"]["mbase"])
            data["ne_branch"][string(i)]=branch_plat(data["ne_branch"][string(i)], hoa["candidates"]["wf"]["ac"][string(k)]["plat"], data["gen"]["1"]["mbase"])
        end
    end

    for (i,k) in enumerate(vcat(candidates[1],candidates[1]))
        if (i<=length(candidates[2]))
            data["branchdc_ne"][string(i)]=branchdc_ne(data["branchdc_ne"][string(i)],hoa["candidates"]["ic"][string(k)]["cableBE"],data["gen"]["1"]["mbase"])
            data["branchdc_ne"][string(i)]=branchdc_on(data["branchdc_ne"][string(i)],hoa["candidates"]["ic"][string(k)],data["gen"]["1"]["mbase"])
        else
            data["branchdc_ne"][string(i)]=branchdc_ne(data["branchdc_ne"][string(i)],hoa["candidates"]["ic"][string(k)]["cableUK"],data["gen"]["1"]["mbase"])
            data["branchdc_ne"][string(i)]=branchdc_on(data["branchdc_ne"][string(i)],hoa["candidates"]["ic"][string(k)],data["gen"]["1"]["mbase"])
        end
    end

    for (i,k) in enumerate(candidates[2])
        data["convdc_ne"][string(i)]=converterdc_ne(data["convdc_ne"][string(i)],hoa["candidates"]["wf"]["dc"][string(k)],data["gen"]["1"]["mbase"])
    end
    size_set=[];for (k,v) in hoa["candidates"]["ic"]; push!(size_set,parse(Float64,k));end
    mx_key=string(maximum(size_set))
    data["convdc"]["1"]=converterdc(data["convdc"]["1"],hoa["candidates"]["ic"][mx_key]["conv_on"],data["gen"]["1"]["mbase"])
    data["convdc"]["2"]=converterdc(data["convdc"]["2"],hoa["candidates"]["ic"][mx_key]["conv_on"],data["gen"]["1"]["mbase"])

    return data
end=#


#


#=function additional_candidatesWFSdc(data,candidate_wfs,wf_data)
    #DC, IC
    wfs=[]
    for (i,dcb) in data["branchdc_ne"]; push!(wfs,deepcopy(dcb));end

    data["branchdc_ne"]=Dict{String,Any}()
    for (i,wf) in enumerate(wfs)
        for j=1:1:length(candidate_wfs);
            wf["source_id"][2]=j+length(candidate_wfs)*(i-1);
            wf["index"]=j+length(candidate_wfs)*(i-1);
            wf["rateA"]=candidate_wfs[j];
            wf["rateB"]=first(wf_data[i]);
            wf["rateC"]=last(wf_data[i]);
            push!(data["branchdc_ne"],string(j+length(candidate_wfs)*(i-1))=>deepcopy(wf));
        end
    end
    return data
end=#

#=function additional_candidatesWFSac(data,candidate_wfs,wf_data)
    #DC
    #IC
    wfs=[]
    for (i,acb) in data["ne_branch"]; push!(wfs,deepcopy(acb));end

    data["ne_branch"]=Dict{String,Any}()
    for (i,wf) in enumerate(wfs)
        for j=1:1:length(candidate_wfs);
            wf["source_id"][2]=j+length(candidate_wfs)*(i-1);
            wf["index"]=j+length(candidate_wfs)*(i-1);
            wf["rate_a"]=candidate_wfs[j];
            wf["rate_b"]=first(wf_data[i]);
            wf["rate_c"]=last(wf_data[i]);
            push!(data["ne_branch"],string(j+length(candidate_wfs)*(i-1))=>deepcopy(wf));
        end
    end
    return data
end=#

#=function additional_candidatesWFS(data,candidates,ic_data)
    #DC
    #IC
    ics=[]
    for (i,dcb) in data["branchdc_ne"]; push!(ics,deepcopy(dcb));end

    data["branchdc_ne"]=Dict{String,Any}()
    for (i,ic) in enumerate(ics)
        for j=1:1:length(candidates);
            ic["source_id"][2]=j+length(candidates)*(i-1);
            ic["index"]=j+length(candidates)*(i-1);
            ic["rateA"]=candidates[j];
            ic["rateB"]=first(ic_data[i]);
            ic["rateC"]=last(ic_data[i]);
            push!(data["branchdc_ne"],string(j+length(candidates)*(i-1))=>deepcopy(ic));
        end
    end
    return data
end=#

#=function additional_candidates(data,candidates)
    #AC
    c2be=deepcopy(data["ne_branch"]["1"]);c2uk=deepcopy(data["ne_branch"]["2"])
    data["ne_branch"]=Dict{String,Any}()
    #WF
    for i=1:1:length(candidates[2]); c2be["source_id"][2]=i;c2be["index"]=i; push!(data["ne_branch"],string(i)=>deepcopy(c2be)); end
    for i=length(candidates[2])+1:1:length(candidates[2])*2; c2uk["source_id"][2]=i;c2uk["index"]=i; push!(data["ne_branch"],string(i)=>deepcopy(c2uk)); end
    #DC
    #WF
    convwf=deepcopy(data["convdc_ne"]["1"])
    data["convdc_ne"]=Dict{String,Any}()
    for i=1:1:length(candidates[2]); convwf["source_id"][2]=i;convwf["index"]=i; push!(data["convdc_ne"],string(i)=>deepcopy(convwf)); end
    #IC
    c2be=deepcopy(data["branchdc_ne"]["1"]);c2uk=deepcopy(data["branchdc_ne"]["2"])
    data["branchdc_ne"]=Dict{String,Any}()
    for i=1:1:length(candidates[1]); c2be["source_id"][2]=i;c2be["index"]=i; push!(data["branchdc_ne"],string(i)=>deepcopy(c2be)); end
    for i=length(candidates[1])+1:1:length(candidates[1])*2; c2uk["source_id"][2]=i;c2uk["index"]=i; push!(data["branchdc_ne"],string(i)=>deepcopy(c2uk)); end
    return data
end=#

#calculates candidate equipment and sorts in struct
#=function hoa_datastruct_candidateIC_w_onshore_converters(hoa,ic_mva,ic_length,candis,f,t)

    #datastructure
    if !(haskey(hoa,"nodes")); push!(hoa,"nodes"=>Dict{String,Any}()); end
    if !(haskey(hoa["nodes"],"dc")); push!(hoa["nodes"],"dc"=>Dict{String,Any}()); end
    if !(haskey(hoa["nodes"]["dc"],string(f))); push!(hoa["nodes"]["dc"],string(f)=>Dict{String,Any}()); end
    if !(haskey(hoa["nodes"]["dc"],string(t))); push!(hoa["nodes"]["dc"],string(t)=>Dict{String,Any}()); end

    #rating of nodes
    if !(haskey(hoa["nodes"]["dc"][string(f)],"mva")); push!(hoa["nodes"]["dc"][string(f)],"mva"=>0.0);end
    if !(haskey(hoa["nodes"]["dc"][string(t)],"mva")); push!(hoa["nodes"]["dc"][string(t)],"mva"=>0.0);end
    hoa["nodes"]["dc"][string(f)]["mva"]=hoa["nodes"]["dc"][string(f)]["mva"]+ic_mva
    hoa["nodes"]["dc"][string(t)]["mva"]=hoa["nodes"]["dc"][string(t)]["mva"]+ic_mva

    #interconnector candidate lines
    if !(haskey(hoa,"candidates")); push!(hoa,"candidates"=>Dict{String,Any}()); end
    if !(haskey(hoa["candidates"],"ic")); push!(hoa["candidates"],"ic"=>Dict{String,Any}()); end
    if !(haskey(hoa["candidates"]["ic"],string(f)*"_"*string(t))); push!(hoa["candidates"]["ic"],string(f)*"_"*string(t)=>Dict{String,Any}()); end
    [push!(hoa["candidates"]["ic"][string(f)*"_"*string(t)],string(p)=>Dict{String,Any}()) for p in candis]

    #ic dc
    [push!(hoa["candidates"]["ic"][string(f)*"_"*string(t)][k],"xfm_on"=>on_xfm(ic_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["ic"][string(f)*"_"*string(t)]]
    [push!(hoa["candidates"]["ic"][string(f)*"_"*string(t)][k],"conv_on"=>on_conv(ic_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["ic"][string(f)*"_"*string(t)]]
    [push!(hoa["candidates"]["ic"][string(f)*"_"*string(t)][k],"cable"=>DC_cbl(ic_mva/2*parse(Float64,k), ic_length)) for (k,p) in hoa["candidates"]["ic"][string(f)*"_"*string(t)]]

    return hoa
end=#

#=n=2;mva=2000
off_xfm(mva).costs.cpx_p+off_xfm(mva).costs.cpx_i+off_conv(mva).costs.cpx+DC_plat(mva).costs.cpx
(1+0.25*(n-2))*0.084*mva + 0.0589*mva+54.895 + 0.04*mva+165
(1+0.25*(n-2))*0.08*mva + 0.1*mva+220
(1+0.25*(n-2))*0.0975*mva + 0.1925*mva#offshore
0.1925*mva#onshore
on_xfm(mva).costs.cpx_p+on_xfm(mva).costs.cpx_i+on_conv(mva).costs.cpx
0.01*mva+0.06*mva+55=#
#=
function on_off_ic(mva,km)
    #ic dc
    cb=DC_cbl(mva, km)
    xf_on=on_xfm(cb.num*cb.elec.mva)
    cn_on=on_conv(cb.num*cb.elec.mva)
    xf_off=off_xfm(cb.num*cb.elec.mva)
    cn_off=off_conv(cb.num*cb.elec.mva)
    dc_plat=DC_plat(cb.num*cb.elec.mva)
    #println("DC connection cost: cb "*string(cb.costs.cpx_i+cb.costs.cpx_p)*" x_off "*string(xf_off.costs.cpx_p+xf_off.costs.cpx_i)*" x_on "*string(xf_on.costs.cpx_p+xf_on.costs.cpx_i)*" c_off "*string(cn_off.costs.cpx)*" c_on "*string(cn_off.costs.cpx)*" plat "*string(dc_plat.costs.cpx))
    #xf=on_xfm(mva)
    #cn=on_conv(mva)
    transformer()
    #################### Original Costs
    #cost=dc_plat.costs.cpx+xf_off.costs.cpx_p+xf_off.costs.cpx_i+cn_off.costs.cpx+xf_on.costs.cpx_p+xf_on.costs.cpx_i+cn_on.costs.cpx+cb.costs.cpx_i+cb.costs.cpx_p
    #################### SEST Costs
    cost=0.1925*(cb.num*cb.elec.mva)*2+((1+0.25*(xf_off.num-2))*0.0975*cb.num*cb.elec.mva)+cb.costs.cpx_i+cb.costs.cpx_p
    println("DC connection cost: total "*string(cost)*" mva "*string(cb.num*cb.elec.mva)*" km "*string(km))
    return cost,(cb.elec.ohm/cb.num)*km,cb.num*cb.elec.mva
    #return cost,(cb.elec.ohm/cb.num)*km,mva
end
=#

#=
function candidateACWF_cost(owpp_mva,bac)
    if (bac["rate_c"]<0)#on-on
        bac["construction_cost"],bdc["r"],bac["rate_a"]=on_on_ic(owpp_mva*bac["rate_a"],bac["rate_b"])
    elseif (bac["rate_c"]>0)#off-off
        bac["construction_cost"],bac["r"],bac["rate_a"]=on_off_wfac(owpp_mva*bac["rate_a"],bac["rate_b"])
    else#on-off
        bac["construction_cost"],bac["r"],bac["rate_a"]=off_off_ic(owpp_mva*bac["rate_a"],bac["rate_b"])
    end
    bac["rate_a"]=bac["rate_c"]=bac["rate_b"]=bac["rate_a"]/bac["mva"]
    return bac
end

function on_off_wfac(mva,km)
    #ic dc
    cb=AC_cbl(mva, km)
    xf_on=on_xfm(cb.num*cb.elec.mva)
    xf_off=off_xfm(cb.num*cb.elec.mva)
    ac_plat=AC_plat(cb.num*cb.elec.mva)
    #xf=on_xfm(mva)
    #cn=on_conv(mva)
    cost=ac_plat.costs.cpx+xf_on.costs.cpx_p+xf_on.costs.cpx_i+xf_off.costs.cpx_p+xf_off.costs.cpx_i+cb.costs.cpx_i+cb.costs.cpx_p
    println("AC connection cost: total "*string(cost)*" mva "*string(cb.num*cb.elec.mva))
    return cost,(cb.elec.ohm/cb.num)*km,cb.num*cb.elec.mva
    #return cost,(cb.elec.ohm/cb.num)*km,mva
end
#cb=DC_cbl(2000, 550)
#cb.costs.cpx_i+cb.costs.cpx_p+(0.1925*cb.num*cb.elec.mva)*2
=#
#=
function on_on_ic(mva,km)
    #ic dc
    cb=DC_cbl(mva, km)
    xf=on_xfm(cb.num*cb.elec.mva)
    cn=on_conv(cb.num*cb.elec.mva)

    #println("DC connection cost: cb "*string(cb.costs.cpx_i+cb.costs.cpx_p)*" x_on "*string(xf.costs.cpx_p+xf.costs.cpx_i)*" c_on "*string(cn.costs.cpx))
    #xf=on_xfm(mva)
    #cn=on_conv(mva)


    #xf=on_xfm(mva)
    #cn=on_conv(mva)
    ##################### Original cost
    #cost=(xf.costs.cpx_p+xf.costs.cpx_i+cn.costs.cpx)*2+cb.costs.cpx_i+cb.costs.cpx_p

    ##################### for SEST
    cost=(0.1925*cb.num*cb.elec.mva)*2+cb.costs.cpx_i+cb.costs.cpx_p
    println("DC connection cost: total "*string(cost)*" mva "*string(cb.num*cb.elec.mva)*" km "*string(km))
    return cost,(cb.elec.ohm/cb.num)*km,cb.num*cb.elec.mva
    #return cost,(cb.elec.ohm/cb.num)*km,mva
end
=#
#=
#calculates candidate equipment and sorts in struct
function hoa_datastruct_candidateICs(ic_mva,owpp_mva,ic_length,owpp_km,candidates)
    #location
    uk_pcc=xy();uk_pcc.x=ic_length;uk_pcc.y=0;
    be_pcc=xy();be_pcc.x=0;be_pcc.y=0;
    owpp_xy=xy();owpp_xy.x=owpp_km;owpp_xy.y=0
    #datastructure
    hoa=Dict{String,Any}();
    push!(hoa,"nodes"=>Dict{String,Any}());
    push!(hoa["nodes"],"ac"=>Dict{String,Any}());
    push!(hoa["nodes"],"dc"=>Dict{String,Any}());

    [push!(hoa["nodes"]["ac"],n=>Dict{String,Any}()) for n in ["1","2","3"]]
    [push!(hoa["nodes"]["dc"],n=>Dict{String,Any}()) for n in ["1","2","3","4"]]
    #position of nodes
    push!(hoa["nodes"]["ac"]["1"],"xy"=>be_pcc);
    push!(hoa["nodes"]["ac"]["2"],"xy"=>uk_pcc);
    push!(hoa["nodes"]["ac"]["3"],"xy"=>owpp_xy);
    #position of nodes
    push!(hoa["nodes"]["dc"]["1"],"xy"=>be_pcc);
    push!(hoa["nodes"]["dc"]["2"],"xy"=>uk_pcc);
    push!(hoa["nodes"]["dc"]["3"],"xy"=>owpp_xy);
    push!(hoa["nodes"]["dc"]["4"],"xy"=>owpp_xy);
    #rating of nodes
    push!(hoa["nodes"]["ac"]["1"],"mva"=>ic_mva);
    push!(hoa["nodes"]["ac"]["2"],"mva"=>ic_mva);
    push!(hoa["nodes"]["ac"]["3"],"mva"=>owpp_mva);
    #rating of nodes
    push!(hoa["nodes"]["dc"]["1"],"mva"=>ic_mva);
    push!(hoa["nodes"]["dc"]["2"],"mva"=>ic_mva);
    push!(hoa["nodes"]["dc"]["3"],"mva"=>owpp_mva);
    push!(hoa["nodes"]["dc"]["4"],"mva"=>owpp_mva);
    #lines
    push!(hoa,"candidates"=>Dict{String,Any}());
    push!(hoa["candidates"],"wf"=>Dict{String,Any}());
    push!(hoa["candidates"],"ic"=>Dict{String,Any}());
    push!(hoa["candidates"]["wf"],"ac"=>Dict{String,Any}());
    push!(hoa["candidates"]["wf"],"dc"=>Dict{String,Any}());

    [push!(hoa["candidates"]["wf"]["ac"],string(p)=>Dict{String,Any}()) for p in candidates[2]]
    [push!(hoa["candidates"]["wf"]["dc"],string(p)=>Dict{String,Any}()) for p in candidates[2]]
    [push!(hoa["candidates"]["ic"],string(p)=>Dict{String,Any}()) for p in candidates[1]]

    #wf-ac
    [push!(hoa["candidates"]["wf"]["ac"][k],"xfm_on"=>on_xfm(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    [push!(hoa["candidates"]["wf"]["ac"][k],"xfm_off"=>off_xfm(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    [push!(hoa["candidates"]["wf"]["ac"][k],"plat"=>AC_plat(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    [push!(hoa["candidates"]["wf"]["ac"][k],"xfm"=>on_xfm(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    [push!(hoa["candidates"]["wf"]["ac"][k],"cableBE"=>AC_cbl(owpp_mva*parse(Float64,k), owpp_km)) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    [push!(hoa["candidates"]["wf"]["ac"][k],"cableUK"=>AC_cbl(owpp_mva*parse(Float64,k), ic_length-owpp_km)) for (k,p) in hoa["candidates"]["wf"]["ac"]]
    #wf-dc
    [push!(hoa["candidates"]["wf"]["dc"][k],"xfm_off"=>off_xfm(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["dc"]]
    [push!(hoa["candidates"]["wf"]["dc"][k],"conv_off"=>off_conv(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["dc"]]
    [push!(hoa["candidates"]["wf"]["dc"][k],"plat"=>DC_plat(owpp_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["wf"]["dc"]]
    [push!(hoa["candidates"]["wf"]["dc"][k],"cableWF"=>DC_cbl(owpp_mva/2*parse(Float64,k), 1)) for (k,p) in hoa["candidates"]["wf"]["dc"]]
    #ic dc
    [push!(hoa["candidates"]["ic"][k],"xfm_on"=>on_xfm(ic_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["ic"]]
    [push!(hoa["candidates"]["ic"][k],"conv_on"=>on_conv(ic_mva*parse(Float64,k))) for (k,p) in hoa["candidates"]["ic"]]
    [push!(hoa["candidates"]["ic"][k],"cableBE"=>DC_cbl(ic_mva/2*parse(Float64,k), owpp_km)) for (k,p) in hoa["candidates"]["ic"]]
    [push!(hoa["candidates"]["ic"][k],"cableUK"=>DC_cbl(ic_mva/2*parse(Float64,k), ic_length-owpp_km)) for (k,p) in hoa["candidates"]["ic"]]

    return hoa
end
=#
