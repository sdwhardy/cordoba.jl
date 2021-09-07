
function gen(data,mva)
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
end

function branchdc_ne(b_ne,dc_can,mbase)
    b_ne["rateA"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["rateB"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["rateC"]=(dc_can.num*dc_can.elec.mva)/mbase
    b_ne["r"]=dc_can.elec.ohm

    #NOTE change
    b_ne["cost"]=dc_can.costs.cpx_p+dc_can.costs.cpx_i
    #b_ne["cost"]=0
    #######
    return b_ne
end

function ne_branch(ne_b,ac_can,mbase)
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
end

function branch_xfo(ne_b,ne_xfo,mbase)
    #ne_b["construction_cost"]=ne_b["construction_cost"]+(ne_xfo.costs.cpx_p+ne_xfo.costs.cpx_i)/mbase

    #NOTE change
    ne_b["construction_cost"]=ne_b["construction_cost"]+ne_xfo.costs.cpx_p+ne_xfo.costs.cpx_i
    #ne_b["construction_cost"]=ne_b["construction_cost"]+0
    ##########
    return ne_b
end

function branch_plat(ne_b,ne_plat,mbase)
    #ne_b["construction_cost"]=ne_b["construction_cost"]+(ne_plat.costs.cpx)/mbase

    #NOTE change
    ne_b["construction_cost"]=ne_b["construction_cost"]+ne_plat.costs.cpx
    #ne_b["construction_cost"]=ne_b["construction_cost"]+0
    ##########
    return ne_b
end

function pm_dict_candidateICs(data,hoa,candidates)
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
end

function additional_candidates(data,candidates)
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
end

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
