

function admm_4_AjAwAgAuAo_intialize(mn_data_nw, fixed_variables, wfz, genz)
    for (key,nw) in mn_data_nw
        push!(fixed_variables,key=>Dict{String,Any}())
        push!(fixed_variables[key],"baseMVA" => nw["baseMVA"])
        #initialize all gen values in dictionary to zero
        push!(fixed_variables[key],"gen" => Dict{String, Any}())
        for (key_g,g) in nw["gen"]
            int_g=parse(Int64,key_g)
            push!(fixed_variables[key]["gen"],key_g => Dict{String, Any}())
            push!(fixed_variables[key]["gen"][key_g],"qg" => 0.00)
            if (issubset(int_g, first.(genz)) && int_g < minimum(first.(wfz)))#conventional gens
                push!(fixed_variables[key]["gen"][key_g],"pg" => 0.00)
            elseif (issubset(int_g, first.(genz)) && int_g > maximum(first.(wfz)))#consumers
                push!(fixed_variables[key]["gen"][key_g],"pg" => 0.00)
            elseif (issubset(parse(Int64,key_g), first.(wfz)))#wind farms
                push!(fixed_variables[key]["gen"][key_g],"wf_pacmax" => 0.00)
                push!(fixed_variables[key]["gen"][key_g],"pg" => 0.00)
            end
        end
        #initialize all bus values in dictionary to zero
        push!(fixed_variables[key],"bus" => Dict{String, Any}())
        for (key_b,b) in nw["bus"]
            push!(fixed_variables[key]["bus"],key_b => Dict{String, Any}())
            push!(fixed_variables[key]["bus"][key_b],"va" => 0.00)
            push!(fixed_variables[key]["bus"][key_b],"lambda" => 1.0)
            push!(fixed_variables[key]["bus"][key_b],"imbalance" => 0.00)
        end
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"branchdc_ne" => Dict{String, Any}())
        for (key_br,br) in nw["branchdc_ne"]
            push!(fixed_variables[key]["branchdc_ne"],key_br => Dict{String, Any}())
            push!(fixed_variables[key]["branchdc_ne"][key_br],"isbuilt" => 0.00)
            push!(fixed_variables[key]["branchdc_ne"][key_br],"pt" => 0.00)
            push!(fixed_variables[key]["branchdc_ne"][key_br],"pf" => 0.00)
            push!(fixed_variables[key]["branchdc_ne"][key_br],"fbusdc" => br["fbusdc"])
            push!(fixed_variables[key]["branchdc_ne"][key_br],"tbusdc" => br["tbusdc"])
            push!(fixed_variables[key]["branchdc_ne"][key_br],"rateA" => br["rateA"])
        end

        push!(fixed_variables[key],"branchdc" => Dict{String, Any}())
        for (key_br,br) in nw["branchdc"]
            push!(fixed_variables[key]["branchdc"],key_br => Dict{String, Any}())
            push!(fixed_variables[key]["branchdc"][key_br],"p_rateA" => 0.00)
            push!(fixed_variables[key]["branchdc"][key_br],"pt" => 0.00)
            push!(fixed_variables[key]["branchdc"][key_br],"pf" => 0.00)
            push!(fixed_variables[key]["branchdc"][key_br],"fbusdc" => br["fbusdc"])
            push!(fixed_variables[key]["branchdc"][key_br],"tbusdc" => br["tbusdc"])
        end
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"storage" => Dict{String, Any}())
        for (key_s,s) in nw["storage"]
            push!(fixed_variables[key]["storage"],key_s => Dict{String, Any}())
            push!(fixed_variables[key]["storage"][key_s],"qs" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"e_absmax" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"e_abs" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"qsc" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"se" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"sc" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"sd" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"ps" => 0.00)
        end
        #if needed add convdc_ne here
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"convdc" => Dict{String, Any}())
        for (key_c,c) in nw["convdc"]
            push!(fixed_variables[key]["convdc"],key_c => Dict{String, Any}())
            push!(fixed_variables[key]["convdc"][key_c],"vaf" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv_pr_fr" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv_dc" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"p_pacmax" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"vac" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv_ac" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv_tf_to" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv_tf_fr" => 0.00)
        end
    end
    return fixed_variables
end


################################################### Power Balance equation ########################
################################## OBZ #######################
# Constraint template: Power balance constraint including candidate storage
# this is the function 2 of the power balance constraint for a OBZ secanario in cluding storage.
# All nodes can have unique shadow prices
#Power balance constraint including candidate storage
function constraint_power_balance_acne_dcne_strg_hm_admm(pm::_PM.AbstractDCPModel, n::Int, is::Set{Int64}, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
    #lambda = pm.setting["fixed_variables"][string(n)]["bus"][string(first(is))]["lam_kcl_r"]
    lambda = pm.setting["fixed_variables"][string(n)]["bus"][string(first(is))]["lambda"]
    p = _PM.var(pm, n, :p)
    pg = _PM.var(pm, n, :pg)
    pconv_grid_ac_ne = _PM.var(pm, n, :pconv_tf_fr_ne)
    pconv_grid_ac = _PM.var(pm, n, :pconv_tf_fr)
    pconv_ac = _PM.var(pm, n, :pconv_ac)
    pconv_ac_ne = _PM.var(pm, n, :pconv_ac_ne)
    p_ne = _PM.var(pm, n, :p_ne)
    ps   = _PM.var(pm, n, :ps)
    #ps_ne   = _PM.var(pm, n, :ps_ne)
    v = 1
    _beta=1
    nodal_balance=0.0
    if !(isempty(bus_arcs))
        nodal_balance+=sum(p[a] for a in bus_arcs);end
    if !(isempty(bus_arcs_ne))
        nodal_balance+=sum(p_ne[a] for a in bus_arcs_ne);end
    if !(isempty(bus_convs_ac))
        nodal_balance+=sum(pconv_grid_ac[c] for c in bus_convs_ac);end
    if !(isempty(bus_convs_ac_ne))
        nodal_balance+=sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne);end
    if !(isempty(bus_gens))
        nodal_balance-=sum(pg[g] for g in bus_gens);end
    if !(isempty(bus_storage))
        nodal_balance+=sum(ps[s] for s in bus_storage);end
    if !(isempty(bus_loads))
        nodal_balance+=sum(pd[d] for d in bus_loads);end
    if !(isempty(bus_shunts))
        nodal_balance+=sum(gs[s] for s in bus_shunts)*v^2;end
    for i in is
        _PM.sol(pm, n, :bus, i)[:imbalance] = nodal_balance
    end

    return lambda*nodal_balance+_beta*(nodal_balance^2)
end
#Power balance constraint including candidate storage
function constraint_power_balance_acne_dcne_strg_admm(pm::_PM.AbstractDCPModel, n::Int, i::Int, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_gens, bus_convs_ac, bus_convs_ac_ne, bus_loads, bus_shunts, bus_storage, bus_storage_ne, pd, qd, gs, bs)
    #lambda = pm.setting["fixed_variables"][string(n)]["bus"][string(i)]["lam_kcl_r"]
    lambda = pm.setting["fixed_variables"][string(n)]["bus"][string(i)]["lambda"]
    p = _PM.var(pm, n, :p)
    pg = _PM.var(pm, n, :pg)
    pconv_grid_ac_ne = _PM.var(pm, n, :pconv_tf_fr_ne)
    pconv_grid_ac = _PM.var(pm, n, :pconv_tf_fr)
    pconv_ac = _PM.var(pm, n, :pconv_ac)
    pconv_ac_ne = _PM.var(pm, n, :pconv_ac_ne)
    p_ne = _PM.var(pm, n, :p_ne)
    ps   = _PM.var(pm, n, :ps)
    #ps_ne   = _PM.var(pm, n, :ps_ne)
    v = 1
    _beta=1
    nodal_balance=0.0
    nodal_balance_l2=0.0
    if !(isempty(bus_arcs))
        nodal_balance_l2+=sum(p[a]^2 for a in bus_arcs)
        nodal_balance+=sum(p[a] for a in bus_arcs)
        end
    if !(isempty(bus_arcs_ne))
        nodal_balance_l2+=sum(p[a]^2 for a in bus_arcs_ne)
        nodal_balance+=sum(p_ne[a] for a in bus_arcs_ne)
        end
    if !(isempty(bus_convs_ac))
        nodal_balance+=sum(pconv_grid_ac[c] for c in bus_convs_ac)
        nodal_balance_l2+=sum(pconv_grid_ac[c]^2 for c in bus_convs_ac)
    end
    if !(isempty(bus_convs_ac_ne))
        nodal_balance+=sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)
        nodal_balance_l2+=sum(pconv_grid_ac_ne[c]^2 for c in bus_convs_ac_ne)
    end

    if !(isempty(bus_gens))
        nodal_balance-=sum(pg[g] for g in bus_gens)
        nodal_balance_l2+=sum(pg[g]^2 for g in bus_gens)
        #=for g in bus_gens
        if (g <= maximum(first.(pm.setting["wfz"])))
            nodal_balance_l2+=pg[g]^2;
        elseif (g > maximum(first.(pm.setting["wfz"])))
            nodal_balance_l2+=-1*pg[g]^2;
        end;end=#
    end
    if !(isempty(bus_storage))
        nodal_balance+=sum(ps[s] for s in bus_storage)
        nodal_balance_l2+=sum(ps[s]^2 for s in bus_storage)
    end
    if !(isempty(bus_loads))
        nodal_balance+=sum(pd[d] for d in bus_loads)
        nodal_balance_l2+=sum(pd[d]^2 for d in bus_loads)
    end
    if !(isempty(bus_shunts))
        nodal_balance+=sum(gs[s] for s in bus_shunts)*v^2
        nodal_balance_l2+=sum(gs[s]^2 for s in bus_shunts)*v^2
    end
    #cstr=JuMP.@constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_ne[a] for a in bus_arcs_ne) + sum(pconv_grid_ac[c] for c in bus_convs_ac) + sum(pconv_grid_ac_ne[c] for c in bus_convs_ac_ne)  == sum(pg[g] for g in bus_gens) - sum(ps[s] for s in bus_storage) -sum(ps_ne[s] for s in bus_storage_ne) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2)
    #if (pm.setting["dual_update"]==false)
    #println("nodal balnce: "*string(nodal_balance))
    #println("nodal balnce_l2: "*string(nodal_balance_l2))

    _PM.sol(pm, n, :bus, i)[:imbalance] = nodal_balance
    return lambda*nodal_balance+_beta*(nodal_balance^2)
    #else
    #    return nodal_balance
    #end
    #return _beta*(nodal_balance^2)
end
 ############################## updating value of lambda #######################

 #n="1";nw=s["fixed_variables"][n]
 function dual_variable_update(fixed_variables)
     beta=1
     residuals=[]
     for (n,nw) in fixed_variables
         for (i_b,b) in nw["bus"]
             lambda0=deepcopy(b["lambda"])
             b["lambda"]=lambda0+beta*b["imbalance"]
             push!(residuals,abs(b["lambda"]-lambda0))
         end
     end
     return fixed_variables, maximum(residuals)
 end


#options:
###Cables
#update_branchdc_ne(n,nw,fixed_variables) #updates fixed variable dictionary's binary cables using convex cable data - sets convex data back to zero
#update_branch_ne(key,nw,fixed_variables)  #update binary cable dictionary
#update_branchdc(key,nw,fixed_variables   #update convex cable dictionary

#update_genz(key,nw,fixed_variables,wfz)  #updates conventional and wind generation
#update_imbalance(key,nw,fixed_variables)  #updates recorded imbalance for each node
#update_storage(key,nw,fixed_variables) #update all storage values in dictionary
#update_dcconv(key,nw,fixed_variables) #update all convdc values in dictionary
 function admm_4_AjAwAgAuAo_update(rez, fixed_variables, agent, wfz)
     for (n,nw) in rez
         if (agent=="Ao")#Only convex cables are used
             fixed_variables=update_branchdc(n,nw,fixed_variables)
             fixed_variables=update_dcconv(n,nw,fixed_variables)
             fixed_variables=update_imbalance(n,nw,fixed_variables)
         elseif (agent=="Aob" || agent=="Aobx")#binary cables/binary + convex
             fixed_variables=update_branch_ne(n,nw,fixed_variables)
             fixed_variables=update_dcconv(n,nw,fixed_variables)
             fixed_variables=update_imbalance(n,nw,fixed_variables)
         elseif (agent=="Aox")#convex cables finding binary cables
             fixed_variables=update_branchdc_ne(n,nw,fixed_variables)
             fixed_variables=update_dcconv(n,nw,fixed_variables)
             fixed_variables=update_imbalance(n,nw,fixed_variables)

         elseif (agent=="Au" || agent=="Ag" || agent=="Aw")#consumpttion, conventional, wind
             fixed_variables=update_genz(n,nw,fixed_variables,wfz)
             fixed_variables=update_imbalance(n,nw,fixed_variables)
         elseif (agent=="Aj")#storage
             fixed_variables=update_storage(n,nw,fixed_variables)
             fixed_variables=update_imbalance(n,nw,fixed_variables)
         #=elseif (agent=="Aobx")#binary cables found from convex cables
             fixed_variables=update_branch_ne(n,nw,fixed_variables)
             fixed_variables=update_dcconv(n,nw,fixed_variables)
                 fixed_variables=update_imbalance(n,nw,fixed_variables)
         elseif (agent=="Ag")#conventional generation
             fixed_variables=update_genz(n,nw,fixed_variables,wfz)
             fixed_variables=update_imbalance(n,nw,fixed_variables)
         elseif (agent=="Aw")#wind power plants
             fixed_variables=update_genz(n,nw,fixed_variables,wfz)
             fixed_variables=update_imbalance(n,nw,fixed_variables)=#
         end
     end
     return fixed_variables
 end

######### Cables
#update convex cable dictionary
function update_branchdc(key,nw,fixed_variables)
    for (key_br,br) in nw["branchdc"]
        fixed_variables[key]["branchdc"][key_br]["p_rateA"] = nw["branchdc"][key_br]["p_rateA"]
        fixed_variables[key]["branchdc"][key_br]["pt"] = nw["branchdc"][key_br]["pt"]
        fixed_variables[key]["branchdc"][key_br]["pf"] = nw["branchdc"][key_br]["pf"]
    end
    return fixed_variables
end

#update binary cable dictionary
function update_branch_ne(key,nw,fixed_variables)
    #update all storage values in dictionary
    for (key_br,br) in nw["branchdc_ne"]
        fixed_variables[key]["branchdc_ne"][key_br]["isbuilt"] = nw["branchdc_ne"][key_br]["isbuilt"]
        if (round(Int64,nw["branchdc_ne"][key_br]["isbuilt"])==1)
            fixed_variables[key]["branchdc_ne"][key_br]["pt"] = nw["branchdc_ne"][key_br]["pt"]
            fixed_variables[key]["branchdc_ne"][key_br]["pf"] = nw["branchdc_ne"][key_br]["pf"]
            println("branch: "*key_br);
        else
            fixed_variables[key]["branchdc_ne"][key_br]["pt"] = 0.0
            fixed_variables[key]["branchdc_ne"][key_br]["pf"] = 0.0;end
    end
    return fixed_variables
end

#updates fixed variable dictionary's binary cables using convex cable data - sets convex data back to zero
function update_branchdc_ne(n,nw,fixed_variables)
    for (i,br) in fixed_variables[n]["branchdc_ne"]
     br["isbuilt"]=0
     br["pt"]=0.0
     br["pf"]=0.0
    end
    for (key_br,br) in nw["branchdc"]
     if (br["p_rateA"]>1)
         matching_br_ne=[]
         for (key_br_ne,br_ne) in fixed_variables[n]["branchdc_ne"]
             if (fixed_variables[n]["branchdc"][key_br]["fbusdc"] ==br_ne["fbusdc"] && fixed_variables[n]["branchdc"][key_br]["tbusdc"]==br_ne["tbusdc"])
                 if ((br_ne["rateA"]-br["p_rateA"])>=0)
                 if ((2*br["p_rateA"]>br_ne["rateA"]))
                     push!(matching_br_ne,(key_br_ne,abs(br_ne["rateA"]-br["p_rateA"])))
                 end;end
             end
         end
         if (length(matching_br_ne)!=0)
             arg=argmin(last.(matching_br_ne))
             fixed_variables[n]["branchdc_ne"][string(first.(matching_br_ne)[arg])]["isbuilt"]=1
             fixed_variables[n]["branchdc_ne"][string(first.(matching_br_ne)[arg])]["pt"]=br["pt"]
             fixed_variables[n]["branchdc_ne"][string(first.(matching_br_ne)[arg])]["pf"]=br["pf"]
         end
     end
     fixed_variables[n]["branchdc"][key_br]["pt"]=0
     fixed_variables[n]["branchdc"][key_br]["pf"]=0
     fixed_variables[n]["branchdc"][key_br]["p_rateA"]=0
    end
    return fixed_variables
end


 #update from realxed set of cables (depricated)
  #=function update_dcbranch_con(key,nw,fixed_variables)
      #update all storage values in dictionary
      eps=0.5
      for (key_br,br) in nw["branchdc_ne"]
          #fixed_variables[key]["branchdc_ne"][key_br]["isbuilt"] = nw["branchdc_ne"][key_br]["isbuilt"]
          if (nw["branchdc_ne"][key_br]["isbuilt"]>eps)
              fixed_variables[key]["branchdc_ne"][key_br]["isbuilt"] = 1
              fixed_variables[key]["branchdc_ne"][key_br]["pt"] = nw["branchdc_ne"][key_br]["pt"]
              fixed_variables[key]["branchdc_ne"][key_br]["pf"] = nw["branchdc_ne"][key_br]["pf"]
          else
              fixed_variables[key]["branchdc_ne"][key_br]["isbuilt"] = 0
              fixed_variables[key]["branchdc_ne"][key_br]["pt"] = 0
              fixed_variables[key]["branchdc_ne"][key_br]["pf"] = 0
          end
      end
      return fixed_variables
  end=#

  #updates conventional and wind generation
   function update_genz(key,nw,fixed_variables,wfz)
       for (key_g,g) in nw["gen"]
           fixed_variables[key]["gen"][key_g]["qg"] = nw["gen"][key_g]["qg"]
           fixed_variables[key]["gen"][key_g]["pg"] = nw["gen"][key_g]["pg"]
           if (issubset(parse(Int64,key_g), first.(wfz)))
               fixed_variables[key]["gen"][key_g]["wf_pacmax"] = nw["gen"][key_g]["wf_pacmax"]
           end
       end
       return fixed_variables
   end

 #updates recorded imbalance for each node
  function update_imbalance(key,nw,fixed_variables)
      #update all storage values in dictionary
      for (key_b,b) in nw["bus"]
          fixed_variables[key]["bus"][key_b]["imbalance"] = nw["bus"][key_b]["imbalance"]
      end
      return fixed_variables
  end

#update all storage values in dictionary
 function update_storage(key,nw,fixed_variables)
     for (key_s,s) in nw["storage"]
         fixed_variables[key]["storage"][key_s]["qs"] = nw["storage"][key_s]["qs"]
         fixed_variables[key]["storage"][key_s]["e_absmax"] = nw["storage"][key_s]["e_absmax"]
         fixed_variables[key]["storage"][key_s]["e_abs"] = nw["storage"][key_s]["e_abs"]
         fixed_variables[key]["storage"][key_s]["qsc"] = nw["storage"][key_s]["qsc"]
         fixed_variables[key]["storage"][key_s]["se"] = nw["storage"][key_s]["se"]
         fixed_variables[key]["storage"][key_s]["sc"] = nw["storage"][key_s]["sc"]
         fixed_variables[key]["storage"][key_s]["sd"] = nw["storage"][key_s]["sd"]
         fixed_variables[key]["storage"][key_s]["ps"] = nw["storage"][key_s]["ps"]
     end
     return fixed_variables
 end

#update all convdc values in dictionary
 function update_dcconv(key,nw,fixed_variables)
     for (key_c,c) in nw["convdc"]
         fixed_variables[key]["convdc"][key_c]["vaf"] = nw["convdc"][key_c]["vafilt"]
         fixed_variables[key]["convdc"][key_c]["pconv_pr_fr"] = nw["convdc"][key_c]["ppr_fr"]
         fixed_variables[key]["convdc"][key_c]["pconv_dc"] = nw["convdc"][key_c]["pdc"]
         fixed_variables[key]["convdc"][key_c]["p_pacmax"] = nw["convdc"][key_c]["p_pacmax"]
         fixed_variables[key]["convdc"][key_c]["vac"] = nw["convdc"][key_c]["vaconv"]
         fixed_variables[key]["convdc"][key_c]["pconv_ac"] = nw["convdc"][key_c]["pconv"]
         fixed_variables[key]["convdc"][key_c]["pconv_tf_to"] = nw["convdc"][key_c]["ptf_to"]
         fixed_variables[key]["convdc"][key_c]["pconv_tf_fr"] = nw["convdc"][key_c]["pgrid"]
     end
     return fixed_variables
 end

 #updates the voltage angle (depricated)
 #=
 function update_bus(key,nw,fixed_variables)
     #update all convdc values in dictionary
     #update all branchdc_ne values in dictionary
     for (key_b,b) in nw["bus"]
         fixed_variables[key]["bus"][key_b]["va"] = nw["bus"][key_b]["va"]
         #push!(fixed_variables[key]["bus"][key_b],"lambda" => 1)
     end
     return fixed_variables
 end
 =#
################################################################################
################# fixing variables in JuMP model ############################
#options:
#fix_storage(pm, n)#set storage vars
#fix_genz(pm, n)#set conv genz
#fix_loadz(pm, n)#set loads
#fix_wind(pm, n)#set wind

#cables
#fix_ne_branches_PGonly(pm, n)
#fix_ne_branches_BINonly(pm, n)
#fix_ne_branches_back2zero(pm, n)
#fix_branchesdc_back2zero(pm, n)
#fix_branchesdc(pm, n)#set branch vars
function fix_variables(pm, n)
    if (pm.setting["agent"]=="Ao" || pm.setting["agent"]=="Aox")
        fix_storage(pm, n)#set storage vars
        fix_genz(pm, n)#set gen vars
        fix_wind(pm, n)#set wind vars
        fix_loadz(pm, n)#set loads
        fix_ne_branches_back2zero(pm, n)
        #undo_relax=JuMP.relax_integrality(pm.model)
    elseif (pm.setting["agent"]=="Aob")
        fix_storage(pm, n)#set storage vars
        fix_genz(pm, n)#set gen vars
        fix_wind(pm, n)#set wind vars
        fix_loadz(pm, n)#set loads
        fix_branchesdc_back2zero(pm, n)#set branch vars
    elseif (pm.setting["agent"]=="Aobx")
        #fix_storage(pm, n)#set storage vars
        #fix_genz(pm, n)#set gen vars
        #fix_wind(pm, n)#set wind vars
        #fix_loadz(pm, n)#set loads
        fix_ne_branches_BINonly(pm, n)#set branch vars
        fix_branchesdc_back2zero(pm, n)#set branch vars
    elseif (pm.setting["agent"]=="fixed_cables_cons")
        #fix_storage(pm, n)#set storage vars
        #fix_genz(pm, n)#set gen vars
        #fix_wind(pm, n)#set wind vars
        #fix_loadz(pm, n)#set loads
        fix_ne_branches_BINonly(pm, n)#set branch vars
        fix_branchesdc_back2zero(pm, n)#set branch vars
        fix_convdc(pm,n)#set convdc vars
    elseif (pm.setting["agent"]=="Aw")
        fix_storage(pm, n)#set storage vars
        fix_genz(pm, n)#set gen vars
        fix_loadz(pm, n)#set loads
        fix_ne_branches_BINonly(pm, n)#set branch vars
        fix_ne_branches_PGonly(pm, n)#set branch vars
        fix_branchesdc(pm, n)#set branch vars
        fix_convdc(pm,n)#set convdc vars
    elseif (pm.setting["agent"]=="Aj")
        fix_wind(pm, n)#set wind vars
        fix_genz(pm, n)#set gen vars
        fix_loadz(pm, n)#set loads
        fix_ne_branches_BINonly(pm, n)#set branch vars
        fix_ne_branches_PGonly(pm, n)#set branch vars
        fix_branchesdc(pm, n)#set branch vars
        fix_convdc(pm,n)#set convdc vars
    elseif (pm.setting["agent"]=="Ag")
        fix_storage(pm, n)#set storage vars
        fix_wind(pm, n)#set wind vars
        fix_loadz(pm, n)#set loads
        fix_ne_branches_BINonly(pm, n)#set branch vars
        fix_ne_branches_PGonly(pm, n)#set branch vars
        fix_branchesdc(pm, n)#set branch vars
        fix_convdc(pm,n)#set convdc vars
    elseif (pm.setting["agent"]=="Au")
        fix_storage(pm, n)#set storage vars
        fix_wind(pm, n)#set wind vars
        fix_genz(pm, n)#set gen vars
        fix_ne_branches_BINonly(pm, n)#set branch vars
        fix_ne_branches_PGonly(pm, n)#set branch vars
        fix_branchesdc(pm, n)#set branch vars
        fix_convdc(pm,n)#set convdc vars
    #=elseif (pm.setting["agent"]=="Aox")
        fix_storage(pm, n)#set storage vars
        fix_genz(pm, n)#set gen vars
        fix_wind(pm, n)#set wind vars
        fix_loadz(pm, n)#set loads
        fix_ne_branches_back2zero(pm, n)
        #fix_binaries(pm, n)#set binary branch vars=#
    end
end

function fix_storage(pm, n)#set storage vars
    for (i,s) in pm.setting["fixed_variables"][string(n)]["storage"]
        for (kv,v) in s
            storage = _PM.var(pm, n, Symbol(kv), parse(Int64,i))
            JuMP.fix(storage,v,force=true)
        end
    end
end

function fix_genz(pm, n)#set conv genz

    for (i,g) in pm.setting["fixed_variables"][string(n)]["gen"]
        i_int=parse(Int64,i)
        if (i_int<minimum(first.(pm.setting["wfz"])))
            #println("genz: "*i)
            for (kv,v) in g
                gen = _PM.var(pm, n, Symbol(kv), i_int)
                JuMP.fix(gen,v,force=true)
            end
        end
    end
end

function fix_loadz(pm, n)#set loads

    for (i,g) in pm.setting["fixed_variables"][string(n)]["gen"]
        i_int=parse(Int64,i)
        if (i_int>maximum(first.(pm.setting["wfz"])))
            #println("loadz: "*i)
            for (kv,v) in g
                #println(kv*" value "*string(v))
                gen = _PM.var(pm, n, Symbol(kv), i_int)
                JuMP.fix(gen,v,force=true)
            end
        end
    end
end

function fix_wind(pm, n)#set wind

    for (i,g) in pm.setting["fixed_variables"][string(n)]["gen"]
        i_int=parse(Int64,i)
        if (issubset(i_int, first.(pm.setting["wfz"])))
            for (kv,v) in g
                gen = _PM.var(pm, n, Symbol(kv), parse(Int64,i))
                JuMP.fix(gen,v,force=true)
            end
        end
    end
end

function fix_ne_branches_PGonly(pm, n)
    if haskey(_PM.var(pm, n), :p_dcgrid_ne)
        p_dcgrid_ne = _PM.var(pm, n, :p_dcgrid_ne)
        if !isempty(p_dcgrid_ne)
            brchs=[]
            #for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_ne):arcs_dcgrid_from_ne
            for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_from_ne)
                #if !(issubset([l],brchs))
                if (pm.setting["fixed_variables"][string(n)]["branchdc_ne"][string(l)]["isbuilt"]==1)
                    JuMP.fix(p_dcgrid_ne[(l,i,j)],pm.setting["fixed_variables"][string(n)]["branchdc_ne"][string(l)]["pf"],force=true)
                end
            end
            for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_to_ne)
                if (pm.setting["fixed_variables"][string(n)]["branchdc_ne"][string(l)]["isbuilt"]==1)
                    JuMP.fix(p_dcgrid_ne[(l,i,j)],pm.setting["fixed_variables"][string(n)]["branchdc_ne"][string(l)]["pt"],force=true)
                end
            end
        end
    end
end

function fix_ne_branches_BINonly(pm, n)
    #set branches
    if haskey(_PM.ref(pm, n), :branchdc_ne)
        branchdc_ne = _PM.ref(pm, n, :branchdc_ne)

        if !isempty(branchdc_ne)
            for (i,br) in pm.setting["fixed_variables"][string(n)]["branchdc_ne"]
                for (kv,v) in br
                    if (kv=="isbuilt")
                        brnch=_PM.var(pm, n, :branchdc_ne, parse(Int64,i))
                        JuMP.fix(brnch,round(Int64,v),force=true)
                    end
                end
            end
        end
    end
end

function fix_ne_branches_back2zero(pm, n)
    if haskey(_PM.ref(pm, n), :branchdc_ne)
        branchdc_ne = _PM.ref(pm, n, :branchdc_ne)

        if !isempty(branchdc_ne)
            for (i,br) in pm.setting["fixed_variables"][string(n)]["branchdc_ne"]
                for (kv,v) in br
                    if (kv=="isbuilt")
                        brnch=_PM.var(pm, n, :branchdc_ne, parse(Int64,i))
                        JuMP.fix(brnch,0,force=true)
                    end
                end
            end
        end
    end
    if haskey(_PM.var(pm, n), :p_dcgrid_ne)
        p_dcgrid_ne = _PM.var(pm, n, :p_dcgrid_ne)
        if !isempty(p_dcgrid_ne)
            brchs=[]
            for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_from_ne)
                if (pm.setting["fixed_variables"][string(n)]["branchdc_ne"][string(l)]["isbuilt"]==1)
                    JuMP.fix(p_dcgrid_ne[(l,i,j)],0.0,force=true)
                end
            end
            for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_to_ne)
                if (pm.setting["fixed_variables"][string(n)]["branchdc_ne"][string(l)]["isbuilt"]==1)
                    JuMP.fix(p_dcgrid_ne[(l,i,j)],0.0,force=true)
                end
            end
        end
    end
end

function fix_branchesdc_back2zero(pm, n)
    #set branches
    if haskey(_PM.ref(pm, n), :branchdc)
        branchdc = _PM.ref(pm, n, :branchdc)

        if !isempty(branchdc)
            for (i,br) in pm.setting["fixed_variables"][string(n)]["branchdc"]
                brnch=_PM.var(pm, n, :p_rateA, parse(Int64,i))
                JuMP.fix(brnch,0.0,force=true)
            end
        end
    end
    if haskey(_PM.var(pm, n), :p_dcgrid)
        p_dcgrid = _PM.var(pm, n, :p_dcgrid)
        if !isempty(p_dcgrid)
            brchs=[]
            for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_from)
                    JuMP.fix(p_dcgrid[(l,i,j)],0.0,force=true)
            end
            for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_to)
                    JuMP.fix(p_dcgrid[(l,i,j)],0.0,force=true)
            end
        end
    end
end

function fix_branchesdc(pm, n)#set branch vars
    #set convdc
    if haskey(_PM.var(pm, n), :p_dcgrid)
        p_dcgrid = _PM.var(pm, n, :p_dcgrid)
        if !isempty(p_dcgrid)
            brchs=[]
            for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_from)
                JuMP.fix(p_dcgrid[(l,i,j)],pm.setting["fixed_variables"][string(n)]["branchdc"][string(l)]["pf"],force=true)
            end
            for (l,i,j) in _PM.ref(pm, n, :arcs_dcgrid_to)
                JuMP.fix(p_dcgrid[(l,i,j)],pm.setting["fixed_variables"][string(n)]["branchdc"][string(l)]["pt"],force=true)
            end
        end
    end
end

#=
["vaf","pconv_pr_fr","pconv_dc","p_pacmax","vac","pconv_ac","pconv_tf_to","pconv_tf_fr"]
=#
#"vaf", "p_pacmax",    "pconv_tf_to","vac","pconv_tf_fr","pconv_pr_fr","pconv_dc","pconv_ac"
function fix_convdc(pm,n)
    #set convdc
    for (i,cv) in pm.setting["fixed_variables"][string(n)]["convdc"]
        for (kv,v) in cv
            #if !(issubset([kv],[]))
                conv = _PM.var(pm, n, Symbol(kv), parse(Int64,i))
                JuMP.fix(conv,v,force=true)
            #end
        end
    end
end


#= Depricated
function fix_bus(pm, n)
    #set storagebus vars
    for (i,bu) in pm.setting["fixed_variables"][string(n)]["bus"]
        for (kv,v) in bu
            if (kv!="lambda" && kv!="imbalance")
                bus = _PM.var(pm, n, Symbol(kv), parse(Int64,i))
                JuMP.fix(bus,v,force=true)
            end
        end
    end
end=#
################################### Depricated #################################
#=
#n="1";nw=s["fixed_variables"][n]
function dual_variable_update(fixed_variables)
    beta=1
    residuals=[]
    for (n,nw) in fixed_variables
        for (i_b,b) in nw["bus"]
            lambda0=deepcopy(b["lambda"])
            b["lambda"]=lambda0+beta*b["imbalance"]
            #println(string(i_b)*" "*string(b["lambda"]))
            push!(residuals,abs(b["lambda"]-lambda0))
        end
    end
    #println(residuals)
    return fixed_variables, maximum(residuals)
end


function admm_4_AjAwAgAuAo_update(rez, fixed_variables, agent, wfz)
    for (n,nw) in rez
        if (agent=="Ao")
            fixed_variables=update_genz(n,nw,fixed_variables,wfz)
            #fixed_variables=update_bus(n,nw,fixed_variables)
            fixed_variables=update_dcbranch_con(n,nw,fixed_variables)
            fixed_variables=update_dcconv(n,nw,fixed_variables)
            fixed_variables=update_imbalance(n,nw,fixed_variables)
        elseif (agent=="Aob" || agent=="Aox")
            fixed_variables=update_genz(n,nw,fixed_variables,wfz)
            #fixed_variables=update_bus(n,nw,fixed_variables)
            fixed_variables=update_branch_ne(n,nw,fixed_variables)
            fixed_variables=update_dcconv(n,nw,fixed_variables)
            fixed_variables=update_imbalance(n,nw,fixed_variables)
        elseif (agent=="Au")
            fixed_variables=update_genz(n,nw,fixed_variables,wfz)
            fixed_variables=update_imbalance(n,nw,fixed_variables)
        elseif (agent=="Ag")
            fixed_variables=update_genz(n,nw,fixed_variables,wfz)
            fixed_variables=update_imbalance(n,nw,fixed_variables)
        elseif (agent=="Aw")
            fixed_variables=update_genz(n,nw,fixed_variables,wfz)
            fixed_variables=update_imbalance(n,nw,fixed_variables)
        elseif (agent=="Aj")
            fixed_variables=update_storage(n,nw,fixed_variables)
            fixed_variables=update_imbalance(n,nw,fixed_variables)
        end
    end
    return fixed_variables
end


function update_genz(key,nw,fixed_variables,wfz)
    for (key_g,g) in nw["gen"]
        fixed_variables[key]["gen"][key_g]["qg"] = nw["gen"][key_g]["qg"]
        fixed_variables[key]["gen"][key_g]["pg"] = nw["gen"][key_g]["pg"]
        if (issubset(parse(Int64,key_g), first.(wfz)))
            fixed_variables[key]["gen"][key_g]["wf_pacmax"] = 40#nw["gen"][key_g]["wf_pacmax"]
        end
    end
    return fixed_variables
end

function update_bus(key,nw,fixed_variables)
    #update all storage values in dictionary
    for (key_b,b) in nw["bus"]
        fixed_variables[key]["bus"][key_b]["va"] = nw["bus"][key_b]["va"]
    end
    return fixed_variables
end

function update_imbalance(key,nw,fixed_variables)
    #update all storage values in dictionary
    for (key_b,b) in nw["bus"]
        fixed_variables[key]["bus"][key_b]["imbalance"] = nw["bus"][key_b]["imbalance"]
    end
    return fixed_variables
end

function update_dcbranch_con(key,nw,fixed_variables)
    #update all storage values in dictionary
    eps=0.5
    for (key_br,br) in nw["branchdc_ne"]
        #fixed_variables[key]["branchdc_ne"][key_br]["isbuilt"] = nw["branchdc_ne"][key_br]["isbuilt"]
        if (nw["branchdc_ne"][key_br]["isbuilt"]>eps)
            println("cable "*string(key_br)*" built!!!!")
            fixed_variables[key]["branchdc_ne"][key_br]["isbuilt"] = 1
            fixed_variables[key]["branchdc_ne"][key_br]["pt"] = nw["branchdc_ne"][key_br]["pt"]
            fixed_variables[key]["branchdc_ne"][key_br]["pf"] = nw["branchdc_ne"][key_br]["pf"]
        else
            fixed_variables[key]["branchdc_ne"][key_br]["isbuilt"] = 0
            fixed_variables[key]["branchdc_ne"][key_br]["pt"] = 0
            fixed_variables[key]["branchdc_ne"][key_br]["pf"] = 0
        end
        #if (round(Int64,nw["branchdc_ne"][key_br]["isbuilt"])>0)
        #    println("a cable exists@!!!!!!!!!")

        #else
        #    fixed_variables[key]["branchdc_ne"][key_br]["pt"] = nw["branchdc_ne"][key_br]["pt"]
        #    fixed_variables[key]["branchdc_ne"][key_br]["pf"] = nw["branchdc_ne"][key_br]["pf"]
        #end
    end
    return fixed_variables
end

function update_branch_ne(key,nw,fixed_variables)
    #update all storage values in dictionary
    for (key_br,br) in nw["branchdc_ne"]
        fixed_variables[key]["branchdc_ne"][key_br]["isbuilt"] = nw["branchdc_ne"][key_br]["isbuilt"]
        if (round(Int64,nw["branchdc_ne"][key_br]["isbuilt"])==1)
            fixed_variables[key]["branchdc_ne"][key_br]["pt"] = nw["branchdc_ne"][key_br]["pt"]
            fixed_variables[key]["branchdc_ne"][key_br]["pf"] = nw["branchdc_ne"][key_br]["pf"]
            println("branch: "*key_br);
        else
            fixed_variables[key]["branchdc_ne"][key_br]["pt"] = 0.0
            fixed_variables[key]["branchdc_ne"][key_br]["pf"] = 0.0;end
    end
    return fixed_variables
end

function update_storage(key,nw,fixed_variables)
    #update all storage values in dictionary
    for (key_s,s) in nw["storage"]
        fixed_variables[key]["storage"][key_s]["qs"] = nw["storage"][key_s]["qs"]
        fixed_variables[key]["storage"][key_s]["e_absmax"] = nw["storage"][key_s]["e_absmax"]
        fixed_variables[key]["storage"][key_s]["e_abs"] = nw["storage"][key_s]["e_abs"]
        fixed_variables[key]["storage"][key_s]["qsc"] = nw["storage"][key_s]["qsc"]
        fixed_variables[key]["storage"][key_s]["se"] = nw["storage"][key_s]["se"]
        fixed_variables[key]["storage"][key_s]["sc"] = nw["storage"][key_s]["sc"]
        fixed_variables[key]["storage"][key_s]["sd"] = nw["storage"][key_s]["sd"]
        fixed_variables[key]["storage"][key_s]["ps"] = nw["storage"][key_s]["ps"]
    end
    return fixed_variables
end

function update_dcconv(key,nw,fixed_variables)
    #update all convdc values in dictionary
    for (key_c,c) in nw["convdc"]
        fixed_variables[key]["convdc"][key_c]["vaf"] = nw["convdc"][key_c]["vafilt"]
        fixed_variables[key]["convdc"][key_c]["pconv_pr_fr"] = nw["convdc"][key_c]["ppr_fr"]
        fixed_variables[key]["convdc"][key_c]["pconv_dc"] = nw["convdc"][key_c]["pdc"]
        fixed_variables[key]["convdc"][key_c]["p_pacmax"] = nw["convdc"][key_c]["p_pacmax"]
        fixed_variables[key]["convdc"][key_c]["vac"] = nw["convdc"][key_c]["vaconv"]
        fixed_variables[key]["convdc"][key_c]["pconv_ac"] = nw["convdc"][key_c]["pconv"]
        fixed_variables[key]["convdc"][key_c]["pconv_tf_to"] = nw["convdc"][key_c]["ptf_to"]
        fixed_variables[key]["convdc"][key_c]["pconv_tf_fr"] = nw["convdc"][key_c]["pgrid"]
    end
    return fixed_variables
end

function update_bus(key,nw,fixed_variables)
    #update all convdc values in dictionary
    #update all branchdc_ne values in dictionary
    for (key_b,b) in nw["bus"]
        fixed_variables[key]["bus"][key_b]["va"] = nw["bus"][key_b]["va"]
        #push!(fixed_variables[key]["bus"][key_b],"lambda" => 1)
    end
    return fixed_variables
end

=#

#=function admm_4_AjAwAgAuAo_intialize(mn_data_nw, fixed_variables)
    for (key,nw) in mn_data_nw
        push!(fixed_variables,key=>Dict{String,Any}())
        push!(fixed_variables[key],"baseMVA" => nw["baseMVA"])
        #initialize all gen values in dictionary to zero
        push!(fixed_variables[key],"gen" => Dict{String, Any}())
        for (key_g,g) in nw["gen"]
            push!(fixed_variables[key]["gen"],key_g => Dict{String, Any}())
            push!(fixed_variables[key]["gen"][key_g],"qg" => 0.00)
            push!(fixed_variables[key]["gen"][key_g],"pg" => 0.00)
            push!(fixed_variables[key]["gen"][key_g],"wf_pacmax" => 0.00)
        end
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"branchdc_ne" => Dict{String, Any}())
        for (key_br,br) in nw["branchdc_ne"]
            push!(fixed_variables[key]["branchdc_ne"],key_br => Dict{String, Any}())
            push!(fixed_variables[key]["branchdc_ne"][key_br],"isbuilt" => 0.00)
            push!(fixed_variables[key]["branchdc_ne"][key_br],"pt" => 0.00)
            push!(fixed_variables[key]["branchdc_ne"][key_br],"pf" => 0.00)
        end
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"storage" => Dict{String, Any}())
        for (key_s,s) in nw["storage"]
            push!(fixed_variables[key]["storage"],key_s => Dict{String, Any}())
            #push!(fixed_variables[key]["storage"][key_s],"qs" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"e_absmax" => 0.00)
            #push!(fixed_variables[key]["storage"][key_s],"e_abs" => 0.00)
            #push!(fixed_variables[key]["storage"][key_s],"qsc" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"se" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"sc" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"sd" => 0.00)
            push!(fixed_variables[key]["storage"][key_s],"ps" => 0.00)
        end
        #if needed add convdc_ne here
        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"convdc" => Dict{String, Any}())
        for (key_c,c) in nw["convdc"]
            push!(fixed_variables[key]["convdc"],key_c => Dict{String, Any}())
            push!(fixed_variables[key]["convdc"][key_c],"vafilt" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"ppr_fr" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pdc" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"p_pacmax" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"vaconv" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"pconv" => 0.00)
            push!(fixed_variables[key]["convdc"][key_c],"ptf_to" => 0.00)
        end

        #initialize all branchdc_ne values in dictionary to zero
        push!(fixed_variables[key],"bus" => Dict{String, Any}())
        for (key_b,b) in nw["bus"]
            push!(fixed_variables[key]["bus"],key_b => Dict{String, Any}())
            push!(fixed_variables[key]["bus"][key_b],"va" => 0.00)
            push!(fixed_variables[key]["bus"][key_b],"lam_kcl_i" => NaN)
            push!(fixed_variables[key]["bus"][key_b],"vm" => 0.00)
            push!(fixed_variables[key]["bus"][key_b],"lam_kcl_r" => 0.00)
        end
    end
    return fixed_variables
end=#
