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
        nodal_balance_l2+=sum(p_ne[a]^2 for a in bus_arcs_ne)
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

