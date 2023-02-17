""

function optimize_model!(aim::_IM.AbstractInfrastructureModel; relax_integrality=false, optimizer=nothing, solution_processors=[])
    start_time = time()

    if relax_integrality
        JuMP.relax_integrality(aim.model)
    end

    if JuMP.mode(aim.model) != JuMP.DIRECT && optimizer !== nothing
        if JuMP.backend(aim.model).optimizer === nothing
            JuMP.set_optimizer(aim.model, optimizer)
        else
            Memento.warn(_PM._LOGGER, "Model already contains optimizer, cannot use optimizer specified in `optimize_model!`")
        end
    end

    if JuMP.mode(aim.model) != JuMP.DIRECT && JuMP.backend(aim.model).optimizer === nothing
        Memento.error(_PM._LOGGER, "No optimizer specified in `optimize_model!` or the given JuMP model.")
    end

    _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(aim.model)
    #_, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(aim.model)
  
    try
        solve_time = JuMP.solve_time(aim.model)
    catch
        Memento.warn(_PM._LOGGER, "The given optimizer does not provide the SolveTime() attribute, falling back on @timed.  This is not a rigorous timing value.");
    end
    
    Memento.debug(_PM._LOGGER, "JuMP model optimize time: $(time() - start_time)")

    start_time = time()
    result = _IM.build_result(aim, solve_time; solution_processors=solution_processors)
    Memento.debug(_PM._LOGGER, "solution build time: $(time() - start_time)")

    aim.solution = result["solution"]

    return result
end

################# Important trouble shooting code do not erase !!!!!!!!!!!!!!!!!!!!!!
#=
if JuMP.termination_status(aim.model) == _MOI.INFEASIBLE_OR_UNBOUNDED
        JuMP.@assert JuMP.termination_status(aim.model) == _MOI.INFEASIBLE
        JuMP.compute_conflict!(aim.model)
        contypes=JuMP.list_of_constraint_types(aim.model)
        for contype in contypes
        as=_MOI.get.(aim.model, _MOI.ConstraintConflictStatus(),JuMP.all_constraints(aim.model, JuMP.VariableRef, _MOI.EqualTo{Float64}))
        as=_MOI.get.(aim.model, _MOI.ConstraintConflictStatus(),JuMP.all_constraints(aim.model, JuMP.VariableRef, _MOI.GreaterThan{Float64}))
        as=_MOI.get.(aim.model, _MOI.ConstraintConflictStatus(),JuMP.all_constraints(aim.model, JuMP.VariableRef, _MOI.LessThan{Float64}))
        for (i,a) in enumerate(as)
        if (a!=_MOI.NOT_IN_CONFLICT)   
            println(i, " ", a);end
        end 
        JuMP.all_constraints(aim.model, JuMP.VariableRef, _MOI.LessThan{Float64})[84343]
        JuMP.all_constraints(aim.model, JuMP.VariableRef, _MOI.GreaterThan{Float64})[13959]
    end
=#
########################################################################################
