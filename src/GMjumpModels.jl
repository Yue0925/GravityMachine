# ==============================================================================
# Modele JuMP pour calculer la relaxation linéaire du 2SPA sur un objectif donne, activant eventuellement une ϵ-contrainte

function computeLinearRelax2SPA(  nbvar::Int,
                                  nbctr::Int,
                                  A::Array{Int,2},
                                  c1::Array{Int,1},
                                  c2::Array{Int,1},
                                  epsilon,
                                  obj::Int
)
  model = Model(CPLEX.Optimizer) ; set_silent(model)
  @variable(model, 0.0 <= x[1:nbvar] <= 1.0 )
  @constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
  if obj == 1
    @objective(model, Min, sum((c1[i])*x[i] for i in 1:nbvar))
    @constraint(model, sum((c2[i])*x[i] for i in 1:nbvar) <= epsilon)
  else
    @objective(model, Min, sum((c2[i])*x[i] for i in 1:nbvar))
    @constraint(model, sum((c1[i])*x[i] for i in 1:nbvar) <= epsilon)
  end
  optimize!(model)
  return objective_value(model), value.(x)
end


function computeLinearRelax( model::JuMP.Model,  x::Array{JuMP.VariableRef}, 
  c::Matrix{Float64},
  epsilon,
  obj::Int
)

  cons = Vector{ConstraintRef}()
  if obj == 1
    set_objective(model, MOI.MIN_SENSE, c[1, 1] + c[1, 2:end]'*x  )
    if epsilon == -1 
       nothing 
    else
      con = @constraint(model, c[2, 1] + c[2, 2:end]'*x  <= epsilon) ; push!(cons, con)
    end

  else
    set_objective(model, MOI.MIN_SENSE, c[2, 1] + c[2, 2:end]'*x )
    if epsilon == -1 
      nothing 
    else
      con = @constraint(model, c[1, 1] + c[1, 2:end]'*x  <= epsilon) ; push!(cons, con)
    end
  end

  optimize!(model)
  v_obj, v_x = objective_value(model), value.(x)

  for con in cons
    JuMP.delete( model, con) ; JuMP.unregister( model, :con)
  end
  return v_obj, v_x
end