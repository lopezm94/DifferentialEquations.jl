"""
FEMSolution

Holds the data for the solution to a finite element problem.

### Fields

* `fem_mesh::FEMmesh`: The finite element mesh the problem was solved on.
* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueknown::Bool`: Boolean flag for if the true solution is given.
* `u_analytic::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `appxTrue::Bool`: Boolean flag for if u_analytic was an approximation.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
is specified in the solver.
* `t::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
is specified in the solver.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.

"""
type FEMSolution <: DESolution
  fem_mesh::FEMmesh
  u#::Array{Number}
  trueknown::Bool
  u_analytic::AbstractArrayOrVoid
  errors#::Dict{String,Float64}
  appxTrue::Bool
  timeseries#::GrowableArray
  t::AbstractArrayOrVoid
  prob::DEProblem
  save_timeseries::Bool
  function FEMSolution(fem_mesh::FEMmesh,u,u_analytic,sol,Du,timeSeries,t,prob;save_timeseries=true)
    errors = Dict(:L2=>getL2error(fem_mesh,sol,u),:H1=>getH1error(fem_mesh,Du,u),
                  :l∞=> maximum(abs(u-u_analytic)), :l2=> norm(u-u_analytic,2))
    return(new(fem_mesh,u,true,u_analytic,errors,false,timeSeries,t,prob,true))
  end
  FEMSolution(fem_mesh,u,u_analytic,sol,Du,prob) = FEMSolution(fem_mesh::FEMmesh,u,u_analytic,sol,Du,[],[],prob,save_timeseries=false)
  function FEMSolution(fem_mesh::FEMmesh,u::AbstractArray,prob)
    return(FEMSolution(fem_mesh,u,[],[],prob,save_timeseries=false))
  end
  function FEMSolution(fem_mesh::FEMmesh,u::AbstractArray,timeseries,t,prob;save_timeseries=true)
    return(new(fem_mesh,u,false,nothing,Dict{String,Float64},false,timeseries,t,prob,save_timeseries))
  end
end

"""
SDESolution

Holds the data for the solution to a SDE problem.

### Fields

* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueknown::Bool`: Boolean flag for if the true solution is given.
* `u_analytic::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
is specified in the solver.
* `t::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
is specified in the solver.
* `Ws`: All of the W's in the solution. Only saved if `save_timeseries=true` is specified
in the solver.
* `timeseries_analytic`: If `save_timeseries=true`, saves the solution at each save point.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.
* `appxTrue::Bool`: Boolean flag for if u_analytic was an approximation.

"""
type SDESolution <: DESolution
  u#::AbstractArrayOrNumber
  trueknown::Bool
  u_analytic#::AbstractArrayOrNumber
  errors#::Dict{}
  timeseries::AbstractArrayOrVoid
  t::AbstractArrayOrVoid
  Δt::AbstractArrayOrVoid
  Ws::AbstractArrayOrVoid
  timeseries_analytic::AbstractArrayOrVoid
  appxTrue::Bool
  save_timeseries::Bool
  maxStackSize::Int
  W
  function SDESolution(u::Union{AbstractArray,Number};timeseries=nothing,timeseries_analytic=nothing,t=nothing,Δt=nothing,Ws=nothing,maxStackSize=nothing,W=nothing)
    save_timeseries = timeseries == nothing
    trueknown = false
    return(new(u,trueknown,nothing,Dict(),timeseries,t,Δt,Ws,timeseries_analytic,false,save_timeseries,maxStackSize,W))
  end
  function SDESolution(u,u_analytic;timeseries=nothing,timeseries_analytic=nothing,t=nothing,Δt=nothing,Ws=nothing,maxStackSize=nothing,W=nothing)
    save_timeseries = timeseries != nothing
    trueknown = true
    errors = Dict(:final=>mean(abs(u-u_analytic)))
    if save_timeseries
      errors = Dict(:final=>mean(abs(u-u_analytic)),:l∞=>maximum(abs(timeseries-timeseries_analytic)),:l2=>sqrt(mean((timeseries-timeseries_analytic).^2)))
    end
    return(new(u,trueknown,u_analytic,errors,timeseries,t,Δt,Ws,timeseries_analytic,false,save_timeseries,maxStackSize,W))
  end
  #Required to convert pmap results
  SDESolution(a::Any) = new(a.u,a.trueknown,a.u_analytic,a.errors,a.timeseries,a.t,a.Δt,a.Ws,a.timeseries_analytic,a.appxTrue,a.save_timeseries,a.maxStackSize,a.W)
end

"""
ODESolution

Holds the data for the solution to an ODE problem.

### Fields

* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueknown::Bool`: Boolean flag for if the true solution is given.
* `u_analytic::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
is specified in the solver.
* `t::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
is specified in the solver.
* `timeseries_analytic`: If `save_timeseries=true`, saves the solution at each timestep.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.
* `appxTrue::Bool`: Boolean flag for if u_analytic was an approximation.

"""
type ODESolution <: DESolution
  u#::AbstractArrayOrNumber
  trueknown::Bool
  u_analytic#::AbstractArrayOrNumber
  errors#::Dict{}
  timeseries::AbstractArrayOrVoid
  t::AbstractArrayOrVoid
  timeseries_analytic::AbstractArrayOrVoid
  appxTrue::Bool
  save_timeseries::Bool
  function ODESolution(u;timeseries=nothing,timeseries_analytic=nothing,t=nothing)
    save_timeseries = timeseries == nothing
    trueknown = false
    return(new(u,trueknown,nothing,Dict(),timeseries,t,timeseries_analytic,false,save_timeseries))
  end
  function ODESolution(u,u_analytic;timeseries=nothing,timeseries_analytic=nothing,t=nothing)
    save_timeseries = timeseries != nothing
    trueknown = true
    errors = Dict(:final=>mean(abs(u-u_analytic)))
    if save_timeseries
      errors = Dict(:final=>mean(abs(u-u_analytic)),:l∞=>maximum(abs(timeseries-timeseries_analytic)),:l2=>sqrt(mean((timeseries-timeseries_analytic).^2)))
    end
    return(new(u,trueknown,u_analytic,errors,timeseries,t,timeseries_analytic,false,save_timeseries))
  end
end

"""
StokesSolution

Holds the data for the solution to a Stokes problem.

### Fields

* u
* v
* p
* u_analytic
* vTrue
* pTrue
* mesh
* trueknown
* errors
* converrors

"""
type StokesSolution <: DESolution
  u
  v
  p
  u_analytic
  vTrue
  pTrue
  mesh::FDMMesh
  trueknown::Bool
  errors
  converrors
  StokesSolution(u,v,p,u_analytic,vTrue,pTrue,mesh,trueknown;errors=nothing,converrors=nothing) = new(u,v,p,u_analytic,vTrue,pTrue,mesh,trueknown,errors,converrors)
end

"""
appxTrue!(res,res2)

Adds the solution from res2 to the FEMSolution object res.
Useful to add a quasi-true solution when none is known by
computing once at a very small time/space step and taking
that solution as the "true" solution
"""
function appxTrue!(res::FEMSolution,res2::FEMSolution)
  res.u_analytic = res2.u
  res.errors = Dict(:l∞=>maximum(abs(res.u-res.u_analytic)),:l2=>norm(res.u-res.u_analytic,2))
  res.appxTrue = true
end

"""
S = FEMSolutionTS(timeseries::GrowableArray,numvars::Int)
S[i][j] => Variable i at time j.
"""
function FEMSolutionTS(timeseries::GrowableArray,numvars::Int)
  G = GrowableArray(timeseries[1][:,1])
  for j = 2:length(timeseries)
    push!(G,timeseries[j][:,1])
  end
  ts = GrowableArray(G)
  if numvars > 1
    for i=2:numvars
      G = GrowableArray(timeseries[1][:,i])
      for j = 2:length(timeseries)
        push!(G,timeseries[j][:,i])
      end
      push!(ts,G)
    end
  end
  return(ts)
end

Base.length(sol::DESolution) = length(sol.t)
Base.size(sol::DESolution) = (length(sol.t),size(sol.u))
Base.endof(sol::DESolution) = length(sol)
Base.getindex(sol::DESolution,i::Int) = sol.timeseries[i,..]
Base.getindex(sol::DESolution,i::Int,I::Int...) = sol.timeseries[i,I...]
Base.getindex(sol::DESolution,::Colon) = sol.timeseries

function print(io::IO, sol::DESolution)
  if sol.trueknown
    str="Analytical solution is known."
  else
    str="No analytical solution is known."
  end
  println(io,"$(typeof(sol)) with $(length(sol)) timesteps. $str")
  println(io,"u: $(sol.u)")
  sol.trueknown && println(io,"errors: $(sol.errors)")
  sol.t!=[] && println(io,"t: $(sol.t)")
  sol.timeseries!=[] && println(io,"timeseries: $(sol.timeseries)")
  sol.trueknown && sol.timeseries_analytic!=[] && println(io,"timeseries_analytic: $(sol.timeseries_analytic)")
  nothing
end

function show(io::IO,sol::DESolution)
  print(io,"$(typeof(sol)), $(length(sol)) timesteps, final value $(sol.u)")
end
