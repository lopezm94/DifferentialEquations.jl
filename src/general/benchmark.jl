BenchmarkSimulation
  solutions::Vector{Array{DESolution}}
  elapsed_times
  algs::Vector{Symbol}
  errors
  N
  auxdata
  benchmark_axis
  function BenchmarkSimulation(solutions::Vector{Array{DESolution}},elapsed_times,algs,benchmark_axis;auxdata=nothing)
    N = size(solutions[1],1)
    errors = Dict() #Should add type information
    for solArray in solutions
      for k in keys(solArray[1].errors)
        errors[k] = reshape(Float64[sol.errors[k] for sol in solArray],size(solArray)...)
      end
    end
    return(new(solutions,,elapsed_times,algs,errors,N,auxdata,benchmark_axis))
  end
  end
  #function ConvergenceSimulation(solutions::Array{SDESolution},convergence_axis;auxdata=nothing)
  #ConvergenceSimulation(convert(Array{DESolution},solutions),convergence_axis;auxdata=auxdata)
end


function benchmark_algs(Δts::AbstractArray,prob::ODEProblem;tspan=[0,1],save_timeseries=true,adaptive=false,kwargs...)

end
