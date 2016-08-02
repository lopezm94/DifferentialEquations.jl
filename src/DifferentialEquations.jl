__precompile__()

module DifferentialEquations

using IterativeSolvers, Parameters, Plots, GenericSVD, SIUnits,
      EllipsisNotation, GrowableArrays, ChunkedArrays
import Base: length, size, getindex, endof, show, print

"PdeProblem: Defines differential equation problems via its internal functions"
abstract DEProblem
"PdeSolution: Wrapper for the objects obtained from a solver"
abstract DESolution
"Mesh: An abstract type which holds a (node,elem) pair and other information for a mesh"
abstract Mesh
"Tableau: Holds the information for a Runge-Kutta Tableau"
abstract Tableau
typealias String AbstractString
typealias KW Dict{Symbol,Any}
AbstractArrayOrVoid = Union{AbstractArray,Void}
NumberOrVoid = Union{Number,Void}
FunctionOrVoid = Union{Function,Void}

#Constants

const TEST_FLOPS_CUTOFF = 1e10
const initialized_backends = Set{Symbol}()

include("general/backends.jl")
include("general/misc_utils.jl")
include("fem/fem_heatintegrators.jl")
include("fem/meshes.jl")
include("fem/fem_assembly.jl")
include("fem/fem_boundary.jl")
include("fem/fem_error.jl")
include("general/problems.jl")
include("general/solutions.jl")
include("general/stochastic_utils.jl")
include("general/convergence.jl")
include("premades/premade_problems.jl")
include("premades/premade_meshes.jl")
include("fem/fem_solve.jl")
include("fdm/stokes_solve.jl")
include("sde/sde_solve.jl")
include("sde/sde_integrators.jl")
include("ode/ode_tableaus.jl")
include("ode/ode_constants.jl")
include("ode/ode_integrators.jl")
include("ode/ode_solve.jl")
include("general/plotrecipes.jl")
include("sde/sde_tableaus.jl")

#Types
export DEProblem, DESolution, HeatProblem, PoissonProblem, FEMSolution, Mesh,
       ConvergenceSimulation, FEMmesh, SimpleMesh, SDEProblem, StokesProblem,
       SDESolution, ODESolution, ODEProblem, FDMMesh, ExplicitRK, MonteCarloSimulation

#SDE Example Problems
export linearSDEExample, cubicSDEExample, waveSDEExample, additiveSDEExample,
       multiDimAdditiveSDEExample, twoDimlinearSDEExample, oval2ModelExample,
       lorenzAttractorSDEExample

#ODE Example Problems
export twoDimlinearODEExample, linearODEExample, lorenzAttractorODEExample

#FEM Example Problems
export  heatProblemExample_moving, heatProblemExample_diffuse, heatProblemExample_pure,
        poissonProblemExample_wave, poissonProblemExample_noisyWave, heatProblemExample_birthdeath,
        poissonProblemExample_birthdeath, heatProblemExample_stochasticbirthdeath,
        homogeneousStokesExample, dirichletzeroStokesExample, poissonProblemExample_birthdeathsystem,
        poissonProblemExample_birthdeathinteractingsystem,heatProblemExample_birthdeathinteractingsystem,
        heatProblemExample_birthdeathsystem,heatProblemExample_grayscott,heatProblemExample_diffusionconstants,heatProblemExample_gierermeinhardt

#Example Meshes
export  meshExample_bunny, meshExample_flowpastcylindermesh, meshExample_lakemesh,
        meshExample_Lshapemesh, meshExample_Lshapeunstructure, meshExample_oilpump,
        meshExample_wavymesh, meshExample_wavyperturbmesh

#Plot Functions
export  plot, animate

#General Functions
export conv_ests, appxTrue!, accumarray, solve, test_convergence, monteCarloSim

#FEM Functions
export  assemblematrix, findboundary, setboundary, findbdtype, getL2error, quadpts, getH1error,
        gradu, gradbasis, quadfbasis, fem_squaremesh, CFLμ, CFLν,
        meshgrid, notime_squaremesh, parabolic_squaremesh, quadpts1

#Tableus
export constructRalston, constructHuen, constructRKF, constructBogakiShampine,
       constructCashKarp, constructDormandPrince, constructRKF8, constructDormandPrince8,
       constructMSRI1,constructFeagin10, constructFeagin12, constructFeagin14

#Misc Tools
export quadfbasis2, CG2, numparameters, checkSRIOrder, checkSRAOrder,
       constructSRIW1, constructSRA1, def

end # module
