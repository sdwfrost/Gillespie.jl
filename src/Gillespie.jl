module Gillespie

using Distributions
using DataFrames
using QuadGK
using Roots

export
    ssa,
    gillespie,
    jensen,
    truejump,
    ssa_data,
    pfsample,
    SSAArgs,
    SSAStats,
    SSAResult

include("SSA.jl")

"
This function performs stochastic simulation, currently using either the Doob-Gillespie algorithm or Jensen's method. It takes the following arguments:

- **x0** : a `Vector` of `Int64`, representing the initial states of the system.
- **F** : a `Function` or a callable type, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system. In the case of time-varying rates (for algorithms `:jensen` and `:tjm`), there should be a third argument, the time of the system.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`)

There are several named arguments:

- **algo**: the algorithm to use (`Symbol`, either `:gillespie` (default), ':jensen', or ':tjm').
- **max_rate**: the maximum rate (`Float64`, for Jensen's method only).
- **thin**: (`Bool`) whether to thin jumps for Jensens method (default: `true`).

"
function ssa(x0::Vector{Int64},F::Base.Callable,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64; algo=:gillespie, max_rate::Float64=0.0, thin::Bool=true)
  @assert algo in [:gillespie,:jensen,:tjm] "Available algorithms are :gillespie, :jensen, and :tjm"
  if algo == :gillespie
    return gillespie(x0,F,nu,parms,tf)
  end
  if algo == :jensen
    return jensen(x0,F,nu,parms,tf,max_rate,thin)
  end
  if algo == :tjm
    return tjm(x0,F,nu,parms,tf)
  end
end

end # module
