module Gillespie

using Distributions
using DataFrames
using Compat: UTF8String, view

export
    ssa,
    ssa_data,
    pfsample,
    SSAArgs,
    SSAStats,
    SSAResult

include("SSA.jl")

function ssa(x0::Vector{Int64},F::Base.Callable,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64; algo=:gillespie, max_rate::Float64=nothing, all_jumps=false)
  @assert algo in [:gillespie,:jensen] "Available algorithms are :gillespie and :jensen"
  if algo == :gillespie
    return gillespie(x0,F,nu,parms,tf)
  end
  if algo == :jensen
    return jensen(x0,F,nu,parms,tf,max_rate)
  end
end

end # module
