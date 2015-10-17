module Gillespie

using Distributions
using StatsBase
using DataFrames
VERSION < v"0.4-dev" && using Lexicon

export
	ssa,
        ssa_data,
	SSAArgs,
	SSAStats,
	SSAResult

include("SSA.jl")

end # module
