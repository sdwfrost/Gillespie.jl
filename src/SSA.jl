@doc doc"""
  A type storing the status at the end of a call to `ssa`.
  """ ->
type SSAStats
  termination_status::ASCIIString
  nsteps::Int64
end

@doc doc"""
  A type storing the call to `ssa`.
  """ ->
type SSAArgs
  x0::Vector{Int64}
  F::Function
  nu::Matrix{Int64}
  parms::Vector{Float64}
  tf::Float64
end

@doc doc"""
  This type stores the output of `ssa`, and comprises of:

      - **time** : a `Vector` of `Float64`, containing the times of simulated events.
      - **data** : a `Matrix` of `Int64`, containing the simulated states.
      - **stats** : an instance of `SSAStats`.
      - **args** : arguments passed to `ssa`.
  """ ->
type SSAResult
  time::Vector{Float64}
  data::Matrix{Int64}
  stats::SSAStats
  args::SSAArgs
end

@doc doc"""
  This function performs Gillespie's stochastic simulation algorithm. It takes the following arguments:

      - **x0** : a `Vector` of `Int64`, representing the initial states of the system.
      - **F** : a `Function`, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system.
      - **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
      - **parms** : a `Vector` of `Float64` representing the parameters of the system.
      - **tf** : the final simulation time (`Float64`)
  """ ->
function ssa(x0::Vector{Int64},F::Function,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64)
  # Args
  args = SSAArgs(x0,F,nu,parms,tf)
  # Set up time array
  ta = Array(Float64,0)
  t = 0.0
  push!(ta,t)
  # Set up initial x
  x = reshape(x0,1,length(x0))
  xa = deepcopy(x)
  # Main loop
  termination_status = "finaltime"
  nsteps = 0
  while t <= tf
    pf = F(x,parms)
    pf = WeightVec(convert(Array{Float64,1},pf))
    # Update time
    sumpf = sum(pf)
    if sumpf == 0.0
        termination_status = "zeroprop"
        break
    end
    dt = rand(Exponential(1/sumpf))
    t = t + dt
    push!(ta,t)
    # Update event
    ev = sample(pf)
    deltax = nu[ev,:]
    x = x .+ deltax
    xa = vcat(xa,x)
    # update nsteps
    nsteps = nsteps + 1
  end
  stats = SSAStats(termination_status,nsteps)
  result = SSAResult(ta,xa,stats,args)
  return(result)
end

@doc doc"""
  This takes a single argument of type `SSAResult` and returns a `DataFrame`.
  """ ->
function ssa_data(s::SSAResult)
  df = hcat(DataFrame(time=s.time),convert(DataFrame,s.data))
  df
end

