"A type storing the status at the end of a call to `ssa`."
type SSAStats
  termination_status::ASCIIString
  nsteps::Int64
end

"A type storing the call to `ssa`."
type SSAArgs
  x0::Vector{Int64}
  F::Any
  nu::Matrix{Int64}
  parms::Vector{Float64}
  tf::Float64
end

"
This type stores the output of `ssa`, and comprises of:

- **time** : a `Vector` of `Float64`, containing the times of simulated events.
- **data** : a `Matrix` of `Int64`, containing the simulated states.
- **stats** : an instance of `SSAStats`.
- **args** : arguments passed to `ssa`.

"
type SSAResult
  time::Vector{Float64}
  data::Matrix{Int64}
  stats::SSAStats
  args::SSAArgs
end

"
This function is a substitute for `StatsBase.sample(wv::WeightVec)`, which avoids recomputing the sum and size of the weight vector, as well as a type conversion of the propensity vector. It takes the following arguments:

- **w** : an `Array{Float64,1}`, representing propensity function weights.
- **s** : the sum of `w`.
- **n** : the length of `w`.

"
function pfsample(w::Array{Float64,1},s::Float64,n::Int64)
    t = rand() * s
    i = 1
    cw = w[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += w[i]
    end
    return i
end

"
This function performs Gillespie's stochastic simulation algorithm. It takes the following arguments:

- **x0** : a `Vector` of `Int64`, representing the initial states of the system.
- **F** : a `Function` or a callable type, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`)

"
function ssa(x0::Vector{Int64},F::Function,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64)
  # Args
  args = SSAArgs(x0,F,nu,parms,tf)
  # Set up time array
  ta = Array(Float64,0)
  t = 0.0
  push!(ta,t)
  # Set up initial x
  nstates = length(x0)
  x = copy(x0)
  x = reshape(x,1,nstates)
  xa = vec(x)
  # Number of propensity functions
  numpf = size(nu,1)
  # Main loop
  termination_status = "finaltime"
  nsteps = 0
  while t <= tf
    pf = F(x,parms)
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
    ev = pfsample(pf,sumpf,numpf)
    deltax = nu[ev,:]
    for i in 1:nstates
      @inbounds x[1,i] = x[1,i]+deltax[i]
    end
    for xx in x
      push!(xa,xx)
    end
    # update nsteps
    nsteps = nsteps + 1
  end
  stats = SSAStats(termination_status,nsteps)
  xar = transpose(reshape(xa,length(x),nsteps+1))
  result = SSAResult(ta,xar,stats,args)
  return(result)
end

function ssa{F}(x0::Vector{Int64},::Type{F},nu::Matrix{Int64},parms::Vector{Float64},tf::Float64)
  # Args
  args = SSAArgs(x0,F,nu,parms,tf)
  # Set up time array
  ta = Array(Float64,0)
  t = 0.0
  push!(ta,t)
  # Set up initial x
  nstates = length(x0)
  x = copy(x0)
  x = reshape(x,1,nstates)
  xa = vec(x)
  # Number of propensity functions
  numpf = size(nu,1)
  # Main loop
  termination_status = "finaltime"
  nsteps = 0
  while t <= tf
    pf = F(x,parms)
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
    ev = pfsample(pf,sumpf,numpf)
    deltax = nu[ev,:]
    for i in 1:nstates
      @inbounds x[1,i] = x[1,i]+deltax[i]
    end
    for xx in x
      push!(xa,xx)
    end
    # update nsteps
    nsteps = nsteps + 1
  end
  stats = SSAStats(termination_status,nsteps)
  xar = transpose(reshape(xa,length(x),nsteps+1))
  result = SSAResult(ta,xar,stats,args)
  return(result)
end

"This takes a single argument of type `SSAResult` and returns a `DataFrame`."
function ssa_data(s::SSAResult)
  df = hcat(DataFrame(time=s.time),convert(DataFrame,s.data))
  df
end
