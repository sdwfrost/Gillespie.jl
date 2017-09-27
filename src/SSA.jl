"A type storing the status at the end of a call to `ssa`:

- **termination_status** : whether the simulation stops at the final time (`finaltime`) or early due to zero propensity function (`zeroprop`)
- **nsteps** : the number of steps taken during the simulation.

"
type SSAStats
    termination_status::String
    nsteps::Int64
end

"A type storing the call to `ssa`:

- **x0** : a `Vector` of `Int64`, representing the initial states of the system.
- **F** : a `Function` or a callable type, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`).
- **alg** : the algorithm used (`Symbol`, either `:gillespie`, `jensen`, or `tjc`).
- **tvc** : whether rates are time varying.
"
type SSAArgs
    x0::Vector{Int64}
    F::Any
    nu::Matrix{Int64}
    parms::Vector{Float64}
    tf::Float64
    alg::Symbol
    tvc::Bool
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
- **tf** : the final simulation time (`Float64`).
"
function gillespie(x0::Vector{Int64},F::Base.Callable,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64)
    # Args
    args = SSAArgs(x0,F,nu,parms,tf,:gillespie,false)
    # Set up time array
    ta = Vector{Float64}(0)
    t = 0.0
    push!(ta,t)
    # Set up initial x
    nstates = length(x0)
    x = x0'
    xa = copy(x0)
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
        t += dt
        push!(ta,t)
        # Update event
        ev = pfsample(pf,sumpf,numpf)
        deltax = view(nu,ev,:)
        for i in 1:nstates
            @inbounds x[1,i] += deltax[i]
        end
        for xx in x
            push!(xa,xx)
        end
        # update nsteps
        nsteps += 1
    end
    stats = SSAStats(termination_status,nsteps)
    xar = transpose(reshape(xa,length(x),nsteps+1))
    return SSAResult(ta,xar,stats,args)
end

"
This function performs the true jump method for piecewise deterministic Markov processes. It takes the following arguments:

- **x0** : a `Vector` of `Int64`, representing the initial states of the system.
- **F** : a `Function` or a callable type, which itself takes three arguments; x, a `Vector` of `Int64` representing the states, parms, a `Vector` of `Float64` representing the parameters of the system, and t, a `Float64` representing the time of the system.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`).
"
function truejump(x0::Vector{Int64},F::Base.Callable,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64)
    # Args
    args = SSAArgs(x0,F,nu,parms,tf,:tjm,true)
    # Set up time array
    ta = Vector{Float64}(0)
    t = 0.0
    push!(ta,t)
    # Set up initial x
    nstates = length(x0)
    x = x0'
    xa = copy(x0)
    # Number of propensity functions
    numpf = size(nu,1)
    # Main loop
    termination_status = "finaltime"
    nsteps = 0
    while t <= tf
        ds = rand(Exponential(1.0))
        f = (u)->(quadgk((u)->sum(F(x,parms,u)),t,u)[1]-ds)
        newt = fzero(f,t)
        if newt>tf
          break
        end
        t=newt
        pf = F(x,parms,t)
        # Update time
        sumpf = sum(pf)
        if sumpf == 0.0
            termination_status = "zeroprop"
            break
        end
        push!(ta,t)
        # Update event
        ev = pfsample(pf,sumpf,numpf)
        deltax = view(nu,ev,:)
        for i in 1:nstates
            @inbounds x[1,i] += deltax[i]
        end
        for xx in x
            push!(xa,xx)
        end
        # update nsteps
        nsteps += 1
    end
    stats = SSAStats(termination_status,nsteps)
    xar = transpose(reshape(xa,length(x),nsteps+1))
    return SSAResult(ta,xar,stats,args)
end

"
This function performs stochastic simulation using thinning/uniformization/Jensen's method, returning only the thinned jumps. It takes the following arguments:

- **x0** : a `Vector` of `Int64`, representing the initial states of the system.
- **F** : a `Function` or a callable type, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system. In the case of time-varying systems, a third argument, a `Float64` representing the time of the system should be added
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`).
- **max_rate**: the maximum rate (`Float64`).
"
function jensen(x0::Vector{Int64},F::Base.Callable,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64,max_rate::Float64,thin::Bool=true)
    if thin==false
      return jensen_alljumps(x0::Vector{Int64},F::Base.Callable,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64,max_rate::Float64)
    end
    tvc=true
    try
      F(x0,parms,0.0)
    catch
      tvc=false
    end
    # Args
    args = SSAArgs(x0,F,nu,parms,tf,:jensen,tvc)
    # Set up time array
    ta = Vector{Float64}(0)
    t = 0.0
    push!(ta,t)
    # Set up initial x
    nstates = length(x0)
    x = x0'
    xa = copy(x0)
    # Number of propensity functions; one for no event
    numpf = size(nu,1)+1
    # Main loop
    termination_status = "finaltime"
    nsteps = 0
    while t <= tf
        dt = rand(Exponential(1/max_rate))
        t += dt
        if tvc
          pf = F(x,parms,t)
        else
          pf = F(x,parms)
        end
        # Update time
        sumpf = sum(pf)
        if sumpf == 0.0
            termination_status = "zeroprop"
            break
        end
        if sumpf > max_rate
            termination_status = "upper_bound_exceeded"
            break
        end
        # Update event
        ev = pfsample([pf; max_rate-sumpf],max_rate,numpf+1)
        if ev < numpf
          deltax = view(nu,ev,:)
          for i in 1:nstates
              @inbounds x[1,i] += deltax[i]
          end
          for xx in x
            push!(xa,xx)
          end
          push!(ta,t)
          # update nsteps
          nsteps += 1
        end
    end
    stats = SSAStats(termination_status,nsteps)
    xar = transpose(reshape(xa,length(x),nsteps+1))
    return SSAResult(ta,xar,stats,args)
end

"
This function performs stochastic simulation using thinning/uniformization/Jensen's method, returning all the jumps, both real and 'virtual'. It takes the following arguments:

- **x0** : a `Vector` of `Int64`, representing the initial states of the system.
- **F** : a `Function` or a callable type, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system. In the case of time-varying systems, a third argument, a `Float64` representing the time of the system should be added
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`).
- **max_rate**: the maximum rate (`Float64`).
"
function jensen_alljumps(x0::Vector{Int64},F::Base.Callable,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64,max_rate::Float64)
    # Args
    tvc=true
    try
      F(x0,parms,0.0)
    catch
      tvc=false
    end
    # Args
    args = SSAArgs(x0,F,nu,parms,tf,:jensen,tvc)
    # Set up time array
    ta = Vector{Float64}(0)
    t = 0.0
    push!(ta,t)
    while t < tf
      dt = rand(Exponential(1/max_rate))
      t += dt
      push!(ta,t)
    end
    nsteps=length(ta)-1
    # Set up initial x
    nstates = length(x0)
    x = x0'
    xa = Array{Int64,1}((nsteps+1)*nstates)
    xa[1:nstates] = x
    # Number of propensity functions; one for no event
    numpf = size(nu,1)+1
    # Main loop
    termination_status = "finaltime"
    k=1 # step counter
    while k <= nsteps
        if tvc
          t=ta[k]
          pf=F(x,parms,t)
        else
          pf = F(x,parms)
        end
        sumpf = sum(pf)
        if sumpf == 0.0
            termination_status = "zeroprop"
            break
        end
        if sumpf > max_rate
            termination_status = "upper_bound_exceeded"
            break
        end
        # Update event
        ev = pfsample([pf; max_rate-sumpf],max_rate,numpf+1)
        if ev < numpf # true jump
          deltax = view(nu,ev,:)
          for i in 1:nstates
              @inbounds xa[k*nstates+i] = xa[(k-1)*nstates+i]+deltax[i]
          end
        else
          for i in 1:nstates
              @inbounds xa[k*nstates+i] = xa[(k-1)*nstates+i]
          end
        end
        k +=1
    end
    stats = SSAStats(termination_status,nsteps)
    xar = transpose(reshape(xa,length(x),nsteps+1))
    return SSAResult(ta,xar,stats,args)
end

"This takes a single argument of type `SSAResult` and returns a `DataFrame`."
function ssa_data(s::SSAResult)
    hcat(DataFrame(time=s.time),convert(DataFrame,s.data))
end
