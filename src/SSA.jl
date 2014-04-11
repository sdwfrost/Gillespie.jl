type SSAStats
  termination_status::ASCIIString
  nsteps::Int64
end

type SSAArgs
  x0::Vector{Int64}
  F::Function
  nu::Matrix{Int64}
  parms::Vector{Float64}
  tf::Float64
end

type SSAResult
  time::Vector{Float64}
  data::Matrix{Int64}
  stats::SSAStats
  args::SSAArgs
end

function ssa(x0::Vector{Int64},F::Function,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64)
  # Args
  args = SSAArgs(x0,F,nu,parms,tf)
  # Set up time array
  ta=Array(Float64,0)
  t=0.0
  push!(ta,t)
  # Set up initial x
  x=reshape(x0,1,length(x0))
  xa=deepcopy(x)
  # Main loop
  termination_status = "finaltime"
  nsteps = 0
  while t<= tf
    pf = F(x,parms)
    pf = WeightVec(convert(Array{Float64,1},pf))
    # Update time
    sumpf = sum(pf)
    if sumpf==0.0
        termination_status = "zeroprop"
        break
    end
    dt = rand(Exponential(1/sumpf))
    t = t + dt
    push!(ta,t)
    # Update event
    ev = sample(pf)
    deltax = nu[ev,:]
    x =x .+ deltax
    xa=vcat(xa,x)
    # update nsteps
    nsteps = nsteps + 1
  end
  stats = SSAStats(termination_status,nsteps)
  result = SSAResult(ta,xa,stats,args)
  return(result)
end

function ssa_data(s::SSAResult)
  df = cbind(DataFrame(time=s.time),convert(DataFrame,s.data))
  df
end
