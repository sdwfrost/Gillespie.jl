using DataFrames
using Distributions
using Gillespie
using BenchmarkTools

function sir(beta,gamma,N,S0,I0,R0,tf)
    t = 0
    S = S0
    I = I0
    R = R0
    ta=Vector{Float64}(undef,0)
    Sa=Vector{Float64}(undef,0)
    Ia=Vector{Float64}(undef,0)
    Ra=Vector{Float64}(undef,0)
    while t < tf
        push!(ta,t)
        push!(Sa,S)
        push!(Ia,I)
        push!(Ra,R)
        pf1 = beta*S*I
        pf2 = gamma*I
        pf = pf1+pf2
        dt = rand(Exponential(1/pf))
        t = t+dt
        if t>tf
            break
        end
        ru = rand()
        if ru<(pf1/pf)
            S=S-1
            I=I+1
        else
            I=I-1
            R=R+1
        end
    end
    results = DataFrame()
    results[!,:t] = ta
    results[!,:S] = Sa
    results[!,:I] = Ia
    results[!,:R] = Ra
    return(results)
end

function F(x,parms)
  (S,I,R) = x
  (beta,gamma) = parms
  infection = beta*S*I
  recovery = gamma*I
  [infection,recovery]
end

x0 = [9999,1,0]
nu = [[-1 1 0];[0 -1 1]]
parms = [0.1/10000.0,0.05]
tf = 1000.0

ssa(x0,F,nu,parms,tf) # compile
Random.seed!(1234)
@benchmark ssa($x0,$F,$nu,$parms,$tf) samples=1000 seconds=100

function F2(x,parms)
  (S,I,R) = x
  (beta,gamma) = parms
  infection = beta*S*I
  recovery = gamma*I
  [infection,recovery]
end

x0 = SA[9999,1,0]
nu = SA[-1 1 0
         0 -1 1]
parms = SA[0.1/10000.0,0.05]
tf = 1000.0

ssa(x0,F2,nu,parms,tf) # compile
Random.seed!(1234)
@benchmark ssa($x0,$F2,$nu,$parms,$tf) samples=1000 seconds=100

sir(0.1/10000,0.05,10000,9999,1,0,1000) # compile
Random.seed!(1234)
@benchmark sir(0.1/10000,0.05,10000,9999,1,0,1000) samples=1000 seconds=100

using DiffEqBiological

sir_model = @reaction_network rn begin
  0.1/10000.0, s + i --> 2i
  0.05, i --> r
end
sir_prob = DiscreteProblem([9999,1,0],(0.0,tf))
sir_jump_prob = JumpProblem(sir_prob,Direct(),sir_model)
sir_sol = solve(sir_jump_prob,SSAStepper()) # compile
@benchmark solve(sir_jump_prob,SSAStepper()) samples=1000 seconds=100

sir_prob = DiscreteProblem(SA[9999,1,0],(0.0,tf))
sir_jump_prob = JumpProblem(sir_prob,Direct(),sir_model)
sir_sol = solve(sir_jump_prob,SSAStepper()) # compile
@benchmark solve(sir_jump_prob,SSAStepper()) samples=1000 seconds=100
