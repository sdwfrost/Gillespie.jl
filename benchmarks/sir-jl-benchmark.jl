
using DataFrames
using DataArrays
using Distributions
using Gillespie
using BenchmarkTools

function sir(beta,gamma,N,S0,I0,R0,tf)
    t = 0
    S = S0
    I = I0
    R = R0
    ta=DataArray(Float64,0)
    Sa=DataArray(Float64,0)
    Ia=DataArray(Float64,0)
    Ra=DataArray(Float64,0)
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
    results[:t] = ta
    results[:S] = Sa
    results[:I] = Ia
    results[:R] = Ra
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
srand(1234)
@benchmark ssa($x0,$F,$nu,$parms,$tf) samples=1000 seconds=100

immutable G; end
call(::Type{G},x,parms) = F(x,parms)
ssa(x0,G,nu,parms,tf) # compile
srand(1234)
@benchmark ssa($x0,$G,$nu,$parms,$tf) samples=1000 seconds=100

sir(0.1/10000,0.05,10000,9999,1,0,1000) # compile
srand(1234)
@benchmark sir(0.1/10000,0.05,10000,9999,1,0,1000) samples=1000 seconds=100
