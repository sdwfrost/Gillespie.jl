using Gillespie
import Random: seed!

function F_lg(x,parms)
    (N,) = x
    (b,d,K) = parms
    [b*N,(d+(b-d)*N/K)*N]
end

x0 = [10]
nu = reshape([[1];[-1]],2,1)
parms = [2.0,1.0,1000.0]
tf = 15.0
seed!(1234)

result = ssa(x0,F_lg,nu,parms,tf)

data = ssa_data(result)
