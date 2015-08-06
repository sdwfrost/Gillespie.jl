using Gillespie
using Gadfly

function F(x,parms)
  (N) = x
  (b,d,K) = parms
  [b*N,(d+(b-d)*N/K)*N]
end

x0 = [500]
nu = reshape([[1];[-1]],2,1)
parms = [2.0,1.0,1000.0]
tf = 10.0
srand(1234)

result = ssa(x0,F,nu,parms,tf)

data = ssa_data(result)

