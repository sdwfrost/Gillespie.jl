using Gillespie
using Gadfly
using Debug

function F(x,parms)
  (Y1,Y2) = x
  (c1,c2,c3) = parms
  [c1*Y1,c2*Y1*Y2,c3*Y2]
end

x0 = [1000,1000]
nu = [[1 0];[-1 1];[-1 0]]
parms = [10.0,0.01,10.0]
tf = 2.0
srand(1234)

@debug result = ssa(x0,F,nu,parms,tf)

data = ssa_data(result)

F(x0,parms)
