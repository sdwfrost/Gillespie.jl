using Gillespie
using Compat

function F_sir2(x,parms)
    (S,I,R) = x
    (beta,gamma) = parms
    infection = beta*S*I
    recovery = gamma*I
    [infection,recovery]
end

x0 = [999,1,0]
nu = [[-1 1 0];[0 -1 1]]
parms = [0.1/1000.0,0.01]
tf = 250.0
srand(1234)
immutable G_sir2; end
@compat (::Type{G_sir2})(x,parms) = F_sir2(x,parms)

result = ssa(x0,G_sir2,nu,parms,tf)

data = ssa_data(result)

