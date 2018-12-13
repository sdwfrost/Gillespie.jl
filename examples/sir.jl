using Gillespie
import Random: seed!

function F_sir(x,parms)
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
seed!(1234)

result = ssa(x0,F_sir,nu,parms,tf)

data = ssa_data(result)
