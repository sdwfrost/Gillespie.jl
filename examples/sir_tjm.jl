using Gillespie

function F_sir(x,parms,t)
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

result = tjm(x0,F_sir,nu,parms,tf)

data = ssa_data(result)
