using Gillespie
using Gadfly

function F(x,parms)
  (S,I) = x
  (beta,mu) = parms
  infection = beta*S*I
  recovery = mu*I
  [infection,recovery]
end

x0 = [999,1,0]
nu = [[-1 1 0];[0 -1 1]]
parms = [0.1/1000.0,0.05]
tf = 250.0

result = ssa(x0,F,nu,parms,tf)

data = ssa_data(result)

p=plot(data,x="time",y="x2",Geom.step,Guide.xlabel("Time"), Guide.ylabel("I"), Guide.title("SIR"))
draw(PNG("sir.png", 4inch, 3inch), p)
