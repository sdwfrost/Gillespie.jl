using Gillespie
using Gadfly

function F(x,parms)
  (S1,S2,S3) = x
  (c1,c2,c3,c4) = parms
  [c1*S1,c2*S1*S1,c3*S2,c4*S2]
end

x0 = [10000,0,0]
nu = [[-1 0 0];[-2 1 0];[2 -1 0];[0 -1 1]]
parms = [1.0,0.002,0.5,0.04]
tf = 10.0
srand(1234)

result = ssa(x0,F,nu,parms,tf)

data = ssa_data(result)

p=plot(data,
  layer(x="time",y="x1",Geom.step,Theme(default_color=colorant"red")),
  layer(x="time",y="x2",Geom.step,Theme(default_color=colorant"blue")),
  layer(x="time",y="x3",Geom.step,Theme(default_color=colorant"green")),
  Guide.xlabel("Time"),
  Guide.ylabel("Number"),
  Guide.manual_color_key("Population",
                            ["S1", "S2", "S3"],
                            ["red", "blue", "green"]),
  Guide.title("Decaying dimer"))
