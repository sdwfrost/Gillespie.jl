using Gillespie
import Random: seed!

function F_l(x,parms)
    (Y1,Y2) = x
    (c1,c2,c3) = parms
    [c1*Y1,c2*Y1*Y2,c3*Y2]
end

x0 = [1000,1000]
nu = [[1 0];[-1 1];[0 -1]]
parms = [10.0,0.01,10.0]
tf = 2.0
seed!(1234)

gillespie_result = ssa(x0,F_l,nu,parms,tf)
gillespie_data = ssa_data(gillespie_result)

jensen_result = ssa(x0,F_l,nu,parms,tf,algo=:jensen,max_rate=100000.0,thin=true)
jensen_data = ssa_data(jensen_result)

jensen_alljumps_result = ssa(x0,F_l,nu,parms,tf,algo=:jensen,max_rate=100000.0,thin=false)
jensen_alljumps_data = ssa_data(jensen_alljumps_result)
