using Gillespie
using Base.Test

include("../examples/lotka.jl")
@test isequal(data[end,2],741)
@test isequal(data[end,3],1048)

include("../examples/logistic_growth.jl")
@test isequal(data[end,2],1008)

include("../examples/decaying_dimer.jl")
@test isequal(data[end,2],312)
@test isequal(data[end,3],714)
@test isequal(data[end,4],856)

include("../examples/sir.jl")
@test isequal(data[end,2],1)
@test isequal(data[end,3],181)
@test isequal(data[end,4],818)

include("../examples/sir2.jl")
@test isequal(data[end,2],1)
@test isequal(data[end,3],181)
@test isequal(data[end,4],818)
