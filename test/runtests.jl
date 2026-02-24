using Gillespie
using Test
using Random
using DataFrames

@testset "Gillespie.jl" begin

    @testset "pfsample" begin
        Random.seed!(42)
        w = [0.0, 0.0, 1.0]
        @test pfsample(w, 1.0, 3) == 3
        w = [1.0, 0.0, 0.0]
        @test pfsample(w, 1.0, 3) == 1
    end

    @testset "SSAResult construction" begin
        stats = SSAStats("finaltime", 5)
        @test stats.termination_status == "finaltime"
        @test stats.nsteps == 5
        args = SSAArgs([10, 0], identity, [[-1 1]], [0.1], 10.0, :gillespie, false)
        @test args.tf == 10.0
        @test args.alg == :gillespie
    end

    @testset "gillespie algorithm" begin
        # Simple pure-death process: one species, one reaction
        function F_death(x, parms)
            (N,) = x
            (d,) = parms
            [d * N]
        end
        x0 = [100]
        nu = reshape([-1], 1, 1)
        parms = [0.1]
        tf = 100.0
        Random.seed!(1234)
        result = gillespie(x0, F_death, nu, parms, tf)
        @test result.stats.termination_status in ["finaltime", "zeroprop"]
        @test length(result.time) == size(result.data, 1)
        @test result.data[1, 1] == 100
        # Population should decrease or stay same
        @test result.data[end, 1] <= 100
    end

    @testset "ssa dispatch" begin
        function F_death2(x, parms)
            (N,) = x
            (r,) = parms
            [r * N]
        end
        x0 = [10]
        nu = reshape([-1], 1, 1)
        parms = [0.5]
        tf = 100.0

        Random.seed!(1234)
        result_g = ssa(x0, F_death2, nu, parms, tf, algo=:gillespie)
        @test result_g.stats.termination_status in ["finaltime", "zeroprop"]
    end

    @testset "ssa_data returns DataFrame" begin
        function F_test(x, parms)
            (S, I) = x
            (beta,) = parms
            [beta * S * I]
        end
        x0 = [99, 1]
        nu = [-1 1]
        parms = [0.01]
        tf = 10.0
        Random.seed!(1234)
        result = ssa(x0, F_test, nu, parms, tf)
        df = ssa_data(result)
        @test df isa DataFrame
        @test names(df) == ["time", "x1", "x2"]
        @test size(df, 1) == length(result.time)
        @test df.time == result.time
    end

    @testset "Example: Lotka-Volterra" begin
        include("../examples/lotka.jl")
        @test gillespie_data isa DataFrame
        @test size(gillespie_data, 2) == 3  # time + 2 species
        @test gillespie_data[1, 2] == 1000
        @test gillespie_data[1, 3] == 1000
        # Populations should remain non-negative
        @test all(gillespie_data[:, 2] .>= 0)
        @test all(gillespie_data[:, 3] .>= 0)
        # Jensen results should also be valid DataFrames
        @test jensen_data isa DataFrame
        @test jensen_alljumps_data isa DataFrame
    end

    @testset "Example: Logistic growth" begin
        include("../examples/logistic_growth.jl")
        @test data isa DataFrame
        @test data[1, 2] == 10  # initial population
        # Population should grow toward carrying capacity K=1000
        @test data[end, 2] > 100
        @test all(data[:, 2] .>= 0)
    end

    @testset "Example: Decaying dimer" begin
        include("../examples/decaying_dimer.jl")
        @test data isa DataFrame
        @test size(data, 2) == 4  # time + 3 species
        @test data[1, 2] == 10000
        @test data[1, 3] == 0
        @test data[1, 4] == 0
        # All species non-negative
        @test all(data[:, 2] .>= 0)
        @test all(data[:, 3] .>= 0)
        @test all(data[:, 4] .>= 0)
    end

    @testset "Example: SIR" begin
        include("../examples/sir.jl")
        @test data isa DataFrame
        @test size(data, 2) == 4  # time + S, I, R
        # Conservation law: S + I + R = 1000
        for row in eachrow(data)
            @test row[2] + row[3] + row[4] == 1000
        end
        # S should be non-increasing (infections only decrease S)
        @test issorted(data[:, 2], rev=true)
    end

    @testset "Example: SIR2" begin
        include("../examples/sir2.jl")
        @test data isa DataFrame
        @test size(data, 2) == 4
        for row in eachrow(data)
            @test row[2] + row[3] + row[4] == 1000
        end
    end

end
