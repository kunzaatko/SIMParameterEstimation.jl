using SIMParameterEstimation
using Test
using Aqua

@testset "SIMParameterEstimation.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SIMParameterEstimation)
    end
    # Write your tests here.
end
