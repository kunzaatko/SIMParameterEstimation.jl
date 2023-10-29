using SIMParameterEstimation
using Test
using Aqua

@testset "SIMParameterEstimation.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SIMParameterEstimation;
            # TODO: Find what are the ambiguities and link the repos / issues <26-10-23> 
            ambiguities=(; broken=true)
        )
    end
    # Write your tests here.
end
