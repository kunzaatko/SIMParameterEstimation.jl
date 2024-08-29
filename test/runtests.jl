using SIMParameterEstimation
using Test
using Aqua

@testset "SIMParameterEstimation.jl" begin
    if haskey(ENV, "RUNTESTS_FULL") || haskey(ENV, "GITHUB_ACTIONS")
        @testset "Code quality (Aqua.jl)" begin
            Aqua.test_all(
                SIMParameterEstimation;
                ambiguities=(; exclude=VERSION >= v"1.11" ? [checkindex, checkbounds] : []),
                unbound_args=(; broken=true)
            )
        end
    else
        @info "Skipping Aqua.jl quality tests. For a full run set `ENV[\"RUNTESTS_FULL\"]=true`."
    end
    @testset "utils.jl" begin
        using LinearAlgebra
        using SIMParameterEstimation: mixin_matrix, separation_matrix, mixin_components, separate_components
        @test mixin_matrix((0.5, 0.8, 1.0)) isa Matrix{<:Complex}
        @test mixin_matrix(0.0, 3) == mixin_matrix((0.0, 2π / 3, 4π / 3))
        @test mixin_matrix(0.5, 5) == mixin_matrix(0.5, (1, 1, 1, 1, 1))

        @test separation_matrix(0.0, 3) isa Matrix{<:Complex}
        @test separation_matrix((0.5, 0.8, 1.0)) isa Matrix{<:Complex}
        @test separation_matrix((0.5, 0.8, 1.0), (0.3, 0.4, 0.2)) isa Matrix{<:Complex}
        @test_broken separation_matrix(0.5, (1, 1, 1, 1, 1)) isa Matrix{<:Complex}

        M = mixin_matrix(0.0, 3)
        M_inv = separation_matrix(0.0, 3)
        comps = randn(10, 10, 3, 1)
        @test mixin_components(comps[:, :, :], M) isa Array{<:Complex,4}
        @test mixin_components(comps, M) == mixin_components(comps[:, :, :], M)
        raw = mixin_components(comps, M)
        @test separate_components(raw, M_inv) == separate_components(raw[:, :, :], M_inv)
        @test comps ≈ separate_components(raw, M_inv)

        @test separation_matrix((1., 2., 3.), (2., 2., 2.)) * transpose([1 1 1; exp(im) exp(2im) exp(3im); exp(-im) exp(-2im) exp(-3im)]) ≈ I(3)
    end
end
