# TODO: Add performance regression tests <06-09-24> https://docs.juliahub.com/General/RegressionTests/stable/,
# https://juliaci.github.io/BenchmarkTools.jl/dev/manual/

using SIMParameterEstimation
using Test, Documenter, Aqua, CompatHelperLocal

const run_all = isempty(ARGS) ? true : false

skip = Dict{String,Bool}(
    "compat" => !(VERSION >= v"1.9"), # NOTE: `CompatHelperLocal` only compatible with later Julia version <28-02-25> 
    "aqua" => !haskey(ENV, "GITHUB_ACTIONS") && !haskey(ENV, "RUNTESTS_FULL"),
    "doctests" => !haskey(ENV, "RUNTESTS_FULL") && !(haskey(ENV, "RUNNER_OS") && ENV["RUNNER_OS"] == "Linux"),
    "ambiguities" => true # FIX: Fix the ambiguities <24-04-25> 
)

function should_test(arg::String)::Bool
    global run_all
    if run_all
        return !get(skip, arg, false)
    elseif arg in ARGS
        return true
    end
    return false
end

macro cond_testset(name, block)
    quote
        if should_test($name)
            @testset $name begin
                esc($block)
            end
        end
    end
end

@testset "SIMParameterEstimation.jl" begin
    @testset "Code quality" begin
        @cond_testset "aqua" begin
            Aqua.test_all(
                SIMParameterEstimation;
                ambiguities=false,
            )
        end

        @cond_testset "ambiguities" begin
            aqua_ambiguities = false
            if aqua_ambiguities
                Agua.test_ambiguities(SIMParameterEstimation)
            else
                @test length(Test.detect_ambiguities(SIMParameterEstimation)) == 0
            end
        end

        @cond_testset "compat" begin
            @test CompatHelperLocal.check(SIMParameterEstimation; checktest=false)
        end
    end

    @cond_testset "doctests" begin
        # NOTE: Better than doc-testing in `make.jl` because, I can track the coverage
        DocMeta.setdocmeta!(SIMParameterEstimation, :DocTestSetup, :(
                include(joinpath(@__DIR__, "doctestsetup.jl"));
                using Logging;
                # NOTE: Not necessary in `docs/make.jl`. `@warn` should work there <19-12-24>
                Logging.disable_logging(Logging.Warn)
            ); recursive=true)
        !haskey(ENV, "FIX_DOCTESTS") && @info "You can fix doctests by setting `ENV[\"FIX_DOCTESTS\"] = true`."
        doctest(SIMParameterEstimation; fix=ifelse(haskey(ENV, "FIX_DOCTESTS"), true, false))
    end

    @cond_testset "xcorr" begin
        include("xcorr.jl")
    end
end
