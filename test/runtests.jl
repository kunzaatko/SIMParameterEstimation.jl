# TODO: Add performance regression tests <06-09-24> https://docs.juliahub.com/General/RegressionTests/stable/,
# https://juliaci.github.io/BenchmarkTools.jl/dev/manual/
using SIMParameterEstimation
using Test
using Aqua

include("helpers.jl")

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

    @testset_skip "Passes... Testing other" "utils.jl" begin
        # @testset "utils.jl" begin
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

        @test separation_matrix((1.0, 2.0, 3.0), (2.0, 2.0, 2.0)) * transpose([1 1 1; exp(im) exp(2im) exp(3im); exp(-im) exp(-2im) exp(-3im)]) ≈ I(3)
    end

    @testset "xcorr.jl" begin
        using OffsetArrays: Origin, centered
        using SIMParameterEstimation: WindowOverlap, TrueOverlap, Global

        @testset_skip "Passes... Testing other" "SummingSet" begin
            # @testset "SummingSet" begin
            img = ones(1024, 1024)
            img_shifted = Origin(-300, -200)(ones(1024, 1024))

            using SIMParameterEstimation: Window

            @test Window(img) isa Window
            w_img = Window(img)
            @test Base.Indices(w_img) == (1:1024, 1:1024)

            @test Window(img_shifted) isa Window
            w_img_shifted = Window(img_shifted)
            @test Base.Indices(w_img_shifted) == (-300:(1024-301), -200:(1024-201))

            @test w_img .- CartesianIndex(301, 201) == w_img_shifted
            @test w_img .+ CartesianIndex(-301, -201) == w_img_shifted

            # FIX: Same but hash not same... The inner Cartesian Indices are not of the same type... Should be reported
            # <02-09-24> 
            @test_broken hash(w_img .- CartesianIndex(301, 201)) == hash(w_img_shifted)

            @test intersect(w_img, w_img_shifted) == Window(CartesianIndices((1:(1024-301), 1:(1024-201))))
            @test length(w_img) == 1024 * 1024
            @test length(w_img_shifted) == 1024 * 1024

            using IntegralArrays
            I_img = IntegralArray(img)

            @test sum(img) == I_img[w_img]

            # SummingSets
            @test WindowOverlap(img, img_shifted) isa WindowOverlap
            @test_throws MethodError WindowOverlap(ones(3, 3, 3), img) # NOTE: Accept only same dimensions  

            @test TrueOverlap(img, img_shifted) isa TrueOverlap
            @test TrueOverlap(img .== 1, img_shifted .== 1) isa TrueOverlap
            @test TrueOverlap(img .== 1, img_shifted .== 1) == TrueOverlap(img, img_shifted)
            @test_throws MethodError TrueOverlap(img, rand(Bool, 3, 3, 3)) # NOTE: Accept only same dimensions

            @test Global(img, img_shifted) isa Global
        end

        # @testset_skip "Passes... Testing other" "xcorr helpers" begin
        @testset "xcorr helpers" begin
            using SIMParameterEstimation: Direct, Transform, means, energies, means!, energies!

            @testset "SummingSet $set" for set in (WindowOverlap, TrueOverlap)
                # Testing correct output types
                real_input = centered(rand(30, 30))
                @test means(real_input, real_input, set)[1] isa AbstractArray{<:Real}
                @test energies(real_input, real_input, set)[1] isa AbstractArray{<:Real}

                using ImageFiltering: reflect

                s_m, t_m = means(real_input, real_input, set)
                @test s_m.parent[1:(end-1), 1:(end-1)] ≈ reflect(t_m).parent[2:end, 2:end]
                s_e, t_e = energies(real_input, real_input, set)
                @test s_e.parent[1:(end-1), 1:(end-1)] ≈ reflect(t_e).parent[2:end, 2:end]

                real_input_odd = centered(rand(31, 31))
                s_m_odd, t_m_odd = means(real_input_odd, real_input_odd, set)
                @test s_m_odd.parent ≈ reflect(t_m_odd).parent
                s_e_odd, t_e_odd = energies(real_input_odd, real_input_odd, set)
                @test s_e_odd.parent ≈ reflect(t_e_odd).parent

                complex_input = centered(rand(ComplexF32, 30, 30))
                @test means(complex_input, complex_input, set)[1] isa AbstractArray{<:Complex}
                @test energies(complex_input, complex_input, set)[1] isa AbstractArray{<:Real}

                s_m, t_m = means(complex_input, complex_input, set)
                @test s_m.parent[1:(end-1), 1:(end-1)] ≈ reflect(t_m).parent[2:end, 2:end]
                s_e, t_e = energies(complex_input, complex_input, set)
                @test s_e.parent[1:(end-1), 1:(end-1)] ≈ reflect(t_e).parent[2:end, 2:end]

                complex_input_odd = centered(rand(ComplexF32, 31, 31))
                s_m_odd, t_m_odd = means(complex_input_odd, complex_input_odd, set)
                @test s_m_odd.parent ≈ reflect(t_m_odd).parent
                s_e_odd, t_e_odd = energies(complex_input_odd, complex_input_odd, set)
                @test s_e_odd.parent ≈ reflect(t_e_odd).parent
            end

            @testset "Equivalence with $(eltype(input)) $(isodd(size(input, 1)) ? "odd" : "even")" for input in (rand(30, 30), rand(ComplexF32, 30, 30), rand(31, 31), rand(ComplexF32, 31, 31))
                input = centered(input)
                # NOTE: Support over full window
                @test (means(input, input, WindowOverlap) .≈ means(input, input, TrueOverlap)) |> all
                @test (energies(input, input, WindowOverlap) .≈ energies(input, input, TrueOverlap)) |> all
                @test (energies!(similar(input), similar(input), input, input, TrueOverlap, Direct()) .≈ energies!(similar(input), similar(input), input, input, TrueOverlap, Transform())) |> all
            end
        end
    end

    @testset "xcorr" begin
        using SIMParameterEstimation: xcorr

        # @testset_skip "Passes... Testing other" "signature" begin
        @testset "signature" begin
            # NOTE: Do not accept Real signal with Complex template <02-09-24> 
            @test_throws ArgumentError xcorr(rand(Float32, 30, 30), rand(ComplexF32, 30, 30))
        end

        # @testset_skip "Passes... Testing other" "Real linearly dependent" begin
        @testset "Real linearly dependent" begin
            # NOTE: Linearly dependent signals <02-09-24> 
            signal_rand = centered(rand(30, 30))
            template_rand_lin = rand() .* signal_rand

            @test xcorr(signal_rand, template_rand_lin, WindowOverlap)[0, 0] ≈ 1
            @test xcorr(signal_rand, -1 .* template_rand_lin, WindowOverlap)[0, 0] ≈ -1

            @test xcorr(signal_rand, template_rand_lin, TrueOverlap)[0, 0] ≈ 1
            @test xcorr(signal_rand, -1 .* template_rand_lin, TrueOverlap)[0, 0] ≈ -1

            # NOTE: Support over full window
            @test xcorr(signal_rand, template_rand_lin, WindowOverlap) ≈ xcorr(signal_rand, template_rand_lin, TrueOverlap)
        end

        # @testset_skip "Passes... Testing other" "Complex linearly dependent" begin
        @testset "Complex linearly dependent" begin
            signal_rand_complex = centered(rand(ComplexF32, 30, 30))
            template_rand_complex_lin = rand() .* signal_rand_complex

            @test xcorr(signal_rand_complex, template_rand_complex_lin, WindowOverlap)[0, 0] ≈ 1 atol = 1e-5
            @test xcorr(signal_rand_complex, -1 .* template_rand_complex_lin, WindowOverlap)[0, 0] ≈ -1 atol = 1e-5

            @test xcorr(signal_rand_complex, template_rand_complex_lin, TrueOverlap)[0, 0] ≈ 1 atol = 1e-5
            @test xcorr(signal_rand_complex, -1 .* template_rand_complex_lin, TrueOverlap)[0, 0] ≈ -1 atol = 1e-5


            # NOTE: Support over full window
            @test xcorr(signal_rand_complex, template_rand_complex_lin, WindowOverlap) ≈ xcorr(signal_rand_complex, template_rand_complex_lin, TrueOverlap) atol = 1e-5 norm = (ab -> maximum(skipmissing(abs.(ab[1] .- ab[2]))))

            @test xcorr(signal_rand_complex, template_rand_complex_lin, WindowOverlap) isa AbstractMatrix{<:Complex}
            @test xcorr(signal_rand_complex, template_rand_complex_lin, TrueOverlap) isa AbstractMatrix{<:Complex}


            signal_rand_complex_odd = centered(rand(ComplexF32, 31, 31))
            template_rand_complex_odd_lin = rand() .* signal_rand_complex_odd

            @test xcorr(signal_rand_complex_odd, template_rand_complex_odd_lin, WindowOverlap)[0, 0] ≈ 1 atol = 1e-5
            @test xcorr(signal_rand_complex_odd, -1 .* template_rand_complex_odd_lin, WindowOverlap)[0, 0] ≈ -1 atol = 1e-5

            @test xcorr(signal_rand_complex_odd, template_rand_complex_odd_lin, TrueOverlap)[0, 0] ≈ 1 atol = 1e-5
            @test xcorr(signal_rand_complex_odd, -1 .* template_rand_complex_odd_lin, TrueOverlap)[0, 0] ≈ -1 atol = 1e-5

            # NOTE: Support over full window
            @test xcorr(signal_rand_complex_odd, template_rand_complex_odd_lin, WindowOverlap) ≈ xcorr(signal_rand_complex_odd, template_rand_complex_odd_lin, TrueOverlap) atol = 1e-5 norm = (ab -> maximum(skipmissing(abs.(ab[1] .- ab[2]))))
        end

        # @testset_skip "Passes... Testing other" "Constant signals" begin
        @testset "Constant signals" begin
            signal_const = centered(ones(30, 30))
            template_const_lin = signal_const .* 0.5

            @test any(isinf.(xcorr(signal_const, template_const_lin, WindowOverlap))) == false
            @test any(isinf.(xcorr(signal_const, template_const_lin, TrueOverlap))) == false

            # NOTE: Bypassing normalization and gives the two constant signals the "correct" variance <06-09-24> 
            @test all(abs.(xcorr(signal_const, template_const_lin, WindowOverlap)) .== 1)
            @test all(abs.(xcorr(signal_const, template_const_lin, TrueOverlap)) .== 1)

            signal_const_negative = centered(fill(0.6, 30, 30))
            template_const_negative = centered(fill(0.3, 30, 30))

            @test xcorr(signal_const_negative, template_const_negative, WindowOverlap) isa AbstractMatrix
            @test xcorr(signal_const_negative, template_const_negative, TrueOverlap) isa AbstractMatrix

            @test (xcorr(signal_const_negative, template_const_negative, WindowOverlap) .== 1) |> all
            @test (xcorr(signal_const_negative, template_const_negative, TrueOverlap) .== 1) |> all
        end

        # @testset_skip "Passes... Testing other" "No Origin" begin
        @testset "No Origin" begin
            signal_nonorigin = rand(30, 30)
            template_nonorigin_lin = rand() .* signal_nonorigin

            # FIX: Handling of signals that do not contain the origin? i.e. index (0,0) is missing <02-09-24> 
            @test_broken xcorr(signal_nonorigin, template_nonorigin_lin, WindowOverlap) isa AbstractMatrix
            # FIX: This gives an output but maybe it is not correct. The issue is that `imfilter` centres the kernel before
            # filtering. The warning given by `imfilter` is also not wanted. <02-09-24> 
            @test xcorr(signal_nonorigin, template_nonorigin_lin, TrueOverlap) isa AbstractMatrix
        end

        @testset "Limited support" begin
            using ShiftedArrays

            signal_radial = centered(rand(100, 100))
            template_radial_lin = rand() .* ShiftedArray(signal_radial, (10, 15))

            for i in CartesianIndices(signal_radial)
                if hypot(Tuple(i)...) > 50
                    signal_radial[i] = 0
                end
            end

            for i in CartesianIndices(template_radial_lin)
                if hypot(Tuple(i)...) > 20
                    template_radial_lin[i] = 0
                end
            end

            @test xcorr(signal_radial, template_radial_lin, TrueOverlap)[-10, -15] ≈ 1


            signal_elliptical_supp = centered(rand(100, 100))
            template_elliptical_supp_lin = rand() .* ShiftedArray(signal_elliptical_supp, (10, 15))

            x_coef = 1.2
            excentric = (x_coef, sqrt(2 - x_coef^2))

            for i in CartesianIndices(signal_elliptical_supp)
                if hypot((Tuple(i) .* excentric)...) > 50 / (maximum(excentric)^2)
                    signal_elliptical_supp[i] = 0
                end
            end

            for i in CartesianIndices(template_elliptical_supp_lin)
                if hypot((Tuple(i) .* reverse(excentric))...) > 20 / (maximum(excentric)^2)
                    template_elliptical_supp_lin[i] = 0
                end
            end

            @test xcorr(signal_elliptical_supp, template_elliptical_supp_lin, TrueOverlap)[-10, -15] ≈ 1

            signal_elliptical_supp_complex = centered(rand(ComplexF32, 100, 100))
            template_elliptical_supp_complex_lin = rand() .* ShiftedArray(signal_elliptical_supp_complex, (10, 15))

            x_coef = 1.2
            excentric = (x_coef, sqrt(2 - x_coef^2))

            for i in CartesianIndices(signal_elliptical_supp_complex)
                if hypot((Tuple(i) .* excentric)...) > 50 / (maximum(excentric)^2)
                    signal_elliptical_supp_complex[i] = zero(ComplexF32)
                end
            end

            for i in CartesianIndices(template_elliptical_supp_complex_lin)
                if hypot((Tuple(i) .* reverse(excentric))...) > 20 / (maximum(excentric)^2)
                    template_elliptical_supp_complex_lin[i] = zero(ComplexF32)
                end
            end

            @test xcorr(signal_elliptical_supp_complex, template_elliptical_supp_complex_lin, TrueOverlap)[-10, -15] ≈ 1 atol = 1e-5
        end
    end
end
