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

        @test separation_matrix((1.0, 2.0, 3.0), (2.0, 2.0, 2.0)) * transpose([1 1 1; exp(im) exp(2im) exp(3im); exp(-im) exp(-2im) exp(-3im)]) ≈ I(3)
    end

    @testset "xcorr.jl" begin
        using OffsetArrays: Origin, centered

        @testset "SummingSet" begin
            img = ones(1024, 1024)
            img_shifted = Origin(-300, -200)(ones(1024, 1024))

            using SIMParameterEstimation: Window, WindowOverlap, TrueOverlap, Global

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

        @testset "xcorr helpers" begin
            # FIX: The norm_factor should not be run based on the eltype of the denom but based on the eltype of the signal
            # and template... I.e. it should depend on the evaluation of the numerator and not whether the FFT was used or
            # not. <02-09-24> 

            # FIX: The means function should always return the same type Real / Complex as the input... Due to FFT, this is
            # not the case <02-09-24> 
            # Testing correct output types
            using SIMParameterEstimation: Direct, Transform, means, energies, means!, energies!

            @testset "SummingSet $set" for set in (WindowOverlap, TrueOverlap)
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

        # NOTE: Do not accept Real signal with Complex template <02-09-24> 
        @test_throws ArgumentError xcorr(rand(Float32, 30, 30), rand(ComplexF32, 30, 30))

        # NOTE: Linearly dependent signals <02-09-24> 
        signal_rand = centered(rand(30, 30))
        template_rand_lin = rand() .* signal_rand

        @test xcorr(signal_rand, template_rand_lin, WindowOverlap)[0, 0] ≈ 1
        @test xcorr(signal_rand, -1 .* template_rand_lin, WindowOverlap)[0, 0] ≈ -1

        @test xcorr(signal_rand, template_rand_lin, TrueOverlap)[0, 0] ≈ 1
        @test xcorr(signal_rand, -1 .* template_rand_lin, TrueOverlap)[0, 0] ≈ -1

        # NOTE: Support over full window
        @test xcorr(signal_rand, template_rand_lin, WindowOverlap) ≈ xcorr(signal_rand, template_rand_lin, TrueOverlap)

        signal_rand_complex = centered(rand(ComplexF32, 30, 30))
        template_rand_complex_lin = rand() .* signal_rand_complex

        @test xcorr(signal_rand_complex, template_rand_complex_lin, WindowOverlap)[0, 0] ≈ 1
        @test xcorr(signal_rand_complex, -1 .* template_rand_complex_lin, WindowOverlap)[0, 0] ≈ -1

        @test xcorr(signal_rand_complex, template_rand_complex_lin, TrueOverlap)[0, 0] ≈ 1
        @test xcorr(signal_rand_complex, -1 .* template_rand_complex_lin, TrueOverlap)[0, 0] ≈ -1


        # NOTE: Support over full window
        @test xcorr(signal_rand_complex, template_rand_complex_lin, WindowOverlap) ≈ xcorr(signal_rand_complex, template_rand_complex_lin, TrueOverlap)


        signal_rand_complex_odd = centered(rand(ComplexF32, 31, 31))
        template_rand_complex_odd_lin = rand() .* signal_rand_complex_odd

        @test xcorr(signal_rand_complex_odd, template_rand_complex_odd_lin, WindowOverlap)[0, 0] ≈ 1
        @test xcorr(signal_rand_complex_odd, -1 .* template_rand_complex_odd_lin, WindowOverlap)[0, 0] ≈ -1

        @test xcorr(signal_rand_complex_odd, template_rand_complex_odd_lin, TrueOverlap)[0, 0] ≈ 1
        @test xcorr(signal_rand_complex_odd, -1 .* template_rand_complex_odd_lin, TrueOverlap)[0, 0] ≈ -1

        # NOTE: Support over full window
        @test xcorr(signal_rand_complex_odd, template_rand_complex_odd_lin, WindowOverlap) ≈ xcorr(signal_rand_complex_odd, template_rand_complex_odd_lin, TrueOverlap)

        signal_const = centered(ones(30, 30))
        template_const_lin = signal_const .* 0.5

        # FIX: Returns Infs for constant signals <02-09-24> 
        @test_broken any(isinf.(xcorr(signal_const, template_const_lin, WindowOverlap))) == false

        signal_nonorigin = rand(30, 30)
        template_nonorigin_lin = rand() .* signal_nonorigin

        # FIX: Handling of signals that do not contain the origin? i.e. index (0,0) is missing <02-09-24> 
        @test_broken xcorr(signal_nonorigin, template_nonorigin_lin, WindowOverlap) isa AbstractMatrix
        @test xcorr(signal_nonorigin, template_nonorigin_lin, TrueOverlap) isa AbstractMatrix
    end
end
