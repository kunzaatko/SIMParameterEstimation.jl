using FFTW, Interpolations, TransferFunctions, LazyGrids
using ImageFiltering
using ImageFiltering: ImageFiltering as ImFilt
using TransferFunctions: TransferFunction, otf_support
using Unitful: Length
using OffsetArrays: centered, Origin
using PaddedViews
using LinearAlgebra

# TODO: Add edge-taper 
# + MATLAB https://www.mathworks.com/help/images/ref/edgetaper.html 
# + CV https://docs.opencv.org/3.4/d1/dfd/tutorial_motion_deblur_filter.html
# + Python https://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/crosscorr.html
# this should be a <type edge-tapering strategy> <11-06-23> 

# TODO: Probably rename or place algorithms for different parameters in separate modules<10-06-23> 
struct AutoCorrelationThroughInputInterpolation <: PE{Harmonic}
    p_θ_0::Tuple{<:Real,<:Real}
end
const ACTII = AutoCorrelationThroughInputInterpolation

# FIX: This currently does not work because of fft bordereffects <17-07-23> 

function estimate(
    alg::ACTII,
    fft_lr::AbstractArray{<:Complex},
    tf::TransferFunction;
    Δxy::Length,
    interpmode=BSpline(Cubic())
    # solver=HiGHS.Optimizer
)
    C_θ = centered(fftshift(fft_lr)) .* centered(conj(fftshift(otf(tf, fft_lr, Δxy)))) .|> ComplexF32 # FIX: type instability on input
    C_θ_ci = interpolate(conj(C_θ), interpmode) # conjugate, interpolated

    # m = Model(solver)
    # @variable(m, 0 <= p_θ[i=1:2] <= [lastindex(C_θ, 1), lastindex(C_θ, 2)][i])
    function acorr(p_θ)
        xs = intersect(axes(C_θ, 1), splat(range)(extrema(axes(C_θ, 1)) .- ceil(p_θ[1])))
        ys = intersect(axes(C_θ, 2), splat(range)(extrema(axes(C_θ, 2)) .- ceil(p_θ[2])))
        Xs, Ys = ndgrid(xs, ys)
        C_θ_cs = C_θ_ci.(Xs .+ p_θ[1], Ys .+ p_θ[2]) # conjugate, shifted
        C_θ_csd = C_θ_cs .- mean(C_θ_cs) # conjugate, shifted, demeaned
        C_θ_d = C_θ[Int64.(xs), Int64.(ys)] .- mean(C_θ[Int64.(xs), Int64.(ys)]) # demeaned
        ccov = C_θ_d .* C_θ_csd |> sum
        return ccov / (norm(C_θ_d) * norm(C_θ_csd)) |> abs
    end
    # FIX: Does not give autocorrelations between -1 and 1 <15-07-23> 
    acorr(alg.p_θ_0)
    # TODO:
    # @objective(m, Max, ccorr(p_θ))
end

struct AutoCorrelationThroughOutputInterpolation <: PE{Harmonic}
    p_θ_0::Tuple{<:Real,<:Real}
end
const ACTOI = AutoCorrelationThroughOutputInterpolation

function estimate(
    alg::ACTOI,
    fft_lr::AbstractArray{<:Complex},
    tf::TransferFunction;
    Δxy::Length,
    interpmode=BSpline(Cubic()),
    solver=HiGHS.Optimizer
)
    C_θ = centered(fftshift(fft_lr)) .* centered(conj(fftshift(otf(tf, fft_lr, Δxy)))) .|> ComplexF32 # FIX: type instability on input

    C_θ_d = C_θ .- mean(C_θ) # demeaned
    C_θ_dc = conj(C_θ_d)

    acov = imfilter(C_θ_d, C_θ_dc, ImFilt.Fill(0), ImFilt.FIRTiled()) # FIX: FFT does not work since we are using ComplexF32

    C_θ_d_norm = sqrt.(imfilter(C_θ_d .^ 2, centered(Ones(C_θ_d)), ImFilt.Fill(0), ImFilt.FIRTiled()))
    C_θ_dc_norm = sqrt.(imfilter(centered(Ones(C_θ_dc)), centered(C_θ_dc .^ 2), ImFilt.Fill(0), ImFilt.FIRTiled()))

    return acov ./ (C_θ_d_norm .* C_θ_dc_norm)

    # cor(p_θ::Tuple{Integer,Integer}) = begin
    #     C_θ_sc = Origin(C_θ.offsets .- p_θ .+ 1)(conj(C_θ)) # shifted, conjugate
    #     # C_θ_p, C_θ_scp = paddedviews(0, C_θ, C_θ_sc) # padded | shifted, conjugate, padded
    #     xs, ys = intersect(axes(C_θ, 1), axes(C_θ_sc, 1)), intersect(axes(C_θ, 2), axes(C_θ_sc, 2))
    #     sum(C_θ[xs, ys] .* C_θ_sc[xs, ys]) ./ (norm(C_θ[xs, ys]) * norm(C_θ_sc[xs, ys])) |> abs
    # end
    # cor(alg.p_θ_0)

    # TODO:
end
