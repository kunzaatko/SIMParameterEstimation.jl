using FFTW, LinearAlgebra
using LazyGrids
using OffsetArrays
using OffsetArrays: centered
using TransferFunctions
using TransferFunctions: otf_support
using FourierTools
using PyCall

@kwdef struct CrossCorrelationPhaseShift <: PE{Harmonic}
    "Minimum relative phase shift in terms of the relative frequencies in the image (∈ [0,1))"
    Δϕr_min::Real = 0.15
    """
    Minimum value of frequency modulation transfer (see [MTF](@ref`TransferFunctions.mtf`)) that is considered in the
    common region of the base and shifted frequency components in the cross-correlation
    """
    MTF_min::Real = 0.05

    function CrossCorrelationPhaseShift(Δϕr_min::Real, MTF_min::Real)
        one(Δϕr_min) > Δϕr_min >= zero(Δϕr_min) || throw(DomainError(Δϕr_min, "Valid domain for Δϕr_min is [0,1)"))
        one(MTF_min) > MTF_min >= zero(MTF_min) || throw(DomainError(MTF_min, "Valid domain for MTF_min is [0,1)"))
        return new(Δϕr_min, MTF_min)
    end
end
const CCPS = CrossCorrelationPhaseShift

function estimate(
    alg::CCPS,
    fft_lr_1::AbstractArray{<:Complex},
    fft_lr_2::AbstractArray{<:Complex},
    tf::TransferFunction;
    Δxy::Length,
    Δϕ_0::Tuple{<:Real,<:Real}, # "Initial guess for the optimizer"
    box_size=Inf
    # algorithm=LBFGS()
)
    @assert size(fft_lr_1) == size(fft_lr_2)
    otf_lr_1 = otf(tf, fft_lr_1, Δxy)

    function xcorr(Δϕ)
        intersection = support_intersect(tf, -1 .* tuple(Δϕ...), size(fft_lr_1), Δxy; MTF_min=alg.MTF_min, Δr_min=alg.Δϕr_min)

        fft_lr_1_i = zeros(eltype(fft_lr_1), size(fft_lr_1)) # at intersection
        fft_lr_1_i[intersection] = fft_lr_1[intersection] ./ otf_lr_1[intersection] # transfer normalized

        fft_lr_2_sΔϕ = shift(fft_lr_2, -1 .* Δϕ) # shifted by -Δϕ (fourier shift theorem)
        fft_lr_2_i = zeros(eltype(fft_lr_2_sΔϕ), size(fft_lr_2_sΔϕ)) # at intersection
        fft_lr_2_i[intersection] = fft_lr_2_sΔϕ[intersection] ./ otf_lr_1[intersection] # same OTF... Imitating the image formation where `lr_img_1` is at the center

        lr_1_i = ifft(fft_lr_1_i) # spatial from intersection
        lr_2_i = ifft(fft_lr_2_i)

        return -abs(sum(lr_2_i .* conj(lr_1_i)) / sum(abs.(lr_1_i) .^ 2))
    end

    optimize = pyimport("scipy.optimize")
    res = optimize.minimize(
        xcorr,
        Δϕ_0;
        bounds=((Δϕ_0[1] .+ (-box_size, +box_size)),
            (Δϕ_0[2] .+ (-box_size, +box_size)))
    )

    if res["success"]
        return res["x"] |> Tuple
    else
        @error "Optimization did not converge: $(res["message"])"
    end
end
