using TransferFunctions
using FFTW

@kwdef struct CrossCorrelationModulation <: PE{Harmonic}
    """
    Minimum value of frequency modulation transfer (see [MTF](@ref`TransferFunctions.mtf`)) that is considered in the
    common region of the base and shifted frequency components in the cross-correlation
    """
    MTF_min::Real = 0.05
end
const CCM = CrossCorrelationModulation

function estimate(
    alg::CCM,
    fft_lr_1,
    fft_lr_2,
    tf::TransferFunction;
    Δxy::Length,
    Δϕ::Tuple{<:Real,<:Real}
)
    @assert size(fft_lr_1) == size(fft_lr_2)
    otf_lr_1 = otf(tf, fft_lr_1, Δxy)

    intersection = support_intersect(tf, Δϕ, size(fft_lr_1), Δxy; MTF_min=alg.MTF_min)

    fft_lr_1_i = zeros(eltype(fft_lr_1), size(fft_lr_1)) # at intersection
    fft_lr_1_i[intersection] = fft_lr_1[intersection] ./ otf_lr_1[intersection]

    fft_lr_2_sΔϕ = shift(fft_lr_2, -1 .* Δϕ) # shifted by -Δϕ (fourier shift theorem)
    fft_lr_2_i = zeros(eltype(fft_lr_2), size(fft_lr_2))
    fft_lr_2_i[intersection] = fft_lr_2_sΔϕ[intersection] ./ otf_lr_1[intersection]

    lr_1_i = ifft(fft_lr_1_i) # spatial from intersection
    lr_2_i = ifft(fft_lr_2_i)

    return sum(lr_2_i .* conj(lr_1_i)) / sum(abs.(lr_1_i) .^ 2) |> abs
end

