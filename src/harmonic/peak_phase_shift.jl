using FFTW, LinearAlgebra
using TransferFunctions
using TransferFunctions: otf_support
struct PeakPhaseShift <: PE{Harmonic}
    "Minimum relative phase shift in terms of OTF support radius (∈ [0,1))"
    Δϕr_min::Real
    function PeakPhaseShift(Δϕr_min::Real)
        one(Δϕr_min) > Δϕr_min >= zero(Δϕr_min) || throw(DomainError(Δϕr_min, "Valid domain for Δϕr_min is [0,1)"))
        return new(Δϕr_min)
    end
end
const PPS = PeakPhaseShift

function estimate(
    alg::PPS,
    fft_lr_1::AbstractArray{<:Complex},
    fft_lr_2::AbstractArray{<:Complex},
    tf::TransferFunction;
    Δxy::Length,
    norm_fact=0.15
)
    lr_1 = ifft(fft_lr_1 ./ (abs.(fft_lr_1) .+ norm_fact))  # abs for the modulation to disappear in correlation
    lr_2 = ifft(fft_lr_2 ./ (abs.(fft_lr_2) .+ norm_fact))

    C = abs.(fft(lr_2 .* conj(lr_1))) # correlation in the frequency domain
    C[otf_support(tf, C, Δxy; ρ=alg.Δϕr_min)] .= 0 # Δϕ outside the OTF support would mean noise

    return argmax(centered(ifftshift(C))) |> Tuple
end
