using FFTW, LinearAlgebra
using OffsetArrays
using OffsetArrays: centered
using FourierTools
using FourierTools: shift!
using LazyGrids
using PyCall

struct IPCrossCorrelationPhaseOffset <: PE{Harmonic} end
const IPCCPO = IPCrossCorrelationPhaseOffset

function estimate(
    alg::IPCCPO,
    fft_lr_imgs::AbstractArray{<:Complex};
    Δϕ::Tuple{<:Real,<:Real}
)
    @assert size(fft_lr_imgs, 3) == 3
    components = separate_components(fft_lr_imgs, separation_matrix(range(0, 4π / 3; length=3))) # centred

    components_p = pad_fft_components(components)
    I_0 = real.(ifft(components_p[:, :, 1, :], (1, 2)))
    @sync @distributed for i = 1:3
        # FIX: `shift!` does not work for views...? <19-07-23> 
        components_p[:, :, 2, i] = shift(components_p[:, :, 2, i], -1 .* Δϕ)
    end
    @sync @distributed for i = 1:3
        components_p[:, :, 3, i] = shift(components_p[:, :, 3, i], Δϕ)
    end

    fft_s = dropdims(sum(components_p; dims=3); dims=3)
    I_SIM = real.(ifft(fft_s, (1, 2)))

    Y = dropdims(sum(I_SIM .* I_0; dims=(1, 2)) ./ (sum(I_SIM, dims=(1, 2)) .* sum(I_0, dims=(1, 2))), dims=(1, 2))
    Y .-= mean(Y)

    # TODO: Curve is fit by LS <19-07-23> 
    X = range(0, 4π / 3; length=3)
    optimize = pyimport("scipy.optimize")
    (m, φ), _ = optimize.curve_fit((x, m, φ) -> sin.(x .+ φ) .* m, X, Y, [3std(Y) / √2, X[argmax(Y)]])

    return mod(π - sign(m) * π / 2 - φ, 2π)
end
