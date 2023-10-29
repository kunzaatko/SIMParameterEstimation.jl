using FFTW, LinearAlgebra
struct PeakFrequency <: PE{Harmonic}
end

function estimate(alg::PeakFrequency, fft_lr)
    f_lr_img = abs.(fft_lr)
    rx, ry = size(fft_lr)
end
