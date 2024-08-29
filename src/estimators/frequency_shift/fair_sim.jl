using TransferFunctions
using FourierTools
using TransferFunctions: overlap

function fairsim_fitpeak(
    band0::AbstractArray,
    band1::AbstractArray,
    otf_band0::SampledOTF,
    otf_band1::SampledOTF,  # NOTE: Does not need to have a defined `center`... It is done in the function <28-08-24> 
    initial::Tuple{Real,Real};
    band0_weight_limit::Real=0.15,
    band1_weight_limit::Real=0.15,
    iterations=3, search=5)
    @assert size(band0) == size(band1) "The components must have the same size"

    band1_center = initial
    for i in Base.OneTo(iterations)
        band1_shifted = shift(band1, band1_center)
        band0_weighted, band1_weighted = filter_under_overlap(band0, band1_shifted, otf_band0, otf_band1, band1_center; band0_weight_limit, band1_weight_limit)
        f_band0 = ifft(band0_weighted)
        f_band1 = ifft(band1_weighted)
        # TODO: `conj` and multiply and powers <29-08-24> 
    end

end

# NOTE: band1 is assumed to be already shifted <29-08-24> 
function filter_under_overlap(
    band0::AbstractArray,
    band1::AbstractArray,
    otf_band0::SampledOTF,
    otf_band1::SampledOTF,
    band1_center::Tuple{Real,Real};
    band0_weight_limit::Real,
    band1_weight_limit::Real
)
    otf_band1_shifted = SampledOTF(otf_band1.transfer, otf_band1.Δxy, band1_center)
    common = overlap(otf_band0, otf_band1_shifted, band0; a_1=band0_weight_limit, a_2=band1_weight_limit)
    out_band0, out_band1 = copy(band0), copy(band1)
    out_band0[common.==false] .= 0
    out_band1[common.==false] .= 0
    out_band0[common.==true] ./= otf(otf_band0, band0)[common]
    out_band1[common.==true] ./= otf(otf_band1_shifted, band1)[common]
    return out_band0, out_band1
end
