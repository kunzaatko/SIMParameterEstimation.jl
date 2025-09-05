using ImageFiltering: AbstractBorder, Fill, imfilter!, padarray
using IntervalSets, IntegralArrays, Statistics, ShiftedArrays

# TODO: Add function Covariance `xcov` that does not normalize... This will be useful for the modulation determination
# <06-09-24> 
# TODO: When `xcov` is written, should make some figures to demonstrate why we need to normalize the signal in the supp... <06-09-24> 


"""
    Window(arr::AbstractArray)

Window of indices of the array `arr`, representing a Cartesian product of the contained index intervals. 

Useful for computing the overlap of two signals during cross-correlation.


```jldoctest xcorr
julia> signal = [0 0 0 0 0 0 0 0 0;
                 0 0 0 1 1 1 0 0 0;
                 0 0 1 1 1 1 1 0 0;
                 0 1 1 1 1 1 1 1 0;
                 0 1 1 1 1 1 1 1 0;
                 0 0 1 1 1 1 1 0 0;
                 0 0 0 1 1 1 0 0 0;
                 0 0 0 0 0 0 0 0 0];

julia> signal = OAs.Origin(-4,-4)(signal)
8×9 OffsetArray(::Matrix{Int64}, -4:3, -4:4) with eltype Int64 with indices -4:3×-4:4:
 0  0  0  0  0  0  0  0  0
 0  0  0  1  1  1  0  0  0
 0  0  1  1  1  1  1  0  0
 0  1  1  1  1  1  1  1  0
 0  1  1  1  1  1  1  1  0
 0  0  1  1  1  1  1  0  0
 0  0  0  1  1  1  0  0  0
 0  0  0  0  0  0  0  0  0

julia> w_signal = SIM_PE.Window(signal)
Window{2}(CartesianIndices((OffsetArrays.IdOffsetRange(values=-4:3, indices=-4:3), OffsetArrays.IdOffsetRange(values=-4:4, indices=-4:4))))

julia> template = [0 1 0;
                   1 1 1;
                   0 1 0];

julia> template = OffsetArrays.centered(template)
3×3 OffsetArray(::Matrix{Int64}, -1:1, -1:1) with eltype Int64 with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> w_template = Window(template)
Window{2}(CartesianIndices((OffsetArrays.IdOffsetRange(values=-1:1, indices=-1:1), OffsetArrays.IdOffsetRange(values=-1:1, indices=-1:1))))

julia> @assert length(w_template) == 9

julia> intersect(w_signal, w_template .+ CartesianIndex((3, 3)))
Window{2}(CartesianIndices((2:3, 2:4)))
```
"""
struct Window{N}
    inds::CartesianIndices{N}
end
Window(A::AbstractArray) = Window(CartesianIndices(A))
Base.Indices(w::Window) = Tuple(first(d):last(d) for d in w.inds.indices)
Base.broadcasted(::typeof(+), w::Window{N}, ind::CartesianIndex{N}) where {N} = Window(w.inds .+ ind)
Base.broadcasted(::typeof(-), w::Window{N}, ind::CartesianIndex{N}) where {N} = Window(w.inds .- ind)
Base.intersect(w::Window{N}, ws::Window{N}...) where {N} = Window(intersect(w.inds, getfield.(ws, :inds)...))
Base.length(w::Window{N}) where {N} = length(w.inds)

Base.@propagate_inbounds function Base.getindex(A::IntegralArray{T,N}, w::Window{N}) where {T,N}
    A[Tuple(ClosedInterval(first(d), last(d)) for d in w.inds.indices)...]
end

abstract type SummingSet{N} end
struct WindowOverlap{N} <: SummingSet{N}
    signal_win::Window{N}
    template_win::Window{N}
end
function WindowOverlap(signal::AbstractArray{TS,N}, template::AbstractArray{TT,N}) where {N,TS,TT}
    return WindowOverlap(Window(signal), Window(template))
end

# PERF: Is this OK to leave it an abstract type instead of specializing over the array types? <15-08-24> 
struct TrueOverlap{N} <: SummingSet{N}
    # NOTE: It can not be a BitArray, since we need to accept OffsetArrays <02-09-24> 
    signal_supp::AbstractArray{Bool,N}
    template_supp::AbstractArray{Bool,N}
end
# TrueOverlap(signal_supp::, template_supp::BitArray) = TrueOverlap(signal_supp, template_supp)
TrueOverlap(signal::AbstractArray{ST,N}, template::AbstractArray{TT,N}) where {ST,TT,N} = TrueOverlap(signal .!= zero(ST), template .!= zero(TT))

struct Global{N} <: SummingSet{N} end
Global(::AbstractArray{ST,N}, ::AbstractArray{TT,N}) where {ST,TT,N} = Global{N}()

# TODO: `means!` without the count and `means` could be defined as a single method on the abstract types. A similar
# thing can be done with `overlap_length` <15-08-24> 
overlap_length!(out::AbstractArray, A_win::Window, B_win::Window) = map!(ind -> length(intersect(A_win, B_win .+ ind)), out, A_win.inds)
overlap_length(A_win::Window, B_win::Window) = overlap_length!(similar(A_win.inds, Int), A_win, B_win)
overlap_length!(out::AbstractArray, set::WindowOverlap) = overlap_length!(out, set.signal_win, set.template_win)
overlap_length(set::WindowOverlap) = overlap_length!(similar(set.signal_win.inds, Float32), set)

function means!(
    signal_mean::AbstractArray,
    template_mean::AbstractArray,
    signal::AbstractArray,
    template::AbstractArray,
    set::WindowOverlap
)
    counts = overlap_length(set.signal_win, set.template_win)
    return means!(signal_mean, template_mean, signal, template, counts, set)
end


means!(
    signal_mean,
    template_mean,
    signal,
    template,
    counts,
    set::Type{<:SummingSet}
) = means!(signal_mean, template_mean, signal, template, counts, set(signal, template))

function means!(
    signal_mean::AbstractArray,
    template_mean::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray{TT},
    counts::AbstractArray{Int},
    set::WindowOverlap
) where {ST,TT}
    signal_p = padarray(ST, signal, Fill(zero(ST), Base.Indices(set.template_win)))
    template_p = padarray(TT, template, Fill(zero(TT), Base.Indices(set.signal_win)))
    signal_iA = IntegralArray(signal_p)
    template_iA = IntegralArray(template_p)

    @inbounds @simd for uv in CartesianIndices(signal)
        signal_mean[uv] = signal_iA[set.template_win.+uv]
        template_mean[uv] = template_iA[set.signal_win.-uv]
    end
    @inbounds signal_mean ./= counts
    @inbounds template_mean ./= counts
    return signal_mean, template_mean
end

"""
    mean!(out::AbstractArray, A::AbstractArray, set::SummingSet)

Compute the mean of array `A` for each shift  with the specified `set` that depends on the
template for the cross-correlation.

__Warning__: When the template does not contain the origin (`CartesianIndex(0,...)`) in 
    its axes, it may return `NaN`s wherever the shifted template does not overlap with 
    the signal.
"""
means(signal::AbstractArray, template::AbstractArray, set::Type{<:SummingSet}) = means(signal, template, set(signal, template))
means(signal::AbstractArray, template::AbstractArray, set::WindowOverlap) =
    means!(similar(signal), similar(signal), signal, template, set)

energies!(
    signal_energy,
    template_energy,
    signal,
    template,
    set::Type{<:SummingSet},
    args...
) = energies!(signal_energy, template_energy, signal, template, set(signal, template), args...)

# FIX: This is repetition that can be avoided by using a sum function with the `set` argument. Then the `means!` and
# `powers!` functions can be written using this function. <15-08-24> 
function energies!(
    signal_energy::AbstractArray,
    template_energy::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray{TT},
    set::WindowOverlap
) where {ST,TT}
    signal_energy_p = padarray(real(ST), signal .* conj(signal), Fill(zero(real(ST)), Base.Indices(set.template_win)))
    template_energy_p = padarray(real(TT), template .* conj(template), Fill(zero(real(TT)), Base.Indices(set.signal_win)))
    signal_iA = IntegralArray(signal_energy_p)
    template_iA = IntegralArray(template_energy_p)

    @inbounds @simd for uv in CartesianIndices(signal)
        signal_energy[uv] = signal_iA[set.template_win.+uv]
        template_energy[uv] = template_iA[set.signal_win.-uv]
    end
    return signal_energy, template_energy
end

energies(signal::AbstractArray, template::AbstractArray, set::Type{<:SummingSet}) = energies(signal, template, set(signal, template))
energies(signal::AbstractArray{ST}, template::AbstractArray{TT}, set::WindowOverlap) where {ST,TT} =
    energies!(similar(signal, real(ST)), similar(signal, real(TT)), signal, template, set)

# NOTE: There are quite large numerical errors that lead to the `Int` Fourier transform to give `InexactError`s. This is
# not ideal but it should work correctly as the numerical errors are never as large to tip the rounding. <15-08-24> 
overlap_length!(out::AbstractArray, A_supp::AbstractArray{Bool}, B_supp::AbstractArray{Bool}) = round.(Int, imfilter!(out, A_supp, B_supp, Fill(0)))
overlap_length(A_supp::AbstractArray{Bool}, B_supp::AbstractArray{Bool}) = overlap_length!(similar(A_supp, Float32), A_supp, B_supp)
overlap_length!(out::AbstractArray, set::TrueOverlap) = overlap_length!(out, set.signal_supp, set.template_supp)
overlap_length(set::TrueOverlap) = overlap_length!(similar(set.signal_supp, Float32), set)

function means!(
    signal_mean::AbstractArray,
    template_mean::AbstractArray,
    signal::AbstractArray,
    template::AbstractArray,
    set::TrueOverlap
)
    counts = overlap_length(set.signal_supp, set.template_supp)
    return means!(signal_mean, template_mean, signal, template, counts, set)
end

function means!(
    signal_mean::AbstractArray,
    template_mean::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray{TT},
    counts::AbstractArray{Int},
    set::TrueOverlap
) where {ST,TT}
    imfilter!(signal_mean, signal, ST.(set.template_supp), Fill(zero(ST)))
    imfilter!(template_mean, TT.(set.signal_supp), conj(template), Fill(zero(TT)))
    @inbounds signal_mean ./= counts
    @inbounds template_mean ./= counts
    return signal_mean, template_mean
end

means(signal::AbstractArray, template::AbstractArray, set::TrueOverlap) =
    means!(similar(signal), similar(signal), signal, template, set)

abstract type AbstractTrueOverlapSumAlgorithm end
struct Direct <: AbstractTrueOverlapSumAlgorithm end
struct Transform <: AbstractTrueOverlapSumAlgorithm end

function energies!(
    signal_energy::AbstractArray,
    template_energy::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray{TT},
    set::TrueOverlap,
    ::Transform
) where {ST,TT}
    # TODO: Test whether this template_supp doesn't have to be reflected before, i.e. whether this is the correlation
    # rather than convolution as it is supposed to be. <19-08-24> 
    imfilter!(signal_energy, real(signal .* conj(signal)), real(ST).(set.template_supp), Fill(zero(real(ST))))
    imfilter!(template_energy, real(TT).(set.signal_supp), real(template .* conj(template)), Fill(zero(real(TT))))
    return signal_energy, template_energy
end

# PERF: This can be approached in multiple ways, we can directly sum the intersection of the `signal_supp` and
# `template_supp` or a similar method to the window approach can be used. This would lead to better performance if the
# areas under the intersection are large, hence the border is much smaller in terms of pixel counts.  <15-08-24> 
function energies!(
    signal_energy::AbstractArray,
    template_energy::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray{TT},
    set::TrueOverlap,
    ::Direct
) where {ST,TT}
    signal_supp_p = padarray(Bool, set.signal_supp, Fill(false, set.template_supp))
    template_supp_p = padarray(Bool, set.template_supp, Fill(false, set.signal_supp))
    signal_energy_p = padarray(real(ST), signal .* conj(signal), Fill(zero(ST), set.template_supp))
    template_energy_p = padarray(real(TT), template .* conj(template), Fill(zero(ST), set.signal_supp))
    @inbounds @simd for uv in CartesianIndices(signal)
        signal_energy[uv] = sum(signal_energy_p[signal_supp_p.*ShiftedArray(template_supp_p, Tuple(uv); default=false)])
        template_energy[uv] = sum(template_energy_p[signal_supp_p.*ShiftedArray(template_supp_p, Tuple(-uv); default=false)])
    end
    return signal_energy, template_energy
end

energies!(signal_energy::AbstractArray, template_energy::AbstractArray, signal::AbstractArray, template::AbstractArray, set::TrueOverlap) =
    energies!(signal_energy, template_energy, signal, template, set, Transform())
energies(signal::AbstractArray{ST}, template::AbstractArray{TT}, set::TrueOverlap) where {ST,TT} =
    energies!(similar(signal, real(ST)), similar(signal, real(TT)), signal, template, set)

overlap_length!(out::AbstractArray, ::Global) = fill!(out, length(out))
function means!(
    signal_mean::AbstractArray,
    template_mean::AbstractArray,
    signal::AbstractArray,
    template::AbstractArray,
    ::Global
)
    signal_mean .= mean(signal)
    template_mean .= mean(template)
    return signal_mean, template_mean
end

# TODO: There should be an argument min_overlap such that it gives the user an option to set to `missing` if the value
# were to be too untrustworthy <02-09-24> 
# TODO: It doesn't make any sense to calculate the `xcorr` for shifts where there is no overlap. These should instead be
# set to missing... <22-08-24> 
# TODO: There should also be a `range` argument, which determines on which indices the cross-correlation is done. Even
# better, there should be an indices/shifts argument that determines where to calculate the xcorr <10-08-24> 

# NOTE: Step 1: Create the summing set for the demeaning and normalization
xcorr!(
    out::AbstractArray,
    signal::AbstractArray,
    template::AbstractArray,
    normset::Type{<:SummingSet}=WindowOverlap;
    vargs...
) = xcorr!(out, signal, template, normset(signal, template); vargs...)
function xcorr!(
    out::AbstractArray,
    signal::AbstractArray{ST,N},
    template::AbstractArray{TT,N},
    normset::SummingSet{N};
    variance_eps=max(eps(real(ST)), eps(real(TT))),
    deviation_eps=sqrt(variance_eps) # TODO: Document this <06-09-24> 
) where {ST,TT,N}
    # NOTE: Step 1: Throw when the signal is Real and the template is Complex  
    ST <: Real && TT <: Complex && throw(ArgumentError("""
    If signal is `Real`, template must also be `Real`.\nhint: If this is intensional and you want to perform complex\
    cross-correlation, convert the signal to complex using `complex.(signal)`."""))

    template = conj(template) # FIX: Is this a correct interpretation in the statistical sense?

    counts = overlap_length(normset)
    signal_mean, template_mean = means(signal, template, normset)

    # FIX: Is this the definition that I changed in the developed package? <02-09-24> 
    imfilter!(out, signal, conj(template), Fill(zero(ST)))
    out .-= (counts .* signal_mean .* template_mean) # Demeaning

    signal_energy, template_energy = energies(signal, template, normset)

    # TODO: Should be compacted... Possibly by defining a function to do this <06-09-24> 
    # PERF: If the count is 1, then we already know that the input is constant and the energy of the demeaned input will
    # be 0 <06-09-24> 
    # NOTE: Passes some numerical instabilities of the means algorithms leading to negative energies and therefore
    # errors for the square root <06-09-24> 
    # NOTE: Explicit call of `real` ensures that the denominator is a `AbstractArray{<:Real}` <06-09-24> 
    signal_denom = map(counts, signal_energy, signal_mean) do c, e, m
        sd = c > 1 ? e - (c * real(m * conj(m))) : zero(real(ST))
        sd > deviation_eps ? sd : zero(real(ST))
    end
    template_denom = map(counts, template_energy, template_mean) do c, e, m
        td = c > 1 ? e - (c * real(m * conj(m))) : zero(real(TT))
        td > deviation_eps ? td : zero(real(TT))
    end

    denom = sqrt.(signal_denom .* template_denom)

    # TODO: This should be put into separate function and tested probably <06-09-24> 
    # Normalization
    out = map(out, denom, counts, signal_mean, template_mean) do o, d, c, sm, tm
        if d > deviation_eps  # non-constant signal
            o / d
        elseif c > 0 # constant signal but valid overlap, `sm` and `tm` are the constant signal values
            # TODO: Does this make sense for complex signals? Should be documented anyhow <06-09-24> 
            sign(sm * tm)
        else # no overlap
            missing # TODO: This must be documented <06-09-24> 
        end
    end
    return out
end
xcorr(signal::AbstractArray, template::AbstractArray, normset=WindowOverlap; vargs...) =
    xcorr!(similar(signal), signal, template, normset; vargs...)
